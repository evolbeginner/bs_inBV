#! /bin/env ruby


############################################################
require 'getoptlong'
require 'parallel'
require 'fileutils'
require 'tmpdir'
require 'find'
require 'colorize'

require_relative 'lib/Dir'
require_relative 'lib/processbar'
require_relative 'lib/do_mcmctree.rb'


############################################################
#RSCRIPT="/usr/bin/Rscript" # otherwise cannot work w/ cl007
RSCRIPT='Rscript'

$VERBOSE = nil
DIR = File.dirname(__FILE__)
$VERBOSE = true
LIB_DIR ||= File.join(DIR, 'lib')

IQTREE = 'iqtree'
NW_STATS = 'nw_stats'
NW_TOPOLOGY = 'nw_topology'
MCMCTREE ||= 'mcmctree'
REORDER_NODE = File.join(DIR, 'reorder_node.rb')
FROM_BS_TO_HESSIAN = File.join(DIR, 'from_bs_to_hessian.R')
FIGTREE2NWK ||= File.expand_path(File.join(LIB_DIR, 'figtree2tree.sh'))
EXTRACT_PARAM_FROM_IQTREE = File.join(DIR, "extract_param_from_iqtree_file.rb")

#SOURCE_IQTREE = File.expand_path("~/software/phylo/iqtree/source/v3.0.1/iqtree3")

############################################################
def help()
  STDERR.puts "help message"
  exit 1
end


def make_argu(argu, value)
  return([argu, value].join(' '))
end


def parse_pmsf(model_argu)
  argu = model_argu.split('+').reject{|i|i=~/pmsf/i}.join('+')
  is_pmsf = model_argu =~ /pmsf/i ? true : false
  return([argu, is_pmsf])
end


def processbar_for_bootstrapping(file:, b:)
  thr = Thread.new do |i|
    count = 0
    while true do
      next if not File.exist?(file)
      new_count = `wc -l #{file} | awk '{print $1}'`.chomp.to_i
      sleep(0.2)
      if new_count != count
        count = new_count
        processbar(count, b)
      end
    end
  end
  return(thr)
end


def get_files_under_folder(folder, b)
  files = Array.new
  Find.find(folder) do |path|
    files << path if path =~ /#{b}$/
  end
  return(files.join(' '))
end


def create_inBV(mltree_file:, mcmctree_outdir:, inBV_file:, iqtree_outdir:)
  no_species = `#{NW_STATS} #{mltree_file} | grep '^#leaves:' | awk '{print $2}'`.chomp.to_i  
  `bash -c "echo -e '\n#{no_species}\n' " >#{inBV_file}`

  `#{NW_TOPOLOGY} -bI #{mltree_file} >> #{inBV_file}`
  `echo >> #{inBV_file}`

  `cat #{iqtree_outdir}/ml.bls >> #{inBV_file}`
  `echo >> #{inBV_file}`

  gradient = [%w[0] * (2*no_species-3)].join(' ') # no. of branches equals 2n-3 where n is the no. of species
  `echo #{gradient} >> #{inBV_file}`
  `echo >> #{inBV_file}`

  `echo Hessian >> #{inBV_file}`
  `echo >> #{inBV_file}`
  `cat #{iqtree_outdir}/hessian >> #{inBV_file}`
end


def get_iter_num(mcmctree_ctl)
  no = nil
  in_fh = File.open(mcmctree_ctl, 'r')
  in_fh.each_line do |line|
    line.chomp!
    no = $1.to_i if line =~ /nsample[ ]*=[ ]*(\d+)/
  end
  in_fh.close
  return(no)
end


def split_ali_file(ali_file)
  in_fh = File.open(ali_file)
  ali2lines = Hash.new
  count = -1
  in_fh.each_line do |line|
    line.chomp!
    if line =~ /\d+\s+\d+$/
      count += 1
      ali2lines[count] = Array.new
    end
    ali2lines[count] << line
  end
  in_fh.close
  return(ali2lines)
end


def output_gapped_seq(outdir_split1, ali_file1, lines)
  out_fh = File.open(ali_file1, 'w')
  out_fh_all_gap_seq = File.open(File.join(outdir_split1, 'all_gap_seq.list'), 'w')
  lines.each_with_index do |line, index|
    line =~ /^(\S+)(\s+)(\S+)$/
    seq_name = $1; space = $2; seq = $3
    if seq =~ /^[-]+$/
      out_fh_all_gap_seq.puts [seq_name].join("\t")
      k = 1
      seq.sub!(/-{#{k}}/, 'A'*k)
      line = seq_name + space + seq
    end
    out_fh.puts line
  end
  out_fh.close
  out_fh_all_gap_seq.close
end


############################################################
def determine_model(model_argu, ali_file1, iqtree_outdir, cpu, bootstrap_argu, te_argu, add_argu, is_best_fit)

  if is_best_fit
    # no bootstrap
    `#{IQTREE} -redo -s #{ali_file1} -pre #{iqtree_outdir}/iqtree -nt #{cpu} -quiet #{model_argu} #{te_argu} #{add_argu}`
    # best-fitting model
    model_best = `ruby #{EXTRACT_PARAM_FROM_IQTREE} --iqtree #{iqtree_outdir}/iqtree.iqtree --log #{iqtree_outdir}/iqtree.log | tail -n 1`.chomp
    model_argu_new = ['-m', model_best].join(' ')
  else
    model_argu_new = model_argu
  end

  return(model_argu_new)
end


def setup_split_dirs(outdir_split, count, lines)
  outdir_split1 = File.join(outdir_split, count.to_s)
  mkdir_with_force(outdir_split1)

  ali_file1 = File.join(outdir_split1, 'combined.phy')
  output_gapped_seq(outdir_split1, ali_file1, lines)

  iqtree_outdir   = File.join(outdir_split1, 'iqtree')
  mcmctree_outdir = File.join(outdir_split1, 'mcmctree')
  mkdir_with_force(iqtree_outdir)
  mkdir_with_force(mcmctree_outdir)

  boottree_file = File.join(iqtree_outdir, 'iqtree.boottrees')
  mltree_file   = File.join(iqtree_outdir, 'iqtree.treefile')
  inBV_file     = File.join(mcmctree_outdir, 'in.BV')

  return([outdir_split1, ali_file1, iqtree_outdir, mcmctree_outdir, boottree_file, mltree_file, inBV_file])
end


############################################################
# using finite diff. by specifying -flmix
# given up already
# bc numerical stability not good in -flmix
def generate_bl_outdir(iqtree_outdir, boottree_file, ref_tree_file, ali_file1, model_argu, add_argu)
  bl_outdir = File.join(iqtree_outdir, 'bl_outdir')
  `#{NW_TOPOLOGY} -Ib #{boottree_file} | ruby #{REORDER_NODE} -i - --ref #{ref_tree_file} --bl_outdir #{bl_outdir} --force`
  p "#{NW_TOPOLOGY} -Ib #{boottree_file} | ruby #{REORDER_NODE} -i - --ref #{ref_tree_file} --bl_outdir #{bl_outdir} --force"
  bl_sub_outdirs = Dir.glob(File.join(bl_outdir, '*'))
  Parallel.each(bl_sub_outdirs, in_threads: cpu) do |bl_sub_outdir|
    bl_outdir2s = Array.new
    1.upto(2) do |i|
      bl_outdir2 = File.join(bl_sub_outdir, i.to_s)
      bl_outdir2s = bl_outdir2s.push(bl_outdir2)
      `#{IQTREE} -redo -s #{ali_file1} -pre #{bl_outdir2}/iqtree -nt 1 -quiet #{model_argu} -te #{bl_outdir2}/bl_changed.treefile #{add_argu} -blfix`
    end
    in_fh = File.open(File.join(bl_outdir2s[0], "step_size"), 'r')
    step = in_fh.readline.to_f
    in_fh.close
    lnL1 = `grep BEST #{bl_outdir2s[0]}/iqtree.log | awk '{print $NF}'`.to_f
    lnL2 = `grep BEST #{bl_outdir2s[1]}/iqtree.log | awk '{print $NF}'`.to_f
    out_fh = File.open(File.join(bl_sub_outdir, 'gradient'), 'w')
    out_fh.puts [lnL1, lnL2, step].flatten.join("\t")
    out_fh.puts fd(lnL1, lnL2, step)
    out_fh.close
  end
end


def fd(lnL1, lnL2, step)
  (lnL2-lnL1)/(2*step)
end


############################################################
mcmctree_indir = nil
mcmctree_ctl_file = nil
ali_file = nil
model_argu = '-m LG+G'
bootstrap = 1000
bootstrap_argu = "-b #{bootstrap}"
te_argu = nil
is_pmsf = false
add_argu = '-mwopt'
calib_tree_file = nil
ref_tree_file = nil

outdir = nil
is_force = false
cpu = 1
is_run_mcmctree = false
is_best_fit = false


############################################################
ARGV_COPY = Marshal.load(Marshal.dump(ARGV))

opts = GetoptLong.new(
  ['--mcmctree_indir', GetoptLong::REQUIRED_ARGUMENT],
  ['--mcmctree_ctl', GetoptLong::REQUIRED_ARGUMENT],
  ['--ali', GetoptLong::REQUIRED_ARGUMENT],
  ['--ref', '--ref_tree_file', GetoptLong::REQUIRED_ARGUMENT],
  ['-m', GetoptLong::REQUIRED_ARGUMENT],
  ['-b', GetoptLong::REQUIRED_ARGUMENT],
  ['--pmsf', GetoptLong::NO_ARGUMENT],
  ['--te', GetoptLong::REQUIRED_ARGUMENT],
  ['--calib_tree', '--calibrated_tree', GetoptLong::REQUIRED_ARGUMENT],
  ['--outdir', GetoptLong::REQUIRED_ARGUMENT],
  ['--force', GetoptLong::NO_ARGUMENT],
  ['--cpu', GetoptLong::REQUIRED_ARGUMENT],
  ['--run_mcmctree', GetoptLong::NO_ARGUMENT],
  ['--no_mcmctree', GetoptLong::NO_ARGUMENT],
  ['--add_cmd', '--add_argu', GetoptLong::REQUIRED_ARGUMENT],
  ['--no_mwopt', GetoptLong::NO_ARGUMENT],
  ['--best_fit', '--best-fit', GetoptLong::NO_ARGUMENT],
)

opts.each do |opt, value|
  case opt
    when '--mcmctree_indir'
      mcmctree_indir = value
    when '--mcmctree_ctl'
      mcmctree_ctl_file = value
    when '--ali'
      ali_file = value
    when '--ref', '--ref_tree_file'
      ref_tree_file = value
      te_argu = make_argu('-te', value)
    when '-m'
      model_argu = make_argu('-m', value)
      model_argu, is_pmsf = parse_pmsf(model_argu)
    when '-b'
      bootstrap = value.to_i
      bootstrap_argu = make_argu('-b', value)
    when '--te'
      te_argu = make_argu('-te', value)
    when '--pmsf'
      is_pmsf = true
    when '--calib_tree', '--calibrated_tree'
      calib_tree_file = value
    when '--outdir'
      outdir = value
    when '--force'
      is_force = true
    when '--cpu'
      cpu = value.to_i
    when '--run_mcmctree'
      is_run_mcmctree = true
    when '-h'
      help()
    when '--no_mcmctree'
      is_run_mcmctree = false
    when '--add_cmd', '--add_argu'
      add_argu = value
    when '--no_mwopt'
      add_argu.sub!('-mwopt', '')
    when '--best_fit', '--best-fit'
      is_best_fit = true
  end
end


############################################################
mkdir_with_force(outdir, is_force)

cmd_out_fh = File.open(File.join(outdir, 'cmd'), 'w')
cmd_out_fh.puts [__FILE__, ARGV_COPY].flatten.join(" ")
cmd_out_fh.close

ali2lines = split_ali_file(ali_file)

outdir_split = File.join(outdir, "split")

mkdir_with_force(outdir_split, is_force)

STDERR.puts "CPU:\t#{cpu}".colorize(:red)


############################################################
ali2lines.to_a.reverse.each do |count, lines|
  outdir_split1, ali_file1, iqtree_outdir, mcmctree_outdir, boottree_file, mltree_file, inBV_file = 
    setup_split_dirs(outdir_split, count, lines) 

  STDOUT.puts "Running IQ-Tree ......"

  # to record how many bootstrap trees are built
  thr = processbar_for_bootstrapping(file:boottree_file, b:bootstrap)

  if is_pmsf
    `#{IQTREE} -redo -s #{ali_file1} -pre #{iqtree_outdir}/guide -nt #{cpu} -quiet -m LG4M+G #{te_argu} #{add_argu}`
    add_argu = [add_argu, "-ft #{iqtree_outdir}/guide.treefile"].join(' ')
    #`#{IQTREE} -redo -s #{ali_file1} -pre #{iqtree_outdir}/iqtree -nt #{cpu} -quiet #{model_argu} #{bootstrap_argu} #{te_argu} #{add_argu} -ft #{iqtree_outdir}/guide.treefile`
  end

  model_argu_new = determine_model(model_argu, ali_file1, iqtree_outdir, cpu, bootstrap_argu, te_argu, add_argu, is_best_fit)
  #model_argu_new = model_argu

  STDERR.puts "Note: pmsf not used for seemingly profile-mixture model #{model_argu}. Take long.".colorize(:blue) if model_argu =~ /C[0-9]+/

  `#{IQTREE} -redo -s #{ali_file1} -pre #{iqtree_outdir}/iqtree -nt #{cpu} -quiet #{model_argu_new} #{bootstrap_argu} #{te_argu} #{add_argu}`
  if File.read(File.join(iqtree_outdir, "iqtree.log")).include?("ERROR: No mixture model was specified!")
    STDERR.puts "pmsf is used for seemingly non-mixture model #{model_argu}. Exiting ......".colorize(:blue)
    exit
  end

  Thread.kill(thr) and puts if $? == 0

  `#{NW_TOPOLOGY} -Ib #{boottree_file} | ruby #{REORDER_NODE} -i - --ref #{ref_tree_file} > #{iqtree_outdir}/boot.bls`

  #generate_bl_outdir(iqtree_outdir)

  # for the ml tree (iqtree.treefile)
  `ruby #{REORDER_NODE} -i #{mltree_file} --ref #{ref_tree_file} > #{iqtree_outdir}/ml.bls`

  `#{RSCRIPT} #{FROM_BS_TO_HESSIAN} #{iqtree_outdir}/boot.bls #{iqtree_outdir}/hessian`

  mltree_file = File.join(iqtree_outdir, 'iqtree.treefile')
  create_inBV(mltree_file:mltree_file, mcmctree_outdir:mcmctree_outdir, inBV_file:inBV_file, iqtree_outdir:iqtree_outdir)
end


############################################################
inBV_files = get_files_under_folder(outdir, 'in.BV')

mcmctree_outdir = File.join(outdir, 'mcmctree')
mkdir_with_force(mcmctree_outdir)
`cat #{inBV_files} > #{mcmctree_outdir}/in.BV`

FileUtils.cp(ali_file, File.join(mcmctree_outdir,'combined.phy'))
`cp #{calib_tree_file} #{mcmctree_outdir}/species.trees`

############################################################
# optional

if ! mcmctree_ctl_file.nil?
  prepare_paml_ctl(mcmctree_ctl_file, mcmctree_outdir, {'seqfile'=>'combined.phy', 'treefile'=>'species.trees'})
end

if is_run_mcmctree
  STDOUT.puts "Running MCMCtree ......"
  Dir.chdir(mcmctree_outdir)

  thr = processbar_for_bootstrapping(file:'mcmc.txt', b:get_iter_num(File.basename(mcmctree_ctl_file)))

  `#{MCMCTREE} #{File.basename(mcmctree_ctl_file)} > mcmctree.final`

  if $? == 0
    Thread.kill(thr) and puts
    `bash #{FIGTREE2NWK} -i FigTree.tre > figtree.nwk`
    `grep rategram out* >/dev/null && grep -A1 rategram out* | tail -1 > rate.tre`
    puts "Done!" if $? == 0
  end
end


