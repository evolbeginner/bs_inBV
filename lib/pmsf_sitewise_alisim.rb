#! /usr/bin/env ruby
# frozen_string_literal: true

require "getoptlong"
require "fileutils"
require "tmpdir"
require "parallel"
require "bio"

def usage_and_exit(code = 1)
  $stderr.puts <<~USAGE
    Usage:
      ruby pmsf_sitewise_alisim.rb --tree treefile --sitefreq iqtree.sitefreq [options]

    Required:
      --tree FILE                 ML tree with branch lengths
      --sitefreq|--fs FILE             iqtree.sitefreq file

    Options:
      --nrep N                    Number of replicates (default: 1000)
      --cpu N                     Parallel replicates (default: 4)
      -m STR, --model-core STR    Core model before +F{...} (default: LG+G)
      --precision N               Rounding precision for grouping rows (default: 8)
      --outdir DIR                Output directory (default: param_bs_pmsf)
      --iqtree BIN                IQ-TREE binary (default: iqtree3)
      --force                     If outdir exists, remove and continue
      -h, --help                  Show this help
  USAGE
  exit code
end

opts = {
  iqtree: "iqtree3",
  m: "LG+G",       # PMSF injected via +F{...}
  nrep: 1000,
  cpu: 4,
  outdir: "param_bs_pmsf",
  precision: 8,
  force: false
}
opts[:seed] ||= Random.new_seed

parser = GetoptLong.new(
  ['-t', "--tree", GetoptLong::REQUIRED_ARGUMENT],
  ["--fs", "--sitefreq", GetoptLong::REQUIRED_ARGUMENT],
  ["--nrep", GetoptLong::REQUIRED_ARGUMENT],
  ["--cpu", GetoptLong::REQUIRED_ARGUMENT],
  ["-m", "--model-core", GetoptLong::REQUIRED_ARGUMENT],
  ["--precision", GetoptLong::REQUIRED_ARGUMENT],
  #["--seed", GetoptLong::REQUIRED_ARGUMENT],
  ["--outdir", GetoptLong::REQUIRED_ARGUMENT],
  ["--iqtree", GetoptLong::REQUIRED_ARGUMENT],
  ["--force", GetoptLong::NO_ARGUMENT],
  ["-h", "--help", GetoptLong::NO_ARGUMENT]
)

begin
  parser.each do |opt, val|
    case opt
    when "-t", "--tree"
      opts[:tree] = val
    when "--fs", "--sitefreq"
      opts[:sitefreq] = val
    when "--nrep"
      opts[:nrep] = Integer(val)
    when "--cpu"
      opts[:cpu] = Integer(val)
    when "-m", "--model-core"
      opts[:m] = val
    when "--precision"
      opts[:precision] = Integer(val)
    when "--outdir"
      opts[:outdir] = val
    when "--iqtree"
      opts[:iqtree] = val
    when "--force"
      opts[:force] = true
    when "-h", "--help"
      usage_and_exit(0)
    end
  end
rescue ArgumentError, GetoptLong::InvalidOption, GetoptLong::MissingArgument => e
  $stderr.puts "ERROR: #{e.message}"
  usage_and_exit(1)
end


abort("Missing --tree") unless opts[:tree]
abort("Missing --sitefreq") unless opts[:sitefreq]
abort("Tree not found: #{opts[:tree]}") unless File.file?(opts[:tree])
abort("Sitefreq not found: #{opts[:sitefreq]}") unless File.file?(opts[:sitefreq])
abort("--nrep must be > 0") unless opts[:nrep] > 0
abort("--cpu must be > 0") unless opts[:cpu] > 0
abort("--precision must be >= 0") unless opts[:precision] >= 0

def parse_sitefreq(path, precision:)
  rows = []
  File.foreach(path) do |line|
    line = line.strip
    next if line.empty?

    a = line.split
    abort("Bad sitefreq line: #{line}") unless a.size == 21

    probs = a[1..20].map(&:to_f)
    s = probs.sum
    probs = probs.map { |x| x / s } if (s - 1.0).abs > 1e-8
    key = probs.map { |x| format("%.#{precision}f", x) }.join(",")

    rows << { probs: probs, key: key }
  end
  rows
end

def group_sitefreq(rows)
  groups = {} # key => { probs:, positions:[] }
  rows.each_with_index do |row, i|
    groups[row[:key]] ||= { probs: row[:probs], positions: [] }
    groups[row[:key]][:positions] << i # 0-based site index
  end
  groups.values
end

def read_phylip(path)
  lines = File.readlines(path, chomp: true).reject { |x| x.strip.empty? }
  abort("Empty PHYLIP: #{path}") if lines.empty?

  hdr = lines.shift.split
  abort("Bad PHYLIP header in #{path}") unless hdr.size >= 2
  ntax = hdr[0].to_i
  nchar = hdr[1].to_i

  abort("Bad PHYLIP body in #{path}: expected >= #{ntax} lines, got #{lines.size}") if lines.size < ntax

  taxa = []
  seqs = []

  lines.first(ntax).each do |ln|
    parts = ln.strip.split
    name = parts.shift
    seq  = parts.join
    abort("Bad sequence line: #{ln}") if name.nil? || seq.nil? || seq.empty?
    taxa << name
    seqs << seq
  end

  seqs.each do |s|
    abort("Sequence length mismatch in #{path}: expected #{nchar}, got #{s.length}") unless s.length == nchar
  end

  [taxa, seqs, nchar]
end

def run_cmd!(cmd)
  ok = system(*cmd)
  abort("Command failed:\n#{cmd.join(' ')}") unless ok
end

# outdir behavior: fail if exists, unless --force
if Dir.exist?(opts[:outdir])
  if opts[:force]
    #warn "WARNING: --force set, removing existing directory: #{opts[:outdir]}"
    FileUtils.rm_rf(opts[:outdir])
  else
    abort("Output directory already exists: #{opts[:outdir]}\nUse --force to overwrite.")
  end
end
FileUtils.mkdir_p(opts[:outdir])

rows = parse_sitefreq(opts[:sitefreq], precision: opts[:precision])
nsite = rows.size
groups = group_sitefreq(rows)

puts "Sites: #{nsite}"
puts "Unique PMSF profiles: #{groups.size}"
puts "Replicates: #{opts[:nrep]}"
puts "CPUs: #{opts[:cpu]}"
puts "Model core: #{opts[:m]}"
puts "Outdir: #{opts[:outdir]}"

Parallel.each(1..opts[:nrep], in_processes: opts[:cpu]) do |rep|
  rep_dir = File.join(opts[:outdir], rep.to_s)
  FileUtils.mkdir_p(rep_dir)

  taxa = nil
  final_seqs = nil # array of fixed-length mutable strings

  Dir.mktmpdir("pmsf_rep#{rep}_") do |td|
    groups.each_with_index do |g, gi|
      m = g[:positions].length
      fstr = g[:probs].map { |x| format("%.10f", x) }.join(",")
      model = "#{opts[:m]}+F{#{fstr}}"
      outpre = File.join(td, "grp_#{gi + 1}")
      #seed = opts[:seed] #+ rep * 1_000_003 + gi

      cmd = [
        opts[:iqtree],
        "--alisim", outpre,
        "-t", opts[:tree],
        "-m", model,
        "--length", m.to_s,
        #"--seed", seed.to_s,
        "-nt", "1",
        "-redo",
        "-quiet"
      ]
      p cmd
      run_cmd!(cmd)

      phy = "#{outpre}.phy"
      abort("Missing output: #{phy}") unless File.file?(phy)

      t, seqs_grp, nchar = read_phylip(phy)
      abort("Expected #{m} columns but got #{nchar} in #{phy}") unless nchar == m

      if taxa.nil?
        taxa = t
        final_seqs = Array.new(taxa.size) { "\0" * nsite }
      else
        abort("Taxa mismatch at rep #{rep}, group #{gi + 1}") unless taxa == t
      end

      # scatter group columns back to original site indices
      taxa.size.times do |k|
        s = seqs_grp[k]
        g[:positions].each_with_index do |pos, j|
          final_seqs[k].setbyte(pos, s.getbyte(j))
        end
      end
    end
  end

  # BioRuby sequence objects (light validation/consistency)
  aa_objs = final_seqs.map { |s| Bio::Sequence::AA.new(s) }

  outphy = File.join(rep_dir, "combined.phy")
  File.open(outphy, "w") do |f|
    f.puts "#{taxa.size} #{nsite}"
    taxa.each_with_index do |name, k|
      f.puts format("%-10s %s", name, aa_objs[k].to_s)
    end
  end

  $stderr.puts "done #{rep}/#{opts[:nrep]}"
end

puts "Done: #{opts[:outdir]}/<rep>/combined.phy"

