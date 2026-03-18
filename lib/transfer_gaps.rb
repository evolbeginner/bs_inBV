#! /bin/env ruby
# transfer_gaps.rb
# Usage:
#   ruby transfer_gaps.rb gapped.phy ungapped.phy output.phy

abort "Usage: ruby #{$PROGRAM_NAME} <gapped.phy> <ungapped.phy> <output.phy>" if ARGV.size != 3

gapped_file, ungapped_file, out_file = ARGV

def clean_seq(s)
  # keep only sequence-like characters (letters and common gap/missing symbols)
  s.upcase.gsub(/[^A-Z\-\?\.]/, "")
end

def read_phylip(path)
  lines = File.readlines(path, chomp: true)
  abort "ERROR: Empty file: #{path}" if lines.empty?

  header = lines.shift.strip
  m = header.match(/^(\d+)\s+(\d+)/)
  abort "ERROR: Bad PHYLIP header in #{path}: #{header}" unless m

  n = m[1].to_i
  l = m[2].to_i

  # Supports common sequential PHYLIP (relaxed names):
  # taxon_name<space>SEQUENCE
  taxa = []
  seqs = {}
  lines.each do |line|
    next if line.strip.empty?
    name, seq = line.strip.split(/\s+/, 2)
    next if name.nil? || seq.nil?
    taxa << name
    seqs[name] = clean_seq(seq)
  end

  abort "ERROR: #{path} has #{taxa.size} taxa, expected #{n}" unless taxa.size == n

  seqs.each do |name, seq|
    abort "ERROR: #{path} taxon #{name} length #{seq.length}, expected #{l}" unless seq.length == l
  end

  { n: n, l: l, taxa: taxa, seqs: seqs }
end

g1 = read_phylip(gapped_file)
g2 = read_phylip(ungapped_file)

abort "ERROR: Different number of taxa (#{g1[:n]} vs #{g2[:n]})" unless g1[:n] == g2[:n]
abort "ERROR: Different alignment length (#{g1[:l]} vs #{g2[:l]})" unless g1[:l] == g2[:l]

# ensure all taxa in file1 exist in file2
missing = g1[:taxa].reject { |t| g2[:seqs].key?(t) }
abort "ERROR: Missing taxa in #{ungapped_file}: #{missing.join(', ')}" unless missing.empty?

new_seqs = {}

g1[:taxa].each do |taxon|
  s1 = g1[:seqs][taxon] # has gaps
  s2 = g2[:seqs][taxon] # template to receive gaps

  chars = s2.chars
  s1.each_char.with_index do |c, i|
    chars[i] = "-" if c == "-"
  end
  new_seqs[taxon] = chars.join
end

File.open(out_file, "w") do |f|
  f.puts "#{g1[:n]} #{g1[:l]}"
  g1[:taxa].each do |taxon|
    f.puts "#{taxon} #{new_seqs[taxon]}"
  end
end

puts "Done. Wrote: #{out_file}"
