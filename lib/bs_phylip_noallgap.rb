#! /bin/env ruby
# frozen_string_literal: true

require "getoptlong"
require "fileutils"
require "set"

def usage
  <<~USAGE
    Usage:
      ruby bs_phylip.rb [options] -i input.phy

    Options:
      -b, --bs N         Number of bootstrap replicates (default 1)
      -c, --chars N      Number of characters/sites in each replicate (default: original nsite)
      -o, --outdir DIR   Output directory (default: bs_outdir)
      -s, --seed N       Random seed (optional)
      -h, --help         Show this help
  USAGE
end

# ---------------------------
# Parse options
# ---------------------------
infile = nil
n_bs = 1
n_chars = nil
outdir = "bs_outdir"
seed = nil

opts = GetoptLong.new(
  ["--in",     "-i", GetoptLong::REQUIRED_ARGUMENT],
  ["--bs",     "-b", GetoptLong::REQUIRED_ARGUMENT],
  ["--chars",  "-c", GetoptLong::REQUIRED_ARGUMENT],
  ["--outdir", "-o", GetoptLong::REQUIRED_ARGUMENT],
  ["--seed",   "-s", GetoptLong::REQUIRED_ARGUMENT],
  ["--help",   "-h", GetoptLong::NO_ARGUMENT]
)

begin
  opts.each do |opt, arg|
    case opt
    when "--in"
      infile = arg
    when "--bs"
      n_bs = Integer(arg)
      abort("--bs must be >= 1") if n_bs < 1
    when "--chars"
      n_chars = Integer(arg)
      abort("--chars must be >= 1") if n_chars < 1
    when "--outdir"
      outdir = arg
    when "--seed"
      seed = Integer(arg)
    when "--help"
      puts usage
      exit 0
    end
  end
rescue ArgumentError => e
  abort("option parse error: #{e.message}\n\n#{usage}")
end

abort("missing input.phy\n\n#{usage}") unless infile

rng = seed ? Random.new(seed) : Random.new

# ---------------------------
# Read PHYLIP (simple sequential format)
# ---------------------------
raw = File.readlines(infile, chomp: true).map(&:strip).reject(&:empty?)
abort("empty input file") if raw.empty?

hdr = raw.shift.split
abort("bad PHYLIP header") unless hdr.size >= 2

ntax = Integer(hdr[0]) rescue abort("bad ntax in header")
nsite = Integer(hdr[1]) rescue abort("bad nsite in header")

abort("expected at least #{ntax} sequence lines, got #{raw.size}") if raw.size < ntax

names = []
seqs = []

raw.first(ntax).each_with_index do |line, idx|
  parts = line.split
  abort("bad sequence line #{idx + 2}") if parts.size < 2
  name = parts.shift
  seq = parts.join
  names << name
  seqs << seq
end

abort("sequence length mismatch vs header nsite=#{nsite}") unless seqs.all? { |s| s.length == nsite }

n_chars ||= nsite

# ---------------------------
# Precompute valid columns and coverage
# valid_cols: columns that are not all-gap across taxa
# col_covers[col]: taxa indices with non-gap at col
# ---------------------------
valid_cols = []
col_covers = {}

(0...nsite).each do |col|
  covered = []
  ntax.times do |t|
    covered << t if seqs[t].getbyte(col) != 45 # '-' == 45
  end
  next if covered.empty?
  valid_cols << col
  col_covers[col] = covered
end

abort("no valid columns (all columns are all-gap)") if valid_cols.empty?

# Ensure each taxon has at least one non-gap in valid columns
has_support = Array.new(ntax, false)
valid_cols.each { |col| col_covers[col].each { |t| has_support[t] = true } }
unless has_support.all?
  bad = has_support.each_index.select { |i| !has_support[i] }.map { |i| names[i] }
  abort("impossible input: taxa with no non-gap in any valid column: #{bad.join(', ')}")
end

# ---------------------------
# Helpers
# ---------------------------
def build_mandatory_cover(valid_cols, col_covers, ntax, rng)
  uncovered = Set.new((0...ntax).to_a)
  mandatory = []

  while uncovered.any?
    best_gain = -1
    best_cols = []

    valid_cols.each do |col|
      gain = 0
      col_covers[col].each { |t| gain += 1 if uncovered.include?(t) }
      next if gain <= 0

      if gain > best_gain
        best_gain = gain
        best_cols = [col]
      elsif gain == best_gain
        best_cols << col
      end
    end

    raise "failed to cover all taxa" if best_gain <= 0

    chosen = best_cols.sample(random: rng)
    mandatory << chosen
    col_covers[chosen].each { |t| uncovered.delete(t) }
  end

  mandatory
end

def make_bootstrap_columns_rejection(n_chars, valid_cols, rng)
  Array.new(n_chars) { valid_cols.sample(random: rng) }
end

def make_bootstrap_columns_greedy(n_chars, mandatory_cols, valid_cols, rng)
  cols = mandatory_cols.dup
  cols << valid_cols.sample(random: rng) while cols.length < n_chars
  cols.shuffle!(random: rng)
  cols
end

def build_bootstrap_sequences(seqs, bs_cols)
  ntax = seqs.length
  out = Array.new(ntax) { +"" }
  bs_cols.each do |col|
    ntax.times { |t| out[t] << seqs[t].getbyte(col) }
  end
  out
end

def all_gap?(s)
  s.each_byte.all? { |b| b == 45 } # '-'
end

def all_gap_taxa_indices(bs_seqs)
  bad = []
  bs_seqs.each_with_index { |s, i| bad << i if all_gap?(s) }
  bad
end

# Precompute greedy mandatory cover once (used only when fallback needed)
mandatory_cols = build_mandatory_cover(valid_cols, col_covers, ntax, rng)
if mandatory_cols.length > n_chars
  abort("cannot enforce no-all-gap-taxon with --chars=#{n_chars}; need at least #{mandatory_cols.length}")
end

# ---------------------------
# Generate all replicates
# ---------------------------
(1..n_bs).each do |rep|
  rep_dir = File.join(outdir, rep.to_s)
  FileUtils.mkdir_p(rep_dir)
  outfile = File.join(rep_dir, "combined.phy")

  bs_cols = nil
  bs_seqs = nil
  method_used = nil

  # 1) Try rejection approach up to 2 failed attempts
  rejection_failures = 0
  2.times do
    trial_cols = make_bootstrap_columns_rejection(n_chars, valid_cols, rng)
    trial_seqs = build_bootstrap_sequences(seqs, trial_cols)
    bad = all_gap_taxa_indices(trial_seqs)

    if bad.empty?
      bs_cols = trial_cols
      bs_seqs = trial_seqs
      method_used = :rejection
      break
    else
      rejection_failures += 1
    end
  end

  # 2) Fallback to greedy if rejection failed twice
  if bs_seqs.nil?
    bs_cols = make_bootstrap_columns_greedy(n_chars, mandatory_cols, valid_cols, rng)
    bs_seqs = build_bootstrap_sequences(seqs, bs_cols)
    method_used = :greedy_fallback
  end

  # Final sanity check
  bad = all_gap_taxa_indices(bs_seqs)
  unless bad.empty?
    bad_names = bad.map { |i| names[i] }
    abort("internal error: generated all-gap taxa in replicate #{rep}: #{bad_names.join(', ')}")
  end

  File.open(outfile, "w") do |f|
    f.puts "#{ntax} #{n_chars}"
    ntax.times do |i|
      f.puts format("%-10s %s", names[i], bs_seqs[i])
    end
  end

  if method_used == :greedy_fallback
    warn "replicate #{rep}: rejection failed twice; used greedy fallback"
  end
end

puts "Done. Wrote #{n_bs} replicates under: #{outdir}"

