#! /bin/env ruby

require 'getoptlong'

######################################################
# CLI
######################################################
opts = GetoptLong.new(
  ['--log',    '-l', GetoptLong::REQUIRED_ARGUMENT],
  ['--iqtree', '-i', GetoptLong::REQUIRED_ARGUMENT],
  ['--help',   '-h', GetoptLong::NO_ARGUMENT]
)

log_file    = nil
iqtree_file = nil

opts.each do |opt, arg|
  case opt
  when '--log'
    log_file = arg
  when '--iqtree'
    iqtree_file = arg
  when '--help'
    puts <<~HELP
      Usage:
        ruby extract_param_from_iqtree_file.rb \
          --log haha.log \
          --iqtree haha.iqtree
    HELP
    exit
  end
end

abort("ERROR: missing --log")     unless log_file
abort("ERROR: missing --iqtree") unless iqtree_file
abort("ERROR: file not found: #{log_file}")    unless File.exist?(log_file)
abort("ERROR: file not found: #{iqtree_file}") unless File.exist?(iqtree_file)

######################################################
# Step 1: parse haha.log (alias line, FMIX base model)
######################################################
base_model = nil

File.foreach(log_file) do |line|
  if line =~ /is alias for\s+(.+FMIX\{.+\})/
    full_model = Regexp.last_match(1)
    base_model, = full_model.split('+FMIX{', 2)
    break
  end
end

######################################################
# Step 1b: fallback — parse base model from iqtree
######################################################
if base_model.nil?
  File.foreach(iqtree_file) do |line|
    if line =~ /Model of substitution:\s+(.+)/
      # extract only the exchangeability matrix (e.g. LG from LG+G4)
      base_model = Regexp.last_match(1).split('+').first
      break
    end
  end
end

abort("ERROR: could not determine base substitution model") unless base_model

######################################################
# Step 2: parse haha.iqtree (FMIX mixture table)
######################################################
components = []
rates      = []
weights    = []

in_table = false

File.foreach(iqtree_file) do |line|
  if line =~ /^\s*No\s+Component\s+Rate\s+Weight\s+Parameters/
    in_table = true
    next
  end

  next unless in_table

  if line =~ /^\s*\d+\s+(\S+)\s+([\d.eE+-]+)\s+([\d.eE+-]+)/
    component_full = Regexp.last_match(1)
    rate           = Regexp.last_match(2).to_f
    weight         = Regexp.last_match(3).to_f

    if component_full =~ /(FC\d+pi\d+)$/
      components << Regexp.last_match(1)
      rates      << rate
      weights    << (weight < 1e-3 ? 1e-4 : weight)
    end
  elsif line.strip.empty?
    in_table = false
  end
end

######################################################
# Step 3: rebuild FMIX
######################################################
fmix_model =
  if components.empty?
    base_model
  else
    fmix_entries = components.each_with_index.map do |comp, i|
      format("%s:%.6f:%.6f", comp.sub("FC", "C"), rates[i], weights[i])
    end
    "#{base_model}+FMIX{#{fmix_entries.join(',')}}"
  end

######################################################
# Step 4: parse rate heterogeneity
######################################################
r_k       = nil
r_weights = []
r_rates   = []

i_prop  = nil
g_k     = nil
g_alpha = nil

File.foreach(iqtree_file) do |line|
  # FreeRate
  if line =~ /FreeRate with (\d+) categories/
    r_k = Regexp.last_match(1).to_i
  end

  if line =~ /Site proportion and rates:\s+(.+)/
    Regexp.last_match(1).scan(/\(([\d.eE+-]+),([\d.eE+-]+)\)/) do |w, r|
      r_weights << w.to_f
      r_rates   << r.to_f
    end
  end

  # Invariable sites
  if line =~ /Proportion of invariable sites:\s+([\d.eE+-]+)/
    i_prop = Regexp.last_match(1).to_f
  end

  # Gamma
  if line =~ /Gamma shape alpha:\s+([\d.eE+-]+)/
    g_alpha = Regexp.last_match(1).to_f
  end

  if line =~ /(Invar\+)?Gamma with (\d+) categories/
    g_k = Regexp.last_match(2).to_i
  end
end

######################################################
# Step 5: rebuild rate heterogeneity model (ADDITIVE)
######################################################
rate_model = ""

# +R
if r_k
  abort("ERROR: FreeRate parsing failed") unless r_weights.size == r_k &&
                                                 r_rates.size   == r_k

  r_entries = r_weights.zip(r_rates).flat_map do |w, r|
    [format("%.6f", w), format("%.6f", r)]
  end

  rate_model << "+R#{r_k}{#{r_entries.join(',')}}"
end

# +I
if i_prop
  rate_model << format("+I{%.6f}", i_prop)
end

# +G
if g_alpha && g_k
  rate_model << format("+G%d{%.6f}", g_k, g_alpha)
end

######################################################
# Output
######################################################
puts "Updated IQ-TREE model:"
puts "#{fmix_model}#{rate_model}"

