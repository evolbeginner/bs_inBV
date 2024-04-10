#! /usr/bin/env ruby


###########################################################
is_to_install = false


###########################################################
# check ruby version
puts "Checking Ruby version ......"
if RUBY_VERSION =~ /^1\. | ^2\.[0-6]/x
  puts "Your Ruby version is #{RUBY_VERSION}. Pls install Ruby v2.7+."
  is_to_install = true
end
puts


###########################################################
# check ruby gems
puts "Checking Ruby gems ......"
required_gems = %w[parallel colorize find getoptlong time fileutils tmpdir bio bio-nwk]

required_gems.each do |gem_name|
  if Gem::Specification.find_all_by_name(gem_name).empty?
    puts "#{gem_name} is NOT installed. In bash, try \'gem install #{gem_name}\'"
    is_to_install = true
  else
    #puts "#{gem_name} is installed."
    ;
  end
end

puts


###########################################################
# check other tools
puts "Checking other tools ......"

required_tools = {"iqtree"=>"iqtree",
  "nw_topology"=>"Newick Utilities", "mcmctree"=>"MCMCtree (PAML)", "Rscript"=>"R"}

required_tools.each_pair do |k, v|
  if not system("which #{k} > /dev/null 2>&1")
    puts "#{v} NOT installed! Pls install #{v} and add its path to the environmental variable PATH (e.g., https://askubuntu.com/questions/141718/what-is-the-path-environment-variable-and-how-do-i-add-to-it)!"
    is_to_install = true
  end
end
puts


###########################################################
unless is_to_install
  require "colorize"
  puts "Awesome! You seem to have all dependencies installed.".colorize(:red)
end


