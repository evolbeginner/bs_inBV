#! /bin/env ruby

def render_progress(count, total, width: 40)
  total = 1 if total <= 0
  ratio = [[count.to_f / total, 0.0].max, 1.0].min
  filled = (ratio * width).round
  bar = "#" * filled + "-" * (width - filled)

  if $stdout.tty?
    print "\r[#{bar}] #{format('%6.2f', ratio * 100)}% (#{count}/#{total})"
    $stdout.flush
  else
    puts "#{format('%6.2f', ratio * 100)}% (#{count}/#{total})"
  end
end

def progressbar_for_bootstrapping(file:, total:, interval: 0.2)
  stop = false
  cmd_q = Queue.new

  thr = Thread.new do
    count = 0
    pos = 0

    loop do
      # handle commands
      until cmd_q.empty?
        cmd, payload = cmd_q.pop(true) rescue break
        case cmd
        when :reset
          count = 0
          pos = 0
          if payload[:truncate_file] && File.exist?(file)
            File.truncate(file, 0)
          end
          if $stdout.tty?
            print "\r\e[2K" # clear line
            $stdout.flush
          end
          render_progress(0, total)
        when :stop
          stop = true
        end
      end

      break if stop

      if File.exist?(file)
        File.open(file, "r") do |io|
          begin
            io.seek(pos, IO::SEEK_SET)
          rescue Errno::EINVAL
            pos = 0
            io.seek(0, IO::SEEK_SET)
          end

          io.each_line { count += 1 }
          pos = io.pos
        end

        render_progress(count, total)
      end

      sleep interval
    end
  end

  {
    reset: lambda { |truncate_file: false|
      cmd_q << [:reset, { truncate_file: truncate_file }]
    },
    stop: lambda {
      cmd_q << [:stop, {}]
      thr.join
      render_progress(total, total)
      puts if $stdout.tty?
    }
  }
end



##########################################
# old
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


def processbar(count, total)
  maxlen = 80
  prop = (count.to_f/total)*100
  jing = '#' * (prop*maxlen/100)
  printf "\r%-80s%s", jing, prop.round(2).to_s+'%'
end

