
PROJECT_NAME = "step_names"
LIMIT = 5

require 'rubygems'
require 'active_record'


if ARGV.size < 1
  print "first argument should be password"
  break
end
print ARGV[0]

ActiveRecord::Base.establish_connection( 
:adapter => "mysql", 
:host => "sql.media.mit.edu",
:username => "litchfield",  
:password => ARGV[0],
:database => "todogo_game"
) 

class Steps < ActiveRecord::Base
end


# query to return steps, sorted by their task, based on constraints for variance.
steps = Steps.find_by_sql("select step, count(*) as num from steps
group by step")


last_task = ""
last_taskid = ""
f_out = nil #output file pointer


`rm -rf #{PROJECT_NAME}`  # shell command to clear stories directory
`mkdir #{PROJECT_NAME}`
`mkdir #{PROJECT_NAME}/seq`


index_out = File.open("#{PROJECT_NAME}/#{PROJECT_NAME}.index",'w')
all_stories_out =  File.open("#{PROJECT_NAME}/#{PROJECT_NAME}.all",'w')

$step_names = Hash.new(0)

# finds or adds id of item i in hash h
def hash_to_id(h,i)
	if h[i] == 0
		h[i] = h.size+1
	end
	h[i]
end

def replace_with_numbers(h,seq)
  seq.collect {|x| hash_to_id(h,x)}
end

ct = 0

for s in steps
    if s['num'].to_i < LIMIT
      next
    end 
    s = s['step']
   # print s.methods.join("\n")
    # if ct % 100 == 0
    f_out = File.open("#{PROJECT_NAME}/seq/#{ct}",'w')
    index_out.puts "#{PROJECT_NAME}/seq/#{ct}"
    #  end
    all_stories_out.puts s
    replace_with_numbers($step_names,s.split(" ")).each {|x| f_out.puts x }
    ct += 1
end

k_out = File.open("#{PROJECT_NAME}/#{PROJECT_NAME}.keys",'w')
$step_names.each {|k,v| k_out.puts "#{v}\t#{k}"}






