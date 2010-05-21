
PROJECT_NAME = "goal_names"
LIMIT = 2

require 'rubygems'
require 'active_record'

ActiveRecord::Base.establish_connection( 
:adapter => "postgresql",
:host => "localhost",
:username => "dustin",
:password => "",
:database => "goalnet"
) 

class Goal < ActiveRecord::Base
end


# query to return steps, sorted by their task, based on constraints for variance.
steps = Goal.find_by_sql(["select DISTINCT(goal_str) from how_tos where parsed = 1"]) 
    

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
   s  = s['goal_str'].downcase
   # print s.methods.join("\n")
    # if ct % 100 == 0
    puts s
    f_out = File.open("#{PROJECT_NAME}/seq/#{ct}",'w')
    index_out.puts "#{PROJECT_NAME}/seq/#{ct}"
    #  end
    all_stories_out.puts s
    replace_with_numbers($step_names,s.split(" ")).each {|x| f_out.puts x }
    ct += 1
end

k_out = File.open("#{PROJECT_NAME}/#{PROJECT_NAME}.keys",'w')
$step_names.each {|k,v| k_out.puts "#{v}\t#{k}"}






