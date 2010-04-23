
MAX_VARIANCE = 1

require 'rubygems'
require 'active_record'


if ARGV.size < 1
  print "first argument should be password"
  return
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
steps = Steps.find_by_sql("
select b.task, b.avg, s2.* from tasks t2, steps s2, (
select task, id, count(*) as num, sum(num_steps)/count(*)  as avg,
pow(num_steps-(sum(num_steps)/(count(*))),2) as var_sq 
from
(select t.task, t.id, s.step, count(s.step) as num_steps from tasks t
inner join steps s on t.id = s.taskid
group by s.taskid) a
group by task
order by num desc, var_sq asc
limit 100) b 
where s2.taskid = t2.id
and t2.task = b.task
order by b.task, s2.id
")


last_task = ""
last_taskid = ""
f_out = nil #output file pointer
index_out = File.open("story_sequences_index.txt",'w')
all_stories_out = File.open("all_stories.txt",'w')

$story_words = Hash.new(0)

# finds or adds id of item i in hash h
def hash_to_id(h,i)
	if h[i] == 0
		h[i] = h.size+1
	end
	h[i]
end

def replace_words_with_numbers(seq)
  seq.collect {|x| hash_to_id($story_words,x)}
end

`rm -rf story_sequences`  # shell command to clear stories directory
`mkdir story_sequences`  
ct = 0
for s in steps
	if last_task != s.task
		#print "==============\n#{s.task}\n==============\n"
		last_task = s.task
		last_taskid = s.taskid
    dir = "story_sequences/#{s.task.gsub(" ","_")}"
    `mkdir #{dir}`
    f_out = File.open(dir+"/#{ct}",'w')
	end
	if s.taskid != last_taskid
		last_taskid = s.taskid
    f_out.close()
    ct += 1
		#f_out.puts "# #{s.avg}"
    f_out = File.open(dir+"/#{ct}",'w')
    index_out.puts "#{dir}/#{ct}"

	end
  all_stories_out.puts s.step
	f_out.puts "#{replace_words_with_numbers(s.step.split(" ")).join(" ")}" 
end

k_out = File.open("story_keys.txt",'w')
$story_words.each {|k,v| k_out.puts "#{v}\t#{k}"}

