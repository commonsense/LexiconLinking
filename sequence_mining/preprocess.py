import sys,os
import fileinput

BREAK_AT = int(sys.argv[1]) # 0 if no break
os.system("rm -rf seq")
os.system("mkdir seq")

i = 0
words =[]

def replace_with_ids(seq):
    r = []
    for i in seq:
        if words.count(i) == 0:
            words.append(i)
        r.append(str(words.index(i)))
    return '\n'.join(r)


for line in fileinput.input(["goal_names.txt"]):
    line = replace_with_ids(line.strip().lower().split())
    os.system("touch seq/%s" % (str(i)))
    of = open("seq/%s" % (str(i)),'w')
    of.write(line)
    of.close()
    i += 1
    if BREAK_AT != 0 and i >= BREAK_AT:
        break

kf = open("keys.txt",'w')
for t,w in enumerate(words):
    kf.write("%s %s\n" % (str(t),w))
kf.close()

of = open("inseq.txt",'w')
for j in xrange(0,i):
    of.write("seq/%s\n" % (str(j)))