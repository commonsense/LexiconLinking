
from collections import defaultdict
import os, re

patterns = defaultdict(list)
for i, line in enumerate(open("project.txt",'r').readlines()):
    for j, v in enumerate(line.split()):
        if v == "1":
            patterns[int(j)].append(int(i))



for num, vals in patterns.items():
    print "Pattern %i" % num
    for v in vals:
        os.system("head -n %i goal_names.txt | tail -n 1" % v)
    print "\n"

    
