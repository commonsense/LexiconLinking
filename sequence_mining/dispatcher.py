import os
from collections import defaultdict

PSPAN_NAME = "lexbuild"
PROJECT_NAME = "step_names"



top_k = 10 
maxgap = 3 
minlen = 4
total_runs = 10 

for run in xrange(0,total_runs):
    print "Run number %i " % (run)
    os.system("./%s -P %s -v -K %i -G %i --length-min %i" % (PSPAN_NAME,PROJECT_NAME,top_k,maxgap,minlen))

    all_sequence_matches = []
    # load keys
    keynames = {}
    for line in open('%s/%s.keys' % (PROJECT_NAME,PROJECT_NAME),'r').readlines():
        try:
            k,v = line.strip().split('\t')
            keynames[int(k)]=v
        except:
            print "Skipping line", line
            pass

    # parse sequences results
    sequences = map(lambda x: map(lambda z: int(z), x[2:-3].split(" ] [ ")), open('%s/%s.out' % (PROJECT_NAME,PROJECT_NAME),'r').readlines())
    for seq in sequences:
        print seq

    # project sequence results 
    for root, dirs, files in os.walk('%s/seq' % PROJECT_NAME):
        for j,seq in enumerate(sequences):
            sequence_matches = defaultdict(list)
            for file in files:
                pos_in_seq = 0
                slot_vals = []
                tmp_gap = []
                for i,el in enumerate(map(lambda x: int(x.strip()), open("%s/seq/%s" %(PROJECT_NAME,file),'r').readlines())):
                    if i+(len(seq)-pos_in_seq-1) > len(seq) or pos_in_seq == len(seq):
                        break 
                    if el == seq[pos_in_seq]:
                        pos_in_seq += 1
                        if pos_in_seq != 1 and len(tmp_gap) != 0:
                            slot_vals.append(tmp_gap)
                        tmp_gap = []
                    else:
                        tmp_gap.append(el)
                if pos_in_seq == len(seq):
                    # this is a match
                    if len(slot_vals) > 0:
                        print "Match for seq %i " % (j), slot_vals
                        for i, vals in enumerate(slot_vals):
                            sequence_matches[i].append(vals)
            # save found sequence gap fillers
            all_sequence_matches.append(sequence_matches)

    # find longest filler and write them to file
    max_len = 0
    max_seq = 0
    for seq, matches in enumerate(all_sequence_matches): 
        print "Sequence %i " % seq , ' / '.join(map(lambda x: keynames[x], sequences[seq]))
        for slot_i, vals in matches.items():
            if len(vals) >= max_len:
                max_len = len(vals)
                max_seq = (seq,slot_i)

            for val in vals:
                print "\tSlot %i => %s" % (slot_i,' '.join(map(lambda x: keynames[x],val))), val

    print "Writing replacement"
    rf = open('%s/%s.replace' %(PROJECT_NAME,PROJECT_NAME),'a')
    sqn, sln = max_seq
    items = all_sequence_matches[sqn][sln]
# check to make sure there are only one item in each gap
    for item in items:
        if len(item) > 1:
            print "Error: item list longer than one", len(item), item

    rf.write("#%0.3i%0.2i0%i\t%s\n" % (run,sqn,sln,'\t'.join(map( lambda x: keynames[x[0]], items))))
    items = map(lambda x: str(x[0]), items)
    rf.write("%0.3i%0.2i0%i\t%s\n" % (run,sqn,sln,'\t'.join(items)))

#--dspcacount 10  --length-min 3 -K 2   ; python interpret_stories.py")



