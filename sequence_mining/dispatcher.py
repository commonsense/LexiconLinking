import os
from collections import defaultdict

PSPAN_NAME = "lexbuild"
PROJECT_NAME = "step_names"

top_k =  10 
maxgap = 3 
minlen = 4
total_runs = 100
# TODO:  Only load files that are marked in the project file
# don't depend on support from lexicon for replacing keys.  Do that instead within this program
seen_seq = []
for run in xrange(0,total_runs):
    print "Run number %i " % (run)
    os.system("./%s -P %s -v -K %i -G %i --length-min %i" % (PSPAN_NAME,PROJECT_NAME,top_k,maxgap,minlen))

    all_sequence_matches = []
    # load keys
    keynames = {}
    replaced = {}
    for line in open('%s/%s.keys' % (PROJECT_NAME,PROJECT_NAME),'r').readlines():
        try:
            k,v = line.strip().split('\t')
            keynames[int(k)]=v
        except:
            print "Skipping line", line
            pass
    replacements = open('%s/%s.replace' % (PROJECT_NAME,PROJECT_NAME),'r').readlines()
    for line in replacements:
        if line[0] != '#':
            vals = [int(x) for x in line.strip().split('\t')]
            keynames[vals[0]] = 'seq-%i-%s' % (vals[0],'/'.join(map(lambda x: keynames[x], vals[1:])))
            to_replace = vals[1:]
            got = []
            for v in to_replace: 
                if replaced.has_key(v) and got.count(v) == 0:
                    to_replace.append(replaced[v])
                replaced[v] = vals[0]
                got.append(v)

    # parse sequences results
    sequences = map(lambda x:\
                    map(lambda z: int(z), x[2:-3].split(" ] [ ")), open('%s/%s.out' % (PROJECT_NAME,PROJECT_NAME),'r').readlines())
    sequences = filter(lambda y: seen_seq.count('#'.join(map(lambda x: str(x),y)))== 0,sequences)

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
                for i,el in enumerate(map(lambda x: (replaced.has_key(int(x.strip())) and replaced[int(x.strip())])\
                                          or int(x.strip()), \
                                          open("%s/seq/%s" %(PROJECT_NAME,file),'r').readlines())):
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
                            if sequence_matches[i].count(vals) == 0:
                                sequence_matches[i].append(vals)
            # save found sequence gap fillers
            all_sequence_matches.append(sequence_matches)

    # find longest filler and write them to file
    max_len = 0
    max_seq = 0
    rf = open('%s/%s.replace' %(PROJECT_NAME,PROJECT_NAME),'a')
    for seq, matches in enumerate(all_sequence_matches): 
        print "Sequence %i " % seq , ' / '.join(map(lambda x: keynames[x], sequences[seq]))
        for slot_i, vals in matches.items():
            lv = len(set(map(lambda x: x[0], vals)))
            if lv >= max_len:
                max_len = lv 
                max_seq = (seq,slot_i)
            if lv >= 2:
                seen_seq.append('#'.join(map(lambda x: str(x),sequences[seq])))
                rf.write("#%0.3i%0.2i0%i\t%s\n" % (run,seq,slot_i,'\t'.join(map( lambda x: keynames[x[0]], vals))))
                items = map(lambda x: str(x[0]), vals)
                rf.write("%0.3i%0.2i0%i\t%s\n" % (run,seq,slot_i,'\t'.join(items)))
            for val in vals:
                print "\tSlot %i => %s" % (slot_i,' '.join(map(lambda x: keynames[x],val))), val
    rf.close()
    # make sure there's something here
    if max_len < 2 and minlen == 2:
        print "DONE - goodbye!" 
        break
    if max_len < 2:
        minlen -= 1 # lower minimum length
        continue



#--dspcacount 10  --length-min 3 -K 2   ; python interpret_stories.py")



