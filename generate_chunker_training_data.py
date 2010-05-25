from csc.conceptnet4.models import *
en = Language.get('en')
import cPickle
import nltk
import random
all_plans = cPickle.load(open('plan_jar/plans.pickle','r'))

"""
1 = has first sub event
2 = has last subevent
3 = has prereq
7 = used fork
9 = motivated by goal
10 desires
17 causes desire
18 = causes
19 = has subevent

"""
from start_wrapper import remove_nested_phrases

def generate_training_data():
    """ Creates training data for training a chunk extracting classifier.
    Bootstraps relations from conceptnet.
    """
    RELATION_TYPES = { 'SUBEVENT': {'RELATIONS':"1,2,9,10,19", 'MIN':0},\
                       'PARTOF': {'RELATIONS':"21", 'MIN':4},\
                       'LOCATION': {'RELATIONS':"6", "MIN":3},\
                       'NEXTEVENT': {'RELATIONS':"7,18,17", 'MIN':0},\
                       'PREVEVENT': {'RELATIONS':"3",'MIN':0}}

    outFile = open('training_goals.txt','w')

    for entry in all_plans.keys():
        #normalize goal string
        n_goal_names = set([en.nl.normalize(gn) for gn in map(lambda x: str(x['alt_name']), all_plans[entry])+[entry]])
        print "Goal_names", n_goal_names
        relations_for_goal = {}
        total = 0
        for relation, data in RELATION_TYPES.items():
            concepts = []
            for n_goal in n_goal_names:
                concepts += map(lambda x: x.concept2.text, Assertion.objects.filter(concept2__text=n_goal,language=en,score__gte=data['MIN']).extra(where=["relation_id in (%s)" % data['RELATIONS']]))
            total += len(concepts)
	    print concepts
            relations_for_goal[relation] = set(concepts)
        if total*8 / float(len(all_plans[entry])) < 1: 
            print "Skipping Goal", total, total*10/float(len(all_plans[entry])) # not enough knowledge about it
            continue

        for plan in all_plans[entry]:
            sentences = [nltk.word_tokenize(en.nl.normalize(remove_nested_phrases(' '.join(sent)))) for sent in plan['plan']] 
            tagged_sentences = [nltk.pos_tag(sent) for sent in sentences]
            begin_chunk = {}
            end_chunk = {}
            positions = 0
            offset = 0
            print sentences
            for sent in sentences:
                sent =  ' '.join(sent)
                for relation_type, concepts in relations_for_goal.items():
                    for concept in concepts:
                        if concept in sent:
                            print "Concept in sent", concept, "\\\\", relation_type, "\\\\", sent,"\\\\",  concept.split()[-1]
                            # compute number of spaces between starting and ending position of the sentence 
                            start_pos = len(sent[0:sent.index(concept)-1].split())
                            end_pos = len(sent[0:sent.index(concept)+len(concept)-1].split())-1
                            begin_chunk[offset+start_pos] = relation_type
                            end_chunk[offset+end_pos] = relation_type
                offset += len(sent.split())
            normalized_counter = 0
            for t_sent in tagged_sentences:
                default_chunk = "O" 
                if len(begin_chunk.keys()) < 1: break 
                outFile.write("\n")
                for t_toke in t_sent:
                    if begin_chunk.has_key(normalized_counter):
                        tag = "B-%s" % (begin_chunk[normalized_counter])
                        default_chunk = "I-%s" % (begin_chunk[normalized_counter])
                        if end_chunk.has_key(normalized_counter) and end_chunk[normalized_counter] == begin_chunk[normalized_counter]:
                            default_chunk = "O"
                    elif end_chunk.has_key(normalized_counter):
                        tag = "I-%s" % (end_chunk[normalized_counter])
                        default_chunk = "O"
                    else:
                        tag = default_chunk
                    normalized_counter += 1

                    print "%-20s%-10s%-10s" % (t_toke[0],t_toke[1],tag)
                    outFile.write("%-40s%-10s%-10s\n" % (t_toke[0],t_toke[1],tag))



generate_training_data()

