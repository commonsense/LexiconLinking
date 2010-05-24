#!/bin/env python
#import scipy as sp
from operator import itemgetter
from collections import defaultdict
import pymongo
import copy 
import itertools
import cPickle, re
import networkx as nx
from start_wrapper import *
try:
    from nltk.corpus import verbnet as vn
    from nltk.corpus import wordnet as wn
    from nltk import pos_tag, clean_html
except:
    print "Could not import NLTK"

# load psyco when possible
try:
  print "loading Psyco"
  import psyco
  psyco.full()
except:
  pass


def clean_statement(statement):
  statement = clean_html(statement)
  statement = re.sub(r'<[a-zA-Z\/][^>]*>','',statement.split("|")[0])
  # why don't people put spaces after periods!  
  statement = re.sub(r'[.?!]([A-Z])',lambda x: ".  "+x.group(1),statement)
  # replace spaces and ampersands
  statement = statement.replace("&nbsp;"," ").replace("&amp;"," and ")
  # remove duplicate punctuation marks: ?, ! and .
  return re.sub(r'([.?!])[.?!]+',lambda x: x.group(1),statement)

def load_verb_ocean(relations_and_thresholds):
    verbmap = defaultdict(list)
    for line in open('verb_ocean.txt','r'):
        if line.count("::") == 0: continue
        line,score = line.strip().split("::")
        score = float(score)
        rel_start_indx = line.find("[")
        rel_end_indx = line.find("]")
        verb1 = line[0:rel_start_indx-1]
        verb2 = line[rel_end_indx+2:-1]
        relation = line[rel_start_indx+1:rel_end_indx]
        if relation in relations_and_thresholds.keys() and score >= relations_and_thresholds[relation]: 
            #print relation, verb1, verb2, score
            if verb2 not in verbmap[verb1]: 
                verbmap[verb1].append((verb2,relation,score))
                #verbmap[verb2].append(verb1)
    return verbmap

def enforce_and_resolve_transitivity(verb1,verb2,graph,relation):
    if not nx.is_directed_acyclic_graph(graph):
        print "Adding edge", verb1, verb2, "violates dag of", relation
        if graph.has_edge(verb1,verb2) and graph.has_edge(verb2,verb1):
            score1 = graph.get_edge_data(verb1,verb2)['score']
            score2 = graph.get_edge_data(verb2,verb1)['score']
            print score1
            if score1 > score2:
                print " (%s,%s) > (%s,%s)   [%i,%i]" % (verb1,verb2,verb2,verb1,score1,score2)
                graph.remove_edge(verb2,verb1)
            else:
                print " (%s,%s) > (%s,%s)   [%i,%i]" % (verb2,verb1,verb1,verb2,score1,score2)
                graph.remove_edge(verb1,verb2)
        else:
            # removing edge:
            score1 = graph.get_edge_data(verb1,verb2)['score']
            graph.remove_edge(verb1,verb2)
            graph.add_edge(verb2,verb1,relation=relation,score=score1)
            if not nx.is_directed_acyclic_graph(graph):
                graph.remove_edge(verb2,verb1)
            else:
                print "Hola,", verb2, verb1




def load_verb_ocean_into_graph(relations_and_thresholds):
    graphs = {}
    for relation in relations_and_thresholds.keys():
        graphs[relation] = nx.DiGraph()
    verbmap = defaultdict(list)
    for line in open('verb_ocean.txt','r'):
        if line.count("::") == 0: continue
        # parse the file to extract verbs and relation
        line,score = line.strip().split("::")
        score = float(score)
        rel_start_indx = line.find("[")
        rel_end_indx = line.find("]")
        verb1 = line[0:rel_start_indx-1]
        verb2 = line[rel_end_indx+2:-1]
        relation = line[rel_start_indx+1:rel_end_indx]
        if score < relations_and_thresholds[relation]: continue 
        # enforce relational properties
        is_transitive = False
        is_symmetric = False
        if relation in ['stronger-than','happens-before']:
            is_transitive = True
        if relation in ['similar','is-opposite']:
            is_symmetric = True
        # add edge:
        graphs[relation].add_edge(verb1,verb2,relation=True,score=score)
        if is_transitive:
            enforce_and_resolve_transitivity(verb1,verb2,graphs[relation],relation)
        if is_symmetric:
            graphs[relation].add_edge(verb2,verb1,relation=True,score=score)

            if is_transitive:
                enforce_and_resolve_transitivity(verb2,verb1,graphs[relation],relation)
        if relation in relations_and_thresholds.keys() and score >= relations_and_thresholds[relation]: 
            #print relation, verb1, verb2, score
            if verb2 not in verbmap[verb1]: 
                verbmap[verb1].append((verb2,relation,score))
                #verbmap[verb2].append(verb1)
    return graphs 


def load_pspan_results():
    slot_graph = nx.DiGraph()
    PROJECT = "goal_names"
    keys = {}
    nonterminals = set()
    for line in open('sequence_mining/%s/%s.keys' % (PROJECT,PROJECT),'r'):
        line = line.strip().split("\t")
        nonterminals.add(line[1])
        keys[int(line[0])] = line[1]

    f = open("sequence_mining/%s/%s.new_keys" % (PROJECT,PROJECT),'r')
    weight = 0
    slot_to_id = {} 
    seen = False
    id_to_slots = defaultdict(list) 
    for line in f.readlines():
        line = line.strip().split("\t")
        lhs = int(line[0])
        weight += 1
        slot_vals = []
        for e in line[1:]:
            slot_val = []
            for token in e.split():
                if int(token) == 0: continue
                # remove the terminal symbols that appear in the gramma
                if int(token) in keys and keys[int(token)]:
                    slot_val.append(keys[int(token)])
                elif id_to_slots.has_key(int(token)):
                    slot_val = []
                    continue
            slot_vals.append(slot_val)
        for sv in slot_vals:
            sv_string = ' '.join(sv)
            if len(sv_string) < 4: continue
            slot_to_id[sv_string] = lhs
            id_to_slots[lhs].append(sv_string)
    return slot_to_id 




#slot_to_id = load_pspan_results()
#nx.write_dot(slot_graph,"slot_graph.dot")


from pymongo import Connection
connection = Connection('localhost')
db = connection.sm

def find_other_verb_args(verb,nouns):
    results = db.events.find({"verb": {"$in": [verb]}, "noun": {"$in": nouns}})
    return results


def load_verb_args():
    from time import time
    start_time = time()

#events = db.events.find({"verb": {"$in": ["run"]}, "count": {"$gte": 150}})
#events = db.events.find({"count": {"$gte": 100}})
#events = db.events.find({"count": {"$gte": 100}})
    nouns = db.events.find({"verb": {"$in": ["run"]}})
    ns = list(set([t['noun'] for t in nouns]))
    print ns
    events = db.events.find({"noun": {"$in": ns}, "count": {"$gte": 80}})

#Jevents = db.events.find({"count": {"$gte": 200}})
# extract labels
#verbs = list(set(map(lambda x: x[0], bip.keys())))
#nouns = list(set(map(lambda x: x[1], bip.keys())))
    lattice = dual_taxonomy_builder(events)

    print "elapsed time:", time() - start_time

def get_nouns(seq):
    """ Takes a sequence of (token,pos) tuples and removes
    and concatinates noun phrases"""
    nouns = []
    tmp = []
    for (word,tag) in seq:
        if tag[0] == 'N':
            tmp.append(word)
        else:
            if len(tmp) > 0:
                nouns.append(' '.join(tmp))
            tmp = []
    if len(tmp) > 0:
        nouns.append(' '.join(tmp))
    return nouns


colors = { 'happens-before': 'blue', 'can-result-in': 'pink', 'stronger-than': 'green', 'similar':'red', 'opposite-of':'purple','unk':'orange', 'low-vol':'brown'}
relations_and_thresholds = { 'happens-before': 11, 'can-result-in': 11, 'stronger-than': 10, 'similar':11, 'opposite-of':10,'unk':20, 'low-vol':10}
#relations = load_verb_ocean(relations_and_thresholds)
try:
    graphs = cPickle.load(open('verb_ocean.pickle','r'))
except:
    graphs = load_verb_ocean_into_graph(relations_and_thresholds)
    cPickle.dump(graphs,open('verb_ocean.pickle','wb'))

replaced = {}
def replace_with_id(seq):
    if replaced.has_key(seq):
        return replaced[seq]
    else:
        i = len(replaced.keys())
        replaced[seq]=i
        return i

        


OUTPUT_DIR = "boost_seq"
index_file = open("%s.index" % OUTPUT_DIR, 'w') 




def parse_plan_statements(plan_statements,save=True):
    """ Iterates through a list of plan statements and parses each sentence using 
    the web-based START parser, and optionally when save==True, adds the parsed 
    result to the plan['parsed'] entry and saves the pickle file as 'parse_plans.pickle'
    
    It also keeps track of the verb count, used later for TD-IDF normalization."""
    from pymongo import Connection
    connection = Connection('localhost')
    db = connection.plans
    ct = 0
    verb_ct = defaultdict(int)
    parsed_correctly = 0
    parsed_incorrectly = 0
    # go through each goal statement 
    for goal, plans in goals.items()[1:2]:
        ct += 1
        print "--- On goal #", ct, "of ", len(goals.keys())
        # for each 'how_i_did_it' story, which I call a 'plan'
        parsed_plans = []
        goal_verb_count = defaultdict(int)
        for plan in plans:
            parsed = {} 
            unique_verbs = []
            # for each sentence in the plan (found earlier by a sentence segmenter) 
            for sent in plan['plan']:
                # parse the plan
                parsed_sent = start_parse(' '.join(sent))
                # see if it was parse-able
                if len(parsed_sent.keys()) == 0:
                    parsed_incorrectly += 1
                else:
                    parsed_correctly += 1
                    # merge all the parse entries for the entire plan
                    parsed.update(parsed_sent)
                    # add all words (verbs?) with '+' inside to verb
                    # frequency count
                    for word in parsed_sent.keys():
                        if word.count('+') != 0:
                            unique_verbs.append(word)
                            verb_ct[word] += 1
            parsed_plans.append(parsed)
            for verb in set(unique_verbs):
                goal_verb_count[verb] += 1
        db.goals.insert({'goal':goal, 'num_plans': len(plans), 'parsed': parsed_plans,'verb_counts': goal_verb_count})
        # store the parse dictionary in the plan
    db.meta.insert({'verb_counts_global':verb_ct})
    print "Parsed ", parsed_correctly, " of ", parsed_incorrectly+parsed_correctly, " sentences "
    #db.goals.create_index(['goal','num_plans'], unique=True) 


goals = cPickle.load(open('./plan_jar/plans.cornichon.pickle','r'))
parse_plan_statements(goals)

def salient_plans(): 
    mDur = min(filter(lambda x: x != 0, map(lambda y: y['duration'], goals[mGoal])))
    verb_tfidf = defaultdict(float)  
    verb_tf = defaultdict(int)
    for plan in goals[mGoal]:
        verbs = []
        for sent in plan['plan']:
            tagged_sent = pos_tag(sent)
            verbs += filter(lambda x: x[1][0] == 'V', tagged_sent)
        for verb in set(verbs):
            # only count these once per plan
            verb_tf[verb] +=1
    for verb in verb_tf.keys():
        verb_tfidf[verb] = (float(verb_tf[verb])/sum(verb_tf.values())/(float(verb_ct[verb])/sum(verb_ct.values())))
    for k, v in sorted(verb_tfidf.items(),key=lambda x: x[1],reverse=True):
        print k, v

    print mGoal


def parse_plan_statements2():
    file_num = 0
    categories = 0
    for goal, plans in goals.items():
        goal = clean_statement(goal.lower())
        verb = goal.split()[0]
        # create a graph composed of different verb relationships
        # that is centered around the neighboring nodes for the     
        # verbs in VerbOcean
        goal_graph = nx.MultiDiGraph()
        for relation,graph in graphs.items():
            if relation not in ['stronger-than','similar']: continue
            if graph.has_node(verb):
                goal_graph.add_edges_from(graph.edges(graph.neighbors(verb),data=True),color=colors[relation],label=relation)
        nouns = get_nouns(pos_tag(goal.split())[1:])
        print goal 
        for alt_verb in goal_graph.nodes():
            # output sequence output file
            of = file('%s/%i' % (OUTPUT_DIR,file_num),'w')
            output = [alt_verb]+nouns
            for o in output:
                of.write("%i\n" % (replace_with_id(o)))
            of.close()
            index_file.write("%s/%i %i\n" % (OUTPUT_DIR,file_num,categories))
            file_num +=1
        categories +=1
        if categories == 80: break
        #print plans

    kf = open("%s.keys" % (OUTPUT_DIR), "w")
    for k,v in replaced.items():
        kf.write("%i\t%s\n" % (v,k))

        #print len(plans)

