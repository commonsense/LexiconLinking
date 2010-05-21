#!/bin/env python
import scipy as sp
from nltk.corpus import wordnet as wn
from csc.divisi.util import get_picklecached_thing
from nltk.corpus import verbnet as vn
from operator import itemgetter
from collections import defaultdict
from nltk import pos_tag, clean_html
import pymongo
import copy 
import itertools
import cPickle, re
import networkx as nx

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
            slot_to_id[sv_string] = lhs
            id_to_slots[lhs].append(sv_string)





slot_graph = load_pspan_results()
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


relations_and_thresholds = { 'happens-before': 11, 'can-result-in': 11, 'stronger-than': 10, 'similar':11, 'opposite-of':10,'unk':20, 'low-vol':10}
#relations = load_verb_ocean(relations_and_thresholds)
try:
    graphs = cPickle.load(open('verb_ocean.pickle','r'))
except:
    graphs = load_verb_ocean_into_graph(relations_and_thresholds)
    cPickle.dump(graphs,open('verb_ocean.pickle','wb'))

def parse_plan_statements():
    goals = cPickle.load(open('./plan_jar/plans.cornichon.pickle','r'))

    colors = { 'happens-before': 'blue', 'can-result-in': 'pink', 'stronger-than': 'green', 'similar':'red', 'opposite-of':'purple','unk':'orange', 'low-vol':'brown'}

    for goal, plans in goals.items():
        goal = clean_statement(goal.lower())
        verb = goal.split()[0]
        # create a graph composed of different verb relationships
        # that is centered around the neighboring nodes for the     
        # verbs in VerbOcean
        goal_graph = nx.MultiDiGraph()
        for relation,graph in graphs.items():
            if graph.has_node(verb):
                print "subgraph", graph.edges(graph.neighbors(verb),data=True)
                goal_graph.add_edges_from(graph.edges(graph.neighbors(verb),data=True),color=colors[relation],label=relation)
        nouns = get_nouns(pos_tag(goal.split())[1:])
        print goal,"\n", verb, nouns
        if goal_graph != None: nx.write_dot(goal_graph,"goals.dot")
        x = raw_input("Press Enter to continue...")
        #print plans

        #print len(plans)


