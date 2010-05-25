import nltk

from nltk.corpus.reader.conll import * 
from nltk.classify import MaxentClassifier

from nltk.tag import TaggerI, untag
from nltk.chunk import ChunkParserI, tree2conlltags, conlltags2tree
from features import *

class ClassifierTagger(TaggerI):
	'''Abstracted from "Training Classifier-Based Chunkers" section of
    http://nltk.googlecode.com/svn/trunk/doc/book/ch07.html
	'''
	def __init__(self, feature_extractor, classifier):
		self.feature_extractor = feature_extractor
		self.classifier = classifier

	def tag(self, sent):
		history = []

		for i, word in enumerate(sent):
			featureset = self.feature_extractor(sent, i, history)
			tag = self.classifier.classify(featureset)
			history.append(tag)

		return zip(sent, history)

	@classmethod
	def train(cls, train_sents, feature_extractor, classifier_cls, **kwargs):
		train_set = []

		for tagged_sent in train_sents:
			untagged_sent = untag(tagged_sent)
			history = []

			for i, (word, tag) in enumerate(tagged_sent):
				featureset = feature_extractor(untagged_sent, i, history)
				train_set.append((featureset, tag))
				history.append(tag)

		classifier = classifier_cls.train(train_set, **kwargs)
		return cls(feature_extractor, classifier)




class ClassifierChunker(nltk.chunk.ChunkParserI):
	def __init__(self, train_sents, *args, **kwargs):
		tag_sents = [tree2conlltags(sent) for sent in train_sents]
		train_chunks = [[((w,t),c) for (w,t,c) in sent] for sent in tag_sents]
		self.tagger = ClassifierTagger.train(train_chunks, *args, **kwargs)

	def parse(self, tagged_sent):
		if not tagged_sent: return None
		chunks = self.tagger.tag(tagged_sent)
		return conlltags2tree([(w,t,c) for ((w,t),c) in chunks])

train_sents = ConllChunkCorpusReader('./',['training_goals.txt'],('PARTOF','NEXTEVENT','SUBEVENT','PREVEVENT','PLACE')).chunked_sents()
# featx is one of the feature extractors defined above 
chunker = ClassifierChunker(train_sents, pos, MaxentClassifier, min_lldelta=0.01, max_iter=10)


