import scipy

## an exceedingly simple test case
testNouns = ["dog","cat","rabbit"]
testVerbs = ["run","breathe","meow"]
testWeights = scipy.array([5,3,5],[4,4,4],[0,10,0])
testWeights = scipy.array([[5,3,5],[4,4,4],[0,10,0]])
edge = ("meow","rabbit",5)



class Lattice:

    def __init__(self, weights):
		self.verbs = verbs ## array of strings
		self.nouns = nouns ## array of strings
		self.weights = weights
		## weights: a matrix of verb sense vs noun sense weights: weight = weights[verb][noun]
		## each verb gets a list of weights for nouns
		## default value = None

	## given a verbSense and nounSense, return the corresponding coordinates in self.weights
	def getIndices(self, verbSense, nounSense):
		## check that verb/noun senses are in the dataset
		assert (verbSense in self.verbs), "The verb sense %s is not in the dataset." % (verbSense)
		assert (nounSense in self.nouns), "The noun sense %s is not in the dataset." % (nounSense)
		return scipy.where(self.verbs==verbSense)[0], scipy.where(self.nouns==nounSense)[0]
        ## adds an edge between two senses
	def add_edge(self, verbSense, nounSense, weight):
                ## this function may not be necessary

		## check that verb/noun senses are in the dataset
		assert (verbSense in self.verbs), "The verb sense %s is not in the dataset." % (verbSense)
		assert (nounSense in self.nouns), "The noun sense %s is not in the dataset." % (nounSense)

		## get the indices of the correct weight (this part is almost correct)
		indices = scipy.where(self.verbs==verbSense)[0],scipy.where(self.nouns==nounSense)[0]

		if self.weights[verbIndex][nounIndex] != None: ## already has an edge
			self.weights[verbIndex][nounIndex] = weight
		## else, doesn't already have an edge
		self.weights[verbIndex][nounIndex] = weight
		
	## adds a verb or noun class to the existing structure
	def add_class(self, verb, newWeights = None):
		## ***in***
		## verb: is the new class a verb (or a noun)?
		## newWeights: None is default
		## ***actions***
		## (should only be added if class is not already present, correct?)
		## add either row (verb) or column (noun) to existing structure:
		##		-) add to self.verbs/self.nouns list
		##		-) add row/column to self.weights matrix, with values from newWeights (or 0s/None)

		## how should we assign a "name" for the class??
                ## will be taken care of by verb names -> multiple verb senses implementation

                return add_general(verb, newWeights)

        def add_general(verb, newWeights):
            ## a more general add for both rows and columns
            a,b = (self.weights.shape[0]+1, self.weights.shape[1]) if verb else (self.weights.shape[0], self.weights.shape[1]+1)
            self.weights = scipy.resize(self.weights, (a,b))

            a,b = (self.weights.shape[0], self.weights.shape[1]) if verb else (self.weights.shape[0], self.weights.shape[1])

            if newWeights == None:
                newWeights = scipy.ones(a)
            if verb: self.weights[a-1,:] = newWeights
            else:
	        self.weights[:,b-1] = newWeights
            return self.weights

	## removes a verb or noun class from the existing structure
	def kill_class(self, verbOrNoun, className):
		## ***in***
		## verbOrNoun: is the new class a verb or a noun?
		## className: the "name' of the class
		## ***actions***
		## (should only be done if class is already present)
		## remove associated self.weights row/column
		## remove class name from self.verbs/self.nouns

	## merges 2 or more classes into a new class
	def merge(self, *classes)
		## ***in***
		## MUST be at least 2 arguments passed into "*classes"
		## ***actions***
		## add new class, then kill old classes

		## for now, sum weights
		## will deal with multiple senses/pointers issue later
