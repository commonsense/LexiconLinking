import scipy

## an exceedingly simple test case
testWeights = scipy.array([5,3,5],[4,4,4],[0,10,0])
testWeights = scipy.array([[5,3,5],[4,4,4],[0,10,0]])
edge = ("meow","rabbit",5)

testNounSenses = ["dog","cat","rabbit"]
testVerbSenses = ["run","breathe","meow"]

testVerbSenses = [["run",0],["breathe",1],["meow",2]] ## [[verbSenseName,first_Class_Index,second_Class_Index, third...],...]
testNounSenses = [["dog",0],["cat",1],["rabbit",2]]

## translate to unicode:
def toUincode(l):
    return [[unicode(l[i][0]), l[i][1]] for i in range(len(l))]

testNounSenses = toUnicode(testNounSenses)
testVerbSenses = toUnicode(testVerbSenses)

## AS OF TUESDAY APRIL 27:
## I was looking to see what the incoming datatype is, build around that.
## Then, finish methods, write tests.
## (unicode strings (run build_lattices.py in ~/goalmining/build_lattices)) DONE

## TODO as of 5/12/10:
## rewrite init method to put weights into proper place/use class indices properly.
## finish implementing methods
## (from build_lattices.py:- bipartite:  a list of ((:wjk:verb,noun),weight) tuples)

class Lattice:


    ## IN: list of verbSenses, list of nounSenses, and a [(len(verbSenses))x(len(nounSenses))] array
    ## description of verbSenses/nounSenses:
    ## [[verbName,first_Sense_Index,second_Sense_Index, third...],...]
    def __init__(self, verbSenses, nounSenses, weights):
        self.verbSenseNum = 0
        self.nounSenseNum = 0
        self.weights = weights
        ## weights: a matrix of verb sense vs noun sense weights: weight = weights[verb][noun]
        ## each verb gets a list of weights for nouns
        ## default value = None

        self.verbSenses = verbSenses ## array of strings
        self.nounSenses = nounSenses ## array of strings

        self.verbs = [self.verbSenses[i][0] for i in range(len(self.verbSenses))] ## names (strings)
        self.nouns = [self.nounSenses[i][0] for i in range(len(self.nounSenses))]



	## given a verbSense and nounSense, return the corresponding coordinates in self.weights
    def getIndices(self, verbSense, nounSense):
        ## check that verb/noun senses are in the dataset
        assert (verbSense in self.verbs), "The verb sense %s is not in the dataset." % (verbSense)
        assert (nounSense in self.nouns), "The noun sense %s is not in the dataset." % (nounSense)
        return scipy.where(self.verbs==verbSense)[0], scipy.where(self.nouns==nounSense)[0]

    ## gets unique sense number for verb/noun senses
    def getSenseNumber(self, isVerbSense):
        if isVerbSense:
            self.verbSenseNum += 1
            return str(self.verbSenseNum)
        self.nounSenseNum += 1
        return self.nounSenseNum

    ## adds an edge between two senses
    def add_edge(self, verbName, nounName, weight):
        ## this function may not be necessary

        ## check that verb/noun senses are in the dataset
        assert (verbName in self.verbs), "The verb sense %s is not in the dataset." % (verbName)
        assert (nounName in self.nouns), "The noun sense %s is not in the dataset." % (nounName)

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
