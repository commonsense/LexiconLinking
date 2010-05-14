import scipy
import numpy
from scipy import sparse
import pymongo
from numpy.random import randint

## an exceedingly simple test case
testWeights = scipy.array([[5,3,5],[4,4,4],[0,10,0]])
edge = ("meow","rabbit",5)

testNounSenses = ["dog","cat","rabbit"]
testVerbSenses = ["run","breathe","meow"]

##testVerbSenses = [["run",0],["breathe",1],["meow",2]] ## [[verbSenseName,first_Class_Index,second_Class_Index, third...],...]
##testNounSenses = [["dog",0],["cat",1],["rabbit",2]]

## verb = row, noun = column


## translate to unicode:
def toUnicode(l):
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


    ## IN: list of verbs, list of nouns, and a noun-verb pairs list of format [[noun, verb, weight]...]
    def __init__(self, verbs, nouns, weights):
        self.verbSenseNum = 0
        self.nounSenseNum = 0

        self.verbs = numpy.array(verbs)
        self.nouns = numpy.array(nouns)
        self.weights = weights

        ## given a numpy array of format [verb, noun, num], returns a matrix with that info.
    def createWeights(self, nvList):
        ## (just in case, make the nvList a numpy array)
        nvList = numpy.array(nvList)

        weights = numpy.zeros(nvList.shape, dtype=int)
        for nvPair in nvList:
            vrow = numpy.where(self.verbs == nvPair[0])
            ncol = numpy.where(self.nouns == nvPair[1])
            assert (len(vrow) == 1), str.format("Duplicate verb name if >0, Invalid name passed to createWeights if ==0 : {0}", len(vrow))
            assert (len(ncol) == 1), str.format("Duplicate noun name if >0, Invalid name passed to createWeights if ==0 : {0}", len(ncol))
## for now, no print, just assertion
#             if (len(vrow) != 1):
#                 print str.format("Duplicate verb name if >0, Invalid name passed to createWeights if ==0 : {0}", len(vrow))
#             if (len(ncol) != 1):
#                 print str.format("Duplicate noun name if >0, Invalid name passed to createWeights if ==0 : {0}", len(ncol))
            weights[vrow][ncol] = nvPair[2]
        return weights

    ## given a verb, look in the ordered list of verbs to find the correct index
    ## x = verb or noun (ex: "run")
    ## xlist = self.verbs or self.nouns (ex: self.verbs)
    def safeFindIndex(x, xlist):
        assert (x in xlist), str.format("Invalid verb/noun {0}", x)
        loc = numpy.where(xlist == x)
        assert (len(loc) == 1), str.format("{0} occurences of {1} in verb/noun list", (x, len(loc)))
        return loc

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
        pass
        ## ***in***
        ## verbOrNoun: is the new class a verb or a noun?
        ## className: the "name' of the class
        ## ***actions***
        ## (should only be done if class is already present)
        ## remove associated self.weights row/column
        ## remove class name from self.verbs/self.nouns

    ## merges 2 or more classes into a new class
    def merge(self, *classes):
        pass
        ## ***in***
        ## MUST be at least 2 arguments passed into "*classes"
        ## ***actions***
        ## add new class, then kill old classes

        ## for now, sum weights
        ## will deal with multiple senses/pointers issue later

def test1():
    A = sparse.lil_matrix((10,20))
    A[:,:] = randint(0,300,200)
    rows, cols = A.shape
    fake_nouns = []
    fake_verbs = []
    for i in xrange(0,rows): fake_verbs.append("verb_%i" % i)
    for i in xrange(0,cols): fake_nouns.append("noun_%i" % i)
    L = Lattice(fake_verbs,fake_nouns,A)

#if __name__ == "__main__":
A = test1()

def test2()
    from pymongo import Connection
    connection = Connection('localhost')
    db = connection.sm


## sort by cluster membership
## plot a binary matrix

## n v num
