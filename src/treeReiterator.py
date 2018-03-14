from builtins import str
from builtins import map
from builtins import next
from builtins import range
from builtins import object
import itertools

# this is needed for python 2.7, probably.
# from sets import *

from reiterators import *
from utils import *
from functools import reduce

class treeReiterable(object) :
    """
treeFunc is a function which takes in predecessorData and outputs
an iterable generating tuples
(  ((subTreeData, multiplicity),(...),...),  nodeDecoration)

subTreeData is what will be fed to the next calls to treeFunc on the subtrees.
multiplicity is the number of "so far symmetric" copies of each subTreeData.

nodeDecoration contains data that decorates the root of the (big) tree,
when this choice of children is selected. The most succint, lossless option is
to decorate the root with (predecessorData, counter) where counter
refers to the order in the treeFunc reiterable.
"""
    def __init__(self, treeFunc, predecessorData) :
        self.treeFunc = treeFunc
        self.rootReiterable = treeFunc(predecessorData)

    def reiter(self, status) :
        return treeReiterator(self, status)

    def initialStatus(self) :
        return (self.rootReiterable.initialStatus(), False)

    def __iter__(self) :
        return self.reiter(self.initialStatus())

# predecessorData is parent's data only, not a tuple!
# status is a tree of rootReiter's statuses.
class treeReiterator(object) :
    def __init__(self, tr, status) :
        self.treeFunc = tr.treeFunc

        # think of rootReiter is a pointer to the list of options at the root.
        self.rootReiter = tr.rootReiterable.reiter(status[0])
        self.childrenReiterSet = status[1]
        if self.childrenReiterSet :
            self.childrenReiterable = status[2]
            self.childrenReiter = self.childrenReiterable.reiter(status[3])
            self.childrenDataMult = status[4]
            self.rootDecor = status[5]

    def __iter__(self) :
        return self

    def __next__(self) :
        if not self.childrenReiterSet :
            self.childrenDataMult,self.rootDecor = next(self.rootReiter)
            self.childrenReiterable = imapReiterable(\
                productReiterable(*[\
                    sortedProductReiterable(\
                        treeReiterable(self.treeFunc, cdm[0],), cdm[1])\
                for cdm in self.childrenDataMult]),
                lambda childecorts : reduce(lambda x, y : x + y, childecorts, ()))
            self.childrenReiter = iter(self.childrenReiterable)
            self.childrenReiterSet = True

        while True :
            try :
                return (self.rootDecor, next(self.childrenReiter))
            except StopIteration :
                self.childrenDataMult,self.rootDecor = next(self.rootReiter)
                self.childrenReiterable = imapReiterable(\
                    productReiterable(*[\
                        sortedProductReiterable(\
                            treeReiterable(self.treeFunc, cdm[0]), cdm[1])\
                    for cdm in self.childrenDataMult]),
                    lambda childecorts : reduce(lambda x, y : x + y, childecorts, ()))
                self.childrenReiter = iter(self.childrenReiterable)

    def status(self) :
        # (rootReiter.status(), childrenReiterSet)
        # if childrenReiterSet is False otherwise,
        # (rootReiter.status(), childrenReiterSet, childrenReiterable, childrenReiter.status(),
        # ... childrenDataMult, rootDecor)
        if self.childrenReiterSet :
            return (self.rootReiter.status(), self.childrenReiterSet, \
                    self.childrenReiterable, self.childrenReiter.status(),\
                    self.childrenDataMult, self.rootDecor)
        else :
            return (self.rootReiter.status(), self.childrenReiterSet)

class UnmarkedTreesFunc(object) :
    """
A silly example of a tree reiterator function. Below we generate all *six* isomorphism types of binary trees of depth two, with two types of leaves (yellow or red).

.. code-block:: py

    >>> from treeReiterator import *
    >>> f = UnmarkedTreesFunc([2],2)
    >>> tr = list(treeReiterable(f,0))
    >>> len(tr)
    6
    >>> for t in tr :
    ...     print(stree(t), end='')
    ...     print()
    ...
    Depth 0 Node
       Depth 1 Node
          Yellow Leaf
          Yellow Leaf
       Depth 1 Node
          Yellow Leaf
          Yellow Leaf
    <BLANKLINE>
    Depth 0 Node
       Depth 1 Node
          Yellow Leaf
          Yellow Leaf
       Depth 1 Node
          Yellow Leaf
          Red Leaf
    <BLANKLINE>
    Depth 0 Node
       Depth 1 Node
          Yellow Leaf
          Yellow Leaf
       Depth 1 Node
          Red Leaf
          Red Leaf
    <BLANKLINE>
    Depth 0 Node
       Depth 1 Node
          Yellow Leaf
          Red Leaf
       Depth 1 Node
          Yellow Leaf
          Red Leaf
    <BLANKLINE>
    Depth 0 Node
       Depth 1 Node
          Yellow Leaf
          Red Leaf
       Depth 1 Node
          Red Leaf
          Red Leaf
    <BLANKLINE>
    Depth 0 Node
       Depth 1 Node
          Red Leaf
          Red Leaf
       Depth 1 Node
          Red Leaf
          Red Leaf
    <BLANKLINE>
"""
    def __init__(self, childops, maxdepth) :
        self.childops = childops
        self.maxdepth = maxdepth

    def __call__(self, pd) :
        if pd >= self.maxdepth :
            return seqReiterable([((),"Yellow Leaf"), ((), "Red Leaf")] )
        return seqReiterable([(((pd + 1,i),), "Depth {} Node".format(pd)) \
                              for i in self.childops])

def stree(tr,formatter = lambda x,d : ' '*(d*3) + str(x) + '\n',depth = 0) :
    S = formatter(tr[0],depth)
    for c in tr[1] :
        S += stree(c, formatter, depth + 1)
    return S


def vertexDescriptions(forest) :
    for i in range(len(forest)) :
        yield (i,)
        for rnd in vertexDescriptions(forest[i][1]) :
            yield (i,) + rnd

def vertices(tree) :
    return ((),) + tuple(edges(tree))

def edges(tree) :
    """
an edge of a trees is uniquely determined by any non-root node of the tree.
"""
    return vertexDescriptions(tree[1])

def getVertex(vertex, tree) :
    if not vertex :
        return tree[0]
    return getVertex(vertex[1:], tree[1][vertex[0]])

def getSubtree(vertex, tree) :
    if not vertex :
        return tree
    return getSubtree(vertex[1:], tree[1][vertex[0]])


# it might be good to write a memoize decorator for the iterator class,
# so that the iterator constantly saves its status and output.
# This way eventually the function is "mapped out", and this
# should make running times much better.

# a canonical rep for a mutlisetPartition is one in which the
# elements are ordered from small to large.
# That's why the recursive computation needs to also know about the
# smallest allowed element.
def multiSetPartitions(n,k,smallest = 0) :
    if k <= 0 :
        if k < 0 or n != 0:
            return
        else :
            yield tuple()

    for i0 in range(smallest, n + 1) :
        for m0 in range(1, k + 1) : #range(min(n + 1 / i0,k)) :
            for pp in multiSetPartitions(n-i0*m0,k-m0,i0+1) :
                yield ((i0,m0),) + pp

def strMultiSet(ms) :
    i,m = ms[0]
    s = '{} x {}'.format(i,m)
    for (i,m) in ms[1:] :
        s+= ' + {} x {}'.format(i,m)
    return s

def someSymPartitions(n, ks) :
    for coarsePartition in tuplePartitions(n,len(ks)) :
        for finePartition in itertools.product(*[multiSetPartitions(coarsePartition[i],ks[i]) for \
                                    i in range(len(ks))]) :
            yield finePartition

def dummyPartition(k) :
    yield (((),k),)

def addPartition(generator, n) :
    for oldP in generator :
        ks = [t[-1] for t in oldP]
        for newP in someSymPartitions(n, ks) :
            flatP = ()
            for i in range(len(ks)) :
                for j in range(len(newP[i])) :
                    flatP += ((oldP[i][0] + (newP[i][j][0],), newP[i][j][1]),)
            yield flatP

def strSomeSymPartition(p) :
    S = ''
    for t in p :
        S += strMultiSet(t) + '; '
    return S

def tree_map(f, tr) :
    """
apply f recursively to all the tree nodes
"""
    return (f(tr[0]), tuple([tree_map(f, x) for x in tr[1]]))

def tree_reduce(f, tr) :
    """
f is a function of the form
f(rootData, subTreeOutputsList) -> treeOutput
and this applies it recursively to the tree.
"""
    def trf(x) :
        return tree_reduce(f,x)

    return f(tr[0],tuple(map(trf,tr[1])))

def nedges(tree) :
    """
Compute the number of edges in a tree.
"""
    return tree_reduce(lambda r,so : 1 + intsum(so), tree)

# This means if you run the module from the command line like so
# $ python <module-name>.py -v
# 
# the doctest module will automatically generate some tests from the 
# code in the docstrings of the module, and compare 
# the output to the printout. The -v flag means you will see 
# a report of all the tests (including those that passed, which should be all
# of them!).
#
# NOTE: some of the sub-modules have no tests. As it is, the function
# W2(d,l) in evenOddTrees.py runs almost all of the code (except 
# for some printing methods, etc.) and its output is quite fragile
# - so if the tests there pass it probably means everything's ok.
if __name__ == "__main__":
    print("running docstring tests...")
    print()
    import doctest
    doctest.testmod()
