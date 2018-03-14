from __future__ import division
from __future__ import print_function



from builtins import str
from builtins import zip
from builtins import map
from builtins import range
from past.utils import old_div
from treeReiterator import *
from utils import *
from functools import reduce

# the purpose of this file is to generate all fixed point diagrams,
# specifying a T^2 fixed point component inside an S^1 fixed point component
# of a moduli of disks (see discussion in weighted_fp_contribs
# module documentation)

# Memoize is a decorator class defined in utils.py,
# it means we call the function once for each input,
# and keep a hash of the output. This is useful for
# time-consuming functions which have a modestly sized output
# and which are called recursively with each value
# accessed multiple times.
@Memoize
def fixedPoints(totalRelDeg, totalL, totalK, totalT) :
    return list(fixedPointsGen(totalRelDeg, totalL, totalK, totalT))

def fixedPointsGen(totalRelDeg, totalL, totalK, totalT) :
    """
this is a wrapper function generating all fixed point trees, including the special disc data.
"""

    # if there are bmp's, we must have a ghost disc.
    if totalK > 0 :
        if totalRelDeg % 2 == 1 :
            return

        discDeg = 0

        # totalT counts the dummy marked points; in this version
        # it is always zero.
        tr = treeReiterable(fixedPointTreeFunc,\
                            (old_div(totalRelDeg,2), totalL, totalT, 1))

        for subtree in tr :
            yield discDeg, totalK, subtree
        return

    # if there are no bmp's, choose where disc is headed.
    for mu in [0,2] :
        # we don't allow ghost discs unless there are bmp's, so
        # counting of discDeg starts from 1 or from 2 in this case.
        for discDeg in range(2 - (totalRelDeg)%2, totalRelDeg + 1, 2) :
            for subtree in treeReiterable(\
                fixedPointTreeFunc, \
                (old_div((totalRelDeg - discDeg),2), totalL, totalT, mu)) :
                yield discDeg, 0, subtree
    return


@Memoize
def fixedPointTreeFunc(subtreeData) :
    """
this function can be fed to treeReiterable.
subTreeData takes the form:
(totalD,totalL,totalT,rootMu)
where totalD is the total *absolute* degree :math:`\\in H_2(X)`
(twice the relative degree :math:`\\in H_2(X,L)`) of the subtree,
totalL and totalT=0 are the total number of IMP's and Dummy's respectively,
and rootMu is the image fixed point the contracted component is mapped to.

as specified by the treeReiterable interface, the function returns a
reiterable containing
( ((subtree1 data, mult1),...,(subtreek data, multk)) ,nodeDecor).

where nodeDecor is (subtreeData,nop,(mu1,degs),(mu2,degs)) where
subtreeData is the input subtreeData to the function,
mu2 are the two fp indices other than rootMu, and the accompanying degs
is a tuple describing the edge degrees.

"""
    return seqReiterable(list(genFixedPointTreeOps(subtreeData)))


def genFixedPointTreeOps(subtreeData) :
    totalD,totalL,totalT,rootMu = subtreeData
    counter = 0

    # we renumber the fixed points of CP2 (abstractly, for now) as
    # 0,1 and 2
    # where 0 corresponds to rootMu and 1 and 2, which are (rootMu + i)%3
    # respectively, are the other options. We refer to them as option 1 and
    # option2, and generally 1,2 indexes refer to these options,
    # and then 0 index means the root data
    # (e.g. imps0,imps1,imps2 or dumbs0,dumbs1,dumbs2).

    # we will create two partitions, one for option1 and one for option2.
    # but first we need to know how the (double) degrees divide between:
    # edges1, children1, edges2, children2
    # and how the imps divide between
    # root, children1, children2.

    counter = 0

    for nKids in range(totalD + 1) :
        for nKids1,nKids2 in tuplePartitions(nKids, 2) :
            # ddegsEdges1 is at least nKids1, and
            # ddegsEdges2 is at least nKids2, so we partition
            # totalDoubleDeg - nKids and then put back the nKids.
            for degsEdges1, degsEdges2, degsKids1, degsKids2 in \
                tuplePartitions(totalD - nKids,4) :

                # the first two partitions determine how many of the imps
                # and dummies will go on the subtrees labeled mu + 1 %3 and mu+2%3
                # and how many will remain on the current vertex (this is imps0 and dumb0).
                for (imps0,imps1,imps2) in tuplePartitions(totalL, 3) :
                    for dumb0,dumb1,dumb2 in tuplePartitions(totalT, 3) :
                        for multiPart1 in \
                            addPartition(\
                                addPartition(\
                                    addPartition(\
                                        addPartition(\
                                            dummyPartition(nKids1),\
                                            degsEdges1),\
                                        degsKids1),\
                                    imps1),\
                                dumb1) :
                                # here we add back the one to the ddeg of the edges.
                                multiPart1 = tuple(((x[0][0]+1, ((x[0][1:] + ((rootMu + 1)%3,)), x[1])) \
                                                    for x in multiPart1))
                                for multiPart2 in \
                                    addPartition(\
                                        addPartition(\
                                            addPartition(\
                                                addPartition(\
                                                    dummyPartition(nKids2),\
                                                    degsEdges2),\
                                                degsKids2),\
                                            imps2),\
                                        dumb2) :

                                            multiPart2 = tuple(((x[0][0]+1, ((x[0][1:] + ((rootMu + 2)%3,)), x[1])) \
                                                                for x in multiPart2))

                                            twoParts = multiPart1 + multiPart2
                                            yield (tuple(part[1] for part in twoParts), \
                                                   (subtreeData, counter, rootMu, \
                                                    imps0, dumb0,\
                                                    reduce(lambda x, y : x+y, \
                                                           ((part[0],)*part[1][-1]\
                                                            for part in twoParts),\
                                                           ())))
                                            counter += 1


### PAUSE: PRETTY PRINTING ###
# the next four functions are not invoked from the main body of the program,
# but they help visualize and play with fixed point trees.

def FPformatter(x, depth) :
    (subtreeData, counter, rootMu, imps0, dumb0,degs) = x
    return '{:27}'.format('ops({})[{}]'.format(subtreeData, counter)) + \
           '{:20}'.format(' '*(3*depth) + '<p{}>'.format(['+','0','-'][rootMu]) + \
           '|'*imps0 + 'T'*dumb0) + 'degs: ' + str(degs) + '\n'

def discFormatter(discDeg, totalK) :
    return ' '*24 + '{} )'.format(discDeg) +  '.'*totalK + '\n'


def strFP(xxx_todo_changeme) :
    (discDeg, totalK, subtree) = xxx_todo_changeme
    return discFormatter(discDeg, totalK) + stree(subtree, FPformatter)

def printFP(fp) :
    """
This prints fixed-point diagrams in a palatable format:

>>> fps = fixedPoints(6,7,3,0)
>>> len(fps)
3090
>>> printFP(fps[2500])
                        0 )...
ops((3, 7, 0, 1))[182]     <p0>||||            degs: (1, 1)
ops((1, 0, 0, 2))[0]          <p->             degs: (1,)
ops((0, 0, 0, 1))[0]             <p0>          degs: ()
ops((0, 3, 0, 0))[0]          <p+>|||          degs: ()

The first two lines represent the disk (which in the even degree is also the root of the tree part of the diagram). Note the condition that there are no interior markings on a vertex mapping to :math:`p_0` is not enforced.
The disk has 3 boundary markings in this case. It is connected by two edges to vertices mapping to :math:`p_-` and :math:`p_+`. The degree of the edges is shown at the parent, in this case both edges have absoulte degree :math:`1 \\in H_2(X)` (the map :math:`\\mathbb{Z} = H_2(X) \\to H_2(X,L) = \\mathbb{Z}` from the absolute to the relative degree is multiplication by 2, so the total degree of the diagram is :math:`2 \\cdot 3 = 6`, as it should be).


>>> fps = fixedPoints(5,2,0,0)
>>> len(fps)
102
>>> printFP(fps[35])
                        1 )
ops((2, 2, 0, 0))[17]      <p+>                degs: (1, 1)
ops((0, 1, 0, 1))[0]          <p0>|            degs: ()
ops((0, 1, 0, 2))[0]          <p->|            degs: ()

When the total degree is odd, the diagram is fairly similar except now the root of the tree part of the diagram represents the point where the disk irreducible component is attached to the rest of the tree (there may be a contracted component there, of course). This vertex can map either to :math:`p_+` or to :math:`p_-`. We also see that there must be no boundary markings on the disk component.
"""
    print(strFP(fp), end='')

### RESUME MAIN CODE ###

# the novel thing about our generation scheme -
# every subtree appears in its canonical form. This means
# that we can simply check if two subtrees are equal, not up to
# any permutation, to determine if they can be moved.
# since the root is distinguished, the levels are also determined.
def calcNAutomorphisms(tree) :
    ((subtreeData, counter, rootMu, imps0, dumb0,degs),kids) = tree
    if not kids :
        return 1

    rootData, subtrees = tree

    equivs = Counter(list(zip(degs, kids)))

    return myprod(list(map(factorial, list(equivs.values())))) * \
           myprod(list(map(calcNAutomorphisms, kids)))

def oldCalcNAutomorphisms(tree) :
    if not tree[1] :
        return 1
    equivs = {}
    for t in tree[1] :
        if t in equivs :
            equivs[t] += 1
        else :
            equivs[t] = 1
    return myprod(list(map(factorial, list(equivs.values())))) * \
           myprod(list(map(calcNAutomorphisms, tree[1])))


def fpCombFactor(xxx_todo_changeme1) :
    """
compute the combinatorial factor for the fixed point's "amplitude"
(product of relevant moduli integrals) in the total integral.

combFactor =
#{labeling of tree rep} / #{*unlabeled* tree automorphism} / deg product

note the count of number of labelings ignores possible tree automorphisms,
this is accounted for by the fact we divide by unlabeled tree automorphisms:
that is, automorphisms mapping a vertex to a vertex with the same *number*
of markings, but these markings are unlabeled.

the label-preserving automorphisms combine with the product of the degrees to
determine the degree of the cover of the map from the product of modulis we're
integrating over to the fixed suborbifold with corners of the moduli of maps.

Note: the integral of the open discs' moduli already accounts for the possible
orderings of the boundary marked points.

Unimportant note:
This can probably be simplified according to the symmetric-markings yoga
(indeed you see that the multinomials form a telescopic product).
"""
    (discDeg, totalK, tree) = xxx_todo_changeme1
    def degProdReducer(vertexData, subDegProds) :
        # first element is the degs of the outgoing edges.
        return myprod(tuple(vertexData[-1]) + tuple(subDegProds))

    def LTvecReducer(vertexData, subLTvecs) :
        (subtreeData, counter, rootMu, imps0, dumb0,degs) = vertexData
        return ((imps0,dumb0),) + reduce(lambda x, y : x + y, subLTvecs, ())


    degsProd   = tree_reduce(degProdReducer, tree)
    if discDeg > 0 :
        degsProd *= discDeg

    impsVec, dumbVec = list(zip(*tree_reduce(LTvecReducer, tree)))
    nLabelings = multinomial(*impsVec) * multinomial(*dumbVec)

    nAut       = calcNAutomorphisms(tree)
    return rat(nLabelings, degsProd * nAut)

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
