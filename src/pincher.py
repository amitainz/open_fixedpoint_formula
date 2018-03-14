from __future__ import division
from builtins import zip
from builtins import range
from past.utils import old_div
from fixedPoints import *

def pinch(fp) :
    """
go over each fixed point tree, and wherever an even degree edge between
p+ and p- is encountered, pinch it, thus obtaining the unique totally
pinched representative to the class of T^3 fixed point components correspondingcontained in the same S^1 fp connected component."""
    return fp[:2] + (tree_reduce(pincher, fp[2]),)

def pincher(vertexData,pinchedSubtrees) :
    """
This function pinches at the root, if necessary.
"""
    # the only thing that's going to change in the top vertex data
    # is degs and the counter.
    (tmp, counter, rootMu, imps0, dumb0,degs) = vertexData
    totalD,totalL,totalT,rootMu = tmp
    
    if rootMu == 1 :
        return fixRootCounter((vertexData, pinchedSubtrees))

    newDegs = []
    subTrees = []
    
    for i in range(len(degs)) :
        (stmp, scounter, srootMu, simps0, sdumb0,sdegs) = pinchedSubtrees[i][0]
        stotalD,stotalL,stotalT,srootMu = stmp

        # if the edge is not pinchable, don't change anything...
        if srootMu != 2-rootMu or degs[i] % 2 == 1 :
            newDegs.append(degs[i])
            subTrees.append(pinchedSubtrees[i])
            continue
        
        
        # ... otherwise, pinch it. The new subtree is of the form:
        # p0
        #    pinchedSubtree
        # with the degree half the original degree.
        newDeg = old_div(degs[i], 2)
        newDegs.append(newDeg)
        newTmp = (stotalD + newDeg,stotalL,stotalT, 1)
        vec = (1, 0, 0, (newDeg,))
        newSubtree = fixRootCounter(((newTmp,-1) + vec, (pinchedSubtrees[i],)))
        subTrees.append(newSubtree)
    
    mtree = (vertexData[:-1] + (tuple(newDegs),), tuple(subTrees))
    return fixRootCounter(mtree)


def fixRootCounter(tree) :
    """
recalculate the counter for the root after we pinch.
"""
    (tmp, counter, rootMu, imps0, dumb0,degs) = tree[0]
    totalD,totalL,totalT,rootMu = tmp
    
    # construct a mapping from subtree types (by root data 4-tuple, like tmp)
    # to the sorted list of subtrees of that type.
    subtreesByType = defaultdict(lambda : [])
    for (deg,st) in zip(degs,tree[1]) :
        subtreesByType[(deg,st[0][0])].append(st)
    for k in list(subtreesByType.keys()) :
        subtreesByType[k].sort()
    
    ops = list(genFixedPointTreeOps(tmp))
    
    for newCtr in range(len(ops)) :
        # note deg may have been shuffled...
        if ops[newCtr][1][2:-1] != (rootMu, imps0, dumb0) :
            continue
        
        std = ops[newCtr][0]
        goodDegs = ops[newCtr][1][-1]
        # try to construct an ordered subtrees list, from the subtrees
        # we have, by the prescription defined by std and degs.
        
        orderedSubtrees = []
        degPtr = 0
        for (x,n) in std :
            # the sequence of n degrees is all supposed to be the same d.
            d = goodDegs[degPtr]
            degPtr += n

            if len(subtreesByType[(d,x)]) != n :
                break
            orderedSubtrees += subtreesByType[(d,x)]

        # this is executed only if all the counters agree, i.e. if
        # there are NO mismatches.
        else :
            return ((tmp, newCtr, rootMu, imps0, dumb0,goodDegs), tuple(orderedSubtrees))

    # if there's no newCtr which works, we raise an exception
    else :
        raise ValueError("fixRootCounter: can't find counter for given config.")
    
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
