from builtins import range
from utils import *


@Memoize
def goodDLK_2(d,l,k) :
    """
check for parity (= orientability), stability, and sign of conjugation.
See more about this in documentation.
"""
    if (d == 0) and ((l != 0) or (k != 3)) :
        return False
    return ((2*l + k - (3*d -1)) % 4 in [0,3])


# a partition is just a multiset. We like to keep it ordered
# so that it works well with tree Iterator too.
def dlk_partitions(totalD, totalL, totalK,\
                    minD = 0,minL = 0,minK = 0) :
    """
return all the dlk_partitions in a list.
"""
    partitions = []
##    if goodDLK_2(totalD,totalL,totalK+1) and totalE >= 1:
##        partitions.append((((totalD,totalL,totalK,totalE-1),1),))
    if (totalD,totalL,totalK) == (0,0,0) :
        return [()]
    for d1 in range(minD, totalD +1):
        loD = totalD - d1
        for l1 in range(minL, totalL +1):
            loL = totalL - l1
            for k1 in range(minK, totalK +1):
                loK = totalK - k1
                if not goodDLK_2(d1,l1,k1+1) :
                    continue
                
                rest = dlk_partitions(loD,loL,loK,d1,l1,k1)
                partitions += [updatePartition(r, (d1,l1,k1)) for r in rest]
                # this updating of the lower bound of iterations
                # is because bound is on lexicographical order.
            minK = 0
        minK = 0
        minL = 0
    return partitions
                    
def updatePartition(partition, y) :
    if not partition :
        return ((y,1),)
    x, m = partition[0]
    if x == y :
        return ((y,m+1),) + partition[1:]
    return ((y,1),) + partition
            

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
