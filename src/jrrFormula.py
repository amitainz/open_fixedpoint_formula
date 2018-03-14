from __future__ import division
from builtins import range
from past.utils import old_div
from math import *
from numpy import *
from collections import Counter
from utils import *
from functools import reduce

def addTau0s(tauPows, dummies) :
    if dummies == 0 :
        return tuple(tauPows)
    
    if len(tauPows) > 0 and tauPows[0][0] == 0 :
        return ((0,tauPows[0][1] + dummies),) + tuple(tauPows[1:])
    return ((0,dummies),) + tuple(tauPows)


def psiPows2tauPows(psiPows) :
    return list(Counter(psiPows).items())

def tauPows2psiPows(tauPows) :
    return tuple(sorted(reduce(lambda x, y : x + y, ([(x,)*y for x,y in tauPows]), ())))

def dec1(tup, i) :
    return tup[:i] + (tup[i] - 1,) + tup[i+1:]

# if you want to write this more efficiently with multinomial coefficients,
# note that after decreasing some psiPow with a few tau0's you may get another
# tau0, it's like a chain reaction...

def jrrFormula(tauPows, sigPow) :
    psiPows = tauPows2psiPows(tauPows)
    l = len(psiPows)
    k = sigPow
    if l == 0 and k == 3 :
        return 2
    
    deg = 2*intsum(psiPows)
    if 2*l + k -3 != deg :
        raise ValueError("in OGWformula, there's a dim / deg mismatch.")
    
    if psiPows[0] == 0 :
        if l == 1 and k == 1 :
            return 1
        return intsum([jrrFormula(psiPows2tauPows(dec1(psiPows, i)[1:]),sigPow) \
                    for i in range(1,len(psiPows)) if psiPows[i] > 0])
    
    
    
    val = 2**(old_div((k-1),2))*( \
        old_div(factorial(intsum([2*a for a in psiPows]) - l +1),\
        intprod([doubleFactorial(2*a -1) for a in psiPows])) )
    return val
        
def doubleFactorial(n) :
    if n % 2 != 1 :
        raise ValueError("double factorial is defined for odd numbers only.")
    if n == 1 :
        return 1
    return n*doubleFactorial(n-2)
            
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
