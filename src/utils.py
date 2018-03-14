from __future__ import division
from __future__ import print_function
# utilities
from builtins import map
from builtins import range
from builtins import object
from past.utils import old_div
from math import *
import sympy as sp
from numpy import *
from collections import defaultdict, Counter

# this may be needed for python 2.7 compatibility
# from sets import *

import time
from functools import reduce

u   = sp.Symbol('u')
wp,wm,w0 = sp.Symbol("w+"), sp.Symbol("w-"), sp.Symbol("w0")
eps = sp.Symbol("eps")

rat = sp.Rational

def tup_concat(ls) :
    return reduce(lambda x, y : x +y, ls, ())

def intsum(ls) :
    return reduce(lambda x, y : x +y, ls, 0)

def intprod(ls) :
    return reduce(lambda x, y: x*y, ls, 1)
def multinomial(*ls) :
    return old_div(factorial(intsum(ls)),intprod(list(map(factorial, ls))))

def str_concat(ls):
    return reduce(lambda x,y: x+y, ls, '')

def myprod(ls) :
    return reduce(lambda x , y : x *y, ls, 1)

from fractions import gcd
def mgcd(*ls) :
    return reduce(gcd, ls)

class Memoize(object):
    """
This class is a decorator for functions.
Putting @Memoize in front of a function f, like so:

.. code-block:: py 

    @Memoize
    def f :
       ...

means we build a lookup table for the function, so if we call it again with the same input we do not execute the code but instead use the stored value. It is used to run those parts of the code where a function is called with the same input many times, and the output is not too large (this typically happens inside recursions).
"""
    def __init__(self, f):
        self.f = f
        self.memo = {}
        self.memocount = defaultdict(int)
        self.shout_new_call = False
        
    def __call__(self, *args):
        if not args in self.memo:
            if self.shout_new_call :
                print("new call {}{}...".format(self.f.__name__, args), end=' ')
                s = time.time()
                self.memo[args] = self.f(*args)
                self.memocount[args] = 0
                f = time.time()
                print("finished in {} seconds".format(f-s))
            else :
                self.memo[args] = self.f(*args)
                self.memocount[args] = 0
        self.memocount[args] += 1
        return self.memo[args]

    def getcounters(self) :
        return list(self.memocount.values())

    def clearMemo(self) :
        self.memo = {}
        self.memocount = {}

def tuplePartitions(n,k) :
    if k < 0 :
        return
    if k == 0 :
        if n == 0 :
            yield ()
        if n != 0 :
            return

    for i0 in range(n + 1) :
        for pp in tuplePartitions(n - i0, k - 1) :
            yield (i0,) + pp

class frozenCounter(Counter):
    __slots__ = ('_hash',)
    def __hash__(self):
        rval = getattr(self, '_hash', None)
        if rval is None:
            rval = self._hash = hash(frozenset(iter(self.items())))
        return rval


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
