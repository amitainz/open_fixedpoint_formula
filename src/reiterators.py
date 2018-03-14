# -*- coding: cp1255 -*-
from builtins import object
import itertools

# this may be needed for python 2.7 compatibility
# from sets import *

import inspect


################ REITERATORS ##################

class seqReiterable(object) :
    """
This is the simple example of a reiterable.
Here status is just the index to the list. This example highlights why it's better to have status distinct from the value returned by the iterator - we don't want ot have to invoke .index in reiter.
"""
    def __init__(self, seq) :
        self.dat = seq

    def reiter(self, status) :
        return seqReiterator(self,status)

    def initialStatus(self) :
        return 0

    def __iter__(self) :
        return self.reiter(self.initialStatus())
    
class seqReiterator(object) :
    def __init__(self, reiterable, status) :
        self.dat = reiterable.dat
        self.idx = status
    def __iter__(self) :
        return self
    
    def __next__(self) :
        i = self.idx
        self.idx += 1
        try :
            return self.dat[i]
        except IndexError :
            raise StopIteration
    
    def status(self) :
        return self.idx

class productReiterable(object) :
    def __init__(self, *reiterables) :
        self.n = len(reiterables)
        if self.n == 0 :
            pass
        if self.n > 0 :
            self.firstReiterable = reiterables[0]
            self.restReiterable = productReiterable(*reiterables[1:])
        
    def reiter(self, status) :
        return productReiterator(self, status)

    def initialStatus(self) :
        # in the case of an empty product, the status is either True
        # or False. If it is True, we return () and change it to False,
        # if it is False we raise StopIteration().
        if self.n == 0 :
            return True
        
        if self.n > 0 :
            return {"firstValue" : None, "firstValueSet" : False, \
                    "firstStatus" : self.firstReiterable.initialStatus(),\
                    "restStatus" : self.restReiterable.initialStatus()}

    def __iter__(self) :
        return self.reiter(self.initialStatus())

class productReiterator(object) :
    def __init__(self, pr, status) :
        self.n = pr.n
        if self.n == 0 :
            self.trivialStatus = status
            
        if self.n > 0 :
            for k in ['firstValue', 'firstValueSet'] :
                setattr(self, k, status[k])
            
            self.firstReiter = pr.firstReiterable.reiter(status['firstStatus'])
            # we need to keep the reiterABLE of the rest as well, because it is
            # reset every time we advance the first reiterATOR.
            self.restReiterable = pr.restReiterable
            self.restReiter = self.restReiterable.reiter(status['restStatus'])
            
    def __iter__(self) :
        return self

    def __next__(self) :
        if self.n == 0 :
            if self.trivialStatus :
                self.trivialStatus = False
                return ()
            else :
                raise StopIteration()
        if self.n > 0 :
            if not self.firstValueSet :
                self.firstValue = next(self.firstReiter)
                self.firstValueSet = True
            try :
                restValue = next(self.restReiter)
                return (self.firstValue,) + restValue
            except StopIteration :
                self.firstValue = next(self.firstReiter)
                self.restReiter = iter(self.restReiterable)
                restValue = next(self.restReiter)
                return (self.firstValue,) + restValue
                
                # conceptually, maybe more elegant to put the try in a loop -
                # since if we had allowed "back dependencies" maybe a few are 
                # empty before we find a good one. but in a cartesian product 
                # this doesn't happen. so if the restValue can't be computed 
                # after a reset - we're done.
                
    def status(self) :
        if self.n == 0 :
            return self.trivialStatus
        if self.n > 0 :
            return {"firstValue" : self.firstValue, "firstValueSet" : self.firstValueSet, \
                    "firstStatus" : self.firstReiter.status(),\
                    "restStatus" : self.restReiter.status()}

class suffixReiterable(object) :
    """
from a reiterable and a status, produces a new reiterable which goes over
the old reiterable starting form status. 
"""
    def __init__(self, reiterable, beginStatus) :
        self.reiterable = reiterable
        self.beginStatus = beginStatus

## the reiterators of the suffixReiterable are in fact just reiterators
## of self.reiterable. The only things that changes is the initial status.
    def reiter(self, status) :
        return self.reiterable.reiter(status)

##    this is an example where the status must be understood in context,
##    since the statuses of the suffixReiterator look exactly the same
##    as the statuses of the original reiterator. 
    def initialStatus(self) :
        return self.beginStatus

    def __iter__(self) :
        return self.reiter(self.initialStatus())
            
        

class sortedProductReiterable(object) :
    """
this class produces "canonical representatives" for all the (unsorted) multisets of a given Reiterable.
It produces all sequences (i1,...,ik) of reiterator values with i1 <= i2 <= ... <= i3.
"""
    def __init__(self, reiterable, repeat) :
        self.repeat = repeat
        self.reiterable = reiterable

    def reiter(self, status) :
        return sortedProductReiterator(self, status)
        
    def initialStatus(self) :
        if self.repeat > 0 :
            return (False,)
        else :
            return False
        
    def __iter__(self) :
        return self.reiter(self.initialStatus())

class sortedProductReiterator(object) :
    def __init__(self, spr, status) :
        if spr.repeat < 0 :
            raise ValueError("sortedProductReiterator initialization: repeat is negative!")
        self.repeat = spr.repeat
        self.reiterable = spr.reiterable
        if self.repeat > 0 :
            if status[0] :
                (self.firstValueSet, self.firstValue, \
                 firstStat, \
                 self.restReiterable, \
                 restStat) = status
                self.restReiterator = self.restReiterable.reiter(restStat)
                self.firstReiterator = self.reiterable.reiter(firstStat)
            else :
                self.firstValueSet = False
                self.firstReiterator = iter(self.reiterable)
        if self.repeat == 0 :
            self.done = status
        
               
    def __iter__(self) :
        return self

    def __next__(self) :
        # the problem is we want to keep the firstValue,
        # but have the restReiterator a sortedProduct of suffixes of
        # defined by the status pointing to the space 
        # just *before* this firstValue (because we allow repetitions).
        # this means the first step is a bit awkward.
        
        if self.repeat > 0 :
            if self.firstValueSet :
                try :
                    restVal = next(self.restReiterator)
                    return (self.firstValue,) + restVal

                except StopIteration :
                    firstStatus = self.firstReiterator.status()
            # if firstValue is not set
            else :
                firstStatus = self.reiterable.initialStatus()

            # this code is run if either the restReiterator
            # has finished, or this is the first call
            # to next(). In the former case, we use the status
            # of firstReiterator to define the suffixes for the next run.
            # in the latter case, we use initialStatus() for this.
            self.restReiterable = \
                sortedProductReiterable(\
                    suffixReiterable(self.reiterable, firstStatus),\
                    self.repeat - 1)
            self.restReiterator = iter(self.restReiterable)            
            self.firstValue = next(self.firstReiterator)
            self.firstValueSet = True
            
            restVal = next(self.restReiterator)
            return (self.firstValue,) + restVal
        if self.repeat == 0 :
            if self.done :
                raise StopIteration
            else :
                self.done = True
                return ()   
        
    def status(self) :
        if self.repeat == 0 :
            return self.done
        if self.repeat > 0 :
            if self.firstValueSet :
                return (self.firstValueSet, self.firstValue, \
                        self.firstReiterator.status(), \
                        self.restReiterable,\
                        self.restReiterator.status())
            else :
                return (self.firstValueSet,)

    
class imapReiterable(object) :
    """
build a reiterable by composing a function on the output of a given
reiterable.
"""
    def __init__(self, reiterable, func) :
        self.reiterable = reiterable
        self.func = func
        
    def initialStatus(self) :
        return self.reiterable.initialStatus()
        
    def reiter(self, status) :
        return imapReiterator(self, status)
    
    def __iter__(self) :
        return self.reiter(self.initialStatus())

class imapReiterator(object) :
    def __init__(self, imr, status) :
        self.reiter = imr.reiterable.reiter(status)
        self.func = imr.func

    def __next__(self) :
        return self.func(next(self.reiter))

    def status(self) :
        return self.reiter.status()

    def __iter__(self) :
        return self
    
    
    
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
