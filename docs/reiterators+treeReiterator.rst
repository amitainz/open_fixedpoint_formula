Reiterables and tree generation
===============================

The reiterators module
----------------------

A Reiterable is like an iterator, that supports restarting the iteration from a given position. The main reason you need this is in order to generate multisets (which, in turn, are needed in order to generate trees recursively). This module defines the :py:class:`Reiterable` interface, and implements a few useful constructions of reiterables.  

Motivation
----------

A natural way to represent a multiset of a sequence :math:`S` is by a non-decreasing sequence of elements from the sequence:

.. math::

   i_1 \leq i_2 \cdots \leq i_n

to generate the next multiset, we simply try to increment :math:`i_n`. If we can't we increment :math:`i_{n-1}` and then, and this is crucial - *let* :math:`i_n` *start iterating from the new* :math:`i_{n-1}` *value*. This is easy enough to implement when :math:`S` is a list, but it is often very undesirable to hold the entire sequence :math:`S` in memory. Instead, the :py:class:`Reiterable` interface spells out what we need from our object in order to support taking multisets of it. 


To counterbalance the rather abstract (and possibly vague) definition of the interface below, you might want to look at the implementation of the simple :py:class:`seqReiterable` and :py:class:`seqReiterator`, constructed from a list, and then take a look at the implementation of :py:class:`sortedProductReiterable` (which really should be called "multisetReiterable" - that's the main reason to go through this hassle. Then you can take a look at :py:class:`treeReiterable` (in the next module) which iterates over isomorphism types of rooted labeled trees.


The Reiterable-Reiterator Interface
-----------------------------------

Because we're not saavy programs, the interface is just an agreement we've bound ourselves to, about how all the objects *should* behave. It is not enforced through any sort inheritance syntax or whatever, but if you see a class ending in Reiterable or Reiterator then that's an indication that this class supports the structure indicated below.

If you find the description below too abstract

A Reiterable class has
  - a method :py:obj:`reiter(status)` which accepts :py:obj:`status` (see below) and returns a Reiterator
  - a method :py:obj:`initialStatus()` which returns the initial status used for implementation of :py:obj:`__iter__()` which should just be :py:obj:`reiter(initialStatus())`.

A Reiterator class
 - is the class returned by calling :py:obj:`reiter` on a Reiterable.
 - supports method :py:obj:`getStatus()` which returns a status which can be fed to a reiter call to allow reiteration.
 - supports method :py:obj:`__next__()` which advances status one step.


One should think of :py:class:`Reiterable` as a generlization of a list, or other sequences. Something that represents an ordered set of elements, but is generated "on the fly". The :py:class:`Reiterator` is a kind of pointer to this sequence, though really it can "live on its own" even if the parent :py:class:`Reiterable` is destroyed, since it contains all the information necessary to complete the iteration. There may be more than one reiterator pointing to the same reiterable.

It is best to think of a :py:class:`Reiterator` as a pointer to "the spaces between the sequence elements". Indeed, if a reiterator generates n elements it will generally have n+1 different conditions (or statuses, see below). 

When next is called, the element between the current space and the next space is outputted and the pointer is advanced. If next is called on the last space, stopIteration() is raised - but for consistent results it is important keep this possibility in mind.

:py:obj:`status` is that portion of the information about the iteration that
changes through the iteration. Note there may still be some parameters
which are constant throughout the iteration (think of the endpoints of an
interval, say, as opposed to the actual place in the interval we're
pointng to). The proper way to deal with these constants is to store them as variables of the ReiterABLE, and *not* return them as part of the status. Remember that to produce a reiterator from a status the reiterable must be invoked, so the status is only meaningful in the context of some :py:class:`Reiterable` which must be kept in hand.

.. automodule:: reiterators
   :members:
   :undoc-members:

The treeReiterator module
-------------------------

This implements the :py:class:`treeReiterable` which generates isomorphism types of rooted labeled trees, generated recursively by a function which takes the data at a vertex and generates possibilities for data for the children.

.. automodule:: treeReiterator
   :members:
   :undoc-members:
