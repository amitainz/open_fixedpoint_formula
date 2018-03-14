from __future__ import print_function
from builtins import map
from builtins import range
from weighted_fp_contribs import *
from treeReiterator import *
from utils import *
from dlk_partitions import *

def W2(d,l) :
    """
A wrapper for computing the Welschinger counts as a special case of the open Gromov-Witten invariants. It computes the signed count of degree  :py:obj:`d` real rational curves through a generic configuration of :py:obj:`l` conjugate pairs of complex points and a suitable number :py:obj:`k` (computed from :py:obj:`d` and :py:obj:`l`) of real points.

For example, the following calls show that there is a unique conic through 5 real points, and 40 degree 4 curves through 3 complex pairs and 5 real points.

>>> W2(2,0)
-1

>>> W2(4,3)
40
"""
    if 2*l > 3*d-1 :
        raise ValueError("too many interior constraints!")
    k = 3*d-1-2*l
    return 2**l * OGW(d,l,k)/2

def OGW(d,l,k) :
    """
This is the only function you need to call to compute an equivariant open Gromov-Witten invariant with the fixed point formula. It is really just a wrapper for :py:func:`superPotential`.
"""
    p = int((2*l+k-3*d+1)/2)
    c = factorial(l)*factorial(k)
    return c*superPotential(d,l,k)*u**p

def superPotential(totalD,totalL,totalK) :
    """
:py:func:`superPotential` computes a coefficient of the equivariant superpotential (up to some power of :math:`\sqrt{-1}`).

More precisely, given three non-negative integers

* :py:obj:`totalD` = :math:`d \\in \mathbb{Z} = H_2(\mathbb{C}P^2,\mathbb{R}P^2)`, the total degree

* :py:obj:`totalL` = :math:`l`, the total number of interior markings, and

* :py:obj:`totalK` = :math:`k`, the total number of boundary markings,

:py:func:`superPotential` (:math:`\\beta,l,k`) computes a sympy rational number :math:`c` such that

.. math::

    \\frac{1}{k! \\, l!} OGW(d,l,k) = c \\, u^\\frac{k + 2\\, l - (3 \\, \\beta-1)}{2}.

See the documentation for more details on the notation.
"""
    # check stability and sign of conjugation
    if not goodDLK_2(totalD, totalL,totalK) :
        print("superPotential: moduli is unstable or conjugation sign is trivial, returning a `trivial` zero.")
        return 0

    # check relative orientability of moduli space
    # (also, if this fails there's a mismatch between the parity
    # of the dimension of the moduli space and the degree of the form)
    if (totalD + totalK) % 2 != 1 :
        print("superPotential: moduli is not relatively oriented, returning at `rivial` zero.")
        return 0


    # the contribution of the unique diagram with a single vertex,
    # and no edges, is computed first.
    res = wfp_nou(totalD,totalL,totalK,0)*rat(1,factorial(totalL)*factorial(totalK))

#    print("{}".format(res))

    # We compute the sum of contributions of diagrams with more than one edge.
    #
    # This is the same as sum over diagrams with a highlighted
    # edge*, of the contribution of the underlying diagram divided by the number of
    # edges. Each edge connects an even degree vertex to an odd degree vertex,
    # so if we sever the highlighted edge we obtain an ordered pair of rooted
    # trees.
    #
    # In practice, we go over all ways to write the 3-tuple
    # (totalD,totalL,totalK) as a sum (dodd,lodd,kodd) + (deven,leven,keven),
    # and compute the contribution of all pairs of rooted trees, one with an
    # even-degree root and one with an odd root. The contribution is the product
    # of the two amplitudes, divided by the product of automorphism groups,
    # and the total number of edges.
    #
    # It will probably be more efficient to change the sum of products
    # to a product of sums.
    #
    # * all sums are automorphism-weighted; highlighting an edge
    # means we also reduce the automorphism group.
    for dodd in range(totalD+1):
        for lodd in range(totalL +1):
            for kodd in range(totalK+1):
                deven = totalD - dodd
                leven = totalL - lodd
                keven = totalK - kodd

                # we make sure that both sides of the edge
                # are stable and have a good sign of conjugation.
                if not (goodDLK_2(dodd,lodd,kodd+1) and\
                        goodDLK_2(deven,leven,keven+1)) :
                    continue

                # iterate over the rooted tree on the odd side
                for trodd in treeReiterable(evenOddTreeFunc, \
                                            (dodd,lodd,kodd, 1)) :
                    # iterate over the rooted tree on the even side
                    for treven in treeReiterable(evenOddTreeFunc,\
                                                 (deven,leven,keven, 0)) :
                        # compute the product of amplitudes,
                        # divide by the size of automorphism group of the
                        # highlighted diagram (which is the product of the
                        # automorphism groups of the trees) and by
                        # the total number of edges, including the severed edge.
                        tmp = calcAmplitude(trodd) * calcAmplitude(treven) *\
               rat(1,nedges(trodd) + nedges(treven) - 1) *\
               rat(1,calcEOAutomorphisms(trodd) * calcEOAutomorphisms(treven))

#                        print("{}".format(tmp))

                        res += tmp
    return res



################ EVEN-ODD TREES  ##########
# see more about this in the documentation.

def evenOddTreeFunc(subtreeData) :
    """
This function complies with the :py:class:`treeReiterable` interface (see the documentation there for details). It represents the recursive step in the construction of an even-odd rooted tree.

The input subtreeData is a tuple :py:obj:`(totalD,totalL,totalK,rootDegParity)`,
where :py:obj:`totalD`, :py:obj:`totalL` and :py:obj:`totalK` are the total degree, number of interior markings and number of boundary markings of the tree to be generated, and the :py:obj:`rootDegParity` is the parity of the degree of the root (0 or 1).

As specified by the :py:class:`treeReiterable` interface, the function returns an
iterable generating tuples of tree specifications of the form:
:py:obj:`(  ((subTreeData, multiplicity),(...),...),  nodeDecoration)`.

where

* :py:obj:`subTreeData` has the same format as the input (indeed, it is fed to **evenOddTreeFunc** recursively when generating trees),

* :py:obj:`multiplicity` is a positive integer representing the number of subtrees of the specified type.

* :py:obj:`nodeDecoration` is a nested tuple of integers formatted like so:
  :py:obj:`(totalD,totalL,totalK,rootDegParity), counter, (droot,lroot,kroot))`

Here :py:obj:`(totalD,totalL,totalK,rootDegParity)` is a copy of the initial :py:obj:`subtreeData`, and :py:obj:`counter` is the index we're at in the iterable (so the first tree specification has :py:obj:`counter` = 0, the second has :py:obj:`counter` = 1, etc.). :py:obj:`(droot,lroot,kroot)` represent the degree and number of interior and boundary markings on the root vertex of the tree.
"""
    (totalD,totalL,totalK, parity) = subtreeData
    return seqReiterable(genEvenOddOps(totalD,totalL,totalK, parity))


def genEvenOddOps(totalD,totalL,totalK, rootDegParity) :
    """
This generates a list of tuples of the form
(  ((subTreeData, multiplicity),(...),...),  nodeDecoration)
as explained in the documentation for :py:func:`evenOddTreeFunc`.
"""
    # verifies that the tuple is valid in terms of stability and sign of
    # conjugation.
    if not goodDLK_2(totalD,totalL,totalK+1) :
        # you're not supposed to get here
        raise ValueError("in genEvenOddOps: {} have bad sign of conjugation!".format((totalD,totalL,totalK+1)))

    counter = 0
    ops = []

    # go over the allowed degrees for the root, even or odd as specified
    # by rootDegParity.
    for droot in range(rootDegParity, totalD+1,2) :
        # interior markings on the root
        for lroot in range(totalL+1) :
            # boundary markings on the root. This does *not*
            # include markings coming from edges.
            for kroot in range(totalK +1) :
                # you cannot have any boundary markings when the degree
                # is odd (no fixed points in that case). Though
                # of course you can have special points corresponding
                # to edges, all outgoing.
                if rootDegParity == 1 and kroot > 0 :
                    continue

                # the "resources" left over to
                # be partitioned among the subtrees.
                leftoverD = totalD - droot
                leftoverL = totalL - lroot
                leftoverK = totalK - kroot

                # this goes over all partitions of of the leftover resources
                # into subtrees. A partition is a kind of multiset,
                # a list of ((d,l,k),m) where
                # m is the multiplicity, or number of times each subtree
                # with (d,l,k) appears.
                for dlk_part in dlk_partitions(leftoverD,leftoverL,leftoverK) :
                    # the total number of special points on the root =
                    # (edges going down to the subtrees)  +
                    # + (boundary marked points) +
                    # (the edge going up)
                    ktag = intsum([m for x,m in dlk_part]) + kroot +1
                    # if the root is even degree, all edges are incoming
                    # and ...
                    if (droot % 2 == 0) :
                        # we check that the moduli space is relatively
                        # oriented ...
                        if ktag % 2 == 0 :
                            continue
                        # ... and stable.
                        if droot == 0 and \
                           ((ktag != 3) or lroot > 0) :
                            continue
                    # since trees are even-odd, the parity of the
                    # subtrees is always opposite that of the root.
                    newParity = 1-rootDegParity

                    # cast the generated result in the format
                    # (((subTreeData, multiplicity),(...),...), nodeDecoration)
                    # as dictated by the treeReiterator interface
                    # (see also evenOddTreeFunc(subtreeData) documentation).
                    ops.append((tuple((dlk + (newParity,),mult) \
                                      for dlk,mult in dlk_part),\
                                ((totalD,totalL,totalK,rootDegParity),\
                                 counter, (droot,lroot,kroot))))
                    counter += 1
    return ops


#### TREE PRINTING - Not Used by superPotential ####
# the following pair of functions should help examine and play
# with even odd diagrams, they're not invoked in the
# computation itself.
def evenOddTrees(d,l,k,rp) :
    """
    return a list from the treeReiterable generating even-odd trees.
"""
    return list(treeReiterable(evenOddTreeFunc, (d,l,k,rp)))

def printTrees(d,l,k,rp, n=None) :
    """
A nice print of all the trees generated by treeReiterable.
Each tree is displayed using :py:func:`stree` from the :py:mod:`treeReiterator` module;
so a tree is represented by a block of lines in DFS order,
with each line containing the nodeDecoration of a single vertex,
indented by its depth.

Recall the nodeDecoration for an evenOddTree is:
(totalD,totalL,totalK,rootDegParity), counter, (droot,lroot,kroot))

For example:

.. code-block:: python

    >>> printTrees(3,0,7,0)
    found 13 trees
    --------------
    tmp[0] is:
    ((3, 0, 7, 0), 0, (0, 0, 0))
       ((1, 0, 0, 1), 0, (1, 0, 0))
       ((2, 0, 7, 1), 0, (1, 0, 0))
          ((0, 0, 2, 0), 0, (0, 0, 2))
          ((0, 0, 2, 0), 0, (0, 0, 2))
          ((0, 0, 2, 0), 0, (0, 0, 2))
          ((1, 0, 1, 0), 0, (0, 0, 1))
             ((1, 0, 0, 1), 0, (1, 0, 0))
    amplitude: -2; automorphisms: 48
    <BLANKLINE>
    tmp[1] is:
    ((3, 0, 7, 0), 0, (0, 0, 0))
       ((1, 0, 0, 1), 0, (1, 0, 0))
       ((2, 0, 7, 1), 1, (1, 0, 0))
          ((0, 0, 2, 0), 0, (0, 0, 2))
          ((1, 0, 5, 0), 0, (0, 0, 1))
             ((1, 0, 4, 1), 0, (1, 0, 0))
                ((0, 0, 2, 0), 0, (0, 0, 2))
                ((0, 0, 2, 0), 0, (0, 0, 2))
    amplitude: -2; automorphisms: 16
    <BLANKLINE>
    tmp[2] is:
    ((3, 0, 7, 0), 1, (0, 0, 0))
       ((1, 0, 4, 1), 0, (1, 0, 0))
          ((0, 0, 2, 0), 0, (0, 0, 2))
          ((0, 0, 2, 0), 0, (0, 0, 2))
       ((2, 0, 3, 1), 0, (1, 0, 0))
          ((0, 0, 2, 0), 0, (0, 0, 2))
          ((1, 0, 1, 0), 0, (0, 0, 1))
             ((1, 0, 0, 1), 0, (1, 0, 0))
    amplitude: -2; automorphisms: 16
    <BLANKLINE>
    tmp[3] is:
    ((3, 0, 7, 0), 2, (0, 0, 1))
       ((3, 0, 6, 1), 1, (1, 0, 0))
          ((0, 0, 2, 0), 0, (0, 0, 2))
          ((0, 0, 2, 0), 0, (0, 0, 2))
          ((0, 0, 2, 0), 0, (0, 0, 2))
          ((2, 0, 0, 0), 0, (0, 0, 0))
             ((1, 0, 0, 1), 0, (1, 0, 0))
             ((1, 0, 0, 1), 0, (1, 0, 0))
    amplitude: -2; automorphisms: 96
    <BLANKLINE>
    tmp[4] is:
    ((3, 0, 7, 0), 2, (0, 0, 1))
       ((3, 0, 6, 1), 1, (1, 0, 0))
          ((0, 0, 2, 0), 0, (0, 0, 2))
          ((0, 0, 2, 0), 0, (0, 0, 2))
          ((0, 0, 2, 0), 0, (0, 0, 2))
          ((2, 0, 0, 0), 1, (2, 0, 0))
    amplitude: 1; automorphisms: 48
    <BLANKLINE>
    tmp[5] is:
    ((3, 0, 7, 0), 2, (0, 0, 1))
       ((3, 0, 6, 1), 2, (1, 0, 0))
          ((0, 0, 2, 0), 0, (0, 0, 2))
          ((0, 0, 2, 0), 0, (0, 0, 2))
          ((1, 0, 1, 0), 0, (0, 0, 1))
             ((1, 0, 0, 1), 0, (1, 0, 0))
          ((1, 0, 1, 0), 0, (0, 0, 1))
             ((1, 0, 0, 1), 0, (1, 0, 0))
    amplitude: -2; automorphisms: 16
    <BLANKLINE>
    tmp[6] is:
    ((3, 0, 7, 0), 2, (0, 0, 1))
       ((3, 0, 6, 1), 4, (1, 0, 0))
          ((0, 0, 2, 0), 0, (0, 0, 2))
          ((2, 0, 4, 0), 0, (0, 0, 0))
             ((1, 0, 0, 1), 0, (1, 0, 0))
             ((1, 0, 4, 1), 0, (1, 0, 0))
                ((0, 0, 2, 0), 0, (0, 0, 2))
                ((0, 0, 2, 0), 0, (0, 0, 2))
    amplitude: -2; automorphisms: 16
    <BLANKLINE>
    tmp[7] is:
    ((3, 0, 7, 0), 2, (0, 0, 1))
       ((3, 0, 6, 1), 4, (1, 0, 0))
          ((0, 0, 2, 0), 0, (0, 0, 2))
          ((2, 0, 4, 0), 1, (0, 0, 1))
             ((2, 0, 3, 1), 0, (1, 0, 0))
                ((0, 0, 2, 0), 0, (0, 0, 2))
                ((1, 0, 1, 0), 0, (0, 0, 1))
                   ((1, 0, 0, 1), 0, (1, 0, 0))
    amplitude: -2; automorphisms: 4
    <BLANKLINE>
    tmp[8] is:
    ((3, 0, 7, 0), 2, (0, 0, 1))
       ((3, 0, 6, 1), 4, (1, 0, 0))
          ((0, 0, 2, 0), 0, (0, 0, 2))
          ((2, 0, 4, 0), 3, (2, 0, 4))
    amplitude: 32; automorphisms: 48
    <BLANKLINE>
    tmp[9] is:
    ((3, 0, 7, 0), 2, (0, 0, 1))
       ((3, 0, 6, 1), 5, (1, 0, 0))
          ((1, 0, 1, 0), 0, (0, 0, 1))
             ((1, 0, 0, 1), 0, (1, 0, 0))
          ((1, 0, 5, 0), 0, (0, 0, 1))
             ((1, 0, 4, 1), 0, (1, 0, 0))
                ((0, 0, 2, 0), 0, (0, 0, 2))
                ((0, 0, 2, 0), 0, (0, 0, 2))
    amplitude: -2; automorphisms: 8
    <BLANKLINE>
    tmp[10] is:
    ((3, 0, 7, 0), 2, (0, 0, 1))
       ((3, 0, 6, 1), 6, (3, 0, 0))
          ((0, 0, 2, 0), 0, (0, 0, 2))
          ((0, 0, 2, 0), 0, (0, 0, 2))
          ((0, 0, 2, 0), 0, (0, 0, 2))
    amplitude: -182/3; automorphisms: 48
    <BLANKLINE>
    tmp[11] is:
    ((3, 0, 7, 0), 6, (2, 0, 3))
       ((1, 0, 4, 1), 0, (1, 0, 0))
          ((0, 0, 2, 0), 0, (0, 0, 2))
          ((0, 0, 2, 0), 0, (0, 0, 2))
    amplitude: 32; automorphisms: 48
    <BLANKLINE>
    tmp[12] is:
    ((3, 0, 7, 0), 8, (2, 0, 7))
       ((1, 0, 0, 1), 0, (1, 0, 0))
    amplitude: 6144; automorphisms: 5040
    <BLANKLINE>

We examine this example more closely in the documentation for :py:func:`evenOddTrees.calcEOAutomorphisms`.
"""
    tmp = evenOddTrees(d,l,k,rp)
    N = len(tmp)

    print("found {} trees".format(N))
    print("--------------")
    if n is None :
        n= N

    for i in range(n) :
##            print "{}th tree".format(i)
        tr = tmp[i]
        print("tmp[{}] is:".format(i))
        # stree is defined in treeReiterator.py, it turns
        # tree into a string of lines, each line represents a node
        # indented by depth
        print(stree(tmp[i]), end='')
        print("amplitude: {}; automorphisms: {}".format(calcAmplitude(tr),\
                                                        calcEOAutomorphisms(tr)))
        print()

############ AMPLITUDE COMPUTATION ##############################

def calcEOAutomorphisms(tree) :
    """
Computes the size of the automorphism group of the input :py:obj:`tree`.

We think of :py:obj:`tree` as a rooted tree, whose vertices are decorated by degrees and which has additional "exterior" edges of two distinct types, corresponding to the boundary and the interior markings, but which are otherwise indistinguishable and thus contribute to the symmetries of the tree.

For example:

.. code-block:: python

    >>> printTrees(3,0,7,0) # doctest: +SKIP
    found 13 trees
    --------------
    tmp[0] is:
    ((3, 0, 7, 0), 0, (0, 0, 0))
       ((1, 0, 0, 1), 0, (1, 0, 0))
       ((2, 0, 7, 1), 0, (1, 0, 0))
          ((0, 0, 2, 0), 0, (0, 0, 2))
          ((0, 0, 2, 0), 0, (0, 0, 2))
          ((0, 0, 2, 0), 0, (0, 0, 2))
          ((1, 0, 1, 0), 0, (0, 0, 1))
             ((1, 0, 0, 1), 0, (1, 0, 0))
    amplitude: -2; automorphisms: 48

    ...

We see that the root is connected to two distinct subtrees.
The subtree starting at ((2, 0, 7, 1), 0, (1, 0, 0)) is connected to 4
subtrees, three of which are identical, each containing a pair of boundary markings. Thus the automorphism group is :math:`3! \\times 2^3 = 48` in this case.
"""
    (root,subtrees) = tree

    dr,lr,kr = root[-1]

    # Compute the equivalence classes of the subtrees.
    # Note the sizes of these equivalence classes are smaller
    # than the multiplicites of the subTreeData's
    # as seen from the root vertex. Two subtrees are equivalent
    # only if they are equal "all the way down". Here it is important
    # that we keep the subtrees ordered lexicographically.
    equivs = Counter(subtrees)

    # the automorphism group is the semidirect product
    # of the automorphism groups of the subtrees, with the product of the
    # symmetry groups
    # S_{lr} x S_{kr} x S_{equivs}
    # where S_{equivs} = S_{C1} x S_{C2} x ...
    # is the automorphism of the set of subtrees with its equivalence relation.
    return intprod(list(map(factorial, list(equivs.values())))) * \
           factorial(kr) * factorial(lr) *\
           intprod(list(map(calcEOAutomorphisms, subtrees)))

def calcAmplitude(tree) :
    """
Computes the amplitude of the evenOddTree :py:obj:`tree`, which is just :math:`\\prod_v \\sum_{\\tilde F} A_{\\tilde F}`.

(so this does *not* include the symmetry factor).
"""
    (root,subtrees) = tree

    ((td,tl,tk,rootDegParity),ctr, (rootD,rootL,rootK)) = root

    rootE = len(subtrees) +1
    if rootDegParity == 0 :
        rootContrib = wfp_nou(rootD,rootL,rootK+rootE,0)
    else :
        rootContrib = wfp_nou(rootD,rootL,rootK,rootE)

    return rootContrib*intprod(list(map(calcAmplitude,subtrees)))


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
