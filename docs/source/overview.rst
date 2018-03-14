Overview of the Program
=======================

If you just want to compute the invariants, all you need to do is run 
:py:func:`evenOddTrees.OGW` (or :py:func:`evenOddTrees.W2`, for the Welschinger counts). It is my hope, however, that this extended documentation and the many comments in the code will invite the reader to look under the bonnet and study the fixed point formula itself.

The fixed-point formula can be recast in the following form [#details]_, which can be used as a schematic representation for the entire algorithm.

.. math::
   :label: fpf-concise

   OGW(d,l,k) = \overbrace{\sum_{T \in EO(d,l,k)} C_T \prod_{v \in V(T)}}^{\mbox{evenOddTrees.py}} \underbrace{\sum_{\tilde F} A_{\tilde F}}_{\mbox{weighted\_fp\_contribs.py}}.

Let us explain the main features of this formula; more details can be found in the specific modules and their docs. 

First, we have a sum over all *even-odd diagrams*. An even-odd diagram is built by gluing a pair of even-odd rooted trees, each of which looks something like this [#verbose-tree]_ :

>>> eo = evenOddTrees(3,0,7,0)
>>> len(eo)
13
>>> print(stree(eo[5], lambda sd,d : "  "*d + str(sd[-1]) + "\n"))
(0, 0, 1)
  (1, 0, 0)
    (0, 0, 2)
    (0, 0, 2)
    (0, 0, 1)
      (1, 0, 0)
    (0, 0, 1)
      (1, 0, 0)


Each vertex :math:`v \in V(T)` in the diagram :math:`T` (represented by a single line above) specifies a moduli space of marked stable disk maps [#combinatorics]_. For example, the moduli spaces corresponding to the first and second vertices are :math:`\overline{\mathcal{M}}_{0,k = 3, l = 0}(d = 0)` and :math:`\overline{\mathcal{M}}_{0,k = 0,l=0}(d=1)`, respectively. An edge in  :math:`T` roughly corresponds to a real node where the adjacent disks *could* have been attached. The first vertex has two incident edges - the one going down to the second vertex, and an edge going up (where it will be attached to the other "half" of the even-odd diagram). The second vertex has 5 neighbours - the first vertex as well as the vertices on lines 3,4,5 and 7. Note how the parity of the degrees (the first elements in each 3-tuple) alternates. This is why these diagrams are called even-odd. The enumeration of the diagrams is managed in :py:mod:`evenOddTrees`, and this is probably where one should begin exploring this program. This module also contains the source for the function :py:func:`evenOddTrees.superPotential` which implements the automorphism-weighted sum over the even odd diagrams of the product of vertex contributions. The main functions :py:func:`evenOddTrees.OGW`, :py:func:`evenOddTrees.W2`, are simple wrappers for :py:func:`evenOddTrees.superPotential`, and are also implemented in that module.

The contribution of each vertex :math:`v \in V(T)` is given by a sum :math:`\sum_{\tilde F} A_{\tilde F}`. This sum ranges over *fixed-point trees*. These are also tree-like diagrams, but are otherwise quite different than the even odd diagrams. They look something like this:

>>> fps = fixedPoints(6,7,3,0)
>>> len(fps)
3090
>>> printFP(fps[2500])
                        0 )...
ops((3, 7, 0, 1))[182]     <p0>||||            degs: (1, 1)
ops((1, 0, 0, 2))[0]          <p->             degs: (1,)
ops((0, 0, 0, 1))[0]             <p0>          degs: ()
ops((0, 3, 0, 0))[0]          <p+>|||          degs: ()

These are a variant of the diagrams appearing in the well-known fixed-point formula for closed curves. Edges correspond to positive energy components covering a line, and are labeled by their degree. Vertices map to a fixed point :math:`\{p_+,p_-,p_0\} \subset \mathbb{C}P^2`, as indicated. | is used to indicate an interior marking. Each fixed point diagram determines a fixed-point component for the :math:`S^1` action [#T2-actually]_ on the moduli space of disks corresponding to a single vertex of the even-odd diagram. The computation of :math:`\sum_{\tilde F} A_{\tilde F}` is managed by :py:mod:`weighted_fp_contribs` which is probably the second module you will want to examine.

To understand the program, it is probably best to start by looking at the two main modules :py:mod:`evenOddTrees` and :py:mod:`weighted_fp_contribs` together with their docs (which contain more mathematical detail). The other modules can then be read in context, when they're invoked, to get a more detailed picture of the implementation [#autodoc-problem]_.

Short Description of Modules
----------------------------

Here's a list of all the modules with a short description of what they do. The two main modules which we've already discussed are

* :py:mod:`evenOddTrees` - generate even odd diagrams and sum their contributions.
* :py:mod:`weighted_fp_contribs` - compute the contribution of a fixed point diagram, and sum over all such contributions (to get a single even-odd diagram's vertex factor).


The next two modules describe the geometry of the fixed points.

* :py:mod:`fixedPoints` - generate fixed point diagrams, specifying :math:`T^2` fixed point components inside the :math:`S^1` fixed point components of a given moduli space.
* :py:mod:`pincher` - compute a canonical, "totally pinched" fixed point diagram, representing the class of all fixed point digrams belonging to the same :math:`S^1` fixed-point component.


The following two modules provide a general tool for generating 
tree-like diagrams recursively.

* :py:mod:`reiterators` - sets up a general interface, :py:class:`reiterable`,  for an iterator-type object with an assumed extra functionality allowing you to restart the iteration from a given pointer object, :py:class:`reiterator` (hence the "re" in reiterable). Some simple reiterables are constructed to support enumeration of simple objects.

* :py:mod:`treeReiterator` - a reiterable which generates trees recursively.

The smaller modules below each fill in some specific gap.

* :py:mod:`jrrFormula` - computes the Pandharipande-Solomon-Tessler formula for descendent integrals of disks, used for integrating the inverse Euler over the contracted disk-moduli in the even degree vertices.

* :py:mod:`utils` - some utility functions.

* :py:mod:`dlk_partitions` - generates all the partitions of a given moduli specification (d,l,k), used in generation of even odd diagrams. 


Footnotes
---------

.. [#details] To see this, start from Remark 31 and Definition 26 (c) of `arXiv:1703.02950 <https://arxiv.org/pdf/1703.02950.pdf>`_. By acting first by :math:`Sym(r)` and then by :math:`Sym(\boldsymbol s)` we find that the fixed point formula can be written as a sum over isomorphism types :math:`T` of sorted odd-even trees, followed by a sum over the isomorphism types of fixed-point profiles :math:`\boldsymbol \phi \in \mathcal P (T)` for each :math:`T`. Turning now to the formula in Proposition 32, we see that the contribution of each fixed-point profile can be written as a product of contributions, one factor for each vertex. So we get 
   
   .. math::

      \sum_T \cdots \sum_{\boldsymbol{\phi}} \cdots \prod_{v \in V(T)} \cdots . 

   The fixed-point profiles for each :math:`T` are given as a product over the vertices of :math:`T` of "local choices" (cf. Definition 26 (b)),

   .. math:: 
      
      \boldsymbol \phi = (...,\phi_v,...)
       
   It is not hard to see the factor associated with each vertex in Proposition 32 depends only on these local choices. So we can interchange the inner sum and the product, and write the OGW invariants as 

   .. math::

      \sum_T \cdots \prod_{v \in V(T)} \cdots \sum_{\phi_v} \cdots . 

   This is essentially the formula we use, except that we also "unwind" the recursion defining the constraint correlator (see Proposition 6). This means that instead of summing over fixed point profiles, which just record information about the behaviour of the fixed point configuration near :math:`\mathbb{R}P^2`, we sum over "fixed point diagrams" which specify the entire fixed point configuration.

.. [#verbose-tree] The function :py:func:`treeReiterator.stree` is a general printing method for tree-like objects. If you run it without the additional formatting parameter, like so

   >>> print(stree(eo[5]))
   ((3, 0, 7, 0), 2, (0, 0, 1))
      ((3, 0, 6, 1), 2, (1, 0, 0))
	 ((0, 0, 2, 0), 0, (0, 0, 2))
	 ((0, 0, 2, 0), 0, (0, 0, 2))
	 ((1, 0, 1, 0), 0, (0, 0, 1))
	    ((1, 0, 0, 1), 0, (1, 0, 0))
	 ((1, 0, 1, 0), 0, (0, 0, 1))
	    ((1, 0, 0, 1), 0, (1, 0, 0))
    
   you can see what eo[5] is really made of. Above, we've only shown the right-most 3-tuple of integers in each line, which are the moduli specs for the vertices. They contain all of the relevant information, but this redundant representation allows us to compute more efficiently (e.g. we don't have to sum the degrees of the vertices in a subtree, because the total degree is recorded in the left 4-tuple).

.. [#combinatorics] Actually, it is probably best to think of all of the moduli spaces as having *symmetric markings*. In other words, when we refer to :math:`\overline{\mathcal{M}}_{0,k,l}(d)` we really mean the (homotopy, or stacky) quotient of the usual moduli space of disks by the action of the group :math:`\operatorname{Sym}(k) \times \operatorname{Sym}(l)` permuting the labels. This explains 

  * why :py:func:`evenOddTrees.superPotential` naturally computes :math:`\frac{1}{k!\,l!}` times the open Gromov-Witten invariant

  * why the even-odd and fixed-point diagrams only keep track of the *number* of boundary and interior labels, not the actual subsets (this is, of course, much more efficient)

  * why when we divide by the size of the automorphism group of a diagram, it is always of the diagram labeled by *numbers* and not by *subsets*.

.. [#T2-actually] In fact, each f.p. diagram specifies a :math:`T^2` -fixed -point component inside such an :math:`S^1` fixed-point component - this is discussed in some detail at the documentation for :py:mod:`weighted_fp_contribs`.

.. [#autodoc-problem] In various places we use the :py:class:`utils.Memoize` decorator to create lookup tables for functions. This causes Sphinx's autodoc not to show these functions at all, so you have to look at the source to see the complete picture.

