The evenOddTrees module
=======================
The main functions
------------------

.. automodule:: evenOddTrees
   :members: OGW, superPotential
   :noindex:

An Interlude to discuss Trivial Invariants
------------------------------------------

The first thing the function :py:func:`superPotential` does is check that

.. math::
    :label: good_dlk

    3\, d + 2 \, l + k \geq 3 \mbox{ and } (2\,l + k  + d) \equiv 3 \mod 4

If this condition does not hold, it immediately returns zero. 

Following `arXiv:1703.02950 <https://arxiv.org/pdf/1703.02950.pdf>`_, the equivariant open Gromov-Witten invariants are defined for any 3-tuple of non-negative integers :math:`(d,l,k)`. However, if :math:`k + d` is not odd (so the moduli spaces are not relatively orientable), or if they correspond to unstable disk configurations (which is when the first clause of :eq:`good_dlk` is violated) they are zero essentially *by definition*.

:math:`\mathbb{Z}/2` acts on :math:`\mathbb{C}P^2` fixing :math:`\mathbb{R}P^2` and thus also on the moduli spaces we're integrating over. Since we've decided
to consider here only conjugation-invariant integrands (see discussion in :ref:`conj-invariance`), we can consider the sign of conjugation on :eq:`invariants`. These additional vanishings account for the second clause of :eq:`good_dlk`. 

Let us remark that, unlike the the vanishing of under-determined invariants discussed :ref:`above <vanishing-invts>`, when the sign of conjugation condition is violated we do *not* get interesting relations since the contributions of conjugate fixed point components cancel in pairs. This is why we chose to simply return zero immediately in this case, but we do compute the under-determined invariants' zero.



even-odd tree generation
------------------------

.. automodule:: evenOddTrees
   :members: evenOddTreeFunc, genEvenOddOps
   :noindex:
	     
pretty printing of even odd trees
---------------------------------

.. automodule:: evenOddTrees
   :members: printTrees
   :noindex:

computation of the contribution of an even odd tree
---------------------------------------------------

.. automodule:: evenOddTrees
   :members: calcEOAutomorphisms, calcAmplitude
   :noindex:

All of the functions in this module
-----------------------------------

.. automodule:: evenOddTrees
   :members:
   :undoc-members:

