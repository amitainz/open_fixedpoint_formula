The fixed-point components
==========================

The fixedPoints module
----------------------
This takes care of generation of fixed point diagrams, each of which specifies a :math:`T^2` -fixed component inside each :math:`S^1` -fixed component.

.. automodule:: fixedPoints
   :members:
   :undoc-members:

The pincher module
------------------
This handles the mapping from a :math:`T^2` -fixed component to its :math:`S^1` -fixed component. Each :math:`S^1` component is represented by a "totally-pinched" fixed point diagram, in a unique way.

A fixed-point diagram is *totally-pinched* if there are no edges with even absolute degree :math:`2\,d` connecting :math:`p_+` and :math:`p_-`. Indeed, such an edge can be "pinched" inside the :math:`S^1` -fixed component to two edges of degree :math:`d` covering the line :math:`p_+ p_0` and the line :math:`p_0 p_-`, respectively, connected by a node at :math:`p_0`.

.. automodule:: pincher
   :members:
   :undoc-members:
