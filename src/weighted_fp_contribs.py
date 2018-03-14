from __future__ import division
from __future__ import print_function
from builtins import input
from builtins import str
from builtins import map
from builtins import zip
from builtins import range
from past.utils import old_div
from pincher import *

# may be needed for python 2.7 compatiblity
# from sets import *


# the formula for open descendent invariants of a point,
# see Section 1.7 of https://arxiv.org/abs/1409.2191.
from jrrFormula import *

# implements fixed point diagram generation.
from fixedPoints import *

def wfp_nou(d,l,k,e) :
    """
this is short for "weighted fixed point contribution no :math:`u`".

It's a wrapper for :py:func:`weighted_fpContrib` that removes the formal generator u and returns a rational number. This is the function we actually use.

The inputs are as follows:

* :py:obj:`d` is the relative degree :math:`\\in H_2(X,L)`

* :py:obj:`l` is the number of interior markings.

* :py:obj:`k` is the number of boundary markings and incoming edges.

* :py:obj:`e` is the number of outgoing edges.

If :math:`d` is even we must have :math:`e = 0` and if :math:`d` is odd, we must have :math:`k = 0`.

The output is a rational number.
"""
    return weighted_fpContrib(d,l,k,e)*u**(old_div((3*d-1-k-2*l),2))

def weighted_fpContrib(d,l,k,e = 0) :
    """
This is the main function which computes

.. math::

    \\sum_{A_{\\tilde F}} A_{\\tilde F} = \\sum \\xi_{F(\\tilde F)} \\int_{\\tilde F} \\frac{\\omega|_{\\tilde F}}{e^{T^2}(N_{\\tilde F})}`.

There's a detailed explanation in the introduction to the module documentation.
"""

    # generate all the fixed point diagrams
    # The zero here is an optional
    # argument which is never used in this program (allows for
    # "dummy" interior markings which don't carry
    # a constraint).
    fps = fixedPoints(d, l, k, 0)

    N = len(fps)

    # a buffer to hold the conributions of the fixed points.
    conts = zeros([N], sp.Number)

    # this set records the indices of the T^2 fixed point contributions,
    # which would become singular if we set eps = 0,
    # which corresponds to pulling back along BS^1 -> BT^2 (see
    # detailed discussion in the intro to the module)
    sing = set()

    # go over all of the fixed point components...
    for i in range(N) :

        # compute the integral of \omega|_F / e(N_F)
        # for the i'th fixed point. That is,
        # fact is the integral without the weight prefactor
        # \xi_F.
        fact = fpContrib(fps[i])

        # the degree of the irreducible disk component
        discDeg = fps[i][0]

        # where the disk component maps to.
        rootMu = fps[i][-1][0][0][3]

        # the integral of the propagator
        # picks up a sign according to the boundary orientation,
        # which depends on mu.
        sign = (-1)**(old_div(rootMu,2) + 1)

        # if there are outgoing edges...
        if e != 0 :
            # ... the contribution is the \xi_F x fact.
            conts[i] = (sign*rat(discDeg,2))**e*fact
        else :
            # otherwise, the contribution is just fact.
            conts[i] = fact

        # this is an ugly hack, it tries
        # to set eps = 0. It tests whether we get a singularity
        # by looking for "oo" in the string representation of this
        # evaluation. I couldn't find a better way to do this
        # (maybe things have improved since I last looked).
        # if we get a singularity...
        tmp = str(conts[i].subs(eps,0))
        if tmp.find('oo') != -1 or tmp.find('nan') != -1 :
            # ... we record the singularity's index in the list.
            sing.add(i)
        else :
            # otherwise, we replace conts[i] by the eps = 0 subtitution.
            conts[i] = conts[i].subs(eps,0)

    # we want to group the singular T^2 fixed point contributions
    # by the S^1 fixed point component they belong to.
    # as discussed in the intro, the sum over each group should (and does!)
    # become regular at eps = 0.
    #
    # we do NOT want to sum over all of the fixed point components
    # to avoid getting enormous (memory-clogging enormous,
    # if I remember correctly) denominators coming from
    # summing too many unrelated rational functions.
    p2c = defaultdict(lambda : [])
    for i in sing :
        # pinch is implemented in pincher.py
        # and it returns a unique representative to each
        # S^1 equivalence class of fixed point diagrams.
        # so we just add the index i to the dictionary
        # at this key specifying the S^1 fixed point component.
        p2c[pinch(fps[i])].append(i)

    # Each cell of this array corresponds to one of the groups of
    # singular contributions above...
    sconts = zeros([len(p2c)], sp.Number)
    keys = list(p2c.keys())

    for j in range(len(keys)) :
        # ... and it holds the sum of the singular contributions in the group
        # evaluated at eps = 0 (this evaluation should be okay after we sum).
        sconts[j] = sumRats([conts[i] for i in p2c[keys[j]]], vrs = [eps, u]).subs(eps,0)

    # total sum, of all the singular contributions.
    s1 = intsum(sconts)

    # in case we have an empty sum, we replace False by zero.
    if s1 == False :
        s1 = 0

    # sum of the singular contributions and the contributions that
    # were regular to begin with.
    res = s1 + intsum([conts[i] for i in range(len(conts)) if i not in sing])
    return res

def fpContrib(fp) :
    """
Computes

.. math::

    \\int_{\\tilde F} \\frac{\\omega|_{\\tilde F}}{e^{T^2}(N_{\\tilde F})}

(without the outgoing edge weights).
"""
    # get the specs of the moduli spaces of the contracted components.
    # as well as a description of the (factorized) Euler form.
    openModuli,closedModulis = fpModuliAndEuler(fp)
    res = 1

    # if the disk component is contracted (even degree)...
    if openModuli :
        # ... extract its specs ...
        coeff,dumb0,totalK,omegas = openModuli

        # and compute the relevant open descendent integral.
        res *= coeff * openIntegral(dumb0,totalK,omegas)

    # go over the moduli spaces of the other contracted components...
    for coeff,l,omegas in closedModulis :
        # ... and multiply in the corresponding descendent integral.
        res *= coeff * closedIntegral(l,omegas)

    # return res / Aut(fp) (divide by the size of the automorphism group).
    return fpCombFactor(fp) * res

def closedIntegral(l, omegas) :
    """
Let :math:`e` denote :py:obj:`len(omegas)`.
This function computes the integral of

.. math::

    \\prod_{1 \leq i \leq e} \\frac{1}{1-\\psi_i/\\omega_i} := \\sum_{j_1,...,j_e} \\prod_i \\left(\\frac{\\psi_i}{\\omega_i}\\right)^{j_i}

over :math:`\\overline{\\mathcal{M}}_{0,l + e}`.

If :math:`l + e < 3` this is just the unit class. Otherwise, we apply the string equation (see Exercise 25.2.8 in `MS book <http://www.claymath.org/library/monographs/cmim-1.pdf>`_) and simplify to find this is just

.. math::

    (\\sum_i \\frac{1}{\\omega_i})^{l + e -3}

"""
    dim = max(l + len(omegas)-3,0)
    return sum(o**(-1) for o in omegas)**dim  #* (u**(2*l)

from collections import Counter

def openIntegral(dummies, k, omegas) :
    """
computes the integral over open moduli of discs of

.. math::

    \\prod_{1 \leq i \leq e} \\frac{1}{1-\\psi_i/\\omega_i} := \\sum_{j_1,...,j_e} \\prod_i \\left(\\frac{\\psi_i}{\\omega_i}\\right)^{j_i}

by invoking the Pandharipande-Solomon-Tessler formula (see intro to module for more info).
"""
    nnodes = len(omegas)
    # this is the real dimension of the moduli space.
    l = nnodes + dummies
    dim = 2*l + k - 3
    if dim % 2 == 1 :
        raise ValueError("dimension of moduli space is odd, something's wrong...")

    # complex dimension
    cdim = old_div(dim,2)

    openIntegral = 0

    # psiPows describes how many times each of the nodes' psi appears
    # in the integrand. We know the total power must equal the complex dimension
    # of the moduli space.
    for psiPows in tuplePartitions(cdim, nnodes) :
        coefficient = myprod(o**(-p) for (o,p) in zip(omegas, psiPows))

        # if (p,n) appears in tauPows, it should be interpreted as
        # "n of the psi_i's appear to the power p". which means
        # tau_p^n should be added to the invariant we compute.
        tauPows       = list(Counter(psiPows).items())
        openIntegral += coefficient * jrrFormula(addTau0s(tauPows,dummies),k)

    return openIntegral

# For most of the fixed point formula, it is convenient to think of
# CP^2 more symmetrically as being equipped with a (degenerate) T^3 action,
# each factor acting by complex multiplication on one of the three coordinates.
#
# The T^2 action factors through this T^3 action via a homomorphism
# T^2 -> T^3. Pullback in cohomology along BT^2 -> BT^3 is given by
# the assignment described by the list below.
alpha = [u+eps,0,-u+eps]

# Poincare dual to the real fixed point :math:`p_0` inside
# :math:`\\mathbb{R}P^2`. dumb is always zero in this version.
def PDp0(imps,dumb,mu) :
    return ((alpha[mu] - alpha[0])*(alpha[mu]-alpha[2]))**imps*u**(2*dumb)

# Poincare dual to the complex fixed point :math:`p_+` or :math:`p_-`
# inside :math:`\\mathbb{C}P^2`. dumb is always zero in this version.
def PDppm(imps,dumb,mu) :
    return (rat(1,2) * (alpha[mu] - alpha[0]) * (alpha[mu] - alpha[1]) +
            rat(1,2) * (alpha[mu] - alpha[2]) * (alpha[mu] - alpha[1]))**imps *\
            u**(2*dumb)

def fpModuliAndEuler(fp,integrand_at_F = PDppm) :
    """
Given a fixed point diagram :py:obj:`fp`, this function computes

* the parameter spaces for the vertices in the fixed point diagram (so the product of these spaces is diffeomorphic to the fixed point component specified by :py:obj:`fp`). [#this-is-a-lie]_

* the integrand for this fixed point locus, as a product of factors, one for each of the parameter spaces above.

The name of the function, :py:obj:`fpModuliAndEuler`, is thus a bit of a misnomer, perhaps :py:obj:`fpModuliAndIntegrand` would've been better.

The output is of the form :py:obj:`(openModuli, closedModulis)`.

**The open Moduli space and integrand factor**

:py:obj:`openModuli` is :py:obj:`False` if the disk degree is positive (and odd), so there's no contracted component there and the parameter space is just a point. The relevant "square root of edge" term corresponding to deformations of the map is then appended to the coefficient of the root vertex in the fixed point diagram (in general, each vertex's :py:obj:`coeff` is responsible for the unique edge going up towards the root).

Otherwise, if the total degree is even :py:obj:`openModuli = (coeff,dumb0,totalK,omegas)` is a tuple where :py:obj:`coeff`, and each element of the list :py:obj:`omegas`, is a rational function in :py:obj:`u,eps`. :py:obj:`dumb0` can always be assumed to be zero in this version (it's meant to accomodate "dummy" variables, not used here), and :py:obj:`totalK` is an integer. This tuple specifies

* the parameter space :math:`\overline{\mathcal{M}}_{0,k,l}`

here :math:`k` is :py:obj:`totalK` (this includes both the boundary markings and the incoming edges), :math:`l` should be :py:obj:`len(omegas)` assuming [#this-is-a-lie]_ there are no interior markings on a disk mapping to :math:`p_0`, only markings associated with the complex nodes, which are in bijection with the elements of :py:obj:`omegas`).

* a factor :math:`\\mbox{coeff}\\,\\prod_i \\left(\\psi_i - \\mbox{omegas}[i]\\right)^{-1}` of the integrand.

Here :py:obj:`omegas` is a tuple of "omega classes", see Definition 27.3.1 in the `MS book <http://www.claymath.org/library/monographs/cmim-1.pdf>`_. For each flag :math:`F` in a fixed point diagram, :math:`\omega_F` is the first chern class of the tangent space to the positive energy component corresponding to the edge of the flag, evaluated at the special point where it is attached to the flag's vertex (which in our case, is the contracted component of the disk)

**The closed moduli spaces and integrand factors**

Each element of :py:obj:`closedModulis` is a tuple :py:obj:`(coeff,l,omegas)`  where :py:obj:`coeff`, and each element of the list :py:obj:`omegas`, is a rational function in :py:obj:`u,eps`. :py:obj:`l` is an integer specifying the number of interior markings on the vertex. Let :math:`n = l + e` where :math:`e` is the number of positive energy components (including the disk component, if we're at the root) incident to the vertex.

In case :math:`n \geq 3`, :py:obj:`omegas` contains a list of length :math:`e` of the flags' "omega classes" as discussed above, so :math:`omegas` specifies the factor involving :math:`\psi` classes on the parameter space :math:`\\overline{\\mathcal{M}}_{0,n}`. If :math:`n < 3` then this factor does not appear; :py:obj:`omegas` is empty; and the parameter space is just a point. In either case, :py:obj:`coeff` holds the other contributions to the integrand (except for the :math:`\psi` classes, all classes live in the cohomology ring, pulled back from a point).

See the intro to the module for more on the computation of the inverse Euler.

.. [#this-is-a-lie] Actually, this is a somewhat idealized description of what the function does: in case there are interior markings on a contracted disk component you wouldn't know about them, so you cannot really reconstruct the true parameter space in this case.

   .. code-block:: py

    >>> fps = fixedPoints(6,7,3,0)
    >>> printFP(fps[2500])
                            0 )...
    ops((3, 7, 0, 1))[182]     <p0>||||            degs: (1, 1)
    ops((1, 0, 0, 2))[0]          <p->             degs: (1,)
    ops((0, 0, 0, 1))[0]             <p0>          degs: ()
    ops((0, 3, 0, 0))[0]          <p+>|||          degs: ()

    >>> openModuli,closedModulis = fpModuliAndEuler(fps[2500])
    >>> openModuli
    (0, 0, 3, (-eps + u, -eps - u))
    >>> coeff,dumb0,totalK,omegas = openModuli

   Of course this does not matter, since our assumption on the support of the form carried by the interior markings force these contributions to vanish (indeed, we see that coeff is zero). Still, it's a bit ugly and at some point you may want to rewrite this. As we mentioned elsewhere, throwing out diagrams like fps[2500] earlier will probably also speed things up.

"""
    # parse fp:
    # extract the degree of the disk irreducible component,
    # the number of boundary markings (including incoming edges)
    # and the T^2 fixed point tree.
    (discDeg,totalK,tree) = fp

    # list of the closed moduli
    closedModulis = []
    for v in vertices(tree) :
        # we deal with the root separately
        if v == () :
            continue
        (subtreeData, counter, mu, imps0, dumb0,degs) = getVertex(v, tree)

        parentDat = getVertex(v[:-1],tree)

        # get incoming degree - is listed in the parents degs list...
        inDeg = parentDat[-1][v[-1]]

        # get parent's mu
        parentMu = parentDat[2]

        # number of imp's on the vertex
        l = imps0 + dumb0

        # this is the specialization of the extended form
        # defining the invariants to the fixed point.
        coeff = integrand_at_F(imps0,dumb0,mu)
        omegas = []

        # nmarks includes all half edges coming into the moduli,
        # i.e. imp's and nodes with other edges (including one for the parent)
        nMarks = len(degs) + l + 1


        ## FLAGS TERM ##
        # see the intro to the module for the formula we refer to
        # by "flags term", "vertices term" etc.
        #
        # if there's a contracted moduli here, we add the first flags term
        if nMarks >= 3 :
            omega = old_div((alpha[mu] - alpha[parentMu]), inDeg)

            # first deal with parent
            coeff /= omega
            omegas.append(omega)

            for i in range(len(degs)) :
                outDeg = degs[i]
                kidMu  = getVertex(v + (i,), tree)[2]
                omega = old_div((alpha[mu] - alpha[kidMu]),outDeg)
                coeff /= omega
                omegas.append(omega)

        # rest of flags term minus 1 for the vertices term, so we get len(degs)
        # and not len(degs) +1.
        coeff *= myprod([(alpha[mu] - alpha[nu])**(len(degs)) \
                         for nu in [0,1,2] if nu != mu])

        # if valency is 2, and no imps
        if l == 0 and len(degs) == 1 :
            outDeg = degs[0]
            kidMu  = getVertex(v + (0,), tree)[2]
            omega = old_div((alpha[mu] - alpha[kidMu]),outDeg)

            tmp = (old_div((alpha[mu] - alpha[parentMu]),inDeg) + \
                   old_div((alpha[mu] - alpha[kidMu]),outDeg))
            coeff /= tmp

        # if valency is 1 and there are no imps
        if l == 0 and len(degs) == 0 :
            coeff *= old_div((alpha[mu] - alpha[parentMu]),inDeg)

        ## EDGES TERM ##
        coeff *= old_div(rat((-1)**(inDeg)*(inDeg)**(2*inDeg),factorial(inDeg)**2), (alpha[parentMu]-alpha[mu])**(2*inDeg))

        # an akward way of getting the k in the formula...
        # here i,j,k refer to the formula in the MS book,
        # see the intro to the module.
        k = [k for k in [0,1,2] if ((k != mu) and (k != parentMu))][0]
        for a in range(inDeg + 1) :
            b = inDeg - a
            tmp = (rat(a,inDeg) * alpha[mu] + \
                   rat(b,inDeg) * alpha[parentMu] - \
                   alpha[k])
            coeff /= tmp

        closedModulis.append((coeff,l,tuple(omegas)))


    # ITERSECTION TERM AND FIRST VERTEX CONTRIBUTION

    # if disc is contracted...
    if discDeg == 0 :
        if totalK == 0 :
            raise ValueError("when there's a ghost disc need positive number of bmps.")

        v = ()
        (subtreeData, counter, mu, imps0, dumb0,degs) = getVertex(v, tree)
        omegas = []

        if mu != 1 :
            raise ValueError("expected mu to be 1 on ghost disc!, mu = ".format(mu))


        # ugly code - the part of the integrand that comes from
        # the bmp's is written outside the integrand at F...
        coeff = u**(totalK) * integrand_at_F(imps0, dumb0, mu)

        # we must have a bmp so there's a contracted component there.
        for i in range(len(degs)) :
            outDeg = degs[i]
            kidMu  = getVertex(v + (i,), tree)[2]
            omega = old_div((alpha[mu] - alpha[kidMu]),outDeg)
            coeff /= omega
            omegas.append(omega)

        # flags contribution (both val = 2 and val = 1 are irrelevant
        # since there are marked points)

        # (mu == 1 in this case, but okay...)
        coeff *= myprod([(alpha[mu] - alpha[nu])**(len(degs)) \
                        for nu in [0,1,2] if nu != mu])

        # vertex contribution is just 1/u, since the vertex can only move
        # along the lagrangian RP^2. Here's where a subtle sign may come in...
        coeff /= u

        tD,totalL,totalT,rm = subtreeData

        openModuli = (coeff,dumb0,totalK,tuple(omegas))

    # if disc is not contracted
    if discDeg > 0 :
        openModuli = False
        # in this case the first vertex is similar to any other, except for the
        # edges term (which also contains the elusive sign).
        v = ()
        (subtreeData, counter, mu, imps0, dumb0,degs) = getVertex(v, tree)

        # parentMu, for the sake of omega computation, is the conjugate of mu
        parentMu = 2-mu

        # incoming deg in this case is the disc degree.
        inDeg = discDeg

        # number of imp's on the vertex
        l = imps0 + dumb0

        # pullback term
        coeff = (rat(1,2) * (alpha[mu] - alpha[0]) * (alpha[mu] - alpha[1]) +
                 rat(1,2) * (alpha[mu] - alpha[2]) * (alpha[mu] - alpha[1]))**imps0 *\
                u**(2*dumb0)

        omegas = []

        # nmarks includes all half edges coming into the moduli,
        # i.e. imp's and nodes with other edges (including one for the parent)
        nMarks = len(degs) + l + 1

        ## FLAGS TERM ##
        # if there's a contracted moduli here, we add the first flags term
        if nMarks >= 3 :
            omega = old_div((alpha[mu] - alpha[parentMu]), inDeg)

            # first deal with parent
            coeff /= omega
            omegas.append(omega)

            for i in range(len(degs)) :
                outDeg = degs[i]
                kidMu  = getVertex(v + (i,), tree)[2]
                omega = old_div((alpha[mu] - alpha[kidMu]),outDeg)
                coeff /= omega
                omegas.append(omega)

        # rest of flags term minus 1 for the vertices term, so we get len(degs)
        # and not len(degs) +1.
        coeff *= myprod([(alpha[mu] - alpha[nu])**(len(degs)) \
                         for nu in [0,1,2] if nu != mu])

        # if valency is 2, and no imps
        if l == 0 and len(degs) == 1 :
            outDeg = degs[0]
            kidMu  = getVertex(v + (0,), tree)[2]
            omega = old_div((alpha[mu] - alpha[kidMu]),outDeg)
            coeff /= (old_div((alpha[mu] - alpha[parentMu]),inDeg) + \
                      old_div((alpha[mu] - alpha[kidMu]),outDeg))

        # if valency is 1 and there are no imps
        if l == 0 and len(degs) == 0 :
            coeff *= old_div((alpha[mu] - alpha[parentMu]),inDeg)

        ## EDGES TERM ##
        coeff *= old_div(rat((inDeg)**(inDeg),factorial(inDeg)), (alpha[parentMu]-alpha[mu])**(inDeg))

        # an akward way of getting the k in the formula...
        k = [k for k in [0,1,2] if ((k != mu) and (k != parentMu))][0]
        for a in range(old_div((inDeg + 1),2)) :
            b = inDeg - a
            coeff /= (rat(a,inDeg) * alpha[mu] + rat(b,inDeg) * alpha[parentMu] - alpha[k])
        tD,totalL,totalT,rm = subtreeData
        closedModulis.append((coeff,l,omegas))
    return openModuli, closedModulis


def sumRats(rats, vrs = None) :
    """
Rewrites a sum of rational functions :math:`a_i/b_i` as a single rational function by taking the naive common denominator. E.g.

.. math::

    a_1/b_1 + a_2/b_2 + a_3/b_3 = \\frac{a_1 b_2 b_3 + b_1 a_2 b_3 + b_1 b_2 a_3}{b_1 b_2 b_3}.

If vrs is specified, it also factors out from the numerator and the denominator
the highest power possible of each of the variables vrs. It returns a quotient
of polynomials in this case.
"""
    n = len(rats)
    if n == 0 :
        return sp.Number(0)

    nds = [r.as_numer_denom() for r in rats]
    newnumer = sum([nds[i][0] * myprod([nds[j][1] for j in range(n) if j != i])\
                 for i in range(n)])
    newdenom = myprod([x[1] for x in nds])
    if vrs :
        npoly = newnumer.as_poly(vrs)
        dpoly = newdenom.as_poly(vrs)
        cexps = list(map(min, list(zip(list(map(min,list(zip(*list(npoly.as_dict().keys()))))),\
                              list(map(min,list(zip(*list(dpoly.as_dict().keys())))))))))
        common = myprod([vrs[i] ** cexps[i] for i in range(len(cexps))])

        return old_div((old_div(npoly,common)).as_poly(*vrs), (old_div(dpoly,common)).as_poly(*vrs))
    return old_div(newnumer,newdenom)

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
