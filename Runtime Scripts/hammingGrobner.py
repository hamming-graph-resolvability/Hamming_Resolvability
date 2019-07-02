###############################################################################################
# Main functions for checking resolvability of a set R on Hamming Graph. 
# Author: Lucas Laird, Richard Carter Tillquist
# Last Updated: June 2019
###############################################################################################


###############################################################################################
#                                   Modified SymPy Groebner Class                             #
###############################################################################################

import multiprocessing as mp

from sympy import *
from itertools import combinations

###from __future__ import print_function, division

from sympy.core import (
    S, Basic, Expr, I, Integer, Add, Mul, Dummy, Tuple
)

from sympy.core.mul import _keep_coeff
from sympy.core.symbol import Symbol
from sympy.core.basic import preorder_traversal
from sympy.core.relational import Relational
from sympy.core.sympify import sympify
from sympy.core.decorators import _sympifyit
from sympy.core.function import Derivative

from sympy.logic.boolalg import BooleanAtom

from sympy.polys.polyclasses import DMP

from sympy.polys.polyutils import (
    basic_from_dict,
    _sort_gens,
    _unify_gens,
    _dict_reorder,
    _dict_from_expr,
    _parallel_dict_from_expr,
)

from sympy.polys.rationaltools import together
from sympy.polys.rootisolation import dup_isolate_real_roots_list
from sympy.polys.groebnertools import groebner as _groebner
from sympy.polys.fglmtools import matrix_fglm
from sympy.polys.monomials import Monomial
from sympy.polys.orderings import monomial_key

from sympy.polys.polyerrors import (
    OperationNotSupported, DomainError,
    CoercionFailed, UnificationFailed,
    GeneratorsNeeded, PolynomialError,
    MultivariatePolynomialError,
    ExactQuotientFailed,
    PolificationFailed,
    ComputationFailed,
    GeneratorsError,
)

from sympy.utilities import group, sift, public, filldedent

import sympy.polys
import mpmath
from mpmath.libmp.libhyper import NoConvergence

from sympy.polys.domains import FF, QQ, ZZ
from sympy.polys.constructor import construct_domain

from sympy.polys import polyoptions as options

from sympy.core.compatibility import iterable, range, ordered
from sympy.polys.polytools import *


from sympy.polys.monomials import monomial_mul, monomial_lcm, monomial_divides, term_div
from sympy.polys.orderings import lex
from sympy.polys.polyerrors import DomainError
from sympy.polys.polyconfig import query
from sympy.core.symbol import Dummy
from sympy.core.compatibility import range
from sympy.polys.groebnertools import *

@public
def groebnerbasis(iterative, F, n, *gens, **args):
    return GroebnerBasis(iterative, F, n, *gens, **args)


@public
class GroebnerBasis(Basic):
    """Represents a reduced Groebner basis. """

    def __new__(cls, iterative, F, n, *gens, **args):
        """Compute a reduced Groebner basis for a system of polynomials. """
        options.allowed_flags(args, ['polys', 'method'])
        try:
            polys, opt = parallel_poly_from_expr(F, *gens, **args)
        except PolificationFailed as exc:
            raise ComputationFailed('groebner', len(F), exc)

        from sympy.polys.rings import PolyRing
        ring = PolyRing(opt.gens, opt.domain, opt.order)

        polys = [ring.from_dict(poly.rep.to_dict()) for poly in polys if poly]
        if iterative:
            G = iter_groebner(polys,n,ring, method = opt.method)#Assumes last element is the "new" polynomial
        else:
            G = _groebner(polys, ring, method=opt.method)
        G = [Poly._from_dict(g, opt) for g in G]

        return cls._new(G, opt)

    @classmethod
    def _new(cls, basis, options):
        obj = Basic.__new__(cls)

        obj._basis = tuple(basis)
        obj._options = options

        return obj

    @property
    def args(self):
        return (Tuple(*self._basis), Tuple(*self._options.gens))

    @property
    def exprs(self):
        return [poly.as_expr() for poly in self._basis]

    @property
    def polys(self):
        return list(self._basis)

    @property
    def gens(self):
        return self._options.gens

    @property
    def domain(self):
        return self._options.domain

    @property
    def order(self):
        return self._options.order

    def __len__(self):
        return len(self._basis)

    def __iter__(self):
        if self._options.polys:
            return iter(self.polys)
        else:
            return iter(self.exprs)

    def __getitem__(self, item):
        if self._options.polys:
            basis = self.polys
        else:
            basis = self.exprs

        return basis[item]

    def __hash__(self):
        return hash((self._basis, tuple(self._options.items())))

    def __eq__(self, other):
        if isinstance(other, self.__class__):
            return self._basis == other._basis and self._options == other._options
        elif iterable(other):
            return self.polys == list(other) or self.exprs == list(other)
        else:
            return False

    def __ne__(self, other):
        return not self == other

    @property
    def is_zero_dimensional(self):
        """
        Checks if the ideal generated by a Groebner basis is zero-dimensional.

        The algorithm checks if the set of monomials not divisible by the
        leading monomial of any element of ``F`` is bounded.

        References
        ==========

        David A. Cox, John B. Little, Donal O'Shea. Ideals, Varieties and
        Algorithms, 3rd edition, p. 230

        """
        def single_var(monomial):
            return sum(map(bool, monomial)) == 1

        exponents = Monomial([0]*len(self.gens))
        order = self._options.order

        for poly in self.polys:
            monomial = poly.LM(order=order)

            if single_var(monomial):
                exponents *= monomial

        # If any element of the exponents vector is zero, then there's
        # a variable for which there's no degree bound and the ideal
        # generated by this Groebner basis isn't zero-dimensional.
        return all(exponents)

    def fglm(self, order):
        """
        Convert a Groebner basis from one ordering to another.

        The FGLM algorithm converts reduced Groebner bases of zero-dimensional
        ideals from one ordering to another. This method is often used when it
        is infeasible to compute a Groebner basis with respect to a particular
        ordering directly.

        Examples
        ========

        >>> from sympy.abc import x, y
        >>> from sympy import groebner

        >>> F = [x**2 - 3*y - x + 1, y**2 - 2*x + y - 1]
        >>> G = groebner(F, x, y, order='grlex')

        >>> list(G.fglm('lex'))
        [2*x - y**2 - y + 1, y**4 + 2*y**3 - 3*y**2 - 16*y + 7]
        >>> list(groebner(F, x, y, order='lex'))
        [2*x - y**2 - y + 1, y**4 + 2*y**3 - 3*y**2 - 16*y + 7]

        References
        ==========

        J.C. Faugere, P. Gianni, D. Lazard, T. Mora (1994). Efficient
        Computation of Zero-dimensional Groebner Bases by Change of
        Ordering

        """
        opt = self._options

        src_order = opt.order
        dst_order = monomial_key(order)

        if src_order == dst_order:
            return self

        if not self.is_zero_dimensional:
            raise NotImplementedError("can't convert Groebner bases of ideals with positive dimension")

        polys = list(self._basis)
        domain = opt.domain

        opt = opt.clone(dict(
            domain=domain.get_field(),
            order=dst_order,
        ))

        from sympy.polys.rings import xring
        _ring, _ = xring(opt.gens, opt.domain, src_order)

        for i, poly in enumerate(polys):
            poly = poly.set_domain(opt.domain).rep.to_dict()
            polys[i] = _ring.from_dict(poly)

        G = matrix_fglm(polys, _ring, dst_order)
        G = [Poly._from_dict(dict(g), opt) for g in G]

        if not domain.is_Field:
            G = [g.clear_denoms(convert=True)[1] for g in G]
            opt.domain = domain

        return self._new(G, opt)

    def reduce(self, expr, auto=True):
        """
        Reduces a polynomial modulo a Groebner basis.

        Given a polynomial ``f`` and a set of polynomials ``G = (g_1, ..., g_n)``,
        computes a set of quotients ``q = (q_1, ..., q_n)`` and the remainder ``r``
        such that ``f = q_1*f_1 + ... + q_n*f_n + r``, where ``r`` vanishes or ``r``
        is a completely reduced polynomial with respect to ``G``.

        Examples
        ========

        >>> from sympy import groebner, expand
        >>> from sympy.abc import x, y

        >>> f = 2*x**4 - x**2 + y**3 + y**2
        >>> G = groebner([x**3 - x, y**3 - y])

        >>> G.reduce(f)
        ([2*x, 1], x**2 + y**2 + y)
        >>> Q, r = _

        >>> expand(sum(q*g for q, g in zip(Q, G)) + r)
        2*x**4 - x**2 + y**3 + y**2
        >>> _ == f
        True

        """
        poly = Poly._from_expr(expr, self._options)
        polys = [poly] + list(self._basis)

        opt = self._options
        domain = opt.domain

        retract = False

        if auto and domain.is_Ring and not domain.is_Field:
            opt = opt.clone(dict(domain=domain.get_field()))
            retract = True

        from sympy.polys.rings import xring
        _ring, _ = xring(opt.gens, opt.domain, opt.order)

        for i, poly in enumerate(polys):
            poly = poly.set_domain(opt.domain).rep.to_dict()
            polys[i] = _ring.from_dict(poly)

        Q, r = polys[0].div(polys[1:])

        Q = [Poly._from_dict(dict(q), opt) for q in Q]
        r = Poly._from_dict(dict(r), opt)

        if retract:
            try:
                _Q, _r = [q.to_ring() for q in Q], r.to_ring()
            except CoercionFailed:
                pass
            else:
                Q, r = _Q, _r

        if not opt.polys:
            return [q.as_expr() for q in Q], r.as_expr()
        else:
            return Q, r


    def contains(self, poly):
        """
        Check if ``poly`` belongs the ideal generated by ``self``.

        Examples
        ========

        >>> from sympy import groebner
        >>> from sympy.abc import x, y

        >>> f = 2*x**3 + y**3 + 3*y
        >>> G = groebner([x**2 + y**2 - 1, x*y - 2])

        >>> G.contains(f)
        True
        >>> G.contains(f + 1)
        False

        """
        return self.reduce(poly)[1] == 0

def iter_groebner(seq,n, ring, method=None):
    """
    Computes Groebner basis for a set of polynomials in `K[X]`.

    Wrapper around the (default) improved Buchberger and the other algorithms
    for computing Groebner bases. The choice of algorithm can be changed via
    ``method`` argument or :func:`setup` from :mod:`sympy.polys.polyconfig`,
    where ``method`` can be either ``buchberger`` or ``f5b``.

    """

    domain, orig = ring.domain, None
    if not domain.is_Field or not domain.has_assoc_Field:
        try:
            orig, ring = ring, ring.clone(domain=domain.get_field())
        except DomainError:
            raise DomainError("can't compute a Groebner basis over %s" % domain)
        else:
            seq = [ s.set_ring(ring) for s in seq ]

    G = _incr_buch(seq,n, ring)

    if orig is not None:
        G = [ g.clear_denoms()[1].set_ring(orig) for g in G ]

    return G

def _incr_buch(f,n, ring):
    """
    Incremental Computation of a reduced Grobner basis, given that f[:-1] is a reduced Grobner basis
    and f[-1] is a new polynomial to add into the basis.
    New polynomial f[-1] is assumed to be reduced by the grobner basis f[:-1] already.

    """
    order = ring.order
    domain = ring.domain

    monomial_mul = ring.monomial_mul
    monomial_div = ring.monomial_div
    monomial_lcm = ring.monomial_lcm

    def select(P):
        # normal selection strategy
        # select the pair with minimum LCM(LM(f), LM(g))
        pr = min(P, key=lambda pair: order(monomial_lcm(f[pair[0]].LM, f[pair[1]].LM)))
        return pr

    def normal(g, J):
        h = g.rem([ f[j] for j in J ])

        if not h:
            return None
        else:
            h = h.monic()

            if not h in I:
                I[h] = len(f)
                f.append(h)

            return h.LM, I[h]

    def update(G, B, ih):
        # update G using the set of critical pairs B and h
        # [BW] page 230
        h = f[ih]
        mh = h.LM

        # filter new pairs (h, g), g in G
        C = G.copy()
        D = set()

        while C:
            # select a pair (h, g) by popping an element from C
            ig = C.pop()
            g = f[ig]
            mg = g.LM
            LCMhg = monomial_lcm(mh, mg)

            def lcm_divides(ip):
                # LCM(LM(h), LM(p)) divides LCM(LM(h), LM(g))
                m = monomial_lcm(mh, f[ip].LM)
                return monomial_div(LCMhg, m)

            # HT(h) and HT(g) disjoint: mh*mg == LCMhg
            if monomial_mul(mh, mg) == LCMhg or (
                not any(lcm_divides(ipx) for ipx in C) and
                    not any(lcm_divides(pr[1]) for pr in D)):
                D.add((ih, ig))

        E = set()

        while D:
            # select h, g from D (h the same as above)
            ih, ig = D.pop()
            mg = f[ig].LM
            LCMhg = monomial_lcm(mh, mg)

            if not monomial_mul(mh, mg) == LCMhg:
                E.add((ih, ig))

        # filter old pairs
        B_new = set()

        while B:
            # select g1, g2 from B (-> CP)
            ig1, ig2 = B.pop()
            mg1 = f[ig1].LM
            mg2 = f[ig2].LM
            LCM12 = monomial_lcm(mg1, mg2)

            # if HT(h) does not divide lcm(HT(g1), HT(g2))
            if not monomial_div(LCM12, mh) or \
                monomial_lcm(mg1, mh) == LCM12 or \
                    monomial_lcm(mg2, mh) == LCM12:
                B_new.add((ig1, ig2))

        B_new |= E

        # filter polynomials
        G_new = set()

        while G:
            ig = G.pop()
            mg = f[ig].LM

            if not monomial_div(mg, mh):
                G_new.add(ig)

        G_new.add(ih)

        return G_new, B_new
        # end of update ################################

    if not f:
        return []
    f1 = [func for func in f[:-n]]
    for p in f[-n:]:
        r = p.rem(f1)
        if r != 0:
            f1.append(r)
    f = f1
    I = {} # ip = I[p]; p = f[ip]
    F = set()
    G = set()         # set of indices of intermediate would-be Groebner basis
    CP = set()        # set of pairs of indices of critical pairs

    for i, h in enumerate(f):
        I[h] = i #Setup polynomial-index dictionary
        if i >= len(f)-n:
            F.add(i)
        else:
            G.add(i)

    #####################################
    # algorithm GROEBNERNEWS2 in [BW] page 232
    while F:
        # select p with minimum monomial according to the monomial ordering
        h = min([f[x] for x in F], key=lambda f: order(f.LM))
        ih = I[h]
        F.remove(ih)
        G, CP = update(G, CP, ih)

    # count the number of critical pairs which reduce to zero
    reductions_to_zero = 0

    while CP:
        ig1, ig2 = select(CP)
        CP.remove((ig1, ig2))
        
        h = spoly(f[ig1], f[ig2], ring)
        # ordering divisors is on average more efficient [Cox] page 111
        G1 = sorted(G, key=lambda g: order(f[g].LM))
        ht = normal(h, G1)

        if ht:
            G, CP = update(G, CP, ht[1])
        else:
            reductions_to_zero += 1

    ######################################
    # now G is a Groebner basis; reduce it
    Gr = set()

    for ig in G:
        ht = normal(f[ig], G - set([ig]))

        if ht:
            Gr.add(ht[1])

    Gr = [f[ig] for ig in Gr]

    # order according to the monomial ordering
    Gr = sorted(Gr, key=lambda f: order(f.LM), reverse=True)

    return Gr

import numpy as np
from math import floor
import time
from itertools import product
from itertools import combinations

###############################################################################################
#                                   GENERAL FUNCTIONS                                         #
###############################################################################################

def make_linearEqns(A,z):
    #Takes matrix A and variable list z and converts it into a list of linear equations
    zmat = Matrix(z)
    A = Matrix(A)
    A = A.rref()[0]
    lin_fcns = A*zmat
    lin_fcns = [f for f in lin_fcns if f != 0]
    lin_fcns = lin_fcns[::-1]
    return(list(lin_fcns))

def OneHot(R,k,a,alphabet = None):
    #Converts list of strings to list of one-hot encodings
    if alphabet == None:
        alphabet = [str(i) for i in range(a)]
        #alphabet = ''.join(alphabet)
    encodings = []
    for r in R:
        encoding = np.zeros((a,k))
        for i in range(k):
            for j,letter in enumerate(alphabet):
                if r[i] == letter:
                    encoding[j,i] = 1
        encodings.append(encoding)
    return(encodings)

###############################################################################################
#                                 ALL HAMMING GRAPH FUNCTIONS                                 #
###############################################################################################


def create_groebner(k,a):
    #Constructs known groebner basis based on pattern for P and creates fs and variable vector.
    variable_string = ''
    num_digits = len(str(k*a))
    for i in range(0,k*a):
        count = str(i+1).zfill(num_digits)
        variable = ',z'+count
        variable_string = variable_string +variable
    variable_string = variable_string[1:]
    z = var(variable_string)
    z = tuple(z)
    G0 = []
    for i in range(1,a):
        G0.append(z[i]**3-z[i])
    pairs = combinations(range(1,a),2)
    for pair in pairs:
        i = pair[0]
        j = pair[1]
        G0.append(z[i]**2*z[j]+z[i]*z[j]**2)
    triplets = combinations(range(1,a),3)
    for triple in triplets:
        i = triple[0]
        j = triple[1]
        l = triple[2]
        G0.append(z[i]*z[j]*z[l])
    G = [g for g in list(G0)]
    for i in range(1,k):
        replacements = [(z[j],z[a*i+j]) for j in range(a)]
        for g in list(G0):
            G.append(g.subs(replacements))
    squares = [zi**2 for zi in z]
    sum_squares = sum(squares)
    fs = [sum_squares-2*i for i in range(1,k+1)]
    return(G,fs,z)

def make_matrix(R,k,a):
    #Converts list of one-hot encodings to the linear system
    temp = [r.flatten('F') for r in R]
    for i in range(k):
        added_row = np.zeros(k*a)
        for j in range(a):
            #This is the 2nd condition where sum(zi)_i*a+1^i*a+a = 0
            added_row[i*a+j] = 1
        temp.append(list(added_row))
    return(np.array(temp,dtype = int))

def check_resolving_grobner(R,k,a,alphabet = None, procs=1):
    print('Create Grobner')
    (G,fs,z) = create_groebner(k,a)
    print('One Hot')
    OH_encodedR = OneHot(R,k,a,alphabet = alphabet)
    A = make_matrix(OH_encodedR,k,a)
    print('lin equations')
    lin_fcns = make_linearEqns(A,z)
    n = len(lin_fcns)
    #print("Made Linear Functions")
    print('Grobner basis')
    G = groebnerbasis(True,G+lin_fcns,n,z, order = 'lex')
    #print("Finished Linear Functions")
    if procs>1:
        print('make jobs')
        jobs = [(i, True, list(G)+[f], 1, z) for i,f in enumerate(fs)] 
        print('start pool')
        pool = mp.Pool(processes=procs)
        results = pool.map_async(grobnerbasisPar, jobs)
        results = results.get()
        pool.close()
        pool.join()
        print('check result list')
        for Gi in results:
            if Gi!=[1]: #return False
                return False
#contructive algo, binary search, increment
#test parallel grobner
    else:
        for i,fi in enumerate(fs):
            Gi = groebnerbasis(True,list(G)+[fi],1,z, order = 'lex')
            if list(Gi) != [1]:
                return False
    return True

def grobnerbasisPar(args):
    print('job num', args[0], 'args', args[1:])
    with open('grobner_octamer_check.log', 'a+') as o:
      o.write('\t'.join(map(str, args))+'\n')
    return list(groebnerbasis(*args[1:], order = 'lex'))

###############################################################################################
#                                   HYPERCUBE FUNCTIONS                                       #
###############################################################################################

def hypercube_matrix(R):
    A = np.zeros((len(R),len(R[0])))
    for i,r in enumerate(R):
        for j,letter in enumerate(r):
            if letter == '0':
                A[i,j] = -1
            else:
                A[i,j] = int(letter)
    return(A)

def hypercube_polys(k):
    variable_string = ''
    num_digits = len(str(k))
    for i in range(0,k):
        count = str(i+1).zfill(num_digits)
        variable = ',z'+count
        variable_string = variable_string +variable
    variable_string = variable_string[1:]
    z = var(variable_string)
    P = []
    f = 0
    for i in range(k):
        P.append(z[i]*(z[i]-1)*(z[i]+1))
        f = f + z[i]**2
    fs = []
    for i in range(floor(k/2)):
        fi = f-2*(i+1)
        fs.append(fi)
    return(P,fs,z)

def check_hcube_resolving(R,k):
    #Create matrix A from R
    A = hypercube_matrix(R)
    #Get pre-computed Groebner basis and variables for H_k,2
    G,fs,z = hypercube_polys(k)
    #Get linear functions from A matrix
    lin_fcns = make_linearEqns(A,z)
    n = len(lin_fcns)
    #Get Grobner basis of P and linear functions
    G = groebnerbasis(True,G+lin_fcns,n,order = 'lex')
    for i,fi in enumerate(fs):
        #Compute Grobner basis of G+fi
        Gi = groebnerbasis(True,list(G)+[fi],1,order = 'lex')
        #Solutions iff Gi neq 1, if Gi neq 1 then R is not resolving
        if not (list(Gi) == [1]):
            return False
    return True

###############################################################################################
#                                   BRUTE-FORCE FUNCTIONS                                     #
###############################################################################################

def dist_mat(nodes,R):
    distMat = np.zeros((len(nodes), len(R)))
    for i,n in enumerate(nodes):
        for j,r in enumerate(R):
            distMat[i,j] = sum(1 for c1, c2 in zip(n, r) if c1 != c2)
    return(distMat)

def brute_force_resolve(R,k,a, alphabet = None):
    if alphabet == None:
        temp = [str(i) for i in range(a)]
        alphabet = ''.join(temp)
    H = [''.join(x) for x in product(alphabet,repeat = k)]
    distMat = dist_mat(H,R)
    n = len(H)
    for i in range(n-1):
        for j in range(i+1,n):
            if all(distMat[i,:] == distMat[j,:]):
                return(False)
    return(True)



