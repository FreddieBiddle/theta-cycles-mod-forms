"""
Computation of theta cycles of modular forms using SageMath.
"""

from functools import lru_cache
from sage.all import *

NN = 10000
R = LaurentSeriesRing(QQ, default_prec=NN, names=('q',)); (q,) = R._first_ngens(1)


@lru_cache(maxsize=None)
def modular_form_basis(k0, Q, bd):
    """
    Cached basis computation
    """
    return [[c % Q for c in coeffs(g, bd)] for g in ModularForms(1, k0).q_expansion_basis(bd)]


def coeffs(qexp, prec):
    """
    Gives list of coeffs of q expansion that includes 0-coeffs
    """
    coeffs_list = qexp.list()
    return coeffs_list[:prec] + [0] * max(0, prec - len(coeffs_list))


def theta(f, k, Q, prec=None):
    """
    Ramanujan theta operator
    """
    p = Q.prime_factors()[0]
    # Sturm-type bound from Chen-Kiming-Rasmussen: "On congruences mod p^m between eigenforms and
    # their attached Galois representations"
    bd = int((k * Gamma1(p).index()) / 12) + 1
    if prec is None:
        prec = bd
    return sum(((n * f[n])%Q) * q**n for n in range(len(coeffs(f, prec))))


def phi(Q):
    """
    Euler totient function
    """
    return Q * prod(1 - 1/p for p in Q.prime_factors())


def filt(f, k, Q, prec=None):
    """
    Weight filtration of f modulo Q

    Args:
        f (sage.rings.power_series_poly.PowerSeries_poly): input modular form q-series
        k (int): input weight of f mod Q (need not be f.weight() if e.g. f = theta(g))
        Q (int): power of prime

    Output:
        int: weight filtration of f modulo Q
    """
    p = Q.prime_factors()[0]
    bd = int((k * Gamma1(p).index()) / 12) + 1
    if prec is None:
        prec = bd
    # Get f's q_expansion modulo Q
    f_qexp = sum((f[n] % Q) * q**n for n in range(len(coeffs(f, prec))))
    # if f == 1 + O(q^bd) then f == 0 (mod Q)
    if f_qexp.coefficients() == [1]:
        return 0
    # Initialize list of lower weights for filtration candidates
    weights = []
    phi_Q = phi(Q)
    for k0 in [i for i in range(1, k) if (i-k) % phi_Q == 0][::-1]:
        M_k0_Q = modular_form_basis(k0, Q, prec)
        f_modQ = vector(QQ, coeffs(f_qexp, bd))[:prec]
        V = VectorSpace(QQ, prec)
        M_k0_vecs = V.subspace([V(v) for v in M_k0_Q])
        if f_modQ in M_k0_vecs:
            weights.append(k0)
    if weights == []:
        return k
    elif min(weights) == 2:
        return 0
    else:
        return min(weights)


def theta_cycle(f, k, Q, prec=None):
    """
    Theta cycle of f modulo Q

    Args:
        f (sage.rings.power_series_poly.PowerSeries_poly): input modular form q-series
        k (int): input weight of f mod Q (need not be f.weight() if e.g. f = theta(g))
        Q (int): power of prime

    Output:
        list: weight filtrations of theta^i(f) for 0 < i < phi(Q)+1
    """
    p = Q.prime_factors()[0]
    phi_Q = phi(Q)
    bd = int((k * Gamma1(p).index()) / 12) + 1
    if prec is None:
        prec = bd

    # Determine weight of image of theta(f) modulo Q (see Chen-Kiming)
    if Q.is_prime():
        init_wt = filt(f, k, Q, prec) + p + 1
        weight_step = p + 1
    else:
        init_wt = filt(f, k, Q, prec) + 2 + 2*phi_Q
        weight_step = 2 + 2*phi_Q

    # Initializing theta cycle
    cycle = []
    F = theta(f, k, Q, prec)
    cycle.append(filt(F, init_wt, Q))

    for _ in range(1, phi_Q):
        next_wt = cycle[-1] + weight_step
        F = theta(F, next_wt, Q, prec)
        cycle.append(filt(F, next_wt, Q))

    return cycle
