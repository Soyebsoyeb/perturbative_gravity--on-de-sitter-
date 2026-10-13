"""Symbolic scalar-vector-tensor algebra in Fourier space.

For a symmetric spatial tensor h_ij(k) with k != 0 define khat = k/|k|,
P_ij = delta_ij - khat_i khat_j, and (T = h_ii, L = khat_i khat_j h_ij,
V_i = P_il khat_j h_lj):

    psi   = -(T - L)/4
    E     = -(3 L - T)/(4 k^2)
    F_j   = -i V_j / k
    gamma = Lambda h,    Lambda_{ij,kl} = (1/2)(P_ik P_jl + P_il P_jk - P_ij P_kl)

with the Fourier convention d_i -> i k_i so that

    h_ij = -2 psi delta_ij - 2 k_i k_j E + i (k_i F_j + k_j F_i) + gamma_ij.

``verify_svt_identities`` proves these statements for a fully generic symmetric h
(six independent symbols) and a generic wavevector; no special direction is used.
"""

from __future__ import annotations

import functools

import itertools
import sympy as sp


def symbols_k():
    k1, k2, k3 = sp.symbols("k1 k2 k3", real=True)
    return sp.Matrix([k1, k2, k3])


def generic_symmetric_h():
    h = sp.zeros(3, 3)
    for i in range(3):
        for j in range(i, 3):
            s = sp.Symbol(f"h{i+1}{j+1}")
            h[i, j] = h[j, i] = s
    return h


def projectors(k):
    k2 = (k.T * k)[0]
    P = sp.eye(3) - k * k.T / k2
    return P, k2


def tt_projector_action(P, h):
    """Lambda h = 1/2 (P h P + (P h P)^T - P tr(P h))  for symmetric h  (= P h P - P tr(P h P)/2)."""
    PhP = P * h * P
    return (PhP - P * (P * h).trace() / 2)


def svt_components(h, k):
    P, k2 = projectors(k)
    kn = sp.sqrt(k2)
    kh = k / kn
    T = h.trace()
    L = (kh.T * h * kh)[0]
    V = P * h * kh
    psi = -(T - L) / 4
    E = -(3 * L - T) / (4 * k2)
    F = -sp.I * V / kn
    gam = tt_projector_action(P, h)
    return psi, E, F, gam


def svt_reconstruct(psi, E, F, gam, k):
    h = -2 * psi * sp.eye(3) - 2 * E * (k * k.T) + sp.I * (k * F.T + F * k.T) + gam
    return h


@functools.lru_cache(maxsize=None)
def verify_svt_identities():
    """Return a dict of residuals (all should simplify to zero)."""
    k = symbols_k()
    h = generic_symmetric_h()
    P, k2 = projectors(k)
    psi, E, F, gam = svt_components(h, k)
    kn = sp.sqrt(k2)
    res = {}
    res["P_idempotent"] = (P * P - P).applyfunc(sp.simplify)
    res["P_transverse"] = (P * k).applyfunc(sp.simplify)
    res["gamma_transverse"] = (gam * k).applyfunc(sp.simplify)
    res["gamma_traceless"] = sp.simplify(gam.trace())
    res["gamma_symmetric"] = (gam - gam.T).applyfunc(sp.simplify)
    res["F_transverse"] = sp.simplify((k.T * F)[0])
    rec = svt_reconstruct(psi, E, F, gam, k)
    res["reconstruction"] = (rec - h).applyfunc(sp.simplify)
    # Lambda is a projector: Lambda(Lambda h) = Lambda h
    res["Lambda_idempotent"] = (tt_projector_action(P, gam) - gam).applyfunc(sp.simplify)
    # Lambda annihilates the scalar and vector parts
    scal = -2 * psi * sp.eye(3) - 2 * E * (k * k.T)
    vec = sp.I * (k * F.T + F * k.T)
    res["Lambda_kills_scalar"] = tt_projector_action(P, scal).applyfunc(sp.simplify)
    res["Lambda_kills_vector"] = tt_projector_action(P, vec).applyfunc(sp.simplify)
    return res


def all_zero(res):
    def z(v):
        if isinstance(v, sp.MatrixBase):
            return all(e == 0 for e in v)
        return v == 0
    return {k: z(v) for k, v in res.items()}
