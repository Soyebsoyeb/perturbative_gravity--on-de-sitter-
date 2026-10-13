"""Constraint analysis of the vacuum linearized equations on de Sitter space.

Scalar sector (Fourier mode, k != 0), written with the verified identities of
``scalar_sector``:

    traceless ij :  k_i k_j (Phi_B - Psi_B) = 0      =>  Phi_B = Psi_B
    0i           :  Psi_B' + calH Phi_B     = 0      =>  (a Psi_B)' = 0
    00           :  -2 k^2 Psi_B - 6 calH (Psi_B' + calH Phi_B) = 0   =>  Psi_B = 0

Hence Phi_B = Psi_B = 0: the scalar sector carries no gauge-invariant degree of
freedom in vacuum (it is pure gauge plus constraints).

Vector sector: -(1/2) d_z^2 V = 0 => V = 0 for k != 0.

Tensor sector: no constraint removes gamma_ij; it propagates (2 polarizations).

Counting (4D GR): 10 metric components - 4 constraints - 4 gauge = 2 dynamical
degrees of freedom, the two tensor polarizations.

Zero-momentum (k = 0) modes are excluded above: there d_i d_j drops out, the
constraints no longer fix Psi_B, and the corresponding modes are homogeneous
rescalings/shifts associated with large gauge transformations (boundary
conditions decide whether they are physical).
"""

from __future__ import annotations

import functools

import sympy as sp

from .scalar_sector import newtonian_gauge_dE


@functools.lru_cache(maxsize=None)
def scalar_constraint_solution():
    """Solve the Newtonian-gauge constraints for a plane wave Psi = psi(eta) cos(kz)."""
    d = newtonian_gauge_dE()
    s = d["setup"]
    eta, z = s["eta"], s["z"]
    k = sp.Symbol("k", positive=True)
    Phi_f = sp.Function("Phi")(eta, z)
    Psi_f = sp.Function("Psi")(eta, z)
    phi = sp.Function("phi")(eta)
    ps = sp.Function("ps")(eta)
    subs = {Phi_f: phi * sp.cos(k * z), Psi_f: ps * sp.cos(k * z)}
    dE = d["dE"]

    def sub(e):
        return sp.simplify(e.subs(subs).doit())

    e00 = sub(dE[0, 0]) / sp.cos(k * z)
    e0z = sub(dE[0, 3]) / (-k * sp.sin(k * z))
    ett = sub(dE[1, 1] - dE[3, 3]) / (-k**2 * sp.cos(k * z))
    e00, e0z, ett = (sp.simplify(e) for e in (e00, e0z, ett))
    # solve: ett = phi - ps = 0
    sol_phi = sp.solve(sp.Eq(ett, 0), phi)[0]
    e0z_s = sp.simplify(e0z.subs(phi, sol_phi))
    e00_s = sp.simplify(e00.subs(phi, sol_phi))
    # with 0z: ps' + calH ps = 0 =>  ps' = -calH ps ; insert into 00
    ps_p = sp.solve(sp.Eq(e0z_s, 0), sp.diff(ps, eta))[0]
    e00_final = sp.simplify(e00_s.subs(sp.diff(ps, eta), ps_p))
    return {"e00": e00, "e0z": e0z, "e_traceless": ett, "phi_solution": sol_phi,
            "psi_prime_solution": ps_p, "e00_after_constraints": e00_final, "k": k, "ps": ps}
