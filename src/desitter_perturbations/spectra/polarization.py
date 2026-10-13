"""Polarization tensors and their normalisation.

Convention: gamma_ij(k) = sum_lam e^lam_ij(khat) gamma_lam(k), with real linear
polarizations e^+ and e^x obeying

    e^lam_ij e^lam'_ij = 2 delta^{lam lam'},   khat_i e^lam_ij = 0,   e^lam_ii = 0,
    sum_lam e^lam_ij e^lam_kl = 2 Lambda_{ij,kl}(khat).

The factor 2 in the normalisation is what makes z_T^2 = M_pl^2 a^2 / 2 and the
spectrum Delta_t^2 = 2 H^2/(pi^2 M_pl^2) come out with their standard values.
Circular polarizations e^{R,L} = (e^+ +- i e^x)/sqrt(2) obey e^R_ij e^{R*}_ij = 2.
"""

from __future__ import annotations

import numpy as np

from ..perturbations.tensor_sector import polarization_basis, tt_projector


def circular_basis(k):
    ep, ex = polarization_basis(k)
    eR = (ep + 1j * ex) / np.sqrt(2.0)
    eL = (ep - 1j * ex) / np.sqrt(2.0)
    return eR, eL


def check_normalization(k):
    ep, ex = polarization_basis(k)
    G = np.array([[np.sum(a * b) for b in (ep, ex)] for a in (ep, ex)])
    return float(np.max(np.abs(G - 2 * np.eye(2))))


def check_transverse_traceless(k):
    ep, ex = polarization_basis(k)
    khat = np.asarray(k, float) / np.linalg.norm(k)
    t = max(np.max(np.abs(khat @ e)) for e in (ep, ex))
    tr = max(abs(np.trace(e)) for e in (ep, ex))
    return float(t), float(tr)


def check_completeness(k):
    """max |sum_lam e_ij e_kl - 2 Lambda_{ij,kl}|."""
    ep, ex = polarization_basis(k)
    S = np.einsum("ij,kl->ijkl", ep, ep) + np.einsum("ij,kl->ijkl", ex, ex)
    return float(np.max(np.abs(S - 2 * tt_projector(k))))


def helicity_eigenvalue(k, e):
    """Eigenvalue h of the rotation generator about khat acting on e (R: +2, L: -2)."""
    khat = np.asarray(k, float) / np.linalg.norm(k)
    J = np.array([[0, -khat[2], khat[1]], [khat[2], 0, -khat[0]], [-khat[1], khat[0], 0]])
    # generator acting on a rank-2 tensor: (J e + e J^T) * (i)
    Je = 1j * (J @ e + e @ J.T)
    num = np.sum(np.conj(e) * Je)
    den = np.sum(np.conj(e) * e)
    return complex(num / den)
