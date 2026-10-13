"""Spectral scalar-vector-tensor decomposition on a periodic cubic box.

Domain: [0, L)^3 with periodic boundary conditions, so k = (2 pi / L) n, n in Z^3.
On this domain the projection operators are diagonal in Fourier space and the
decomposition is unique up to the k = 0 mode.

Special modes
-------------
* k = 0: no direction exists; the constant tensor is split as
  h_0 = -2 psi_0 delta + gamma_0 with psi_0 = -tr(h_0)/6 and gamma_0 traceless
  (E = F = 0).  Whether this homogeneous part is physical depends on boundary
  conditions (it is a large gauge/background shift in a closed universe).
* Nyquist modes (even N): the spectral derivative i k at the Nyquist frequency is
  ambiguous; those entries are handled consistently by using the same k grid in
  decomposition and reconstruction, which makes the round trip exact.
* Small |k|: E ~ 1/k^2 and F ~ 1/k amplify errors; see ``noise_amplification``.

Non-periodic domains: the scalar/vector split is a Helmholtz decomposition whose
uniqueness depends on falloff conditions; none of the routines here apply there.
"""

from __future__ import annotations

import numpy as np


def kgrid(N: int, L: float = 2 * np.pi):
    freq = 2 * np.pi * np.fft.fftfreq(N, d=L / N)
    kx, ky, kz = np.meshgrid(freq, freq, freq, indexing="ij")
    return np.stack([kx, ky, kz], axis=-1)  # (N,N,N,3)


def decompose_fourier(hk, kvec):
    """SVT components of hk (shape (...,3,3), complex) at wavevectors kvec (...,3).

    Returns dict(psi, E, F (...,3), gamma (...,3,3)).  k = 0 entries use the homogeneous rule.
    """
    k2 = np.sum(kvec**2, axis=-1)
    zero = k2 == 0
    k2s = np.where(zero, 1.0, k2)
    kn = np.sqrt(k2s)
    kh = kvec / kn[..., None]
    I3 = np.eye(3)
    P = I3 - kh[..., :, None] * kh[..., None, :]
    T = np.einsum("...ii->...", hk)
    Lc = np.einsum("...i,...ij,...j->...", kh, hk, kh)
    V = np.einsum("...il,...l->...i", P, np.einsum("...lj,...j->...l", hk, kh))
    PhP = np.einsum("...ia,...ab,...bj->...ij", P, hk, P)
    trPhP = np.einsum("...ii->...", PhP)
    gam = PhP - 0.5 * P * trPhP[..., None, None]
    psi = -(T - Lc) / 4.0
    E = -(3 * Lc - T) / (4.0 * k2s)
    F = -1j * V / kn[..., None]
    # homogeneous mode
    psi = np.where(zero, -T / 6.0, psi)
    E = np.where(zero, 0.0, E)
    F = np.where(zero[..., None], 0.0, F)
    gam = np.where(zero[..., None, None], hk - (T[..., None, None] / 3.0) * I3, gam)
    return dict(psi=psi, E=E, F=F, gamma=gam)


def reconstruct_fourier(comp, kvec):
    I3 = np.eye(3)
    psi, E, F, gam = comp["psi"], comp["E"], comp["F"], comp["gamma"]
    kk = kvec[..., :, None] * kvec[..., None, :]
    kF = kvec[..., :, None] * F[..., None, :]
    return (-2 * psi[..., None, None] * I3 - 2 * E[..., None, None] * kk
            + 1j * (kF + np.swapaxes(kF, -1, -2)) + gam)


def decompose_periodic(h, L: float = 2 * np.pi):
    """h has shape (3,3,N,N,N) real; returns SVT parts in real space."""
    N = h.shape[-1]
    hk = np.fft.fftn(h, axes=(2, 3, 4))
    hk = np.moveaxis(hk, (0, 1), (-2, -1))  # (N,N,N,3,3)
    kv = kgrid(N, L)
    comp = decompose_fourier(hk, kv)
    return comp, kv


def to_real(comp_k):
    """Inverse FFT of each component (real parts returned)."""
    out = {}
    out["psi"] = np.fft.ifftn(comp_k["psi"]).real
    out["E"] = np.fft.ifftn(comp_k["E"]).real
    out["F"] = np.stack([np.fft.ifftn(comp_k["F"][..., i]).real for i in range(3)])
    g = comp_k["gamma"]
    out["gamma"] = np.stack([[np.fft.ifftn(g[..., i, j]).real for j in range(3)] for i in range(3)])
    return out


def reconstruction_error(h, L: float = 2 * np.pi):
    comp, kv = decompose_periodic(h, L)
    rec_k = reconstruct_fourier(comp, kv)
    rec = np.stack([[np.fft.ifftn(rec_k[..., i, j]).real for j in range(3)] for i in range(3)])
    return float(np.max(np.abs(rec - h)))


def spectral_divergence(gamma_real, L: float = 2 * np.pi):
    """d_j gamma_ij computed spectrally."""
    N = gamma_real.shape[-1]
    kv = kgrid(N, L)
    out = []
    for i in range(3):
        acc = 0
        for j in range(3):
            gk = np.fft.fftn(gamma_real[i, j])
            acc = acc + np.fft.ifftn(1j * kv[..., j] * gk).real
        out.append(acc)
    return np.stack(out)


def finite_difference_divergence(gamma_real, L: float = 2 * np.pi):
    """Second-order centered d_j gamma_ij on the periodic grid."""
    N = gamma_real.shape[-1]
    dx = L / N
    out = []
    for i in range(3):
        acc = 0
        for j in range(3):
            g = gamma_real[i, j]
            acc = acc + (np.roll(g, -1, axis=j) - np.roll(g, 1, axis=j)) / (2 * dx)
        out.append(acc)
    return np.stack(out)


def random_symmetric_field(N: int, seed: int = 0, kmax: int = 3, L: float = 2 * np.pi):
    """Smooth band-limited random real symmetric tensor field (3,3,N,N,N)."""
    rng = np.random.default_rng(seed)
    kv = kgrid(N, L)
    mask = (np.abs(kv[..., 0]) <= kmax) & (np.abs(kv[..., 1]) <= kmax) & (np.abs(kv[..., 2]) <= kmax)
    h = np.zeros((3, 3, N, N, N))
    for i in range(3):
        for j in range(i, 3):
            f = np.fft.fftn(rng.standard_normal((N, N, N))) * mask / N**1.5
            r = np.fft.ifftn(f).real
            h[i, j] = r
            h[j, i] = r
    return h


def tt_part_field(N: int, seed: int = 1, kmax: int = 3, L: float = 2 * np.pi):
    h = random_symmetric_field(N, seed, kmax, L)
    comp, kv = decompose_periodic(h, L)
    g = comp["gamma"].copy()
    g[0, 0, 0] = 0.0  # drop the homogeneous mode
    return np.stack([[np.fft.ifftn(g[..., i, j]).real for j in range(3)] for i in range(3)])


def noise_amplification(N: int = 24, eps: float = 1e-8, seed: int = 5, L: float = 2 * np.pi):
    """Error in recovered E and F induced by perturbing h by random noise of size eps.

    Returns (|k|, errE, errF) flattened over k != 0; expect errE ~ eps/k^2, errF ~ eps/k.
    """
    rng = np.random.default_rng(seed)
    kv = kgrid(N, L)
    hk = rng.standard_normal((N, N, N, 3, 3)) + 0j
    hk = 0.5 * (hk + np.swapaxes(hk, -1, -2))
    base = decompose_fourier(hk, kv)
    noise = rng.standard_normal((N, N, N, 3, 3))
    noise = 0.5 * (noise + np.swapaxes(noise, -1, -2)) * eps
    pert = decompose_fourier(hk + noise, kv)
    kn = np.linalg.norm(kv, axis=-1)
    m = kn > 0
    errE = np.abs(pert["E"] - base["E"])[m]
    errF = np.linalg.norm(pert["F"] - base["F"], axis=-1)[m]
    return kn[m], errE, errF
