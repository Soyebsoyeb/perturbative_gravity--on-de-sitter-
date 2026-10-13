# Scalar-vector-tensor decomposition

## Definition

For the spatial part of the perturbation, write $h_{ij}=a^2\big[-2\psi\delta_{ij}+2\partial_i\partial_jE+\partial_iF_j+\partial_jF_i+\gamma_{ij}\big]$ with $\partial_iF_i=0$, $\partial_i\gamma_{ij}=0$, $\gamma_{ii}=0$. In Fourier space ($\partial_j\to ik_j$, $\hat k=\mathbf k/k$, $P_{ij}=\delta_{ij}-\hat k_i\hat k_j$) define for a symmetric $h_{ij}(\mathbf k)$ (the factor $a^2$ is divided out):

$$
T=h_{ii},\quad L=\hat k_i\hat k_jh_{ij},\quad V_i=P_{il}\hat k_jh_{lj},\quad \Lambda_{ij,kl}=\tfrac12(P_{ik}P_{jl}+P_{il}P_{jk}-P_{ij}P_{kl}).
$$

Then

$$
\psi=-\frac{T-L}{4},\qquad E=-\frac{3L-T}{4k^2},\qquad F_j=-\frac{iV_j}{k},\qquad \gamma_{ij}=\Lambda_{ij,kl}h_{kl}.
$$

Reconstruction: $h_{ij}=-2\psi\delta_{ij}-2k_ik_jE+i(k_iF_j+k_jF_i)+\gamma_{ij}$. (Check: the coefficient of $P_{ij}$ is $-2\psi=(T-L)/2$; the coefficient of $\hat k_i\hat k_j$ is $-2\psi-2k^2E=L$.)

## Proven identities

For a fully generic symmetric $h$ (six independent symbols) and generic wavevector (three symbols) the following are proven by simplification to zero [`svt.*`]: $P^2=P$; $Pk=0$; $\gamma_{ij}k_j=0$; $\gamma_{ii}=0$; $\gamma_{ij}=\gamma_{ji}$; $k_iF_i=0$; exact reconstruction $h=-2\psi\delta-2kkE+i(kF+Fk)+\gamma$; $\Lambda$ is idempotent; $\Lambda$ annihilates the scalar and vector parts.

## Numerical decomposition on a periodic box

`perturbations.svt_numeric` decomposes a real symmetric field on a periodic cube with the FFT.

| question | result | where |
|---|---|---|
| reconstruction error | $\lesssim5\times10^{-17}$ for $N=8$--$24$ (algebraic, exact up to round-off) | `num.svt_roundtrip` |
| transversality of the TT part | spectral divergence $\sim10^{-17}$; second-order finite-difference divergence $\sim2\times10^{-3}$ at $N=16$ and falling as $N^{-2}$ | `num.tt_transverse`, `num.tt_fd_not_transverse`, test `test_finite_difference_divergence_is_second_order` |
| tracelessness | $\sim10^{-18}$ | `num.tt_traceless` |
| small-$\lvert\mathbf k\rvert$ behaviour | error in $E$ grows as $k^{-2.2}$, in $F$ as $k^{-0.96}$ for a perturbation of fixed size | `num.small_k_amplification_*` |
| zero momentum | homogeneous tensor split as $-2\psi_0\delta+\gamma_0$, $E=F=0$ | `test_zero_mode_handling` |

Dependence on resolution: the projection is exact in Fourier space, so the grid size only matters through how the *discrete derivative operator* is defined. A TT field built spectrally is exactly transverse in the spectral sense but not under a second-order stencil. The projections must be used with the derivative operator they were built from.

## Role of boundary conditions

On a periodic box the decomposition is unique except for $\mathbf k=0$. On a non-periodic domain, $E$ and $F_i$ are solutions of Poisson-type problems ($\nabla^2$ must be inverted), and uniqueness depends on falloff or boundary conditions: the kernel of the Laplacian (harmonic functions) leaves ambiguity unless excluded by hypothesis. None of the routines here apply to such domains.

## Physical status of the three parts

It is **not** true that scalar and vector perturbations always propagate as independent gravitational radiation in vacuum. In vacuum de Sitter space the scalar and vector parts are constrained (their gauge-invariant combinations vanish for $k\neq0$, `docs/constraint_equations.md`), and only the tensor part propagates. The decomposition is a kinematic statement about symmetric tensors, not a dynamical one.
