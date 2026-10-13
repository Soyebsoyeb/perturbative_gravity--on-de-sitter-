# Validation report

Checks: 128; passed: 128; failed: 0; of which negative controls: 9.

Environment: python 3.13.16, numpy 2.5.3, scipy 1.18.1, sympy 1.14.0, platform Linux-6.18.44-fc-v114-x86_64-with-glibc2.39, seconds 73.3

Symbolic residuals are 0 (identically zero after simplification) or 1. A *negative control* passes only when the deliberately wrong formula is detected, i.e. its residual exceeds the tolerance.

## background

| id | kind | description | residual | tolerance | result |
|---|---|---|---|---|---|
| `bg.time_relation` | identity | `d eta/dt = 1/a and a(eta(t)) = exp(H t) with eta = -exp(-H t)/H` | 0.000e+00 | 0.0e+00 | PASS |
| `bg.ricci_scalar.cosmic` | identity | `Ricci scalar = 12 H^2 (cosmic slicing)` | 0.000e+00 | 0.0e+00 | PASS |
| `bg.ricci_tensor.cosmic` | identity | `R_mn = 3 H^2 g_mn (cosmic)` | 0.000e+00 | 0.0e+00 | PASS |
| `bg.einstein.cosmic` | identity | `G_mn = -3 H^2 g_mn (cosmic)` | 0.000e+00 | 0.0e+00 | PASS |
| `bg.field_eq.cosmic` | identity | `G_mn + Lambda g_mn = 0 with Lambda = 3H^2 (cosmic)` | 0.000e+00 | 0.0e+00 | PASS |
| `bg.constant_curvature.cosmic` | identity | `R_mnrs = H^2 (g_mr g_ns - g_ms g_nr), all 256 components (cosmic)` | 0.000e+00 | 0.0e+00 | PASS |
| `bg.kretschmann.cosmic` | identity | `R_abcd R^abcd = 24 H^4 (cosmic)` | 0.000e+00 | 0.0e+00 | PASS |
| `bg.ricci_squared.cosmic` | identity | `R_ab R^ab = 36 H^4 (cosmic)` | 0.000e+00 | 0.0e+00 | PASS |
| `bg.ricci_scalar.conformal` | identity | `Ricci scalar = 12 H^2 (conformal slicing)` | 0.000e+00 | 0.0e+00 | PASS |
| `bg.ricci_tensor.conformal` | identity | `R_mn = 3 H^2 g_mn (conformal)` | 0.000e+00 | 0.0e+00 | PASS |
| `bg.einstein.conformal` | identity | `G_mn = -3 H^2 g_mn (conformal)` | 0.000e+00 | 0.0e+00 | PASS |
| `bg.field_eq.conformal` | identity | `G_mn + Lambda g_mn = 0 with Lambda = 3H^2 (conformal)` | 0.000e+00 | 0.0e+00 | PASS |
| `bg.constant_curvature.conformal` | identity | `R_mnrs = H^2 (g_mr g_ns - g_ms g_nr), all 256 components (conformal)` | 0.000e+00 | 0.0e+00 | PASS |
| `bg.kretschmann.conformal` | identity | `R_abcd R^abcd = 24 H^4 (conformal)` | 0.000e+00 | 0.0e+00 | PASS |
| `bg.ricci_squared.conformal` | identity | `R_ab R^ab = 36 H^4 (conformal)` | 0.000e+00 | 0.0e+00 | PASS |
| `bg.sign_control` | negative_control | `Flipped Riemann sign gives R = -12 H^2 (convention test is sensitive)` | 1.000e+00 | 0.0e+00 | PASS |
| `bg.slicing_invariants` | identity | `Curvature scalars agree between cosmic-time and conformal-time slicings` | 0.000e+00 | 0.0e+00 | PASS |
| `bg.embedding` | identity | `Flat slicing embeds in the hyperboloid -X0^2 + X.X + X4^2 = 1/H^2 with induced metric -dt^2 + a^2 dx^2` | 0.000e+00 | 0.0e+00 | PASS |

## metric_expansion

| id | kind | description | residual | tolerance | result |
|---|---|---|---|---|---|
| `exp.inverse_metric` | identity | `g^{mn} = gbar^{mn} - eps h^{mn} + eps^2 h^m_r h^{rn} + O(eps^3) (generic 4D, 10 symbols)` | 0.000e+00 | 0.0e+00 | PASS |
| `exp.inverse_metric.neg_sign` | negative_control | `Wrong sign at second order is detected` | 1.000e+00 | 0.0e+00 | PASS |
| `exp.inverse_metric.neg_index` | negative_control | `Transposed index contraction at second order is detected (non-diagonal background)` | 1.000e+00 | 0.0e+00 | PASS |
| `exp.inverse_metric.nondiagonal` | identity | `Inverse-metric expansion also holds on a non-diagonal background` | 0.000e+00 | 0.0e+00 | PASS |
| `exp.sqrt_det` | identity | `sqrt(-g)/sqrt(-gbar) = 1 + eps h/2 + eps^2 (h^2/8 - h_mn h^mn/4) from the exact determinant` | 0.000e+00 | 0.0e+00 | PASS |
| `exp.sqrt_det.neg_factor` | negative_control | `Factor-of-two error in the h_mn h^mn term is detected` | 1.000e+00 | 0.0e+00 | PASS |
| `exp.two_dim_example` | identity | `Independent explicit 2D example (gbar = diag(-1, a^2))` | 0.000e+00 | 0.0e+00 | PASS |

## linearized_einstein

| id | kind | description | residual | tolerance | result |
|---|---|---|---|---|---|
| `lin.background_zero` | identity | `Background field equation holds in the truncated brute-force pipeline (O(eps^0) part)` | 0.000e+00 | 0.0e+00 | PASS |
| `lin.cov_vs_component` | identity | `Background-covariant dE_mn equals component-level coefficient of eps in G+Lambda g` | 0.000e+00 | 0.0e+00 | PASS |
| `lin.delta_gamma` | identity | `dGamma^l_mn = 1/2 g^{lr}(nabla_m h_rn + nabla_n h_rm - nabla_r h_mn) vs direct derivative of Gamma` | 0.000e+00 | 0.0e+00 | PASS |
| `lin.dE_definition` | identity | `dE_mn = dG_mn + Lambda h_mn` | 0.000e+00 | 0.0e+00 | PASS |
| `lin.cov_vs_component.cosmic` | identity | `Covariant vs component-level linearization in cosmic-time coordinates` | 0.000e+00 | 0.0e+00 | PASS |

## bianchi

| id | kind | description | residual | tolerance | result |
|---|---|---|---|---|---|
| `bianchi.covariant_form` | identity | `gbar^{mr}[nabla_r dG_mn - dGamma^l_rm Gbar_ln - dGamma^l_rn Gbar_ml] = 0` | 0.000e+00 | 0.0e+00 | PASS |
| `bianchi.naive_divergence_nonzero` | negative_control | `Naive nabla^m dG_mn is NOT zero: connection-variation terms are required` | 1.000e+00 | 0.0e+00 | PASS |
| `bianchi.mixed_form` | identity | `Mixed-index linearized Bianchi identity` | 0.000e+00 | 0.0e+00 | PASS |
| `bianchi.field_equation_divergence` | identity | `nabla^m dE_mn = 0 exactly (no connection-variation terms because Ebar = 0)` | 0.000e+00 | 0.0e+00 | PASS |

## gauge

| id | kind | description | residual | tolerance | result |
|---|---|---|---|---|---|
| `gauge.lie_equals_covariant` | identity | `L_xi gbar = nabla_m xi_n + nabla_n xi_m` | 0.000e+00 | 0.0e+00 | PASS |
| `gauge.pure_gauge_solves_linearized` | identity | `Pure-gauge h = -(nabla xi + nabla xi) satisfies dE_mn = 0` | 0.000e+00 | 0.0e+00 | PASS |
| `gauge.field_equation_invariant` | identity | `dE[h + h_pg] = dE[h] (linearized field equation is gauge invariant)` | 0.000e+00 | 0.0e+00 | PASS |
| `gauge.svt_laws` | identity | `SVT transformation laws A,psi,B,E,S,F derived from L_xi gbar` | 0.000e+00 | 0.0e+00 | PASS |
| `gauge.synchronous_residual` | identity | `Synchronous gauge keeps A = 0 under T = C/a: residual gauge freedom exists (a T = const)` | 0.000e+00 | 0.0e+00 | PASS |
| `gauge.newtonian_no_residual` | identity | `Newtonian gauge (B = E = 0) preserved only by T = L = 0 for k != 0` | 0.000e+00 | 0.0e+00 | PASS |

## svt

| id | kind | description | residual | tolerance | result |
|---|---|---|---|---|---|
| `svt.P_idempotent` | identity | `Fourier-space SVT identity: P_idempotent` | 0.000e+00 | 0.0e+00 | PASS |
| `svt.P_transverse` | identity | `Fourier-space SVT identity: P_transverse` | 0.000e+00 | 0.0e+00 | PASS |
| `svt.gamma_transverse` | identity | `Fourier-space SVT identity: gamma_transverse` | 0.000e+00 | 0.0e+00 | PASS |
| `svt.gamma_traceless` | identity | `Fourier-space SVT identity: gamma_traceless` | 0.000e+00 | 0.0e+00 | PASS |
| `svt.gamma_symmetric` | identity | `Fourier-space SVT identity: gamma_symmetric` | 0.000e+00 | 0.0e+00 | PASS |
| `svt.F_transverse` | identity | `Fourier-space SVT identity: F_transverse` | 0.000e+00 | 0.0e+00 | PASS |
| `svt.reconstruction` | identity | `Fourier-space SVT identity: reconstruction` | 0.000e+00 | 0.0e+00 | PASS |
| `svt.Lambda_idempotent` | identity | `Fourier-space SVT identity: Lambda_idempotent` | 0.000e+00 | 0.0e+00 | PASS |
| `svt.Lambda_kills_scalar` | identity | `Fourier-space SVT identity: Lambda_kills_scalar` | 0.000e+00 | 0.0e+00 | PASS |
| `svt.Lambda_kills_vector` | identity | `Fourier-space SVT identity: Lambda_kills_vector` | 0.000e+00 | 0.0e+00 | PASS |

## sectors

| id | kind | description | residual | tolerance | result |
|---|---|---|---|---|---|
| `sector.scalar.00` | identity | `dE component 00 depends on h only through Bardeen potentials` | 0.000e+00 | 0.0e+00 | PASS |
| `sector.scalar.0z` | identity | `dE component 0z depends on h only through Bardeen potentials` | 0.000e+00 | 0.0e+00 | PASS |
| `sector.scalar.xx-zz` | identity | `dE component xx-zz depends on h only through Bardeen potentials` | 0.000e+00 | 0.0e+00 | PASS |
| `sector.vector.0x` | identity | `Vector-sector dE 0x depends only on V = S - F'` | 0.000e+00 | 0.0e+00 | PASS |
| `sector.vector.xz` | identity | `Vector-sector dE xz depends only on V = S - F'` | 0.000e+00 | 0.0e+00 | PASS |
| `sector.tensor.eom` | identity | `TT sector: dE_ij = (1/2)(gamma'' + 2 calH gamma' - nabla^2 gamma), other components vanish` | 0.000e+00 | 0.0e+00 | PASS |
| `sector.tensor.covariant` | identity | `TT ansatz is transverse and traceless, and dE_mn = -(1/2)(box h_mn - 2 H^2 h_mn) (massless spin-2 on de Sitter)` | 0.000e+00 | 0.0e+00 | PASS |

## constraints

| id | kind | description | residual | tolerance | result |
|---|---|---|---|---|---|
| `constraints.scalar_traceless` | identity | `Traceless ij equation gives Phi_B = Psi_B` | 0.000e+00 | 0.0e+00 | PASS |
| `constraints.scalar_vacuum` | identity | `After 0i and traceless constraints the 00 equation is -2 k^2 Psi_B = 0, so Psi_B = 0 for k != 0` | 0.000e+00 | 0.0e+00 | PASS |

## action

| id | kind | description | residual | tolerance | result |
|---|---|---|---|---|---|
| `action.background_lagrangian` | identity | `Background Lagrangian density sqrt(-g)(R - 2 Lambda) = 6/(H^2 eta^4)` | 0.000e+00 | 0.0e+00 | PASS |
| `action.first_order_tensor` | identity | `First-order term of the TT expansion vanishes identically` | 0.000e+00 | 0.0e+00 | PASS |
| `action.route1_vs_target` | identity | `Route 1 (direct expansion): L2 = (a^2/4)(gamma'^2 - (d gamma)^2) + total derivative` | 0.000e+00 | 0.0e+00 | PASS |
| `action.route2_vs_target` | identity | `Route 2 (linearized equations): L2 = -(1/2) sqrt(-gbar) h^mn dE_mn = target + total derivative` | 0.000e+00 | 0.0e+00 | PASS |
| `action.route1_vs_route2` | identity | `Routes 1 and 2 agree up to a total derivative` | 0.000e+00 | 0.0e+00 | PASS |
| `action.normalisation_control` | negative_control | `Negative control: doubling the target normalisation is detected` | 1.000e+00 | 0.0e+00 | PASS |
| `action.tensor_eom` | identity | `E-L equation of the tensor action is gamma'' + 2 calH gamma' + k^2 gamma = 0` | 0.000e+00 | 0.0e+00 | PASS |
| `action.first_order_total_derivative` | identity | `First-order bulk variation (generic diagonal h(eta)) is a total derivative on the background solution` | 0.000e+00 | 0.0e+00 | PASS |
| `action.first_order_wrong_lambda` | negative_control | `Negative control: with Lambda = 2H^2 the first-order variation is not a total derivative` | 1.000e+00 | 0.0e+00 | PASS |

## canonical

| id | kind | description | residual | tolerance | result |
|---|---|---|---|---|---|
| `canon.z_T_squared` | identity | `z_T^2 = M_pl^2 a^2/2 derived from the verified action` | 0.000e+00 | 0.0e+00 | PASS |
| `canon.effective_potential` | identity | `z_T''/z_T = 2/eta^2` | 0.000e+00 | 0.0e+00 | PASS |
| `canon.action_equivalence` | identity | `L[gamma = v/z] equals L_v up to a total derivative` | 0.000e+00 | 0.0e+00 | PASS |
| `canon.mode_equation` | identity | `E-L of the canonical action: v'' + (k^2 - 2/eta^2) v = 0` | 0.000e+00 | 0.0e+00 | PASS |
| `canon.bd_normalisation` | identity | `Wronskian condition fixes \|C\|^2 = H^2/(M_pl^2 k^3)` | 0.000e+00 | 0.0e+00 | PASS |
| `canon.exact_mode_solves_eom` | identity | `gamma = (1 + i k eta) e^{-i k eta} solves gamma'' - (2/eta) gamma' + k^2 gamma = 0` | 0.000e+00 | 0.0e+00 | PASS |
| `canon.positive_frequency` | identity | `sqrt(2k) v_BD e^{i k eta} -> 1 as k eta -> -infinity` | 0.000e+00 | 0.0e+00 | PASS |

## analytic

| id | kind | description | residual | tolerance | result |
|---|---|---|---|---|---|
| `bench.scale_factor_relation` | benchmark | `a(eta(t)) = exp(H t)` | 2.180e-16 | 1.0e-13 | PASS |
| `bench.time_roundtrip` | benchmark | `t -> eta -> t round trip` | 1.110e-16 | 1.0e-12 | PASS |
| `bench.deta_dt` | benchmark | `d eta/dt = 1/a (finite difference)` | 1.667e-07 | 1.0e-05 | PASS |
| `bench.mode_solves_ode` | benchmark | `Exact BD gamma satisfies gamma'' + 2 calH gamma' + k^2 gamma = 0` | 6.664e-16 | 1.0e-11 | PASS |
| `bench.dv_consistency` | benchmark | `Analytic dv_BD/d eta matches centred differences (step ~ 1e-4 of local time scale)` | 1.811e-08 | 1.0e-06 | PASS |
| `bench.v_equals_z_gamma` | benchmark | `\|z_T gamma_BD\| = \|v_BD\| (same state up to a constant phase)` | 6.429e-16 | 1.0e-12 | PASS |
| `bench.subhorizon.x100` | benchmark | `Subhorizon asymptotics \|gamma\|/\|gamma_WKB\| -> 1 with error ~ 1/(2 x^2) (x = 100)` | 5.000e-05 | 1.0e-04 | PASS |
| `bench.subhorizon.x1000` | benchmark | `Subhorizon asymptotics \|gamma\|/\|gamma_WKB\| -> 1 with error ~ 1/(2 x^2) (x = 1000)` | 5.000e-07 | 1.0e-06 | PASS |
| `bench.subhorizon.x10000` | benchmark | `Subhorizon asymptotics \|gamma\|/\|gamma_WKB\| -> 1 with error ~ 1/(2 x^2) (x = 10000)` | 5.000e-09 | 1.0e-08 | PASS |
| `bench.superhorizon.x0.1` | benchmark | `Superhorizon freeze-out \|gamma\| -> H/(M k^1.5) with error x^2/2 (x = 0.1)` | 1.244e-05 | 5.0e-05 | PASS |
| `bench.superhorizon.x0.01` | benchmark | `Superhorizon freeze-out \|gamma\| -> H/(M k^1.5) with error x^2/2 (x = 0.01)` | 1.250e-09 | 5.0e-09 | PASS |
| `bench.superhorizon.x0.001` | benchmark | `Superhorizon freeze-out \|gamma\| -> H/(M k^1.5) with error x^2/2 (x = 0.001)` | 1.249e-13 | 5.1e-13 | PASS |
| `bench.wronskian` | benchmark | `gamma gamma*' - c.c. = 2i H^2 eta^2/M^2` | 2.829e-16 | 1.0e-12 | PASS |
| `bench.canonical_wronskian` | benchmark | `z_T^2 (gamma gamma*' - c.c.) = i` | 4.441e-16 | 1.0e-12 | PASS |
| `bench.power_from_mode` | benchmark | `Delta^2 = 2 k^3\|gamma\|^2/pi^2 equals (2H^2/pi^2 M^2)(1+k^2 eta^2)` | 6.063e-16 | 1.0e-12 | PASS |
| `bench.power_superhorizon` | benchmark | `Delta^2 -> 2 H^2/(pi^2 M_pl^2) for k\|eta\| -> 0` | 1.000e-12 | 1.0e-11 | PASS |
| `bench.power_scaling` | benchmark | `Delta^2 scales as H^2/M_pl^2` | 0.000e+00 | 1.0e-14 | PASS |
| `bench.polarization_norm` | benchmark | `e^lam_ij e^lam'_ij = 2 delta (25 random directions)` | 1.332e-15 | 1.0e-13 | PASS |
| `bench.polarization_completeness` | benchmark | `sum_lam e_ij e_kl = 2 Lambda_{ij,kl}` | 5.551e-16 | 1.0e-13 | PASS |
| `bench.polarization_tt` | benchmark | `k_i e_ij = 0 and e_ii = 0` | 4.441e-16 | 1.0e-13 | PASS |
| `bench.helicity` | benchmark | `Circular polarizations have helicity +2 and -2` | 4.353e-17 | 1.0e-12 | PASS |
| `bench.lattice_commutator` | benchmark | `[gamma(x), pi(y)] = i delta_{xy}/dx on a periodic lattice with BD modes` | 1.744e-15 | 1.0e-12 | PASS |

## numerical

| id | kind | description | residual | tolerance | result |
|---|---|---|---|---|---|
| `num.exact_vs_conformal` | numerical | `Numerical conformal solution vs exact BD (rtol 1e-10)` | 7.605e-10 | 5.0e-09 | PASS |
| `num.exact_vs_cosmic` | numerical | `Numerical cosmic solution vs exact BD (rtol 1e-10)` | 5.533e-10 | 5.0e-09 | PASS |
| `num.exact_vs_efold` | numerical | `Numerical efold solution vs exact BD (rtol 1e-10)` | 5.534e-10 | 5.0e-09 | PASS |
| `num.exact_vs_canonical` | numerical | `Numerical canonical solution vs exact BD (rtol 1e-10)` | 7.857e-10 | 5.0e-09 | PASS |
| `num.cross_cosmic` | numerical | `cosmic vs conformal-time solution, pointwise` | 3.128e-10 | 5.0e-09 | PASS |
| `num.cross_efold` | numerical | `efold vs conformal-time solution, pointwise` | 3.127e-10 | 5.0e-09 | PASS |
| `num.cross_canonical` | numerical | `canonical vs conformal-time solution, pointwise` | 3.193e-10 | 5.0e-09 | PASS |
| `num.wavenumber_independence` | numerical | `Accuracy is independent of k when expressed in x = k/(aH) (k = 0.01, 1, 100)` | 7.767e-10 | 5.0e-09 | PASS |
| `num.tolerance_scaling` | numerical | `Max error <= 1000 x rtol for rtol in {1e-6,1e-8,1e-10}` | 8.413e+00 | 1.0e+03 | PASS |
| `num.tolerance_monotone` | numerical | `Error decreases monotonically with rtol` | 0.000e+00 | 0.0e+00 | PASS |
| `num.rk4_order` | numerical | `Fixed-step RK4 observed order is 4` | 2.989e-03 | 1.0e-01 | PASS |
| `num.no_horizon_spike` | numerical | `Error near horizon crossing (0.3 < x < 3) is no larger than the late-time plateau by more than 2x` | 1.098e+00 | 2.0e+00 | PASS |
| `num.bogoliubov_prediction` | numerical | `Power ratio of leading-WKB initial data equals the exact \|alpha u + beta u*\|^2 prediction` | 1.902e-11 | 1.0e-08 | PASS |
| `num.beta_scaling` | numerical | `\|beta\| = 1/(2 x_start^2) for leading-WKB initial data` | 1.032e-16 | 1.0e-09 | PASS |
| `num.spectrum_vs_exact` | numerical | `Numerical Delta_t^2(k) at fixed eta vs (2H^2/pi^2M^2)(1+k^2 eta^2)` | 9.921e-12 | 1.0e-08 | PASS |
| `num.late_time_freezeout` | numerical | `Delta_t^2 vs exact (1+x^2) for x_end in [1e-4, 3]` | 2.343e-11 | 1.0e-08 | PASS |
| `num.H_M_scaling` | numerical | `Numerical power scales as H^2/M_pl^2 over a grid of (H, M_pl)` | 2.198e-11 | 1.0e-08 | PASS |
| `num.bogoliubov_sweep` | numerical | `Excited-state power vs prediction for 13 phases of beta` | 5.042e-11 | 1.0e-08 | PASS |
| `num.wronskian_conservation` | numerical | `z^2 (gamma gamma*' - c.c.) = i along the numerical solution` | 6.339e-11 | 1.0e-08 | PASS |
| `num.residual_order` | numerical | `Spline-derivative residual converges at ~4th order in output density` | 1.907e-03 | 5.0e-01 | PASS |
| `num.mpmath_crosscheck` | numerical | `Double-precision DOP853 vs 30-digit Taylor integration (mpmath)` | 2.599e-13 | 1.0e-09 | PASS |
| `num.mpmath_vs_exact` | numerical | `30-digit integration vs exact BD solution` | 8.577e-18 | 1.0e-10 | PASS |
| `num.stability_x_start` | numerical | `Exact-BD data give the same late-time power for x_start in {50..400}` | 3.203e-10 | 1.0e-08 | PASS |
| `num.svt_roundtrip` | numerical | `Spectral SVT decomposition/reconstruction round trip (N = 8, 12, 16)` | 4.857e-17 | 1.0e-13 | PASS |
| `num.tt_transverse` | numerical | `Spectral TT part is transverse: d_j gamma_ij = 0` | 2.602e-17 | 1.0e-13 | PASS |
| `num.tt_traceless` | numerical | `Spectral TT part is traceless` | 4.337e-18 | 1.0e-13 | PASS |
| `num.tt_fd_not_transverse` | negative_control | `Finite-difference divergence of the spectral TT field is O(dx^2), not zero (role of the discretisation)` | 2.543e-03 | 1.0e-13 | PASS |
| `num.small_k_amplification_E` | numerical | `Error in E grows like k^-2 at small k` | 1.969e-01 | 5.0e-01 | PASS |
| `num.small_k_amplification_F` | numerical | `Error in F grows like k^-1 at small k` | 4.409e-02 | 5.0e-01 | PASS |
| `num.gauge_invariants` | numerical | `Bardeen potentials are unchanged by gauge shifts` | 3.553e-15 | 1.0e-12 | PASS |
| `num.gauge_dependents` | negative_control | `A, psi, B, E do change under the same shifts (control)` | 4.000e-02 | 1.0e-06 | PASS |
