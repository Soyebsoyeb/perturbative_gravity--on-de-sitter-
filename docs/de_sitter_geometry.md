# De Sitter geometry

## Time coordinates

Flat slicing in cosmic time has $a(t)=e^{Ht}$. Conformal time is defined by $d\eta=dt/a(t)$, so

$$
\eta=\int e^{-Ht}dt=-\frac{e^{-Ht}}{H},
$$

where the integration constant is chosen so that $\eta\to0^-$ as $t\to\infty$. Inverting, $e^{Ht}=-1/(H\eta)$, i.e. $a(\eta)=-1/(H\eta)$ with $\eta<0$. The conformal Hubble parameter is

$$
\mathcal H=\frac{a'}{a}=-\frac1\eta=aH,\qquad \mathcal H'=\frac{1}{\eta^2}=\mathcal H^2,\qquad \frac{a''}{a}=\frac{2}{\eta^2}.
$$

The comoving Hubble radius is $1/(aH)=-\eta$; a comoving wavenumber $k$ is inside the horizon when $x\equiv k/(aH)=-k\eta>1$ and outside when $x<1$. At $x=1$ the physical wavelength is $2\pi/H$.

Verified: [`bg.time_relation`], [`bench.scale_factor_relation`], [`bench.deta_dt`], and `tests/test_background.py`.

## Curvature

Both metrics are plugged into the same engine (`symbolic.curvature.Geometry`), which computes the Christoffel symbols, the full Riemann tensor (all 256 components), the Ricci tensor and scalar, and the curvature invariants, with the sign conventions of `docs/conventions.md`. The results, identical in the two slicings, are

| quantity | value | record |
|---|---|---|
| $R_{\mu\nu\rho\sigma}$ | $H^2(g_{\mu\rho}g_{\nu\sigma}-g_{\mu\sigma}g_{\nu\rho})$ | `bg.constant_curvature.*` |
| $R_{\mu\nu}$ | $3H^2g_{\mu\nu}$ | `bg.ricci_tensor.*` |
| $R$ | $12H^2$ | `bg.ricci_scalar.*` |
| $G_{\mu\nu}$ | $-3H^2g_{\mu\nu}$ | `bg.einstein.*` |
| $R_{\mu\nu}R^{\mu\nu}$ | $36H^4$ | `bg.ricci_squared.*` |
| $R_{\mu\nu\rho\sigma}R^{\mu\nu\rho\sigma}$ | $24H^4$ | `bg.kretschmann.*` |

$G_{\mu\nu}+\Lambda g_{\mu\nu}=0$ with $\Lambda=3H^2$ [`bg.field_eq.*`]. The two slicings give equal scalars [`bg.slicing_invariants`]. For constant curvature, $R_{\mu\nu}=(n-1)H^2g_{\mu\nu}$ in $n=4$ dimensions and the Kretschmann scalar is $2n(n-1)H^4=24H^4$, consistent with the table.

## Embedding

Flat de Sitter space is the portion of the hyperboloid $-X_0^2+\mathbf X^2+X_4^2=1/H^2$ with

$$
X_0=\frac{\sinh Ht}{H}+\frac H2e^{Ht}r^2,\quad X_i=e^{Ht}x_i,\quad X_4=\frac{\cosh Ht}{H}-\frac H2e^{Ht}r^2 .
$$

The pull-back of the five-dimensional Minkowski metric is $-dt^2+e^{2Ht}d\mathbf x^2$ and the hyperboloid relation holds identically [`bg.embedding`]. This is an independent confirmation that the space has constant positive curvature $H^2$, since a hyperboloid of radius $1/H$ in Minkowski space does.

## What the flat slicing does not cover

Flat slicing covers half of the global manifold ($\eta<0$ only). Statements about the full de Sitter space (global modes, the Euclidean sphere, state definitions beyond the Bunch-Davies prescription) are outside this page.
