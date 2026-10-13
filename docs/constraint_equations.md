# Constraint equations

All results below follow from the exactly verified expressions for $\delta E_{\mu\nu}$ in each sector (`perturbations/scalar_sector.py`, `vector_sector.py`, `tensor_sector.py`).

## Scalar sector

In a general gauge with fields depending on $(\eta,z)$, the engine's $\delta E_{\mu\nu}$ satisfies exactly

$$
\delta E_{00}=2\nabla^2\Psi_B-6\mathcal H(\Psi_B'+\mathcal H\Phi_B),\qquad
\delta E_{0i}=2\partial_i(\Psi_B'+\mathcal H\Phi_B),\qquad
\delta E_{ij}\big|_{i\neq j}\propto\partial_i\partial_j(\Phi_B-\Psi_B),
$$

(the last for the traceless part, verified as $\delta E_{xx}-\delta E_{zz}=\partial_z^2(\Phi_B-\Psi_B)$ for dependence on $z$ only) [`sector.scalar.*`]. Every component depends on $(A,B,\psi,E)$ only through $\Phi_B,\Psi_B$, as required by gauge invariance.

Vacuum equations for a Fourier mode with $k\neq0$:

1. traceless $ij$: $k_ik_j(\Phi_B-\Psi_B)=0\Rightarrow\Phi_B=\Psi_B$;
2. $0i$: $\Psi_B'+\mathcal H\Phi_B=0\Rightarrow(a\Psi_B)'=0$;
3. $00$: $-2k^2\Psi_B-6\mathcal H(\Psi_B'+\mathcal H\Phi_B)=0$, which after step 2 is $-2k^2\Psi_B=0$.

Hence $\Phi_B=\Psi_B=0$: the scalar sector has **no** gauge-invariant degree of freedom in vacuum de Sitter space. Steps 1 and 3 are checked symbolically with a plane wave [`constraints.scalar_traceless`, `constraints.scalar_vacuum`]. (In the presence of matter the same equations relate the potentials to sources; none is included here.)

## Vector sector

With $h_{0x}=a^2S_x$, $h_{xz}=a^2\partial_zF_x$ and $V=S-F'$:

$$
\delta E_{0x}=-\tfrac12\partial_z^2V,\qquad \delta E_{xz}=-\tfrac12\partial_z(V'+2\mathcal HV)
$$

[`sector.vector.*`]. For $k\neq0$, $\delta E_{0x}=0$ gives $V=0$.

## Tensor sector

$$
\delta E_{ij}=\tfrac12(\gamma_{ij}''+2\mathcal H\gamma_{ij}'-\nabla^2\gamma_{ij}),\qquad\text{other components}=0,
$$

[`sector.tensor.eom`]: no constraint acts on $\gamma_{ij}$; both polarizations propagate.

## Counting

$10$ metric components $-\;4$ (constraints $\delta E_{0\mu}$, which fix $\Phi_B$-type combinations) $-\;4$ (gauge parameters) $=2$ propagating degrees of freedom, the two helicities $\pm2$. This is the same count as the Hamiltonian analysis (ADM: 12 canonical variables, 4 first-class constraints each removing two phase-space dimensions, leaving 4 phase-space dimensions per point).

## Excluded: $k=0$

The derivation divides by $k^2$. Homogeneous modes ($k=0$) correspond to spatially constant rescalings and shifts, whose status depends on boundary conditions (compact slices, large gauge transformations). They are not claimed to vanish.
