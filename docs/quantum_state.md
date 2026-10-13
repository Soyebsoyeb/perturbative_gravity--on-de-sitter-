# Quantum State, Wronskian and Commutation Relations

Status: free-field canonical quantization of the tensor sector. Standard and established; here implemented and verified.

## Mode expansion

$$
\gamma_\lambda(\eta,\mathbf x)=\int\frac{d^3k}{(2\pi)^3}\Big[\gamma_k(\eta)\,a_{\mathbf k}+\gamma_k^*(\eta)\,a_{-\mathbf k}^\dagger\Big]e^{i\mathbf k\cdot\mathbf x},
\qquad
[a_{\mathbf k},a_{\mathbf k'}^\dagger]=(2\pi)^3\delta^{3}(\mathbf k-\mathbf k').
$$

## Wronskian normalization

The canonical commutator $[\gamma(\mathbf x),\pi(\mathbf x')]=i\,\delta^3(\mathbf x-\mathbf x')$ holds if and only if

$$
z_T^2\big(\gamma_k\gamma_k^{*\prime}-\gamma_k^*\gamma_k'\big)=i
$$

for all $\eta$. Because the mode equation is real and second order, the left-hand side is conserved; the code checks that its drift is $6.3\times10^{-11}$ over the whole integration range (`wronskian.py`).

## Bunch-Davies mode

The exact solution with positive-frequency behaviour $v\to e^{-ik\eta}/\sqrt{2k}$ as $\eta\to-\infty$ is

$$
v_k=\frac{e^{-ik\eta}}{\sqrt{2k}}\Big(1-\frac{i}{k\eta}\Big),\qquad \gamma_k=\frac{v_k}{z_T}.
$$

The code verifies that it solves the mode equation exactly, is correctly normalized, and matches the asymptotic form.

## Other vacua

A Bogoliubov rotation $v=\alpha\,v_{\rm BD}+\beta\,v_{\rm BD}^*$ with $\lvert\alpha\rvert^2-\lvert\beta\rvert^2=1$ gives another Wronskian-normalized state. `vacuum_states.bogoliubov` implements it, and the late-time power is changed by the factor $\lvert\alpha-\beta\rvert^2$, which is compared with the numerical evolution to $1.9\times10^{-11}$.

## Lattice commutator test

On a periodic lattice the equal-time commutator $[\gamma_{\mathbf k},\pi_{\mathbf k'}^\dagger]$ is evaluated from the mode functions (`commutation_relations.lattice_commutator`). The zero mode has no oscillator interpretation and is normalized through the Wronskian with coefficient $1$, as documented in the module.

## What is not claimed

Only the free graviton two-point function is treated. No interacting theory, loop correction or non-linear gauge fixing is attempted.
