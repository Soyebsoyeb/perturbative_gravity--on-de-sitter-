"""desitter_perturbations: perturbative gravity on de Sitter spacetime.

Subpackages
-----------
symbolic        SymPy tensor calculus: metric expansion, curvature, linearized
                Einstein equations, Bianchi identities, gauge transformations,
                second-order action.
perturbations   Scalar, vector and tensor sectors, constraints, gauge invariants.
dynamics        Tensor-mode evolution: analytic solutions and numerical solvers.
spectra         Tensor power spectrum, polarization tensors, parameter sweeps.
quantization    Canonical variables, mode normalization, Wronskian, vacuum states.
visualization   Plotting helpers (matplotlib).
validation      Machine-readable verification of every identity and benchmark.
"""

__version__ = "1.0.0"

from . import conventions  # noqa: F401
from . import background  # noqa: F401
