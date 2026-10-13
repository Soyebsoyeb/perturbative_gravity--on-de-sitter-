"""Verification record type shared by all checks."""

from __future__ import annotations

from dataclasses import dataclass, asdict, field
from typing import Optional


@dataclass
class Record:
    id: str
    category: str
    description: str
    residual: float
    tolerance: float
    passed: bool
    method: str = ""
    assumptions: str = ""
    kind: str = "identity"      # 'identity' (exact), 'benchmark' (analytic), 'numerical', 'negative_control'
    relative_error: Optional[float] = None
    absolute_error: Optional[float] = None
    notes: str = ""

    def as_dict(self):
        return asdict(self)


def record(id, category, description, residual, tolerance, method="", assumptions="", kind="identity",
           relative_error=None, absolute_error=None, notes="", expect_fail=False):
    """Build a record.  ``expect_fail`` marks a negative control: the check passes only if
    the residual EXCEEDS the tolerance (the deliberately wrong formula is detected)."""
    residual = float(residual)
    ok = residual <= tolerance
    if expect_fail:
        ok = residual > tolerance
        kind = "negative_control"
    return Record(id=id, category=category, description=description, residual=residual,
                  tolerance=float(tolerance), passed=bool(ok), method=method, assumptions=assumptions,
                  kind=kind, relative_error=relative_error, absolute_error=absolute_error, notes=notes)


def zero_residual(obj, simplify=None):
    """0.0 if every sympy entry inside ``obj`` (expression, Matrix, list, dict, nested) is
    identically zero after simplification, else 1.0."""
    import sympy as sp
    simplify = simplify or sp.simplify

    def flat(o):
        if isinstance(o, sp.MatrixBase):
            for e in o:
                yield e
        elif isinstance(o, dict):
            for e in o.values():
                yield from flat(e)
        elif isinstance(o, (list, tuple)):
            for e in o:
                yield from flat(e)
        else:
            yield o

    return 0.0 if all(simplify(sp.sympify(e)) == 0 for e in flat(obj)) else 1.0
