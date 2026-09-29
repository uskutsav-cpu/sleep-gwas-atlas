from __future__ import annotations

import importlib.util
from pathlib import Path


SCRIPT = Path(__file__).resolve().parents[1] / "validate_partial_placo_pvalues.py"
SPEC = importlib.util.spec_from_file_location("validate_partial_placo_pvalues", SCRIPT)
MODULE = importlib.util.module_from_spec(SPEC)
assert SPEC and SPEC.loader
SPEC.loader.exec_module(MODULE)


def test_exact_values_match():
    assert MODULE.rel_error(1e-10, 1e-10) == 0
    assert MODULE.is_match(1e-10, 1e-10)


def test_relative_tolerance_is_scale_aware():
    assert MODULE.is_match(1e-10, 1.000000005e-10)
    assert not MODULE.is_match(1e-10, 1.0001e-10)


def test_nonfinite_values_never_pass():
    assert not MODULE.is_match(float("nan"), 1e-10)
    assert not MODULE.is_match(1e-10, float("inf"))
