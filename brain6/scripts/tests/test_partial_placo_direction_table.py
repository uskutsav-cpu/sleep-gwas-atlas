from __future__ import annotations

import importlib.util
from pathlib import Path

SCRIPT = Path(__file__).resolve().parents[1] / "build_partial_placo_direction_table.py"
SPEC = importlib.util.spec_from_file_location("partial_placo_direction_table", SCRIPT)
module = importlib.util.module_from_spec(SPEC)
assert SPEC and SPEC.loader
SPEC.loader.exec_module(module)


def test_same_effect_direction_is_concordant():
    assert module.direction_relation(2.0, 0.5) == "CONCORDANT"
    assert module.direction_relation(-2.0, -0.5) == "CONCORDANT"


def test_opposite_effect_directions_are_opposing():
    assert module.direction_relation(2.0, -0.5) == "OPPOSING"
    assert module.direction_relation(-2.0, 0.5) == "OPPOSING"


def test_zero_z_is_not_called_concordant_or_opposing():
    assert module.direction_relation(0.0, 1.0) == "ZERO_Z_UNDEFINED"
    assert module.direction_relation(1.0, 0.0) == "ZERO_Z_UNDEFINED"
