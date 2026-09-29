from __future__ import annotations

import importlib.util
from pathlib import Path

import numpy as np

SCRIPT = Path(__file__).resolve().parents[1] / "audit_placo_ld_primary_block_psd_v1.py"
SPEC = importlib.util.spec_from_file_location("placo_ld_primary_block_psd", SCRIPT)
module = importlib.util.module_from_spec(SPEC)
assert SPEC and SPEC.loader
SPEC.loader.exec_module(module)


def test_unit_diagonal_replacement_can_break_psd_but_normalized_gram_is_psd():
    factors = np.array([[1.01, 1.01]], dtype=np.float32)
    gram, lava_style, normalized = module.block_matrices(
        factors, np.array([1.0], dtype=np.float32)
    )
    assert np.min(np.linalg.eigvalsh(gram)) >= -1e-12
    assert np.min(np.linalg.eigvalsh(lava_style)) < 0
    assert np.min(np.linalg.eigvalsh(normalized)) >= -1e-12
    assert np.allclose(np.diag(normalized), 1.0)


def test_diagonal_rescaling_does_not_mutate_input_factors():
    factors = np.array([[1.0, 0.5], [0.0, 0.5]], dtype=np.float32)
    original = factors.copy()
    module.block_matrices(factors, np.array([2.0, 1.0], dtype=np.float32))
    assert np.array_equal(factors, original)
