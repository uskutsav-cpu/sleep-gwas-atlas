from __future__ import annotations

import importlib.util
import struct
from pathlib import Path

import numpy as np

SCRIPT = Path(__file__).resolve().parents[1] / "audit_placo_ld_factor_normalization_v1.py"
SPEC = importlib.util.spec_from_file_location("placo_ld_factor_normalization", SCRIPT)
module = importlib.util.module_from_spec(SPEC)
assert SPEC and SPEC.loader
SPEC.loader.exec_module(module)


def test_format_14_header_offsets_are_decoded_and_checked(tmp_path: Path) -> None:
    fields = [21775, 19, 2, 3, 4, 14]
    fields.extend([76, 0, 84, 0, 1, 1] + [0] * 7)
    path = tmp_path / "synthetic.bcor"
    path.write_bytes(struct.pack("<19I", *fields) + bytes(28))
    layout = module.bcor_layout(path)
    assert layout["n_snps"] == 2
    assert layout["n_blocks"] == 1
    assert layout["n_entries"] == 1
    assert layout["block_index_offset"] == 76
    assert layout["entry_index_offset"] == 84


def test_bad_format_or_inconsistent_size_fails_closed(tmp_path: Path) -> None:
    fields = [21775, 19, 2, 3, 4, 11]
    fields.extend([76, 0, 84, 0, 1, 1] + [0] * 7)
    path = tmp_path / "wrong-format.bcor"
    path.write_bytes(struct.pack("<19I", *fields) + bytes(28))
    try:
        module.bcor_layout(path)
    except ValueError as error:
        assert "unexpected LAVA BCOR" in str(error)
    else:
        raise AssertionError("unsupported storage format must be rejected")


def test_unit_diagonal_scaling_is_separate_from_raw_factor_product() -> None:
    factors = np.array([[1.0, 0.5], [0.0, 0.5]], dtype=np.float32)
    raw, diagonal_i, diagonal_j, diagnostic = module.factor_correlation(
        factors, np.array([1.0, 1.0], dtype=np.float32), 0, 1
    )
    assert raw == 0.5
    assert diagonal_i == 1.0
    assert diagonal_j == 0.5
    assert diagnostic == 0.5 / np.sqrt(0.5)
    assert raw != diagnostic
