import sys
from pathlib import Path

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from build_placo_factor_normalized_sensitivity_v1 import (
    clip_roundoff_correlation_boundary,
    edge_chromosome,
    normalized_correlation,
)


def test_factor_gram_normalization_produces_unit_diagonal_psd_correlations():
    factors = np.array([[2.0, 0.0, 1.0], [0.0, 1.0, 1.0]])
    eigenvalues = np.array([3.0, 0.5])
    corr = normalized_correlation(factors, eigenvalues)

    assert np.allclose(np.diag(corr), 1.0, atol=1e-12)
    assert np.allclose(corr, corr.T, atol=1e-12)
    assert np.max(np.abs(corr)) <= 1.0 + 1e-12
    assert np.linalg.eigvalsh(corr)[0] >= -1e-12
    # Normalization differs from the unscaled factor product when its diagonal
    # is not one, which is precisely the diagnostic under examination.
    assert not np.isclose(corr[0, 2], float(factors[:, 0] @ (eigenvalues * factors[:, 2])))


def test_factor_gram_normalization_preserves_signed_perfect_ld():
    corr = normalized_correlation(np.array([[1.0, -1.0]]), np.array([2.0]))
    assert np.isclose(corr[0, 1], -1.0, atol=1e-12)


def test_factor_gram_normalization_resolves_only_floating_point_boundary_excess():
    rng = np.random.default_rng(2)
    vector = rng.normal(size=148).astype(np.float32)
    factors = np.column_stack([vector, vector])
    eigenvalues = np.exp(rng.normal(size=148)).astype(np.float32)
    corr, corrected = normalized_correlation(factors, eigenvalues, report_roundoff=True)

    # Whether this exact input lands one ULP outside the bound depends on the
    # NumPy/BLAS build; either outcome must remain a unit correlation to
    # floating-point precision and never exceed the mathematical boundary.
    assert corrected in {0, 2}
    assert np.array_equal(np.diag(corr), np.ones(2))
    assert np.isclose(corr[0, 1], 1.0, rtol=0.0, atol=2e-16)
    assert np.max(np.abs(corr)) <= 1.0
    assert np.linalg.eigvalsh(corr)[0] >= -1e-12


def test_roundoff_boundary_clipping_is_deterministic_and_fails_closed():
    above = np.nextafter(1.0, np.inf)
    below = np.nextafter(-1.0, -np.inf)
    positive, positive_count = clip_roundoff_correlation_boundary(
        np.array([[1.0, above], [above, 1.0]])
    )
    negative, negative_count = clip_roundoff_correlation_boundary(
        np.array([[1.0, below], [below, 1.0]])
    )
    assert positive_count == negative_count == 2
    assert positive[0, 1] == 1.0 and negative[0, 1] == -1.0
    assert np.linalg.eigvalsh(positive)[0] >= -1e-12
    assert np.linalg.eigvalsh(negative)[0] >= -1e-12
    with pytest.raises(ValueError, match="materially exceeds"):
        clip_roundoff_correlation_boundary(np.array([[1.0, 1.0 + 1e-6],
                                                     [1.0 + 1e-6, 1.0]]))


def test_factor_gram_normalization_rejects_negative_eigenvalues():
    with pytest.raises(ValueError, match="negative eigenvalues"):
        normalized_correlation(np.eye(2), np.array([1.0, -0.1]))


def test_factor_gram_normalization_rejects_bad_dimensions():
    with pytest.raises(ValueError, match="dimensions disagree"):
        normalized_correlation(np.eye(2), np.array([1.0]))


def test_edge_chromosome_uses_precomputed_shared_snp_index():
    index = {"rs1": {1}, "rs2": {1}, "rs3": {2}}
    assert edge_chromosome("rs1", "rs2", index) == 1
    assert edge_chromosome("rs1", "rs3", index) is None


def test_edge_chromosome_rejects_ambiguous_snp_ids():
    index = {"rs1": {1, 2}, "rs2": {1, 2}}
    with pytest.raises(ValueError, match="ambiguous chromosome"):
        edge_chromosome("rs1", "rs2", index)
