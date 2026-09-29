import pytest

from validate_power_optimized_lava_overlap_v1 import check_matrix, log_intercept


def test_check_matrix_accepts_symmetric_positive_definite_overlap(tmp_path):
    matrix = tmp_path / "overlap.txt"
    matrix.write_text("sleep disorder\nsleep 1.0362 -0.0075\ndisorder -0.0075 1.0793\n")
    pair = {
        "trait1": "sleep",
        "trait2": "disorder",
        "sleep_intercept": 1.0362,
        "cross_trait_intercept": -0.0075,
        "disorder_intercept": 1.0793,
    }
    check_matrix(matrix, pair)


def test_check_matrix_rejects_non_positive_definite_overlap(tmp_path):
    matrix = tmp_path / "overlap.txt"
    matrix.write_text("sleep disorder\nsleep 1 -2\ndisorder -2 1\n")
    pair = {
        "trait1": "sleep",
        "trait2": "disorder",
        "sleep_intercept": 1,
        "cross_trait_intercept": -2,
        "disorder_intercept": 1,
    }
    with pytest.raises(ValueError, match="positive definite"):
        check_matrix(matrix, pair)


def test_log_intercept_reads_ldsc3_genetic_covariance_section(tmp_path):
    log = tmp_path / "rg.log"
    log.write_text(
        "Genetic Covariance\n------------------\n"
        "Mean Chi^2: 1.2\nIntercept: -0.0075 (0.0073)\n"
        "Genetic Correlation\n--------------------\n"
    )
    assert log_intercept(log, "Cross trait intercept") == (-0.0075, 0.0073)
