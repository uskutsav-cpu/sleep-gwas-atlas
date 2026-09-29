import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from build_partial_placo_signed_effect_ld import allele_sign, ld_direction_consistency


def test_allele_sign_orients_same_swapped_and_complemented_effects():
    assert allele_sign("A", "C", "A", "C") == 1
    assert allele_sign("A", "C", "C", "A") == -1
    assert allele_sign("A", "C", "T", "G") == 1
    assert allele_sign("A", "C", "G", "T") == -1
    assert allele_sign("A", "C", "A", "G") is None


def test_ld_sign_consistency_is_descriptive_and_fails_closed():
    assert ld_direction_consistency(0.2, 0.4, 0.7) == "CONSISTENT_WITH_LD_SIGN"
    assert ld_direction_consistency(0.2, -0.4, 0.7) == "OPPOSITE_TO_LD_SIGN"
    assert ld_direction_consistency(0.2, 0.4, 1.0001) == "NOT_ASSESSED_INVALID_LD_OR_EFFECT"
    assert ld_direction_consistency(0.0, 0.4, 0.7) == "NOT_ASSESSED_ZERO_EFFECT_OR_LD"
