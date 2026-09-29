from __future__ import annotations

import importlib.util
from pathlib import Path

SCRIPT = Path(__file__).resolve().parents[1] / "audit_lava_pair_context_univariates.py"
SPEC = importlib.util.spec_from_file_location("pair_context_audit", SCRIPT)
assert SPEC is not None and SPEC.loader is not None
MODULE = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(MODULE)


def test_pair_context_conflict_and_canonical_difference_are_separate():
    contexts = [{"p": "0.0742928"}, {"p": "0.00398049"}]
    canonical = {"p": "0.0878159", "status": "TESTED"}

    result = MODULE.summarize_context(contexts, canonical)

    assert result["n_contexts"] == 2
    assert result["contexts_disagree"] is True
    assert abs(result["within_context_p_range"] - (0.0742928 - 0.00398049)) < 1e-15
    assert result["differs_from_canonical"] is True
    assert abs(result["max_abs_context_minus_canonical"] - (0.0878159 - 0.00398049)) < 1e-15


def test_pair_context_agreement_with_missing_canonical_p_is_not_canonical_validation():
    contexts = [{"p": "0.2"}, {"p": "0.2"}]
    result = MODULE.summarize_context(contexts, {"p": "NA", "status": "NOT_RUN"})

    assert result["contexts_disagree"] is False
    assert result["canonical_p"] is None
    assert result["differs_from_canonical"] is False
    assert result["max_abs_context_minus_canonical"] == 0.0
