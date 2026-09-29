from __future__ import annotations

import importlib.util
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[3]
SCRIPT = ROOT / "brain6/scripts/retry_lava_roundoff_environmental_locus.py"
spec = importlib.util.spec_from_file_location("lava_roundoff_retry", SCRIPT)
retry = importlib.util.module_from_spec(spec)
assert spec and spec.loader
spec.loader.exec_module(retry)


def test_only_receipt_failure_from_unwritable_temp_dir_is_retryable():
    assert retry.environmental_failure_reason({
        "status": "FAILED", "reason": "No write permission for directory: /tmp/Rtmp",
    })
    assert not retry.environmental_failure_reason({"status": "FAILED", "reason": "LAVA error"})
    assert not retry.environmental_failure_reason({
        "status": "UNIVARIATE_UNDERPOWERED", "reason": "No write permission for directory: /tmp",
    })


@pytest.mark.parametrize("value", ["", "../retry", "has space", "a"*65, ".retry"])
def test_retry_id_cannot_escape_isolated_output_namespace(value):
    with pytest.raises(ValueError, match="retry id"):
        retry.validate_retry_id(value)


def test_retry_id_is_stable_and_path_safe():
    assert retry.validate_retry_id("env-fix-20260924-001") == "env-fix-20260924-001"


def test_failed_retry_attempts_are_read_from_the_shared_retry_attempt_root(tmp_path):
    retry_dir = tmp_path / "targeted_retries" / "env-fix-1"
    attempts = tmp_path / "targeted_retries" / "failed_attempts"
    attempts.mkdir(parents=True)
    wanted = attempts / "740-01-0001"
    other_locus = attempts / "741-01-0002"
    wanted.mkdir()
    other_locus.mkdir()

    assert retry.failed_attempt_paths(retry_dir, "740") == [wanted]
