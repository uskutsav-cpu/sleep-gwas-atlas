#!/usr/bin/env python3
"""Seal one full pleioFDR MAT result before Snakemake evicts the temporary file."""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import re
from pathlib import Path


PAIR_ID = re.compile(r"[a-z0-9_]+__[a-z0-9_]+$")


def fail(message: str) -> None:
    raise SystemExit(f"ERROR: {message}")


def sha256(path: Path) -> str:
    if not path.is_file() or path.stat().st_size == 0:
        fail(f"missing non-empty checksum-bound artifact: {path}")
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def is_sha256(value: object) -> bool:
    return isinstance(value, str) and re.fullmatch(r"[0-9a-f]{64}", value) is not None


def one_row(path: Path) -> dict[str, str]:
    with path.open(encoding="utf-8", newline="") as handle:
        values = list(csv.DictReader(handle, delimiter="\t"))
    if len(values) != 1:
        fail(f"expected one row in {path}")
    return values[0]


def atomic_json(path: Path, value: object) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(json.dumps(value, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    temporary.replace(path)


def build_receipt(root: Path, pair_id: str) -> tuple[Path, dict[str, object]]:
    task_path = root / f"results/pleiotropy/conjfdr_tasks/{pair_id}.tsv"
    lock_path = root / f"results/pleiotropy/conjfdr_tasks/{pair_id}.lock.tsv"
    completion_path = root / f"results/pleiotropy/conjfdr/{pair_id}/atlas_completion.tsv"
    task, lock, completion = one_row(task_path), one_row(lock_path), one_row(completion_path)
    if (
        lock.get("schema_version") != "sleep-atlas-conjfdr-task.1"
        or lock.get("pair_id") != pair_id
        or lock.get("task_sha256") != sha256(task_path)
        or completion.get("pair_id") != pair_id
        or completion.get("analysis_status") != "CONJFDR_COMPLETE"
        or completion.get("task_sha256") != sha256(task_path)
    ):
        fail(f"conjunction-FDR task/completion identity drifted for {pair_id}")
    result_mat = root / task.get("result_mat", "")
    result_sha256 = sha256(result_mat)
    if completion.get("result_mat_sha256") != result_sha256:
        fail(f"conjunction-FDR MAT differs from completion for {pair_id}")
    receipt_path = result_mat.with_name("result_mat.receipt.json")
    receipt = {
        "schema_version": "sleep-atlas-conjfdr-result-mat.1",
        "analysis_id": task["analysis_id"],
        "pair_id": pair_id,
        "task": str(task_path.relative_to(root)),
        "task_sha256": sha256(task_path),
        "task_lock": str(lock_path.relative_to(root)),
        "task_lock_sha256": sha256(lock_path),
        "completion": str(completion_path.relative_to(root)),
        "completion_sha256": sha256(completion_path),
        "result_mat": str(result_mat.relative_to(root)),
        "result_mat_bytes": result_mat.stat().st_size,
        "result_mat_sha256": result_sha256,
        "sealer_sha256": sha256(Path(__file__)),
        "retention": "TEMPORARY_RESULT_MAT_MAY_BE_EVICTED_AFTER_THIS_RECEIPT",
    }
    return receipt_path, receipt


def validate_receipt(root: Path, pair_id: str, require_mat: bool) -> dict[str, object]:
    receipt_path = root / f"results/pleiotropy/conjfdr/{pair_id}/result_mat.receipt.json"
    try:
        observed = json.loads(receipt_path.read_text(encoding="utf-8"))
    except (OSError, UnicodeError, json.JSONDecodeError) as exc:
        fail(f"conjunction-FDR MAT receipt is unreadable: {exc}")
    if not isinstance(observed, dict):
        fail("conjunction-FDR MAT receipt must be a JSON object")
    expected_path, expected = build_receipt(root, pair_id) if require_mat else (receipt_path, None)
    if expected is not None:
        if expected_path != receipt_path or observed != expected:
            fail(f"conjunction-FDR MAT receipt drifted for {pair_id}")
        return observed
    task_path = root / f"results/pleiotropy/conjfdr_tasks/{pair_id}.tsv"
    lock_path = root / f"results/pleiotropy/conjfdr_tasks/{pair_id}.lock.tsv"
    completion_path = root / f"results/pleiotropy/conjfdr/{pair_id}/atlas_completion.tsv"
    task, lock, completion = one_row(task_path), one_row(lock_path), one_row(completion_path)
    result_mat = root / task.get("result_mat", "")
    if (
        observed.get("schema_version") != "sleep-atlas-conjfdr-result-mat.1"
        or observed.get("pair_id") != pair_id
        or observed.get("analysis_id") != task.get("analysis_id")
        or task.get("pair_id") != pair_id
        or observed.get("task") != str(task_path.relative_to(root))
        or observed.get("task_sha256") != sha256(task_path)
        or observed.get("task_lock") != str(lock_path.relative_to(root))
        or observed.get("task_lock_sha256") != sha256(lock_path)
        or observed.get("completion") != str(completion_path.relative_to(root))
        or observed.get("completion_sha256") != sha256(completion_path)
        or observed.get("result_mat") != str(result_mat.relative_to(root))
        or not is_sha256(observed.get("result_mat_sha256"))
        or observed.get("result_mat_sha256") != completion.get("result_mat_sha256")
        or type(observed.get("result_mat_bytes")) is not int
        or observed["result_mat_bytes"] <= 0
        or observed.get("sealer_sha256") != sha256(Path(__file__))
        or observed.get("retention")
        != "TEMPORARY_RESULT_MAT_MAY_BE_EVICTED_AFTER_THIS_RECEIPT"
        or lock.get("schema_version") != "sleep-atlas-conjfdr-task.1"
        or lock.get("pair_id") != pair_id
        or lock.get("task_sha256") != observed.get("task_sha256")
        or completion.get("pair_id") != pair_id
        or completion.get("analysis_status") != "CONJFDR_COMPLETE"
        or completion.get("task_sha256") != observed.get("task_sha256")
    ):
        fail(f"conjunction-FDR MAT receipt chain drifted for {pair_id}")
    if result_mat.exists() and (
        not result_mat.is_file()
        or result_mat.stat().st_size != observed["result_mat_bytes"]
        or sha256(result_mat) != observed["result_mat_sha256"]
    ):
        fail(f"retained conjunction-FDR MAT drifted for {pair_id}")
    return observed


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("pair_id")
    parser.add_argument("--root", default=".")
    parser.add_argument("--validate-only", action="store_true")
    args = parser.parse_args()
    if not PAIR_ID.fullmatch(args.pair_id):
        fail("invalid pair_id")
    root = Path(args.root).resolve()
    if args.validate_only:
        validate_receipt(root, args.pair_id, require_mat=False)
        print(f"CONJFDR_RESULT_MAT_RECEIPT_OK pair={args.pair_id}")
        return 0
    receipt_path, receipt = build_receipt(root, args.pair_id)
    if receipt_path.exists():
        fail(f"immutable conjunction-FDR MAT receipt already exists: {receipt_path}")
    atomic_json(receipt_path, receipt)
    print(
        f"CONJFDR_RESULT_MAT_SEALED pair={args.pair_id} "
        f"bytes={receipt['result_mat_bytes']}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
