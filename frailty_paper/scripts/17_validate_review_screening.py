#!/usr/bin/env python3
"""Read-only structural and decision-consistency checks for review screening files."""
from __future__ import annotations

import argparse
import csv
from pathlib import Path

DECISIONS = {"include", "exclude", "unclear"}
TITLE_FIELDS = {
    "screening_id", "title_abstract_decision_r1", "reviewer_1",
    "title_abstract_decision_r2", "reviewer_2", "adjudication", "adjudicator",
}
FULLTEXT_FIELDS = {
    "screening_id", "fulltext_location", "fulltext_decision_r1", "reviewer_1",
    "fulltext_decision_r2", "reviewer_2", "adjudication", "adjudicator",
}


def read_tsv(path: Path, required: set[str]) -> tuple[list[dict[str, str]], list[str]]:
    try:
        with path.open(encoding="utf-8-sig", newline="") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            missing = required - set(reader.fieldnames or [])
            if missing:
                return [], [f"{path}: missing columns: {', '.join(sorted(missing))}"]
            return list(reader), []
    except (OSError, UnicodeError, csv.Error) as exc:
        return [], [f"{path}: cannot read TSV: {exc}"]


def validate_decisions(rows: list[dict[str, str]], path: Path, prefix: str) -> list[str]:
    errors: list[str] = []
    seen: set[str] = set()
    for line, row in enumerate(rows, start=2):
        sid = row.get("screening_id", "").strip()
        if not sid:
            errors.append(f"{path}:{line}: empty screening_id")
        elif sid in seen:
            errors.append(f"{path}:{line}: duplicate screening_id {sid}")
        seen.add(sid)
        d1, d2 = row.get(f"{prefix}_r1", "").strip().lower(), row.get(f"{prefix}_r2", "").strip().lower()
        adjudication = row.get("adjudication", "").strip().lower()
        r1, r2, adjudicator = (row.get("reviewer_1", "").strip(),
                               row.get("reviewer_2", "").strip(),
                               row.get("adjudicator", "").strip())
        for column, decision in ((f"{prefix}_r1", d1), (f"{prefix}_r2", d2), ("adjudication", adjudication)):
            if decision and decision not in DECISIONS:
                errors.append(f"{path}:{line}: invalid {column} value {decision!r}")
        if bool(d1) != bool(d2):
            errors.append(f"{path}:{line}: both independent decisions are required before resolving {sid}")
        if d1 and (not r1 or not r2):
            errors.append(f"{path}:{line}: reviewer identities are required for {sid}")
        if d1 and d2 and r1 and r1 == r2:
            errors.append(f"{path}:{line}: independent decisions name the same reviewer for {sid}")
        if d1 and d2 and d1 != d2 and not adjudication:
            errors.append(f"{path}:{line}: disagreement for {sid} requires adjudication")
        if adjudication and not adjudicator:
            errors.append(f"{path}:{line}: adjudicator identity is required for {sid}")
    return errors


def validate_review(review_dir: Path) -> tuple[dict[str, int], list[str]]:
    counts: dict[str, int] = {}
    errors: list[str] = []
    title_path = review_dir / "screening/title_abstract_queue.tsv"
    fulltext_path = review_dir / "screening/fulltext_queue.tsv"
    exclusion_path = review_dir / "screening/fulltext_exclusions.tsv"
    title_rows, e = read_tsv(title_path, TITLE_FIELDS); errors.extend(e)
    full_rows, e = read_tsv(fulltext_path, FULLTEXT_FIELDS); errors.extend(e)
    if not e:
        errors.extend(validate_decisions(title_rows, title_path, "title_abstract_decision"))
    if not any(str(fulltext_path) in x for x in errors):
        errors.extend(validate_decisions(full_rows, fulltext_path, "fulltext_decision"))
        for line, row in enumerate(full_rows, start=2):
            if (row.get("fulltext_decision_r1", "").strip() or row.get("fulltext_decision_r2", "").strip()) and not row.get("fulltext_location", "").strip():
                errors.append(f"{fulltext_path}:{line}: full-text location required once full-text decisions are entered")
    title_resolved: dict[str, str] = {}
    for row in title_rows:
        sid = row.get("screening_id", "").strip()
        d1 = row.get("title_abstract_decision_r1", "").strip().lower()
        d2 = row.get("title_abstract_decision_r2", "").strip().lower()
        title_resolved[sid] = row.get("adjudication", "").strip().lower() or (d1 if d1 == d2 else "")
    full_ids = {row.get("screening_id", "").strip() for row in full_rows}
    for row in full_rows:
        sid = row.get("screening_id", "").strip()
        if sid not in title_resolved:
            errors.append(f"{fulltext_path}: full-text record has no title/abstract screening record: {sid}")
        elif title_resolved[sid] != "include":
            errors.append(f"{fulltext_path}: full-text record lacks a resolved title/abstract include decision: {sid}")
    for sid, decision in title_resolved.items():
        if decision == "include" and sid not in full_ids:
            errors.append(f"{title_path}: title/abstract include is missing from the full-text queue: {sid}")
    counts["title_abstract_records"] = len(title_rows)
    counts["title_abstract_decided"] = sum(bool(r.get("title_abstract_decision_r1", "").strip() and r.get("title_abstract_decision_r2", "").strip()) for r in title_rows)
    counts["fulltext_records"] = len(full_rows)
    counts["fulltext_decided"] = sum(bool(r.get("fulltext_decision_r1", "").strip() and r.get("fulltext_decision_r2", "").strip()) for r in full_rows)
    exclusion_rows, e = read_tsv(exclusion_path, {"screening_id", "primary_exclusion_reason", "reviewer"})
    errors.extend(e)
    exclusion_ids: set[str] = set()
    full_by_id = {r.get("screening_id", "").strip(): r for r in full_rows}
    for line, row in enumerate(exclusion_rows, start=2):
        sid = row.get("screening_id", "").strip()
        if not sid or sid in exclusion_ids:
            errors.append(f"{exclusion_path}:{line}: empty or duplicate screening_id {sid!r}")
        exclusion_ids.add(sid)
        if not row.get("primary_exclusion_reason", "").strip():
            errors.append(f"{exclusion_path}:{line}: primary exclusion reason required for {sid}")
        full = full_by_id.get(sid)
        if full is None:
            errors.append(f"{exclusion_path}:{line}: exclusion record has no full-text queue record: {sid}")
        else:
            d1 = full.get("fulltext_decision_r1", "").strip().lower()
            d2 = full.get("fulltext_decision_r2", "").strip().lower()
            resolved = full.get("adjudication", "").strip().lower() or (d1 if d1 == d2 else "")
            if resolved != "exclude":
                errors.append(f"{exclusion_path}:{line}: exclusion record lacks a resolved exclude decision: {sid}")
    counts["fulltext_exclusions"] = len(exclusion_rows)
    return counts, errors


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--review-dir", default="frailty_paper/review")
    args = parser.parse_args()
    counts, errors = validate_review(Path(args.review_dir).resolve())
    for key, value in counts.items():
        print(f"{key}={value}")
    if errors:
        for error in errors:
            print(f"ERROR: {error}")
        print(f"REVIEW_SCREENING_INVALID errors={len(errors)}")
        return 1
    print("REVIEW_SCREENING_OK")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
