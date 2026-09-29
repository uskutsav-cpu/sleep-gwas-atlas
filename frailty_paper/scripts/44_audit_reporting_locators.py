#!/usr/bin/env python3
"""Check that reporting-checklist source locators resolve to current files/lines."""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import re
from pathlib import Path


TABLES = (
    "prisma.tsv",
    "prisma_abstract.tsv",
    "prisma_s.tsv",
    "strega.tsv",
)
LOCATOR_RE = re.compile(
    r"(?P<path>(?:frailty_paper/)?(?:paper/)?(?:review/)?[A-Za-z0-9_./-]+\.md)"
    r":L(?P<start>[1-9][0-9]*)(?:-L(?P<end>[1-9][0-9]*))?"
)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def resolve_source(repo: Path, raw_path: str) -> Path:
    if raw_path == "manuscript.md":
        return repo / "frailty_paper/paper/manuscript.md"
    if raw_path.startswith("review/"):
        return repo / "frailty_paper" / raw_path
    return repo / raw_path


def audit(repo: Path) -> dict:
    repo = repo.resolve()
    check_dir = repo / "frailty_paper/paper/reporting_checklists"
    errors: list[str] = []
    table_summaries = []
    source_paths: set[Path] = set()
    total_rows = 0
    total_locators = 0

    for filename in TABLES:
        table_path = check_dir / filename
        if not table_path.is_file():
            errors.append(f"missing checklist table: {table_path.relative_to(repo)}")
            continue
        source_paths.add(table_path)
        with table_path.open(newline="", encoding="utf-8") as stream:
            reader = csv.DictReader(stream, delimiter="\t")
            if not reader.fieldnames or "manuscript_page_line" not in reader.fieldnames:
                errors.append(f"missing manuscript_page_line column: {table_path.relative_to(repo)}")
                continue
            rows = list(reader)
        table_locator_count = 0
        total_rows += len(rows)
        for row_number, row in enumerate(rows, start=2):
            locator = (row.get("manuscript_page_line") or "").strip()
            item = row.get("item_id", "")
            if not locator:
                errors.append(f"{filename}:{row_number} item {item}: empty locator")
                continue
            matches = list(LOCATOR_RE.finditer(locator))
            if not matches:
                if not any(token in locator.upper() for token in ("NOT PRESENT", "NOT_DRAFTED", "NOT DRAFTED")):
                    errors.append(
                        f"{filename}:{row_number} item {item}: locator has no line reference or explicit absence status"
                    )
                continue
            for match in matches:
                raw_path = match.group("path")
                source = resolve_source(repo, raw_path).resolve()
                start = int(match.group("start"))
                end = int(match.group("end") or start)
                try:
                    source.relative_to(repo)
                except ValueError:
                    errors.append(f"{filename}:{row_number} item {item}: source escapes repository: {raw_path}")
                    continue
                if not source.is_file():
                    errors.append(f"{filename}:{row_number} item {item}: missing source file: {raw_path}")
                    continue
                source_paths.add(source)
                line_count = len(source.read_text(encoding="utf-8").splitlines())
                if end < start or end > line_count:
                    errors.append(
                        f"{filename}:{row_number} item {item}: {raw_path}:L{start}-L{end} "
                        f"outside current {line_count}-line source"
                    )
                    continue
                table_locator_count += 1
                total_locators += 1
        table_summaries.append({"file": table_path.relative_to(repo).as_posix(), "rows": len(rows), "line_locators": table_locator_count})

    input_hashes = {
        path.relative_to(repo).as_posix(): sha256(path)
        for path in sorted(source_paths)
        if path.is_file()
    }
    return {
        "schema_version": "reporting_checklist_locator_audit.v1",
        "repo": str(repo),
        "result": "PASS" if not errors else "FAIL",
        "checklist_rows": total_rows,
        "line_locators_checked": total_locators,
        "tables": table_summaries,
        "inputs_sha256": input_hashes,
        "errors": errors,
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repo", type=Path, default=Path(__file__).resolve().parents[2])
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    result = audit(args.repo)
    output = args.output or args.repo / "frailty_paper/analysis/reporting_checklist_locator_audit_2026-09-23.json"
    output = output.resolve()
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(
        f"REPORTING_LOCATORS_{result['result']} rows={result['checklist_rows']} "
        f"line_references={result['line_locators_checked']} errors={len(result['errors'])} output={output}"
    )
    if result["errors"]:
        for error in result["errors"]:
            print(f"ERROR: {error}")
        raise SystemExit(1)


if __name__ == "__main__":
    main()
