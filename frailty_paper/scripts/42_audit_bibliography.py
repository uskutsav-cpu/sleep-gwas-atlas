#!/usr/bin/env python3
"""Check Pandoc citation keys against a BibTeX file and emit hash-linked JSON."""
from __future__ import annotations

import argparse
import hashlib
import json
import re
import sys
from pathlib import Path

BIB_ENTRY = re.compile(r"@([A-Za-z]+)\s*[({]\s*([^,\s]+)\s*,")
CITATION_CLUSTER = re.compile(r"\[([^\[\]]*?@[^\[\]]*?)\]")
CITATION_KEY = re.compile(r"-?@([A-Za-z0-9_:.+/-]+)")
NON_ENTRIES = {"comment", "preamble", "string"}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def parse_bib_keys(text: str) -> list[str]:
    return [
        key
        for entry_type, key in BIB_ENTRY.findall(text)
        if entry_type.lower() not in NON_ENTRIES
    ]


def parse_citation_keys(text: str) -> list[str]:
    return [
        key
        for cluster in CITATION_CLUSTER.findall(text)
        for key in CITATION_KEY.findall(cluster)
    ]


def audit(manuscript: Path, references: Path, script: Path) -> dict[str, object]:
    manuscript_text = manuscript.read_text(encoding="utf-8")
    references_text = references.read_text(encoding="utf-8")
    bib_keys = parse_bib_keys(references_text)
    citation_keys = parse_citation_keys(manuscript_text)
    unique_bib = set(bib_keys)
    unique_citations = set(citation_keys)
    duplicates = sorted(key for key in unique_bib if bib_keys.count(key) > 1)
    missing = sorted(unique_citations - unique_bib)
    unused = sorted(unique_bib - unique_citations)
    errors = []
    if duplicates:
        errors.append("duplicate bibliography keys")
    if missing:
        errors.append("citations missing bibliography entries")
    if unused:
        errors.append("uncited bibliography entries")
    return {
        "schema_version": "frailty_bibliography_audit.v1",
        "status": "PASS" if not errors else "FAIL",
        "manuscript": str(manuscript),
        "references": str(references),
        "script": str(script),
        "command": "python frailty_paper/scripts/42_audit_bibliography.py --repo <repo>",
        "inputs_sha256": {
            str(manuscript): sha256(manuscript),
            str(references): sha256(references),
        },
        "script_sha256": sha256(script),
        "citation_uses": len(citation_keys),
        "unique_citations": sorted(unique_citations),
        "bibliography_entries": len(bib_keys),
        "bibliography_keys": sorted(unique_bib),
        "duplicate_keys": duplicates,
        "missing_keys": missing,
        "uncited_keys": unused,
        "errors": errors,
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repo", type=Path, default=Path(__file__).resolve().parents[2])
    parser.add_argument("--manuscript", default="frailty_paper/paper/manuscript.md")
    parser.add_argument("--references", default="frailty_paper/paper/references.bib")
    parser.add_argument("--output", help="Optional JSON report path, relative to repository")
    args = parser.parse_args()
    repo = args.repo.resolve()
    manuscript = (repo / args.manuscript).resolve()
    references = (repo / args.references).resolve()
    script = Path(__file__).resolve()
    result = audit(manuscript, references, script)
    rendered = json.dumps(result, indent=2) + "\n"
    if args.output:
        output = (repo / args.output).resolve()
        output.parent.mkdir(parents=True, exist_ok=True)
        output.write_text(rendered, encoding="utf-8")
    sys.stdout.write(rendered)
    return 0 if result["status"] == "PASS" else 1


if __name__ == "__main__":
    raise SystemExit(main())
