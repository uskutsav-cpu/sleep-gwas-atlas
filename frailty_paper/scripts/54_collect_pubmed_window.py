#!/usr/bin/env python3
"""Collect a bounded PubMed date window into a separate resumable cache."""
from __future__ import annotations

import argparse
import csv
import hashlib
import importlib.util
import json
import os
import sys
from datetime import date, datetime, timezone
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
DEFAULT_CONFIG = ROOT / "frailty_paper/config/pubmed_queries.json"
ACTIVE_CACHE = (ROOT / "frailty_paper/review/pubmed").resolve()
COLLECTOR_PATH = ROOT / "frailty_paper/scripts/01_pubmed_search.py"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def parse_date(value: str) -> date:
    try:
        return date.fromisoformat(value)
    except ValueError as exc:
        raise argparse.ArgumentTypeError("dates must use YYYY-MM-DD") from exc


def date_window(start: date, end: date) -> tuple[date, date]:
    if start > end:
        raise ValueError("start date must be on or before end date")
    return start, end


def load_collector():
    spec = importlib.util.spec_from_file_location("frailty_pubmed_collector", COLLECTOR_PATH)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"could not load existing collector: {COLLECTOR_PATH}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def collect(config_path: Path, outdir: Path, start: date, end: date,
            email: str = "", api_key: str = "") -> dict:
    date_window(start, end)
    output = outdir.resolve()
    active_cache = ACTIVE_CACHE.resolve()
    if output == active_cache or active_cache in output.parents:
        raise ValueError("bounded-window collection must use a separate cache, not the active full-history cache")
    config = json.loads(config_path.read_text(encoding="utf-8"))
    queries = config.get("queries")
    if not isinstance(queries, dict) or not queries:
        raise ValueError("query config must contain a nonempty 'queries' mapping")
    output.mkdir(parents=True, exist_ok=True)
    manifest_path = output / "window_manifest.json"
    if manifest_path.exists():
        existing = json.loads(manifest_path.read_text(encoding="utf-8"))
        if (existing.get("start_date"), existing.get("end_date")) != (start.isoformat(), end.isoformat()):
            raise ValueError("output cache already belongs to a different date window")

    collector = load_collector()
    import requests

    session = requests.Session()
    session.headers.update({"User-Agent": "sleep-frailty-systematic-review/1.1"})
    rows = []
    logs = []
    for query_id, base_query in queries.items():
        qdir = output / query_id
        query_path = qdir / "query.txt"
        term = collector.dated_term(base_query, start, end)
        if query_path.exists() and query_path.read_text(encoding="utf-8").strip() != term:
            raise ValueError(f"existing incremental query differs for {query_id}; use a new output directory")
        qdir.mkdir(parents=True, exist_ok=True)
        query_path.write_text(term + "\n", encoding="utf-8")
        segments = list(collector.split_segments(session, base_query, start, end, email, api_key))
        fetched = 0
        for seg_start, seg_end, count in segments:
            fetched += collector.fetch_segment(
                session, query_id, base_query, seg_start, seg_end, qdir,
                email, api_key, rows,
            )
        if fetched != sum(count for _, _, count in segments):
            raise RuntimeError(f"retrieved count differs from ESearch for {query_id}")
        logs.append({
            "query_id": query_id,
            "count": fetched,
            "fetched": fetched,
            "segments": len(segments),
            "start_date": start.isoformat(),
            "end_date": end.isoformat(),
            "query": term,
        })

    with (output / "search_log.tsv").open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(logs[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(logs)
    result = {
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "start_date": start.isoformat(),
        "end_date": end.isoformat(),
        "config_path": str(config_path.resolve()),
        "config_sha256": sha256(config_path),
        "collector_path": str(COLLECTOR_PATH),
        "collector_sha256": sha256(COLLECTOR_PATH),
        "queries": logs,
        "records_returned_across_queries": sum(row["fetched"] for row in logs),
        "integration_status": "SEPARATE_CACHE_NOT_MERGED",
        "command": "frailty_paper/scripts/54_collect_pubmed_window.py",
    }
    manifest_path.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(result, indent=2))
    return result


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", type=Path, default=DEFAULT_CONFIG)
    parser.add_argument("--outdir", type=Path, required=True)
    parser.add_argument("--start-date", type=parse_date, required=True)
    parser.add_argument("--end-date", type=parse_date, required=True)
    parser.add_argument("--email", default=os.environ.get("NCBI_EMAIL", ""))
    parser.add_argument("--api-key", default=os.environ.get("NCBI_API_KEY", ""))
    args = parser.parse_args(argv)
    try:
        collect(args.config, args.outdir, args.start_date, args.end_date, args.email, args.api_key)
    except (ValueError, OSError, RuntimeError) as exc:
        parser.error(str(exc))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
