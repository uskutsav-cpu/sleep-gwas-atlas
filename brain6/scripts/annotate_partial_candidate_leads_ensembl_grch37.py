#!/usr/bin/env python3
"""Add descriptive Ensembl GRCh37 gene context to gated partial PLACO leads.

This performs positional annotation only. It does not fine-map signals, test
colocalization, or prioritize causal genes. The source candidate regions must
already have a hash-valid NOT_TIERED pre-annotation gate.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import time
from datetime import datetime, timezone
from pathlib import Path
from typing import Any
from urllib.error import HTTPError, URLError
from urllib.request import Request, urlopen

ROOT = Path(__file__).resolve().parents[2]
GATE_PATH = ROOT / "brain6/results/annotation/partial_candidate_region_gate_v1.tsv"
GATE_PROVENANCE = ROOT / "brain6/results/annotation/partial_candidate_region_gate_v1.provenance.json"
VARIANT_PATH = ROOT / "brain6/results/loci/placo_factor_normalized_sensitivity_v2/placo_candidate_variants.tsv"
SENSITIVITY_PROVENANCE = ROOT / "brain6/results/loci/placo_factor_normalized_sensitivity_v2/provenance.json"
OUT_DIR = ROOT / "brain6/results/annotation"
RAW_PATH = OUT_DIR / "partial_candidate_lead_ensembl_grch37_responses_v1.json"
OUT_PATH = OUT_DIR / "partial_candidate_lead_nearest_gene_grch37_v1.tsv"
PROVENANCE_PATH = OUT_DIR / "partial_candidate_lead_nearest_gene_grch37_v1.provenance.json"
BASE_URL = "https://grch37.rest.ensembl.org/overlap/region/human"
WINDOW_BP = 500_000
FIELDS = [
    "candidate_region_id", "pair_id", "lead_variant", "CHR", "BP",
    "gene_id", "gene_symbol", "biotype", "gene_start", "gene_end",
    "strand", "tss_bp", "distance_to_gene_body_bp", "distance_to_tss_bp",
    "annotation_method", "evidence_tier", "interpretation",
]


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as stream:
        return list(csv.DictReader(stream, delimiter="\t"))


def fetch(url: str, attempts: int = 4) -> tuple[int, bytes]:
    last_error: Exception | None = None
    for attempt in range(attempts):
        try:
            request = Request(url, headers={
                "Accept": "application/json",
                "User-Agent": "Brain6-descriptive-annotation/1.0",
            })
            with urlopen(request, timeout=45) as response:
                return response.status, response.read()
        except (HTTPError, URLError, TimeoutError) as exc:
            last_error = exc
            if isinstance(exc, HTTPError) and exc.code < 500 and exc.code != 429:
                raise RuntimeError(f"Ensembl request failed with HTTP {exc.code}: {url}") from exc
            if attempt + 1 < attempts:
                time.sleep(2 ** attempt)
    raise RuntimeError(f"Ensembl request failed after {attempts} attempts: {url}") from last_error


def build(refresh: bool = False) -> dict[str, Any]:
    gate_prov = json.loads(GATE_PROVENANCE.read_text(encoding="utf-8"))
    if (gate_prov.get("status") != "PASS_NOT_TIERED_PREREQUISITES_BLOCKED" or
            gate_prov.get("scope", {}).get("functional_annotations_used_for_tier_decision") is not False or
            gate_prov.get("output", {}).get("sha256") != sha256(GATE_PATH)):
        raise ValueError("Annotation gate is missing, stale, or not explicitly NOT_TIERED")
    sensitivity = json.loads(SENSITIVITY_PROVENANCE.read_text(encoding="utf-8"))
    if sensitivity.get("status") != "DIAGNOSTIC_SENSITIVITY_NOT_PROMOTED":
        raise ValueError("Partial candidate source is not the unpromoted factor-normalized sensitivity")

    gates = read_tsv(GATE_PATH)
    variants = read_tsv(VARIANT_PATH)
    variant_by_pair_snp: dict[tuple[str, str], dict[str, str]] = {}
    for row in variants:
        variant_by_pair_snp[(row["pair_id"], row["SNP"])] = row

    leads: list[dict[str, str]] = []
    for gate in gates:
        for snp in gate["lead_variants"].split(";"):
            variant = variant_by_pair_snp.get((gate["pair_id"], snp))
            if variant is None or variant.get("candidate_status") != "LEAD":
                raise ValueError(f"Missing source-bound lead coordinate: {gate['pair_id']}/{snp}")
            leads.append({**gate, "lead_variant": snp, "CHR": variant["CHR"], "BP": variant["BP"]})
    keyset = {(r["pair_id"], r["lead_variant"]) for r in leads}
    if len(leads) != 21 or len(keyset) != len(leads):
        raise ValueError("Expected 21 unique PLACO lead signals from the 19 gated regions")

    cached: dict[tuple[str, str], dict[str, Any]] = {}
    if RAW_PATH.is_file() and not refresh:
        old = json.loads(RAW_PATH.read_text(encoding="utf-8"))
        for record in old.get("records", []):
            key = record.get("pair_id"), record.get("lead_variant")
            raw_body = record.get("response_body")
            if (key in cached or not isinstance(raw_body, str) or
                    hashlib.sha256(raw_body.encode("utf-8")).hexdigest() != record.get("response_sha256")):
                raise ValueError("Cached Ensembl response snapshot is malformed or hash-invalid; use --refresh")
            cached[key] = record
        if cached and set(cached) != keyset:
            raise ValueError("Cached Ensembl responses do not match the current lead set; use --refresh")
    use_cache = bool(cached) and not refresh
    retrieved_at = (old.get("retrieved_at_utc") if use_cache else
                    datetime.now(timezone.utc).isoformat(timespec="seconds"))
    raw_records: list[dict[str, Any]] = []
    output: list[dict[str, str]] = []
    for index, lead in enumerate(leads):
        chromosome, position = int(lead["CHR"]), int(lead["BP"])
        start, end = max(1, position - WINDOW_BP), position + WINDOW_BP
        url = f"{BASE_URL}/{chromosome}:{start}-{end}?feature=gene"
        cached_record = cached.get((lead["pair_id"], lead["lead_variant"])) if use_cache else None
        if cached_record is None:
            status, body = fetch(url)
            body_text = body.decode("utf-8")
        else:
            if cached_record.get("query_url") != url:
                raise ValueError(f"Cached Ensembl query URL changed: {lead['lead_variant']}")
            status = int(cached_record["http_status"])
            body_text = cached_record["response_body"]
            body = body_text.encode("utf-8")
        try:
            features = json.loads(body_text)
        except json.JSONDecodeError as exc:
            raise ValueError(f"Non-JSON Ensembl payload for {lead['lead_variant']}") from exc
        if not isinstance(features, list):
            raise ValueError(f"Unexpected Ensembl response for {lead['lead_variant']}")
        if any(item.get("assembly_name") != "GRCh37" for item in features):
            raise ValueError("Ensembl response contains a non-GRCh37 gene feature")
        normalized = []
        for item in features:
            gene_start, gene_end = int(item["start"]), int(item["end"])
            tss = gene_start if int(item["strand"]) == 1 else gene_end
            if position < gene_start:
                body_distance = gene_start - position
            elif position > gene_end:
                body_distance = position - gene_end
            else:
                body_distance = 0
            normalized.append({
                "gene_id": item.get("gene_id", item.get("id", "")),
                "gene_symbol": item.get("external_name", ""),
                "biotype": item.get("biotype", ""),
                "gene_start": gene_start, "gene_end": gene_end,
                "strand": int(item["strand"]), "tss_bp": tss,
                "distance_to_gene_body_bp": body_distance,
                "distance_to_tss_bp": abs(position - tss),
                "assembly_name": item["assembly_name"],
            })
        if not normalized:
            raise ValueError(f"No Ensembl gene features within {WINDOW_BP} bp of {lead['lead_variant']}")
        normalized.sort(key=lambda g: (g["distance_to_gene_body_bp"], g["distance_to_tss_bp"], g["gene_id"]))
        nearest = normalized[0]
        output.append({
            "candidate_region_id": lead["candidate_region_id"],
            "pair_id": lead["pair_id"], "lead_variant": lead["lead_variant"],
            "CHR": str(chromosome), "BP": str(position),
            "gene_id": nearest["gene_id"], "gene_symbol": nearest["gene_symbol"],
            "biotype": nearest["biotype"], "gene_start": str(nearest["gene_start"]),
            "gene_end": str(nearest["gene_end"]), "strand": str(nearest["strand"]),
            "tss_bp": str(nearest["tss_bp"]),
            "distance_to_gene_body_bp": str(nearest["distance_to_gene_body_bp"]),
            "distance_to_tss_bp": str(nearest["distance_to_tss_bp"]),
            "annotation_method": f"Nearest Ensembl gene within +/-{WINDOW_BP:,} bp; GRCh37 archive",
            "evidence_tier": "NOT_TIERED",
            "interpretation": "Positional context only; not a causal-gene assignment or fine-mapping result",
        })
        raw_records.append({
            "pair_id": lead["pair_id"], "candidate_region_id": lead["candidate_region_id"],
            "lead_variant": lead["lead_variant"], "chr": chromosome, "bp": position,
            "query_url": url, "http_status": status,
            "response_sha256": hashlib.sha256(body).hexdigest(),
            "response_body": body_text,
            "feature_count": len(features),
        })
        if index + 1 < len(leads):
            time.sleep(0.1)

    OUT_DIR.mkdir(parents=True, exist_ok=True)
    raw_document = {
        "retrieved_at_utc": retrieved_at,
        "endpoint": "https://grch37.rest.ensembl.org",
        "assembly": "GRCh37",
        "scope": "Nearest-gene positional context for partial PLACO leads; no causal prioritization",
        "records": raw_records,
    }
    RAW_PATH.write_text(json.dumps(raw_document, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    with OUT_PATH.open("w", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=FIELDS, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(output)

    provenance = {
        "schema_version": 1,
        "analysis_id": "brain6_partial_candidate_lead_nearest_gene_ensembl_grch37_v1",
        "status": "PASS_DESCRIPTIVE_POSITIONAL_ANNOTATION_NOT_GENE_PRIORITIZATION",
        "scope": {
            "gated_regions": len(gates), "lead_signals": len(output),
            "unique_nearest_gene_rows": len(output), "assembly": "GRCh37",
            "annotation_window_bp_each_side": WINDOW_BP,
            "fine_mapping_performed": False, "colocalization_performed": False,
            "causal_gene_claims": False, "evidence_tier_changed": False,
            "all_gate_rows_remain_not_tiered": all(r["evidence_tier"] == "NOT_TIERED" for r in output),
        },
        "retrieved_at_utc": retrieved_at,
        "sources": {
            "api": "https://grch37.rest.ensembl.org/overlap/region/human",
            "api_method": "GET feature=gene; nearest feature chosen by distance to gene body, then TSS, then Ensembl ID",
            "cache_used": use_cache,
            "gate": {"path": str(GATE_PATH.relative_to(ROOT)), "sha256": sha256(GATE_PATH)},
            "gate_provenance": {"path": str(GATE_PROVENANCE.relative_to(ROOT)), "sha256": sha256(GATE_PROVENANCE)},
            "candidate_variants": {"path": str(VARIANT_PATH.relative_to(ROOT)), "sha256": sha256(VARIANT_PATH)},
            "candidate_sensitivity_provenance": {
                "path": str(SENSITIVITY_PROVENANCE.relative_to(ROOT)),
                "sha256": sha256(SENSITIVITY_PROVENANCE),
            },
        },
        "raw_response_snapshot": {"path": str(RAW_PATH.relative_to(ROOT)), "sha256": sha256(RAW_PATH)},
        "builder": {"path": str(Path(__file__).resolve().relative_to(ROOT)),
                    "sha256": sha256(Path(__file__).resolve())},
        "output": {"path": str(OUT_PATH.relative_to(ROOT)), "sha256": sha256(OUT_PATH),
                   "rows": len(output)},
    }
    PROVENANCE_PATH.write_text(json.dumps(provenance, indent=2, sort_keys=True) + "\n",
                               encoding="utf-8")
    return provenance


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--refresh", action="store_true",
                        help="fetch updated API responses instead of replaying the hash-verified snapshot")
    args = parser.parse_args()
    print(json.dumps(build(refresh=args.refresh), indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
