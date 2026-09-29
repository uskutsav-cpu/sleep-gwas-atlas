#!/usr/bin/env python3
"""Fail-closed current-input audit for optional exploratory native methods."""

from __future__ import annotations

import csv
import hashlib
import json
import shutil
import subprocess
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "brain6/results/exploratory_five_track_v3"
LOCKED_INPUTS = ROOT / "brain6/manifests/locked_dense_input_audit.tsv"
SOURCE_CARDS = ROOT / "extensions/brain6/work/overnight-v03/source-cards"
TEMPLATE = ROOT / "extensions/brain6/configs/followup.template.json"
MOLECULAR_TEMPLATE = ROOT / "extensions/brain6/configs/molecular_followup.template.json"
R = Path("/Volumes/Extreme SSD/brain6-work/r-env-4.3.3-short/bin/Rscript")


def sha(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for block in iter(lambda: f.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def read(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as f:
        return list(csv.DictReader(f, delimiter="\t"))


def main() -> None:
    destination = OUT / "method_readiness_20260927.json"
    report = OUT / "method_readiness_20260927.md"
    if destination.exists() or report.exists():
        raise FileExistsError("Versioned method-readiness artifact already exists")
    dense = read(LOCKED_INPUTS)
    if len(dense) != 7 or len({r["trait_id"] for r in dense}) != 7:
        raise ValueError("Expected seven locked dense GWAS inputs")
    verified = []
    for record in dense:
        path = Path(record["dense_path"])
        if not path.is_file() or path.stat().st_size != int(record["dense_bytes"]) or sha(path) != record["expected_sha256"]:
            raise ValueError(f"Dense GWAS changed or missing: {record['trait_id']}")
        card_path = Path(record["source_card"])
        if sha(card_path) != record["source_card_sha256"]:
            raise ValueError(f"Source card changed: {record['trait_id']}")
        card = json.loads(card_path.read_text())
        verified.append({"trait": record["trait_id"], "path": str(path),
                         "sha256": record["expected_sha256"],
                         "n_semantics": card.get("n_semantics"),
                         "effect_scale": card.get("effect_scale"),
                         "genome_build": card.get("genome_build")})
    if not R.is_file():
        raise ValueError("Pinned R runtime missing")
    check = subprocess.run([str(R), "-e",
                            'cat("susieR=",requireNamespace("susieR",quietly=TRUE)," coloc=",requireNamespace("coloc",quietly=TRUE),"\\n",sep="")'],
                           check=True, capture_output=True, text=True)
    native_status = check.stdout.strip().splitlines()[-1]
    if "susieR=" not in native_status or "coloc=" not in native_status:
        raise ValueError("Could not parse native R package inventory")
    followup = json.loads(TEMPLATE.read_text())
    molecular = json.loads(MOLECULAR_TEMPLATE.read_text())
    reviewed = followup.get("reviewed") is True and molecular.get("reviewed") is True
    if reviewed:
        raise ValueError("Templates unexpectedly claim reviewed method settings")
    qtl_roots = [ROOT / "brain6/data/qtl", ROOT / "data/qtl",
                 Path("/Volumes/Extreme SSD/brain6-work/qtl")]
    available_qtl = [str(p) for p in qtl_roots if p.is_dir() and any(p.iterdir())]
    rows = [
        {"method": "five_track_candidate_reconciliation", "class": "A", "status": "COMPLETE_EXPLORATORY",
         "reason": "All 25 pair candidates and 20 geographic groups are receipt bound."},
        {"method": "independent_ld_and_signed_direction", "class": "A/B", "status": "COMPLETE_EXPLORATORY_WITH_EXCEPTIONS",
         "reason": "Five cross-pair genotype LD edges and 2,686 signed-LD candidate statuses are versioned; six raw R values fail range QC."},
        {"method": "locus_specific_replication", "class": "A", "status": "NOT_ESTABLISHED",
         "reason": "Existing FinnGen ADHD result is global pair-level directional replication only."},
        {"method": "susie_fine_mapping", "class": "B", "status": "NOT_RUN_METHOD_PREREQUISITES_UNREVIEWED",
         "reason": "Pinned R lacks susieR; reviewed trait-specific N/case-fraction and full regional LD contract absent; UKB raw LD has documented out-of-range entries."},
        {"method": "trait_trait_coloc", "class": "B", "status": "NOT_RUN_METHOD_PREREQUISITES_UNREVIEWED",
         "reason": "Pinned R lacks coloc; reviewed regional dataset, sample-size, prior and overlap settings absent; multi-signal path also requires valid fine-mapping."},
        {"method": "eqtl_sqtl_coloc", "class": "B", "status": "NOT_RUN_NO_VERIFIED_QTL_INPUT",
         "reason": "No source-verified full-cis QTL dataset or reviewed molecular follow-up configuration is locally bound."},
        {"method": "regulatory_cell_type_pathway", "class": "B", "status": "DESCRIPTIVE_COORDINATE_CONTEXT_ONLY",
         "reason": "B Ensembl lead-level regulatory lookup was cached; no weighted functional input, tissue panel/background or multiplicity family is reviewed."},
        {"method": "lava_local_rg_and_final_tiers", "class": "C", "status": "BLOCKED_LAVA",
         "reason": "Canonical v3 failed frozen family QC: 3,720 NOT_RUN versus maximum 873."},
    ]
    payload = {
        "analysis_id": "brain6_exploratory_method_readiness_v1",
        "analysis_label": "EXPLORATORY", "dense_inputs_verified": verified,
        "native_R_inventory": native_status,
        "plink_executable": shutil.which("plink"),
        "plink2_executable": shutil.which("plink2"),
        "qtl_directories_with_files": available_qtl,
        "methods": rows,
        "source_sha256": {
            str(LOCKED_INPUTS.relative_to(ROOT)): sha(LOCKED_INPUTS),
            str(TEMPLATE.relative_to(ROOT)): sha(TEMPLATE),
            str(MOLECULAR_TEMPLATE.relative_to(ROOT)): sha(MOLECULAR_TEMPLATE),
            "pinned_Rscript": sha(R),
            **{f"source_card_{r['trait']}": sha(SOURCE_CARDS / (r["trait"] + ".source.json")) for r in verified},
        },
        "script_sha256": sha(Path(__file__)),
    }
    destination.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    lines = ["# Exploratory downstream method readiness — 2026-09-27", "",
             "All seven mounted dense GWAS files match the locked byte sizes and SHA-256 values.",
             "The source cards retain different `N` semantics (`total` or `effective`).",
             "This does not itself validate a trait-specific SuSiE/coloc sample-size model.",
             "", "| Method | Class | Current state | Evidence/limitation |",
             "|---|:---:|---|---|"]
    for item in rows:
        lines.append(f"| {item['method']} | {item['class']} | {item['status']} | {item['reason']} |")
    lines.extend(["", f"Pinned native runtime: `{native_status}`; PLINK: `{shutil.which('plink')}`; PLINK2: `{shutil.which('plink2')}`.",
                  f"Locally bound QTL directories with files: {len(available_qtl)}.",
                  "Missing prerequisites are method blocks, not evidence of a negative biological result.",
                  "No protected LAVA decision or frozen promotion rule was modified.", ""])
    report.write_text("\n".join(lines))
    receipt = {"analysis_id": payload["analysis_id"], "analysis_label": "EXPLORATORY",
               "source_sha256": {str(destination.relative_to(ROOT)): sha(destination)},
               "script_sha256": sha(Path(__file__)),
               "output_sha256": {report.name: sha(report)}}
    (OUT / "method_readiness_20260927.provenance.json").write_text(
        json.dumps(receipt, indent=2, sort_keys=True) + "\n")
    print(json.dumps({"dense_inputs_verified": len(verified),
                      "native_R_inventory": native_status, "qtl_sources": len(available_qtl)},
                     sort_keys=True))


if __name__ == "__main__":
    main()
