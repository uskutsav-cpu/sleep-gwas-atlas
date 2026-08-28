#!/usr/bin/env python3
"""Render a pre-result adversarial review without pretending blocked findings exist."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import subprocess
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--gates", type=Path, default=Path("discovery_extension/results/extension_acceptance_gates.tsv"))
    parser.add_argument("--checklist", type=Path, default=Path("discovery_extension/results/adversarial_review_checklist.tsv"))
    parser.add_argument("--preflight", type=Path, default=Path("discovery_extension/provenance/acquisition_preflight.json"))
    parser.add_argument("--panel-lock", type=Path, default=Path("discovery_extension/config/extension_panel.lock.json"))
    parser.add_argument("--out", type=Path, default=Path("discovery_extension/adversarial_review.md"))
    parser.add_argument("--provenance-out", type=Path, default=Path("discovery_extension/provenance/adversarial_review.json"))
    args = parser.parse_args()

    core = subprocess.run(
        ["python3", "discovery_extension/scripts/00_verify_core_checkpoint.py"],
        check=False, text=True, capture_output=True,
    )
    core_detail = (core.stdout + core.stderr).strip()
    gates = read_tsv(args.gates)
    checklist = read_tsv(args.checklist)
    preflight = json.loads(args.preflight.read_text(encoding="utf-8"))
    panel_lock = json.loads(args.panel_lock.read_text(encoding="utf-8"))
    if len(gates) != 17 or not checklist:
        raise SystemExit("ERROR: acceptance gates/checklist must be generated before adversarial review")
    if core.returncode != 0:
        verdict = "BLOCKING CORE CHECKPOINT FAILURE"
    elif preflight["status"] != "READY":
        verdict = "PRE-RESULT DESIGN REVIEW COMPLETE; FINDINGS-LEVEL REVIEW BLOCKED BY ACQUISITION"
    else:
        verdict = "PRE-RESULT DESIGN REVIEW COMPLETE; FINDINGS-LEVEL REVIEW PENDING ANALYSIS"

    findings = {
        "extension h2": Path("discovery_extension/results/ldsc/extension_trait_readiness.tsv"),
        "global rg": Path("discovery_extension/results/ldsc/extension_rg_matrix.tsv"),
        "pair novelty": Path("discovery_extension/results/novelty/extension_novelty_audit.tsv"),
        "replication": Path("discovery_extension/results/replication/replication_results.tsv"),
        "local architecture": Path("discovery_extension/results/local/local_rg_results.tsv"),
        "pleiotropic loci": Path("discovery_extension/results/pleiotropy/novel_shared_loci.tsv"),
        "fine-mapping/colocalization": Path("discovery_extension/results/fine_mapping/fine_mapping_colocalization.tsv"),
        "mechanistic synthesis": Path("discovery_extension/results/mechanism/mechanistic_synthesis.tsv"),
    }
    absent = [name for name, path in findings.items() if not path.is_file()]
    present = [name for name, path in findings.items() if path.is_file()]
    gate_status_counts = Counter(row["status"] for row in gates)
    risk_status_counts = Counter(row["current_status"] for row in checklist)
    lines = [
        "# Discovery extension adversarial review",
        "",
        f"Generated: {datetime.now(timezone.utc).strftime('%Y-%m-%dT%H:%M:%SZ')}",
        "",
        f"Verdict: **{verdict}**",
        "",
        "## Scope and non-result boundary",
        "",
        f"The immutable core verifier reports: `{core_detail}`",
        "",
        f"The extension panel remains the pre-result lock of {panel_lock['extension_trait_count']} traits and "
        f"{panel_lock['planned_raw_rg_test_count']} planned sleep-by-extension tests. This review challenges the "
        "extension protocol, software readiness, provenance, and current blocking conditions. It does not review "
        "biological findings that have not been generated.",
        "",
        f"Present real downstream finding classes: {', '.join(present) if present else 'none'}.",
        "",
        f"Absent finding classes: {', '.join(absent)}.",
        "",
        "## Blocking finding",
        "",
        f"The exact compressed source family is {preflight['compressed_source_gib']} GiB and the locked "
        f"{preflight['safety_factor']} safety factor requires {preflight['required_free_gib']} GiB free. "
        f"The recorded preflight found {preflight['available_free_gib']} GiB free, a "
        f"{preflight['shortfall_gib']} GiB shortfall. No bulk download was started. The current volume cannot "
        "support the required retained full-resolution inputs, so substituting HapMap3-only data or a result-selected "
        "subset would violate the fine-mapping/colocalization and selection contracts.",
        "",
        "No threshold has been weakened to make this blocker disappear.",
        "",
        "## Challenge ledger",
        "",
        "| Risk | Required control | Current status | Evidence |",
        "|---|---|---|---|",
    ]
    for row in checklist:
        values = [row["risk"], row["required_control"], row["current_status"], row["current_evidence"]]
        escaped = [value.replace("|", "\\|").replace("\n", " ") for value in values]
        lines.append(f"| {escaped[0]} | {escaped[1]} | {escaped[2]} | {escaped[3]} |")
    lines.extend([
        "",
        "## Findings-level questions that remain blocked",
        "",
        "- Heritability: do rerun h2 Z scores, intercepts, liability conventions, effective sample sizes, and phenotype definitions support each retained extension trait?",
        "- Global correlation: are cross-trait intercepts, valid-allele SNP overlaps, direction stability, duplicate phenotypes, and the isolated 1,200-test family acceptable?",
        "- Novelty: does every FDR hit survive direct pair-by-pair literature searches, publication-duplicate checks, and same-phenotype/same-direction/same-sleep-context review?",
        "- Replication: are cohorts truly non-overlapping, ancestry/build/effect definitions compatible, effects concordant, and heterogeneity acceptable?",
        "- Local and pleiotropic analyses: are source-verified LD references compatible, local h2 estimable, simulation seeds locked, and horizontal/vertical/shared-factor/sample-overlap alternatives still plausible?",
        "- Fine-mapping and colocalization: are dense alleles and signed LD exactly aligned, LD-summary diagnostics clean, credible sets pure, models converged, and H0-H4 conclusions robust over the locked priors?",
        "- Mechanism: does every proposed edge have an exact release/accession, local snapshot, checksum, primary citation, compatible tissue/cell context, and an explicit noncausal claim boundary?",
        "",
        "## Current adversarial conclusion",
        "",
        "The extension design has explicit controls for the main foreseeable failure modes, including pre-result panel and family locks, unavailable-dataset retention, ancestry/build/effect-allele checks, independent replication, local LD validation, prior sensitivity, and edge-level mechanistic provenance. Those controls have synthetic execution evidence only. They cannot validate the absent real analyses.",
        "",
        "Publication-grade discovery claims remain disallowed until acquisition and every downstream findings-level challenge above are completed without weakening thresholds or changing the locked families.",
        "",
    ])
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text("\n".join(lines), encoding="utf-8")
    provenance = {
        "schema_version": "1.0.0",
        "generated_utc": datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
        "verdict": verdict,
        "core_check": core_detail,
        "gate_status_counts": dict(sorted(gate_status_counts.items())),
        "risk_status_counts": dict(sorted(risk_status_counts.items())),
        "present_finding_classes": present,
        "absent_finding_classes": absent,
        "gates_sha256": sha256(args.gates),
        "checklist_sha256": sha256(args.checklist),
        "preflight_sha256": sha256(args.preflight),
        "panel_lock_sha256": sha256(args.panel_lock),
        "output": str(args.out),
        "output_sha256": sha256(args.out),
    }
    args.provenance_out.parent.mkdir(parents=True, exist_ok=True)
    args.provenance_out.write_text(json.dumps(provenance, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(f"ADVERSARIAL_REVIEW_OK verdict={verdict!r} risks={len(checklist)} absent_findings={len(absent)}")


if __name__ == "__main__":
    main()
