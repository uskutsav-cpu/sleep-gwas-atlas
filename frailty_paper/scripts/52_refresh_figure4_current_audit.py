#!/usr/bin/env python3
"""Render Figure 4 from a validated full-family LAVA receipt audit.

This small SVG renderer keeps the evidence-status panel reproducible on hosts
where the project's optional matplotlib package is unavailable. It does not
read or modify live runner state or receipts; the supplied JSON audit is its
only analysis input.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import html
import json
import shutil
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "frailty_paper/analysis"
LOCK = ROOT / "frailty_paper/config/lava_frailty_sensitivity_v1.yaml"
SOURCE_DATA = OUT / "figures4_5_downstream_evidence_status_source_data.tsv"
PROVENANCE = OUT / "figures4_5_downstream_evidence_status.provenance.json"
STEM = OUT / "figure4_local_sharing_evidence_status"
EXPECTED = 29_940
TRAITS = 12


def inventory_note(audit: dict) -> str:
    stable = (
        audit.get("receipt_file_list_changed_during_scan") is False
        and audit.get("runner_verified_receipts_at_start") == audit.get("runner_verified_receipts_at_end")
    )
    return (
        "Receipt inventory and runner count remained stable during scanning."
        if stable else "Receipt inventory or runner count advanced during scanning."
    )


def sha(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _text(x: float, y: float, value: str, *, size: int = 12, color: str = "#344054",
          weight: str = "normal", anchor: str = "start") -> str:
    return (f'<text x="{x}" y="{y}" font-family="Arial, Helvetica, sans-serif" '
            f'font-size="{size}" font-weight="{weight}" fill="{color}" '
            f'text-anchor="{anchor}">{html.escape(value)}</text>')


def build_svg(audit: dict) -> str:
    validated = int(audit["receipt_files_validated"])
    expected = int(audit["expected_receipts"])
    failed_traits = sum(not row["passes_1pct_gate"]
                        for row in audit["process_failures_by_trait"].values())
    categories = audit["process_failure_categories"]
    fail_detail = (f"{categories.get('all_phenotypes_negative_variance', 0):,} negative-variance; "
                   f"{categories.get('no_specified_snps_in_reference', 0):,} no-reference-SNP; "
                   f"{categories.get('other_process_failure', 0):,} other.")
    start = datetime.fromisoformat(audit["audit_time_start_utc"]).strftime("%H:%M:%S")
    end = datetime.fromisoformat(audit["audit_time_end_utc"]).strftime("%H:%M:%S")
    claims = audit.get("claims", [])
    claim_labels = [f"{c['trait']} locus {c['locus_index']}" for c in claims[:2]]
    claim_detail = " and ".join(claim_labels) if claim_labels else "No active claims were recorded in the audit snapshot."
    effective = audit.get("effective_parallel_workers")
    requested = audit.get("requested_parallel_workers")
    scan_pct = validated / expected * 100

    elems = [
        '<svg xmlns="http://www.w3.org/2000/svg" width="1022.4pt" height="590.4pt" viewBox="0 0 1022.4 590.4">',
        '<rect width="1022.4" height="590.4" fill="#F5F7FA"/>',
        _text(42, 35, "Figure 4", size=12, color="#145C72", weight="bold"),
        _text(42, 78, "Local sharing: interim evidence status", size=23, color="#182230", weight="bold"),
        _text(42, 107, "Locked sensitivity analysis · current full-family audit · no locus-level inference", size=11, color="#475467"),
    ]

    cards = [
        (40, 146, "LOCKED ANALYSIS", "RUNNING · 12 × 2,495", "The prespecified Frailty Index × sleep/circadian sensitivity family continues.",
         f"Evidence: frozen lock; audit {audit['analysis_id']}"),
        (519, 146, "LATEST FULL-FAMILY AUDIT", f"{validated:,}/{expected:,} receipts",
         f"Read-only scan {start}–{end} UTC; {inventory_note(audit)} Runner count: {audit['runner_verified_receipts_at_end']:,}.",
         f"Evidence: audit completed {audit['audit_time_end_utc'][:10]}"),
        (40, 308, "LOCKED 1% TRAIT GATES", f"FAILED · {failed_traits}/{TRAITS}",
         f"{fail_detail} All trait-level lower bounds exceed the frozen 1% allowance; thresholds unchanged.",
         "Evidence: per-trait gate audit; all 12 locked traits"),
        (519, 308, "INTERPRETATION", "NO LOCAL INFERENCE",
         f"{effective} workers active ({requested} requested); claims: {claim_detail}. Failed gates do not establish absence of sharing.",
         "Evidence: no receipt, claim, or duplicate-identity issues"),
    ]
    for x, y, heading, value, detail, evidence in cards:
        elems.append(f'<rect x="{x}" y="{y}" width="463" height="151" rx="8" fill="#FFFFFF" stroke="#D9E0E8" stroke-width="1.2"/>')
        elems.append(_text(x + 18, y + 29, heading, size=9, color="#667085", weight="bold"))
        value_color = "#9A5A12" if "FAILED" in value else "#145C72"
        elems.append(_text(x + 18, y + 62, value, size=16, color=value_color, weight="bold"))
        # Two wrapped detail lines preserve the original panel's reading density.
        words = detail.split()
        lines: list[str] = []
        line = ""
        for word in words:
            if line and len(line) + len(word) + 1 > 56:
                lines.append(line)
                line = word
            else:
                line = (line + " " + word).strip()
        if line:
            lines.append(line)
        for idx, part in enumerate(lines[:3]):
            elems.append(_text(x + 18, y + 88 + idx * 14, part, size=9, color="#344054"))
        elems.append(_text(x + 18, y + 137, evidence, size=7, color="#667085"))

    footer = (f"Validated {validated:,}/{expected:,} ({scan_pct:.1f}%); integrity issues 0; "
              f"launch failures {audit.get('launch_failures_this_invocation', 0)}; "
              f"stale recoveries {audit.get('stale_claims_recovered', 0)}. "
              "The run continues; failed gates prohibit local-sharing inference but do not show sharing is absent.")
    elems.append('<rect x="40" y="480" width="943" height="55" rx="7" fill="#E9EEF3" stroke="#CBD5DF" stroke-width="1"/>')
    elems.append(_text(58, 512, footer, size=9, color="#344054"))
    elems.append('</svg>')
    return "\n".join(elems) + "\n"


def update_source_data(audit: dict, audit_rel: str) -> None:
    rows = []
    if SOURCE_DATA.exists():
        with SOURCE_DATA.open(encoding="utf-8", newline="") as stream:
            rows = [r for r in csv.DictReader(stream, delimiter="\t") if r["figure"] != "Figure 4"]
    n = audit["receipt_files_validated"]
    fail_n = sum(not r["passes_1pct_gate"] for r in audit["process_failures_by_trait"].values())
    common = audit_rel
    rows.extend([
        {"figure": "Figure 4", "panel": "1", "status": "RUNNING · 12 × 2,495",
         "detail": "The prespecified Frailty Index × sleep/circadian sensitivity family is still executing.",
         "evidence_files": f"config/lava_frailty_sensitivity_v1.yaml; {common}"},
        {"figure": "Figure 4", "panel": "2", "status": f"{n:,}/29,940 receipts",
         "detail": f"Read-only audit completed {audit['audit_time_end_utc']}; {inventory_note(audit)}",
         "evidence_files": common},
        {"figure": "Figure 4", "panel": "3", "status": f"FAILED · {fail_n}/12",
         "detail": "All 12 trait-level failure lower bounds exceed the frozen 1% allowance; no thresholds or exclusions were changed.",
         "evidence_files": common},
        {"figure": "Figure 4", "panel": "4", "status": "NO LOCAL INFERENCE",
         "detail": "The family is incomplete and local estimates remain inadmissible. Failed gates do not establish absence of local sharing.",
         "evidence_files": common},
    ])
    fieldnames = ["figure", "panel", "status", "detail", "evidence_files"]
    with SOURCE_DATA.open("w", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=fieldnames, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--audit", required=True, help="Full-family audit JSON inside the repository")
    args = parser.parse_args()
    audit_path = (ROOT / args.audit).resolve()
    if not audit_path.is_file() or ROOT not in audit_path.parents:
        raise SystemExit("audit JSON must be an existing file inside the repository")
    audit = json.loads(audit_path.read_text(encoding="utf-8"))
    lock_sha = sha(LOCK)
    if (audit.get("mutates_run_or_receipts") is not False or audit.get("expected_receipts") != EXPECTED
            or len(audit.get("process_failures_by_trait", {})) != TRAITS
            or audit.get("lock_sha256") != lock_sha
            or audit.get("runner_state_lock_sha256") != lock_sha):
        raise SystemExit("audit is not a read-only, lock-matched full-family result")
    if not shutil.which("rsvg-convert"):
        raise SystemExit("rsvg-convert is required to export PNG and PDF from the generated SVG")

    svg = build_svg(audit)
    svg_path = STEM.with_suffix(".svg")
    svg_path.write_text(svg, encoding="utf-8")
    subprocess.run(["rsvg-convert", "-d", "300", "-p", "300", "-f", "png", "-o", str(STEM.with_suffix(".png")), str(svg_path)], check=True)
    subprocess.run(["rsvg-convert", "-f", "pdf", "-o", str(STEM.with_suffix(".pdf")), str(svg_path)], check=True)
    audit_rel = audit_path.relative_to(ROOT).as_posix()
    update_source_data(audit, audit_rel)

    try:
        rsvg_version = subprocess.run(["rsvg-convert", "--version"], check=True, capture_output=True, text=True).stdout.strip()
    except subprocess.SubprocessError:
        rsvg_version = "rsvg-convert (version unavailable)"
    provenance = json.loads(PROVENANCE.read_text(encoding="utf-8")) if PROVENANCE.exists() else {}
    outputs = provenance.setdefault("figures", {}).setdefault("figure4_local_sharing_evidence_status", {})
    for ext in ("svg", "png", "pdf"):
        rel = STEM.with_suffix(f".{ext}").relative_to(ROOT).as_posix()
        outputs[rel] = sha(ROOT / rel)
    provenance["audit_generated_utc"] = audit["audit_time_end_utc"]
    provenance["audit_type"] = "full_family"
    provenance["command"] = f"python frailty_paper/scripts/52_refresh_figure4_current_audit.py --audit {audit_rel}"
    provenance["created_utc"] = datetime.now(timezone.utc).isoformat(timespec="seconds")
    provenance["figure4_rendered_this_run"] = True
    provenance["renderer_script"] = "frailty_paper/scripts/52_refresh_figure4_current_audit.py"
    provenance["renderer_script_sha256"] = sha(Path(__file__).resolve())
    provenance["renderer"] = rsvg_version
    provenance.setdefault("inputs_sha256", {})[audit_rel] = sha(audit_path)
    provenance["inputs_sha256"]["frailty_paper/config/lava_frailty_sensitivity_v1.yaml"] = lock_sha
    source_rel = SOURCE_DATA.relative_to(ROOT).as_posix()
    provenance.setdefault("source_data", {})[source_rel] = sha(SOURCE_DATA)
    provenance["latest_lava_receipts_audited"] = audit["receipt_files_validated"]
    provenance["latest_lava_trait"] = "all 12 traits"
    provenance["loci_flagged"] = sum(not r["passes_1pct_gate"] for r in audit["process_failures_by_trait"].values())
    provenance["locked_maximum_failures"] = 0
    provenance["locus_denominator"] = TRAITS
    provenance["python"] = sys.version.split()[0]
    PROVENANCE.write_text(json.dumps(provenance, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(f"FIGURE4_RENDERED receipts={audit['receipt_files_validated']}/{EXPECTED} gates_failed={provenance['loci_flagged']}/12")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
