#!/usr/bin/env python3
"""Audit the repository against the complete atlas-v1.0 finish line.

This is deliberately an acceptance checker, not a result generator. A missing
or partial scientific artifact remains BLOCKED; the presence of pipeline code
alone never satisfies a scientific gate.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import subprocess
import sys
from dataclasses import asdict, dataclass
from pathlib import Path


SYNTHETIC_MARKER = "SYNTHETIC SMOKE-TEST OUTPUT - NOT REAL LDSC RESULTS"


@dataclass
class Gate:
    gate: str
    status: str
    evidence: str
    blocker: str


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def real_nonempty(path: Path) -> bool:
    if not path.is_file() or path.stat().st_size == 0:
        return False
    if "_smoketest" in path.parts:
        return False
    try:
        return SYNTHETIC_MARKER not in path.read_text(encoding="utf-8", errors="ignore")
    except OSError:
        return False


def artifact_gate(root: Path, name: str, paths: list[str], purpose: str) -> Gate:
    resolved = [root / path for path in paths]
    missing = [str(path.relative_to(root)) for path in resolved if not real_nonempty(path)]
    if not missing:
        return Gate(name, "PASS", ", ".join(paths), "")
    return Gate(
        name,
        "BLOCKED",
        ", ".join(path for path in paths if path not in missing),
        f"{purpose}; missing real non-empty artifact(s): {', '.join(missing)}",
    )


def panel_gate(root: Path) -> tuple[Gate, list[dict[str, str]]]:
    manifest_path = root / "config/analysis_panel.tsv"
    lock_path = root / "config/analysis_panel.lock.json"
    if not manifest_path.is_file() or not lock_path.is_file():
        return Gate("locked_panel", "BLOCKED", "", "manifest or lock is missing"), []
    rows = read_tsv(manifest_path)
    lock = json.loads(lock_path.read_text(encoding="utf-8"))
    ids = [row.get("trait_id", "") for row in rows]
    domains = [row.get("domain", "") for row in rows]
    ordered_sha = hashlib.sha256("".join(f"{trait_id}\n" for trait_id in ids).encode()).hexdigest()
    errors = []
    if len(rows) != 45 or len(set(ids)) != 45:
        errors.append(f"expected 45 unique traits, found {len(rows)} rows/{len(set(ids))} IDs")
    if domains.count("sleep") != 12:
        errors.append(f"expected 12 sleep traits, found {domains.count('sleep')}")
    if ordered_sha != lock.get("ordered_trait_ids_sha256"):
        errors.append("ordered trait IDs do not match the lock")
    if errors:
        return Gate("locked_panel", "BLOCKED", str(manifest_path.relative_to(root)), "; ".join(errors)), rows
    return Gate(
        "locked_panel",
        "PASS",
        "config/analysis_panel.tsv + config/analysis_panel.lock.json",
        "",
    ), rows


def source_gate(rows: list[dict[str, str]]) -> Gate:
    verified = [row["trait_id"] for row in rows if row.get("source_status") == "SOURCE_VERIFIED"]
    pending = [row["trait_id"] for row in rows if row.get("source_status") != "SOURCE_VERIFIED"]
    if len(verified) == 45:
        return Gate("source_curation", "PASS", "45/45 SOURCE_VERIFIED in analysis_panel.tsv", "")
    return Gate(
        "source_curation",
        "BLOCKED",
        f"{len(verified)}/45 SOURCE_VERIFIED",
        f"{len(pending)} source records remain: {', '.join(pending)}",
    )


def source_schema_gate(root: Path, rows: list[dict[str, str]]) -> Gate:
    relative = "config/gwas_schemas.tsv"
    path = root / relative
    if not path.is_file():
        return Gate("source_schemas", "BLOCKED", "", f"missing schema registry: {relative}")
    schemas = {
        (row.get("source_id", ""), row.get("trait_id", "")): row
        for row in read_tsv(path)
    }
    verified = []
    blocked = []
    for row in rows:
        schema = schemas.get((row.get("source_id", ""), row.get("trait_id", "")))
        if (
            row.get("source_status") == "SOURCE_VERIFIED"
            and schema
            and schema.get("schema_status", "").startswith("SCHEMA_VERIFIED")
        ):
            verified.append(row["trait_id"])
        else:
            blocked.append(row["trait_id"])
    if len(verified) == 45:
        return Gate("source_schemas", "PASS", "45/45 selected source schemas verified", "")
    return Gate(
        "source_schemas",
        "BLOCKED",
        f"{len(verified)}/45 selected source schemas verified",
        f"effect/allele/statistic schemas remain unresolved for {len(blocked)} traits: {', '.join(blocked)}",
    )


def per_trait_files_gate(
    root: Path,
    rows: list[dict[str, str]],
    name: str,
    template: str,
    purpose: str,
) -> Gate:
    missing = []
    for row in rows:
        relative = template.format(trait=row["trait_id"])
        if not real_nonempty(root / relative):
            missing.append(row["trait_id"])
    if not missing and len(rows) == 45:
        return Gate(name, "PASS", f"45/45 {template}", "")
    return Gate(
        name,
        "BLOCKED",
        f"{len(rows) - len(missing)}/45 traits",
        f"{purpose}; missing or non-real: {', '.join(missing)}",
    )


def exact_trait_table_gate(
    root: Path,
    rows: list[dict[str, str]],
    name: str,
    relative: str,
    trait_column: str,
    required_columns: set[str],
) -> Gate:
    path = root / relative
    if not real_nonempty(path):
        return Gate(name, "BLOCKED", "", f"missing real non-empty artifact: {relative}")
    table = read_tsv(path)
    columns = set(table[0]) if table else set()
    missing_columns = sorted(required_columns - columns)
    expected = {row["trait_id"] for row in rows}
    observed = {row.get(trait_column, "") for row in table}
    if missing_columns or observed != expected:
        detail = []
        if missing_columns:
            detail.append(f"missing columns {missing_columns}")
        if observed != expected:
            detail.append(f"trait set differs (observed {len(observed)}, expected 45)")
        return Gate(name, "BLOCKED", relative, "; ".join(detail))
    return Gate(name, "PASS", relative, "")


def rg_gate(root: Path, rows: list[dict[str, str]]) -> Gate:
    relative = "results/tables/rg_matrix.tsv"
    path = root / relative
    if not real_nonempty(path):
        return Gate("sleep_disease_rg", "BLOCKED", "", f"missing real non-empty artifact: {relative}")
    table = read_tsv(path)
    columns = set(table[0]) if table else set()
    required = {"sleep_trait", "disease_trait", "rg", "se", "p", "fdr"}
    if not required.issubset(columns):
        return Gate("sleep_disease_rg", "BLOCKED", relative, f"missing columns: {sorted(required - columns)}")
    sleep = {row["trait_id"] for row in rows if row.get("domain") == "sleep"}
    other = {row["trait_id"] for row in rows if row.get("domain") != "sleep"}
    expected = {(a, b) for a in sleep for b in other}
    observed = {(row["sleep_trait"], row["disease_trait"]) for row in table}
    if len(table) != 396 or observed != expected:
        return Gate(
            "sleep_disease_rg",
            "BLOCKED",
            relative,
            f"expected the exact 396-pair family; found {len(table)} rows/{len(observed)} unique pairs",
        )
    return Gate("sleep_disease_rg", "PASS", relative, "")


def covariance_gate(root: Path) -> Gate:
    relative = "results/tables/ldsc_covariance_45x45.tsv"
    if not real_nonempty(root / relative):
        return Gate(
            "full_covariance", "BLOCKED", "",
            f"estimate the complete 45-trait covariance structure; missing real non-empty artifact: {relative}",
        )
    validator = root / "scripts/26_validate_covariance.py"
    result = subprocess.run(
        [sys.executable, str(validator), "--root", str(root), "--quiet"],
        capture_output=True, text=True,
    )
    if result.returncode:
        detail = (result.stdout + result.stderr).strip().replace("\n", "; ")
        return Gate("full_covariance", "BLOCKED", relative, f"covariance validation failed: {detail}")
    return Gate(
        "full_covariance", "PASS",
        "45x45 S/Rg/I + 1035x1035 V + pair estimates + pinned metadata/diagnostics", "",
    )


def genomic_sem_gate(root: Path) -> Gate:
    paths = [
        "results/tables/genomic_sem_model_fit.tsv",
        "results/tables/genomic_sem_factor_loadings.tsv",
    ]
    missing = [path for path in paths if not real_nonempty(root / path)]
    if missing:
        return Gate(
            "genomic_sem", "BLOCKED", "",
            "complete EFA, validation, CFA, and model comparison; missing real non-empty artifact(s): "
            + ", ".join(missing),
        )
    fits = read_tsv(root / paths[0])
    validated = [row for row in fits if row.get("validation_status") == "VALIDATED"]
    if not validated:
        return Gate(
            "genomic_sem", "BLOCKED", ", ".join(paths),
            "real chromosome-split models exist, but no candidate passed held-out validation",
        )
    return Gate("genomic_sem", "PASS", f"{len(validated)} held-out validated model(s)", "")


def lava_gate(root: Path) -> Gate:
    paths = [
        "results/tables/lava_locus_status.tsv",
        "results/tables/lava_univariate.tsv",
        "results/tables/lava_bivariate.tsv",
    ]
    missing = [path for path in paths if not real_nonempty(root / path)]
    if missing:
        return Gate(
            "lava", "BLOCKED", "",
            "run complete local-univariate and eligible bivariate LAVA; missing real non-empty artifact(s): "
            + ", ".join(missing),
        )
    validator = root / "scripts/34_validate_lava.py"
    result = subprocess.run(
        [sys.executable, str(validator), "--root", str(root), "--quiet"],
        capture_output=True, text=True,
    )
    if result.returncode:
        detail = (result.stdout + result.stderr).strip().replace("\n", "; ")
        return Gate("lava", "BLOCKED", ", ".join(paths), f"LAVA validation failed: {detail}")
    return Gate(
        "lava", "PASS",
        "2,495 loci x 45 univariate family + locally eligible 396-pair bivariate family", "",
    )


def release_gate(root: Path) -> Gate:
    release = root / "releases/atlas-v1.0"
    required = [
        "analysis_panel.tsv",
        "checksums.sha256",
        "source_provenance.tsv",
        "tool_versions.tsv",
        "code_commit.txt",
        "release_manifest.json",
    ]
    missing = [name for name in required if not real_nonempty(release / name)]
    if missing:
        return Gate(
            "atlas_v1_release",
            "BLOCKED",
            str(release.relative_to(root)),
            f"immutable release is incomplete; missing: {', '.join(missing)}",
        )
    return Gate("atlas_v1_release", "PASS", "releases/atlas-v1.0", "")


def build_gates(root: Path) -> list[Gate]:
    panel, rows = panel_gate(root)
    gates = [panel]
    if not rows:
        return gates
    gates.extend(
        [
            source_gate(rows),
            source_schema_gate(root, rows),
            per_trait_files_gate(
                root, rows, "harmonization", "data/harmonized/{trait}.qc.txt", "complete harmonization QC ledgers"
            ),
            per_trait_files_gate(
                root, rows, "ldsc_inputs", "data/munged/{trait}.sumstats.gz", "complete canonical HapMap3 inputs"
            ),
            exact_trait_table_gate(
                root,
                rows,
                "ldsc_h2",
                "results/tables/h2_summary.tsv",
                "trait",
                {"trait", "h2", "se", "z", "intercept", "ratio", "verdict"},
            ),
            rg_gate(root, rows),
        ]
    )
    artifact_specs = [
        ("mixer", ["results/tables/mixer_univariate.tsv", "results/tables/mixer_bivariate.tsv"], "run real univariate then eligible bivariate MiXeR"),
        ("pleiotropic_loci", ["results/atlas/shared_loci.tsv"], "combine PLACO and conjunction-FDR evidence"),
        ("factor_gwas", ["results/tables/factor_gwas_summary.tsv", "results/tables/q_snp.tsv"], "run factor GWAS and Q_SNP"),
        ("fine_mapping", ["results/atlas/variants.tsv"], "fine-map priority loci with signal-specific credible sets and PIPs"),
        ("colocalization", ["results/tables/colocalization.tsv"], "complete trait-trait and molecular-QTL signal-level colocalization"),
        ("molecular_integration", ["results/tables/molecular_evidence.tsv", "results/atlas/genes.tsv"], "integrate TWAS, sQTL, pQTL/PWAS, and convergent gene evidence"),
        ("regulatory_mapping", ["results/atlas/regulatory_elements.tsv"], "map fine-mapped variants through regulatory elements to genes"),
        ("cell_types", ["results/atlas/cell_types.tsv"], "complete multi-method cell-type analyses"),
        ("pathways", ["results/atlas/pathways.tsv"], "complete high-confidence pathway and network analyses"),
        ("causal_inference", ["results/atlas/causal_tests.tsv"], "complete bidirectional MR and sensitivity analyses"),
        (
            "integrated_atlas",
            [
                "results/atlas/traits.tsv",
                "results/atlas/trait_pairs.tsv",
                "results/atlas/loci.tsv",
                "results/atlas/variants.tsv",
                "results/atlas/genes.tsv",
                "results/atlas/regulatory_elements.tsv",
                "results/atlas/cell_types.tsv",
                "results/atlas/pathways.tsv",
                "results/atlas/causal_tests.tsv",
                "results/atlas/edges.tsv",
            ],
            "assemble the canonical traceable evidence atlas",
        ),
        ("robustness", ["results/tables/robustness_summary.tsv"], "complete the predefined major robustness pass"),
    ]
    gates.append(covariance_gate(root))
    gates.append(lava_gate(root))
    gates.extend(artifact_gate(root, name, paths, purpose) for name, paths, purpose in artifact_specs)
    gates.append(genomic_sem_gate(root))
    gates.append(release_gate(root))
    return gates


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", default=".", help="repository root")
    parser.add_argument("--json-out", help="optional path for the machine-readable audit")
    parser.add_argument(
        "--report-only",
        action="store_true",
        help="return zero even when gates are blocked (useful for progress reports)",
    )
    args = parser.parse_args()
    root = Path(args.root).resolve()
    gates = build_gates(root)
    for gate in gates:
        detail = gate.evidence if gate.status == "PASS" else gate.blocker
        print(f"{gate.status:7} {gate.gate:24} {detail}")
    passed = sum(gate.status == "PASS" for gate in gates)
    print(f"\nAcceptance: {passed}/{len(gates)} gates passed")
    if args.json_out:
        out = Path(args.json_out)
        out.parent.mkdir(parents=True, exist_ok=True)
        out.write_text(json.dumps([asdict(gate) for gate in gates], indent=2) + "\n", encoding="utf-8")
    return 0 if args.report_only or passed == len(gates) else 1


if __name__ == "__main__":
    raise SystemExit(main())
