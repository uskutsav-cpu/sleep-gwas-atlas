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
        "results/tables/lava_results.provenance.json",
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


def mixer_gate(root: Path) -> Gate:
    paths = ["results/tables/mixer_univariate.tsv", "results/tables/mixer_bivariate.tsv"]
    missing = [path for path in paths if not real_nonempty(root / path)]
    if missing:
        return Gate(
            "mixer", "BLOCKED", "",
            "run complete 20-replicate univariate MiXeR then every eligible sleep-by-non-sleep pair; "
            "missing real non-empty artifact(s): " + ", ".join(missing),
        )
    validator = root / "scripts/40_validate_mixer.py"
    result = subprocess.run(
        [sys.executable, str(validator), "--root", str(root), "--quiet"],
        capture_output=True, text=True,
    )
    if result.returncode:
        detail = (result.stdout + result.stderr).strip().replace("\n", "; ")
        return Gate("mixer", "BLOCKED", ", ".join(paths), f"MiXeR validation failed: {detail}")
    return Gate(
        "mixer", "PASS",
        "45 univariate 20-replicate models + exact AIC-eligible sleep-by-non-sleep pair family", "",
    )


def pleiotropy_gate(root: Path) -> Gate:
    paths = [
        "results/atlas/shared_loci.tsv",
        "results/atlas/shared_loci.provenance.json",
    ]
    missing = [path for path in paths if not real_nonempty(root / path)]
    if missing:
        return Gate(
            "pleiotropic_loci", "BLOCKED", "",
            "complete all 396 PLACO+/conjunction-FDR scans and cross-method collation; "
            f"missing real non-empty artifact(s): {', '.join(missing)}",
        )
    validator = root / "scripts/50_collate_pleiotropy.py"
    result = subprocess.run(
        [sys.executable, str(validator), "--root", str(root), "--validate-only", "--quiet"],
        capture_output=True, text=True,
    )
    if result.returncode:
        detail = (result.stdout + result.stderr).strip().replace("\n", "; ")
        return Gate("pleiotropic_loci", "BLOCKED", ", ".join(paths), f"pleiotropy validation failed: {detail}")
    return Gate(
        "pleiotropic_loci", "PASS",
        "396 locked PLACO+/conjunction-FDR pair scans + same-block consensus", "",
    )


def fine_mapping_gate(root: Path) -> Gate:
    paths = [
        "results/atlas/loci.tsv", "results/atlas/variants.tsv",
        "results/tables/fine_mapping_credible_sets.tsv",
        "results/tables/fine_mapping_diagnostics.tsv",
        "results/tables/trait_trait_colocalization.tsv",
        "results/atlas/fine_mapping.provenance.json",
    ]
    purpose = "fine-map every primary cross-method shared locus with signal-specific credible sets and PIPs"
    missing = [path for path in paths if not real_nonempty(root / path)]
    if missing:
        return Gate(
            "fine_mapping", "BLOCKED", "", f"{purpose}; missing real non-empty artifact(s): {', '.join(missing)}",
        )
    validator = root / "scripts/60_collate_finemapping.py"
    result = subprocess.run(
        [sys.executable, str(validator), "--root", str(root), "--validate-only", "--quiet"],
        capture_output=True, text=True,
    )
    if result.returncode:
        detail = (result.stdout + result.stderr).strip().replace("\n", "; ")
        return Gate("fine_mapping", "BLOCKED", ", ".join(paths), f"fine-mapping validation failed: {detail}")
    return Gate("fine_mapping", "PASS", "all locked primary loci have converged SuSiE-RSS outputs", "")


def integrated_atlas_gate(root: Path) -> Gate:
    paths = [
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
        "results/tables/interpretation_coverage.tsv",
        "results/atlas/interpretation.provenance.json",
        "results/tables/major_conclusions.tsv",
        "results/atlas/edges.provenance.json",
    ]
    missing = [path for path in paths if not real_nonempty(root / path)]
    if missing:
        return Gate(
            "integrated_atlas", "BLOCKED", "",
            "assemble the canonical traceable evidence atlas; missing real non-empty artifact(s): "
            + ", ".join(missing),
        )
    validator = root / "scripts/52_validate_integrated_atlas.py"
    result = subprocess.run(
        [sys.executable, str(validator), "--root", str(root), "--quiet"],
        capture_output=True, text=True,
    )
    if result.returncode:
        detail = (result.stdout + result.stderr).strip().replace("\n", "; ")
        return Gate("integrated_atlas", "BLOCKED", ", ".join(paths), f"atlas validation failed: {detail}")
    graph_validator = root / "scripts/78_build_atlas_edges.py"
    graph_result = subprocess.run(
        [sys.executable, str(graph_validator), "--root", str(root), "--validate-only", "--quiet"],
        capture_output=True, text=True,
    )
    if graph_result.returncode:
        detail = (graph_result.stdout + graph_result.stderr).strip().replace("\n", "; ")
        return Gate("integrated_atlas", "BLOCKED", ", ".join(paths), f"atlas graph validation failed: {detail}")
    return Gate("integrated_atlas", "PASS", "ten schema-locked, cross-linked canonical tables", "")


def robustness_gate(root: Path) -> Gate:
    paths = ["results/tables/robustness_summary.tsv", "results/tables/robustness.provenance.json"]
    missing = [path for path in paths if not real_nonempty(root / path)]
    if missing:
        return Gate(
            "robustness", "BLOCKED", "",
            "complete the predefined major robustness pass; missing real non-empty artifact(s): " + ", ".join(missing),
        )
    collator = root / "scripts/81_collate_robustness.py"
    collated = subprocess.run(
        [sys.executable, str(collator), "--root", str(root), "--validate-only", "--quiet"],
        capture_output=True, text=True,
    )
    if collated.returncode:
        detail = (collated.stdout + collated.stderr).strip().replace("\n", "; ")
        return Gate("robustness", "BLOCKED", ", ".join(paths), f"robustness collation validation failed: {detail}")
    validator = root / "scripts/53_validate_robustness.py"
    result = subprocess.run(
        [sys.executable, str(validator), "--root", str(root), "--quiet"],
        capture_output=True, text=True,
    )
    if result.returncode:
        detail = (result.stdout + result.stderr).strip().replace("\n", "; ")
        return Gate("robustness", "BLOCKED", ", ".join(paths), f"robustness validation failed: {detail}")
    return Gate("robustness", "PASS", "complete nine-family robustness matrix", "")


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
    validator = root / "scripts/55_validate_release.py"
    result = subprocess.run(
        [sys.executable, str(validator), "--root", str(root), "--quiet"],
        capture_output=True, text=True,
    )
    if result.returncode:
        detail = (result.stdout + result.stderr).strip().replace("\n", "; ")
        return Gate("atlas_v1_release", "BLOCKED", "releases/atlas-v1.0", f"release validation failed: {detail}")
    return Gate("atlas_v1_release", "PASS", "checksum-validated releases/atlas-v1.0", "")


def molecular_gates(root: Path) -> list[Gate]:
    required = [
        "results/tables/molecular_qtl_colocalization.tsv",
        "results/tables/colocalization.tsv",
        "results/tables/twas.tsv",
        "results/tables/molecular_locus_coverage.tsv",
        "results/tables/molecular_evidence.tsv",
        "results/atlas/genes.tsv",
        "results/atlas/molecular.provenance.json",
    ]
    missing = [path for path in required if not real_nonempty(root / path)]
    if missing:
        detail = "missing real non-empty artifact(s): " + ", ".join(missing)
        return [
            Gate("colocalization", "BLOCKED", "", "complete trait-trait and molecular-QTL signal-level colocalization; " + detail),
            Gate("molecular_integration", "BLOCKED", "", "integrate corrected TWAS, sQTL, pQTL/PWAS, and convergent gene evidence; " + detail),
        ]
    validator = root / "scripts/73_collate_molecular.py"
    result = subprocess.run(
        [sys.executable, str(validator), "--root", str(root), "--validate-only", "--quiet"],
        capture_output=True, text=True,
    )
    if result.returncode:
        detail = (result.stdout + result.stderr).strip().replace("\n", "; ")
        return [
            Gate("colocalization", "BLOCKED", "results/tables/colocalization.tsv", f"molecular validation failed: {detail}"),
            Gate("molecular_integration", "BLOCKED", "results/tables/molecular_evidence.tsv", f"molecular validation failed: {detail}"),
        ]
    return [
        Gate("colocalization", "PASS", "checksum-validated trait-trait plus molecular-QTL colocalization", ""),
        Gate("molecular_integration", "PASS", "complete four-modality locus coverage and supported-gene integration", ""),
    ]


def interpretation_gates(root: Path) -> list[Gate]:
    families = [
        ("regulatory_mapping", "regulatory", "results/atlas/regulatory_elements.tsv", "map fine-mapped variants through all locked regulatory layers/domains"),
        ("cell_types", "cell_type", "results/atlas/cell_types.tsv", "complete all five cell-type strategies across four domains"),
        ("pathways", "pathway", "results/atlas/pathways.tsv", "complete all four frozen pathway resource families"),
        ("causal_inference", "causal", "results/atlas/causal_tests.tsv", "complete both directions and all five estimator families for 396 pairs"),
    ]
    shared = [
        "results/tables/interpretation_task_manifest.tsv",
        "results/tables/interpretation_task_manifest.lock.json",
        "results/tables/interpretation_coverage.tsv",
        "results/atlas/interpretation.provenance.json",
    ]
    required = shared + [item[2] for item in families]
    missing = [path for path in required if not real_nonempty(root / path)]
    if missing:
        detail = "missing real non-empty artifact(s): " + ", ".join(missing)
        return [Gate(name, "BLOCKED", "", purpose + "; " + detail) for name, _family, _path, purpose in families]
    validator = root / "scripts/77_collate_interpretation.py"
    result = subprocess.run(
        [sys.executable, str(validator), "--root", str(root), "--validate-only", "--quiet"],
        capture_output=True, text=True,
    )
    if result.returncode:
        detail = (result.stdout + result.stderr).strip().replace("\n", "; ")
        return [Gate(name, "BLOCKED", path, f"interpretation validation failed: {detail}") for name, _family, path, _purpose in families]
    coverage = read_tsv(root / "results/tables/interpretation_coverage.tsv")
    gates = []
    for name, family, path, purpose in families:
        rows = [row for row in coverage if row.get("analysis_family") == family]
        blocked = sum(int(row.get("access_blocked_unit_count", "0")) for row in rows)
        incomplete = [row for row in rows if row.get("coverage_status") != "COMPLETE"]
        if not rows or blocked or incomplete:
            gates.append(Gate(
                name, "BLOCKED", path,
                f"{purpose}; coverage rows={len(rows)} access-blocked units={blocked} incomplete rows={len(incomplete)}",
            ))
        else:
            gates.append(Gate(name, "PASS", f"{len(rows)} checksum-validated terminal {family} tasks", ""))
    return gates


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
        ("factor_gwas", ["results/tables/factor_gwas_summary.tsv", "results/tables/q_snp.tsv"], "run factor GWAS and Q_SNP"),
    ]
    gates.append(covariance_gate(root))
    gates.append(lava_gate(root))
    gates.append(mixer_gate(root))
    gates.append(pleiotropy_gate(root))
    gates.append(fine_mapping_gate(root))
    gates.extend(molecular_gates(root))
    gates.extend(artifact_gate(root, name, paths, purpose) for name, paths, purpose in artifact_specs)
    gates.extend(interpretation_gates(root))
    gates.append(integrated_atlas_gate(root))
    gates.append(robustness_gate(root))
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
