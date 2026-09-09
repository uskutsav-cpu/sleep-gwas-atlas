#!/usr/bin/env python3
"""Publish the validated Genomic SEM null without fabricating factor-GWAS results."""
from __future__ import annotations

import argparse
import csv
import hashlib
import io
import json
from pathlib import Path


FACTOR_FIELDS = [
    "factor_id", "variant_id", "chromosome", "position", "effect_allele",
    "other_allele", "beta", "se", "z", "p_value", "status", "reason",
    "provenance_id",
]
Q_SNP_FIELDS = [
    "factor_id", "variant_id", "q_statistic", "q_df", "q_p_value", "status",
    "reason", "provenance_id",
]
MODEL_FIT_FIELDS = [
    "model_id", "factors", "indicators", "discovery_cross_loaders",
    "discovery_boundary_residuals", "chisq", "df", "p_chisq", "AIC", "BIC",
    "CFI", "SRMR", "negative_observed_residuals", "warning_count",
    "validation_status", "validation_reason",
]
MODEL_SCRIPTS = (
    "scripts/25_genomicsem_covariance.R", "scripts/25_genomicsem_covariance.sh",
    "scripts/27_genomicsem_trait_qc.py", "scripts/28_prepare_chromosome_split.py",
    "scripts/28_chromosome_split_covariance.R", "scripts/29_genomicsem_model.R",
    "scripts/30_finalize_genomicsem.py",
)
INPUTS = (
    "results/tables/genomicsem_trait_inclusion.tsv",
    "results/tables/chromosome_split_provenance.json",
    "results/tables/genomicsem_discovery_odd_metadata.tsv",
    "results/tables/genomicsem_validation_even_metadata.tsv",
    "results/tables/genomic_sem_efa_models.tsv",
    "results/tables/genomic_sem_model_fit.tsv",
    "results/tables/genomic_sem_factor_loadings.tsv",
    "results/tables/genomic_sem_model_syntax.tsv",
    "results/tables/genomic_sem_split_diagnostics.tsv",
)


def fail(message: str) -> None:
    raise SystemExit(f"ERROR: {message}")


def sha256(path: Path) -> str:
    if not path.is_file() or path.stat().st_size == 0:
        fail(f"missing real non-empty Genomic SEM artifact: {path}")
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def read_tsv(path: Path) -> tuple[list[str], list[dict[str, str]]]:
    with path.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        return list(reader.fieldnames or []), list(reader)


def table_text(fields: list[str]) -> str:
    output = io.StringIO(newline="")
    csv.DictWriter(output, fieldnames=fields, delimiter="\t", lineterminator="\n").writeheader()
    return output.getvalue()


def atomic_text(path: Path, payload: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(payload, encoding="utf-8")
    temporary.replace(path)


def script_hashes(root: Path) -> dict[str, str]:
    return {relative: sha256(root / relative) for relative in MODEL_SCRIPTS}


def build(root: Path, factor_relative: str, q_relative: str) -> tuple[str, str, dict[str, object]]:
    input_hashes = {relative: sha256(root / relative) for relative in INPUTS}
    fit_fields, fits = read_tsv(root / "results/tables/genomic_sem_model_fit.tsv")
    if fit_fields != MODEL_FIT_FIELDS or len(fits) != 10:
        fail("Genomic SEM candidate family is not the exact ten-model validation family")
    model_ids = [row["model_id"] for row in fits]
    if len(model_ids) != len(set(model_ids)):
        fail("Genomic SEM candidate model IDs are duplicated")
    validated = [row for row in fits if row["validation_status"] == "VALIDATED"]
    invalid = [row for row in fits if row["validation_status"] == "NOT_VALIDATED"]
    if validated or len(invalid) != len(fits) or any(not row["validation_reason"] for row in invalid):
        fail("null finalizer is valid only when every candidate honestly failed held-out validation")
    _, loadings = read_tsv(root / "results/tables/genomic_sem_factor_loadings.tsv")
    _, syntax = read_tsv(root / "results/tables/genomic_sem_model_syntax.tsv")
    if (
        {row["model_id"] for row in loadings} != set(model_ids)
        or {row["model_id"] for row in syntax} != set(model_ids)
        or any(row["validation_status"] != "NOT_VALIDATED" for row in loadings)
    ):
        fail("Genomic SEM model fit, loading, and syntax families disagree")
    factor_text, q_text = table_text(FACTOR_FIELDS), table_text(Q_SNP_FIELDS)
    reason = "NO_CANDIDATE_PASSED_HELD_OUT_CFI_SRMR_AND_RESIDUAL_ADMISSIBILITY"
    provenance = {
        "schema_version": "atlas-v1.0-genomicsem-terminal.1",
        "terminal_status": "NOT_APPLICABLE_NO_VALIDATED_MODEL",
        "terminal_reason": reason,
        "candidate_model_count": len(fits), "validated_model_count": 0,
        "candidate_model_ids": model_ids,
        "input_sha256": input_hashes,
        "script_sha256": script_hashes(root),
        "outputs": {
            factor_relative: hashlib.sha256(factor_text.encode()).hexdigest(),
            q_relative: hashlib.sha256(q_text.encode()).hexdigest(),
        },
        "claim_limit": "No latent factor was promoted; factor GWAS and Q_SNP were not estimated.",
    }
    return factor_text, q_text, provenance


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", default=".")
    parser.add_argument("--factor-out", default="results/tables/factor_gwas_summary.tsv")
    parser.add_argument("--q-snp-out", default="results/tables/q_snp.tsv")
    parser.add_argument("--provenance-out", default="results/tables/factor_gwas.provenance.json")
    parser.add_argument("--validate-only", action="store_true")
    parser.add_argument("--quiet", action="store_true")
    args = parser.parse_args()
    root = Path(args.root).resolve()
    factor_text, q_text, provenance = build(root, args.factor_out, args.q_snp_out)
    provenance_text = json.dumps(provenance, indent=2, sort_keys=True) + "\n"
    outputs = (
        (root / args.factor_out, factor_text), (root / args.q_snp_out, q_text),
        (root / args.provenance_out, provenance_text),
    )
    if args.validate_only:
        for path, payload in outputs:
            if not path.is_file() or path.read_text(encoding="utf-8") != payload:
                fail(f"Genomic SEM terminal publication drifted: {path}")
    else:
        if any(path.exists() for path, _ in outputs):
            fail("Genomic SEM terminal publication already exists; refusing overwrite")
        for path, payload in outputs:
            atomic_text(path, payload)
    if not args.quiet:
        print("GENOMICSEM_TERMINAL_OK candidates=10 validated=0 factor_gwas=NOT_APPLICABLE q_snp=NOT_APPLICABLE")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
