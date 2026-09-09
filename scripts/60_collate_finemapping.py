#!/usr/bin/env python3
"""Collate every locked primary locus into fine-mapping and atlas tables."""
from __future__ import annotations

import argparse
import csv
import hashlib
import io
import json
import math
from pathlib import Path

import fine_mapping_contract

LOCUS_FIELDS = fine_mapping_contract.CANONICAL_LOCUS_FIELDS
VARIANT_FIELDS = fine_mapping_contract.CANONICAL_VARIANT_FIELDS


def fail(message: str) -> None:
    raise SystemExit(f"ERROR: {message}")


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def read_tsv(path: Path) -> tuple[list[str], list[dict[str, str]]]:
    if not path.is_file() or path.stat().st_size == 0:
        fail(f"missing real non-empty artifact: {path}")
    with path.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        return reader.fieldnames or [], list(reader)


def probability(value: str, field: str, identity: str) -> float:
    try:
        parsed = float(value)
    except ValueError as exc:
        fail(f"invalid {field} for {identity}")
        raise AssertionError from exc
    if not math.isfinite(parsed) or not 0 <= parsed <= 1:
        fail(f"invalid {field} for {identity}")
    return parsed


def table_text(fields: list[str], rows: list[dict[str, object]]) -> str:
    output = io.StringIO(newline="")
    writer = csv.DictWriter(output, fieldnames=fields, delimiter="\t", lineterminator="\n")
    writer.writeheader()
    writer.writerows(rows)
    return output.getvalue()


def atomic_text(path: Path, text: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(text, encoding="utf-8")
    temporary.replace(path)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", default=".")
    parser.add_argument("--manifest", default="results/tables/fine_mapping_locus_manifest.tsv")
    parser.add_argument("--manifest-lock", default="results/tables/fine_mapping_locus_manifest.lock.json")
    parser.add_argument("--policy", default="config/fine_mapping_analysis_policy.json")
    parser.add_argument("--preflight", default="results/tables/fine_mapping_preflight.json")
    parser.add_argument("--run-dir", default="results/fine_mapping/runs")
    parser.add_argument("--loci-out", default="results/atlas/loci.tsv")
    parser.add_argument("--variants-out", default="results/atlas/variants.tsv")
    parser.add_argument("--credible-out", default="results/tables/fine_mapping_credible_sets.tsv")
    parser.add_argument("--diagnostics-out", default="results/tables/fine_mapping_diagnostics.tsv")
    parser.add_argument("--coloc-out", default="results/tables/trait_trait_colocalization.tsv")
    parser.add_argument("--provenance-out", default="results/atlas/fine_mapping.provenance.json")
    parser.add_argument("--validate-only", action="store_true")
    parser.add_argument("--quiet", action="store_true")
    args = parser.parse_args()
    root = Path(args.root).resolve()
    manifest_path = root / args.manifest
    manifest_lock_path = root / args.manifest_lock
    policy_path, policy = fine_mapping_contract.load_policy(root, args.policy)
    preflight_path = root / args.preflight
    fine_mapping_contract.validate_preflight(root, policy_path, policy, preflight_path)
    manifest, lock = fine_mapping_contract.validate_manifest(
        root, manifest_path, manifest_lock_path, policy_path, policy, preflight_path,
    )

    loci = []
    variants = []
    credible_sets = []
    diagnostics = []
    coloc = []
    input_hashes = {}
    for locus in manifest:
        identity = locus["shared_locus_id"]
        run = root / args.run_dir / identity
        paths = {name: run / name for name in (
            "variants.tsv", "credible_sets.tsv", "colocalization.tsv",
            "shared_variant_posteriors.tsv", "diagnostics.tsv", "provenance.json",
        )}
        if any(not path.is_file() or path.stat().st_size == 0 for path in paths.values()):
            fail(f"fine-mapping run is incomplete for {identity}")
        task_path = root / fine_mapping_contract.safe_relative(locus["task_path"], "fine-mapping task")
        task_lock_path = task_path.with_suffix(".lock.json")
        if not task_path.is_file() or not task_lock_path.is_file():
            fail(f"fine-mapping task/lock is absent for {identity}")
        task_lock = json.loads(task_lock_path.read_text(encoding="utf-8"))
        task_fields, task_rows = read_tsv(task_path)
        del task_fields
        if (
            len(task_rows) != 1 or task_rows[0].get("queue_row_id") != locus["comparison_id"]
            or task_lock.get("schema_version") != "atlas-v1.0-finemapping-task.2"
            or task_lock.get("task_sha256") != sha256(task_path)
            or task_lock.get("manifest_sha256") != sha256(manifest_path)
            or task_lock.get("manifest_lock_sha256") != sha256(manifest_lock_path)
            or task_lock.get("policy_sha256") != sha256(policy_path)
            or task_lock.get("preflight_sha256") != sha256(preflight_path)
            or task_lock.get("script_sha256") != fine_mapping_contract.script_hashes(root)
            or task_lock.get("results_accessed_before_lock") is not False
        ):
            fail(f"fine-mapping task/lock provenance drifted for {identity}")
        provenance = json.loads(paths["provenance.json"].read_text(encoding="utf-8"))
        engine = root / policy["shared_engine"]["path"]
        if (
            provenance.get("schema_version") != "atlas-v1.0-finemapping-run.2"
            or provenance.get("comparison_id") != locus["comparison_id"]
            or provenance.get("task_sha256") != sha256(task_path)
            or provenance.get("task_lock_sha256") != sha256(task_lock_path)
            or provenance.get("manifest_sha256") != sha256(manifest_path)
            or provenance.get("manifest_lock_sha256") != sha256(manifest_lock_path)
            or provenance.get("policy_sha256") != sha256(policy_path)
            or provenance.get("preflight_sha256") != sha256(preflight_path)
            or provenance.get("engine_sha256") != sha256(engine)
            or provenance.get("script_sha256") != fine_mapping_contract.script_hashes(root)
        ):
            fail(f"fine-mapping provenance identity differs for {identity}")
        for name in paths:
            if name == "provenance.json":
                continue
            expected = provenance.get("outputs", {}).get(name, {}).get("sha256")
            expected_rows = provenance.get("outputs", {}).get(name, {}).get("rows")
            observed_fields, observed_rows = read_tsv(paths[name])
            if (
                expected != sha256(paths[name]) or expected_rows != len(observed_rows)
                or observed_fields != fine_mapping_contract.ENGINE_OUTPUT_FIELDS[name]
            ):
                fail(f"fine-mapping output differs from provenance for {identity}/{name}")
        _, variant_rows = read_tsv(paths["variants.tsv"])
        _, credible_rows = read_tsv(paths["credible_sets.tsv"])
        _, diagnostic_rows = read_tsv(paths["diagnostics.tsv"])
        _, coloc_rows = read_tsv(paths["colocalization.tsv"])
        _, shared_rows = read_tsv(paths["shared_variant_posteriors.tsv"])
        if not variant_rows or len(diagnostic_rows) != 2 or any(row["model_converged"] != "TRUE" for row in diagnostic_rows):
            fail(f"SuSiE-RSS did not converge for both traits at {identity}")
        if any(row["comparison_id"] != locus["comparison_id"] for row in variant_rows + credible_rows + diagnostic_rows + coloc_rows + shared_rows):
            fail(f"out-of-family fine-mapping row for {identity}")
        by_dataset: dict[str, dict[str, dict[str, str]]] = {}
        for row in variant_rows:
            probability(row["PIP"], "PIP", f"{identity}/{row['SNP']}")
            current = by_dataset.setdefault(row["dataset_id"], {})
            if row["SNP"] in current:
                fail(f"duplicate variant/dataset row for {identity}/{row['SNP']}")
            current[row["SNP"]] = row
        if set(by_dataset) != {locus["sleep_trait"], locus["non_sleep_trait"]}:
            fail(f"fine-mapping datasets differ from locked traits for {identity}")
        sleep = by_dataset[locus["sleep_trait"]]
        non_sleep = by_dataset[locus["non_sleep_trait"]]
        if set(sleep) != set(non_sleep):
            fail(f"fine-mapping trait variant universes differ for {identity}")
        qc_values = sorted({row["diagnostic_status"] for row in diagnostic_rows})
        qc_status = ";".join(qc_values)
        primary_p12 = float(policy["colocalization"]["p12_primary"])
        primary_coloc = [row for row in coloc_rows if row["prior_role"] == "PRIMARY" and float(row["p12"]) == primary_p12]
        if not primary_coloc:
            fail(f"primary coloc-SuSiE prior is absent for {identity}")
        supported = any(
            row["coloc_method"] == "COLOC_SUSIE"
            and row["fine_mapping_qc"] == "PASS"
            and probability(row["PP_H4"], "PP_H4", identity) >= policy["colocalization"]["shared_signal_min_pp_h4"]
            and float(row["PP_H4_over_PP_H3"]) >= policy["colocalization"]["shared_signal_min_pp_h4_over_pp_h3"]
            for row in primary_coloc if row["PP_H4"] not in {"", "NA"} and row["PP_H4_over_PP_H3"] not in {"", "NA"}
        )
        shared_pp: dict[str, float] = {}
        for row in shared_rows:
            if float(row["p12"]) != primary_p12:
                continue
            value = probability(row["SNP_PP_H4"], "SNP_PP_H4", f"{identity}/{row['SNP']}")
            shared_pp[row["SNP"]] = max(value, shared_pp.get(row["SNP"], 0.0))
        provenance_id = "FM:" + sha256(paths["provenance.json"])
        loci.append({
            "locus_id": identity, "chromosome": locus["chromosome"],
            "start_bp": locus["start_bp"], "end_bp": locus["end_bp"],
            "lead_snp": locus["placo_lead_snp"], "sleep_trait": locus["sleep_trait"],
            "non_sleep_trait": locus["non_sleep_trait"], "analysis_tier": locus["analysis_tier"],
            "placo_lead_p": locus["placo_lead_p"], "conjfdr_min_q": locus["conjfdr_lead_fdr"],
            "effect_direction": locus["effect_direction"], "lava_block_id": locus["lava_block_id"],
            "method_support": "PLACO_PLUS;CONJFDR;SUSIE_RSS;COLOC_SUSIE",
            "evidence_level": "TRAIT_TRAIT_SHARED_SIGNAL_SUPPORTED" if supported else "FINE_MAPPED_NO_TRAIT_TRAIT_COLOC_SUPPORT",
            "provenance_id": provenance_id,
        })
        for snp in sleep:
            first, second = sleep[snp], non_sleep[snp]
            identity_fields = ("CHR", "BP", "A1", "A2", "MAF")
            if any(first[field] != second[field] for field in identity_fields):
                fail(f"fine-mapping allele/coordinate drift for {identity}/{snp}")
            variants.append({
                "variant_id": snp, "locus_id": locus["shared_locus_id"], "rsid": snp,
                "chromosome": first["CHR"], "position_bp": first["BP"],
                "effect_allele": first["A1"], "other_allele": first["A2"],
                "minor_allele_frequency": first["MAF"], "sleep_beta": first["BETA"],
                "sleep_se": first["SE"], "non_sleep_beta": second["BETA"],
                "non_sleep_se": second["SE"], "pip_sleep": first["PIP"],
                "pip_non_sleep": second["PIP"], "credible_set_sleep": first["credible_set_ids"],
                "credible_set_non_sleep": second["credible_set_ids"],
                "shared_signal_posterior": format(shared_pp[snp], ".15g") if snp in shared_pp else "NA",
                "fine_mapping_method": "SuSiE-RSS_0.14.2", "ld_reference": policy["reference"]["id"],
                "qc_status": qc_status, "provenance_id": provenance_id,
            })
        credible_sets.extend(credible_rows)
        diagnostics.extend(diagnostic_rows)
        coloc.extend(coloc_rows)
        input_hashes[identity] = {name: sha256(path) for name, path in paths.items()}

    credible_fields = fine_mapping_contract.ENGINE_OUTPUT_FIELDS["credible_sets.tsv"]
    diagnostic_fields = fine_mapping_contract.ENGINE_OUTPUT_FIELDS["diagnostics.tsv"]
    coloc_fields = fine_mapping_contract.ENGINE_OUTPUT_FIELDS["colocalization.tsv"]
    payloads = {
        root / args.loci_out: table_text(LOCUS_FIELDS, loci),
        root / args.variants_out: table_text(VARIANT_FIELDS, variants),
        root / args.credible_out: table_text(credible_fields, credible_sets),
        root / args.diagnostics_out: table_text(diagnostic_fields, diagnostics),
        root / args.coloc_out: table_text(coloc_fields, coloc),
    }
    provenance = {
        "schema_version": "atlas-v1.0-finemapping-canonical.2",
        "analysis_id": policy["analysis_id"], "locus_count": len(loci),
        "variant_count": len(variants), "credible_set_count": len(credible_sets),
        "colocalization_row_count": len(coloc), "manifest_sha256": sha256(manifest_path),
        "manifest_lock_sha256": sha256(manifest_lock_path), "policy_sha256": sha256(policy_path),
        "preflight_sha256": sha256(preflight_path),
        "shared_loci_sha256": lock["shared_loci_sha256"],
        "shared_loci_provenance_sha256": lock["shared_loci_provenance_sha256"],
        "zero_family_not_applicable": len(manifest) == 0,
        "run_input_hashes": input_hashes,
        "script_sha256": fine_mapping_contract.script_hashes(root),
        "outputs": {str(path.relative_to(root)): hashlib.sha256(text.encode()).hexdigest() for path, text in payloads.items()},
        "claim_limit": policy["claim_limit"],
    }
    provenance_text = json.dumps(provenance, indent=2, sort_keys=True) + "\n"
    provenance_path = root / args.provenance_out
    if args.validate_only:
        for path, text in payloads.items():
            if not path.is_file() or path.read_text(encoding="utf-8") != text:
                fail(f"fine-mapping aggregate drifted: {path}")
        if not provenance_path.is_file() or provenance_path.read_text(encoding="utf-8") != provenance_text:
            fail("fine-mapping provenance drifted")
    else:
        if provenance_path.exists() or any(path.exists() for path in payloads):
            fail("immutable canonical fine-mapping publication already exists")
        for path, text in payloads.items():
            atomic_text(path, text)
        atomic_text(provenance_path, provenance_text)
    if not args.quiet:
        print(f"FINEMAPPING_COLLATION_OK loci={len(loci)} variants={len(variants)} coloc_rows={len(coloc)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
