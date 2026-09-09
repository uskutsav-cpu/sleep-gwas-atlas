#!/usr/bin/env python3
"""Run checksum-bound MAGMA v1.10 gene analysis for one locked atlas trait."""
from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
import math
import re
import subprocess
import tempfile
from pathlib import Path


def fail(message: str) -> None:
    raise SystemExit(f"ERROR: {message}")


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def atomic_json(path: Path, value: object) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(json.dumps(value, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    temporary.replace(path)


def read_panel(path: Path, trait_id: str) -> dict[str, str]:
    with path.open(encoding="utf-8", newline="") as handle:
        rows = [row for row in csv.DictReader(handle, delimiter="\t") if row["trait_id"] == trait_id]
    if len(rows) != 1:
        fail(f"trait is absent or duplicated in the locked panel: {trait_id}")
    return rows[0]


def materialize_pvalues(source: Path, output: Path, policy: dict[str, object]) -> dict[str, object]:
    fields = [policy["gwas_snp_field"], policy["gwas_p_field"], policy["gwas_n_field"]]
    floor = float(policy["p_value_floor"])
    rows, floored = 0, 0
    digest = hashlib.sha256()
    with gzip.open(source, "rt", encoding="utf-8", newline="") as input_handle, output.open(
        "w", encoding="utf-8", newline=""
    ) as output_handle:
        reader = csv.DictReader(input_handle, delimiter="\t")
        if reader.fieldnames is None or not set(fields).issubset(reader.fieldnames):
            fail("harmonized GWAS lacks the locked MAGMA SNP/P/N columns")
        header = "\t".join(fields) + "\n"
        output_handle.write(header)
        digest.update(header.encode())
        for line_number, row in enumerate(reader, start=2):
            snp = row[fields[0]]
            try:
                p_value, sample_size = float(row[fields[1]]), float(row[fields[2]])
            except ValueError as exc:
                fail(f"nonnumeric MAGMA P or N at harmonized line {line_number}")
                raise AssertionError from exc
            if not snp or not math.isfinite(p_value) or not 0 <= p_value <= 1:
                fail(f"invalid MAGMA SNP or P at harmonized line {line_number}")
            if not math.isfinite(sample_size) or sample_size <= 0:
                fail(f"invalid MAGMA per-SNP sample size at harmonized line {line_number}")
            if p_value < floor:
                p_value = floor
                floored += 1
            payload = f"{snp}\t{format(p_value, '.12g')}\t{format(sample_size, '.12g')}\n"
            output_handle.write(payload)
            digest.update(payload.encode())
            rows += 1
    if not rows:
        fail("harmonized GWAS has no rows for MAGMA")
    return {"rows": rows, "p_values_floored": floored, "sha256": digest.hexdigest()}


def magma_gene_count(path: Path) -> int:
    header_seen, rows = False, 0
    with path.open(encoding="utf-8") as handle:
        for line in handle:
            fields = line.split()
            if not fields:
                continue
            if fields[0] == "GENE":
                header_seen = True
                continue
            if header_seen and re.fullmatch(r"ENSG[0-9]{11}", fields[0]):
                rows += 1
    if not header_seen:
        fail("MAGMA gene output lacks its GENE header")
    return rows


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("trait_id")
    parser.add_argument("--root", default=".")
    parser.add_argument("--policy", default="config/interpretation_analysis_policy.json")
    parser.add_argument("--out-prefix")
    parser.add_argument("--provenance-out")
    args = parser.parse_args()
    root = Path(args.root).resolve()
    policy_path = root / args.policy
    policy = json.loads(policy_path.read_text(encoding="utf-8"))
    spec = policy["fuma_scrna"]
    trait = read_panel(root / "config/analysis_panel.tsv", args.trait_id)
    manifest_path = root / spec["component_manifest"]
    if not manifest_path.is_file() or sha256(manifest_path) != spec["component_manifest_sha256"]:
        fail("FUMA component manifest differs from its policy pin")
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    components = {row["component_id"]: row for row in manifest["components"]}
    binary = root / spec["magma_binary_path"]
    binary_pin = components["MAGMA_V1_10_ARM64_BINARY"]
    if (
        not binary.is_file() or binary.stat().st_size != binary_pin["bytes"]
        or sha256(binary) != binary_pin["sha256"]
    ):
        fail("MAGMA binary differs from its exact runtime pin")
    reference_prefix = root / spec["reference_prefix"]
    reference_hashes: dict[str, str] = {}
    for pin in manifest["reference_members"]:
        path = reference_prefix.parent / pin["name"]
        if not path.is_file() or path.stat().st_size != pin["bytes"] or sha256(path) != pin["sha256"]:
            fail(f"extracted MAGMA reference differs from its pin: {pin['name']}")
        reference_hashes[pin["name"]] = sha256(path)
    annotation_prefix = root / spec["gene_annotation_prefix"]
    annotation = Path(str(annotation_prefix) + ".genes.annot")
    annotation_provenance_path = root / spec["gene_annotation_provenance_path"]
    if not annotation.is_file() or not annotation_provenance_path.is_file():
        fail("checksum-bound MAGMA gene annotation and provenance are incomplete")
    annotation_provenance = json.loads(annotation_provenance_path.read_text(encoding="utf-8"))
    if (
        annotation_provenance.get("policy_sha256") != sha256(policy_path)
        or annotation_provenance.get("annotation_sha256") != sha256(annotation)
    ):
        fail("MAGMA gene annotation differs from its policy-bound provenance")
    gwas = root / spec["gwas_path_template"].format(trait_id=args.trait_id)
    if not gwas.is_file() or gwas.stat().st_size == 0:
        fail(f"harmonized GWAS is absent for MAGMA trait {args.trait_id}")
    default_prefix = Path(spec["gene_results_dir"]) / args.trait_id / args.trait_id
    prefix = root / (args.out_prefix or default_prefix)
    genes_raw, genes_out, log = (
        Path(str(prefix) + ".genes.raw"), Path(str(prefix) + ".genes.out"), Path(str(prefix) + ".log")
    )
    provenance_path = root / (
        args.provenance_out or str(Path(spec["gene_results_dir"]) / args.trait_id / "provenance.json")
    )
    if any(path.exists() for path in (genes_raw, genes_out, log, provenance_path)):
        fail("MAGMA gene-analysis output already exists; refusing overwrite")
    prefix.parent.mkdir(parents=True, exist_ok=True)
    temporary_prefix = prefix.with_name(prefix.name + ".tmp")
    temporary_outputs = [
        Path(str(temporary_prefix) + suffix) for suffix in (".genes.raw", ".genes.out", ".log")
    ]
    if any(path.exists() for path in temporary_outputs):
        fail("stale MAGMA gene-analysis temporary output exists; inspect it before retrying")
    with tempfile.TemporaryDirectory(prefix=f"sleep_atlas_magma_{args.trait_id}_") as temporary:
        pvalues = Path(temporary) / f"{args.trait_id}.pval.tsv"
        pvalue_details = materialize_pvalues(gwas, pvalues, spec)
        command = [
            str(binary), "--bfile", str(reference_prefix),
            f"synonym-dup={spec['reference_synonym_duplicate_rule']}",
            "--pval", str(pvalues), f"ncol={spec['gwas_n_field']}",
            f"duplicate={spec['gwas_duplicate_rule']}",
            "--gene-annot", str(annotation), "--gene-model", str(spec["gene_model"]),
            "--out", str(temporary_prefix),
        ]
        result = subprocess.run(command, cwd=root, capture_output=True, text=True, check=False)
    if result.returncode or any(not path.is_file() for path in temporary_outputs):
        fail(f"MAGMA gene analysis failed for {args.trait_id}: {(result.stdout + result.stderr).strip()}")
    log_text = temporary_outputs[2].read_text(encoding="utf-8")
    match = re.search(
        r"read\s+([0-9,]+)\s+lines from file, containing valid SNP p-values for\s+([0-9,]+)\s+SNPs in data",
        log_text,
    )
    if match is None:
        fail("MAGMA log lacks the valid SNP-overlap audit line")
    source_lines, matched_snps = (int(value.replace(",", "")) for value in match.groups())
    overlap_fraction = matched_snps / int(pvalue_details["rows"])
    if source_lines < int(pvalue_details["rows"]) or overlap_fraction < float(spec["minimum_reference_overlap_fraction"]):
        fail(
            f"MAGMA GWAS/reference overlap is below policy: matched={matched_snps} "
            f"rows={pvalue_details['rows']} fraction={overlap_fraction:.6g}"
        )
    gene_rows = magma_gene_count(temporary_outputs[1])
    if gene_rows < int(spec["minimum_genes_in_gene_property_model"]):
        fail(f"MAGMA gene analysis produced too few genes: {gene_rows}")
    for temporary_output, output in zip(temporary_outputs, (genes_raw, genes_out, log)):
        temporary_output.replace(output)
    atomic_json(provenance_path, {
        "schema_version": policy["schema_version"], "analysis_id": policy["analysis_id"],
        "trait_id": args.trait_id, "trait_source_id": trait["source_id"],
        "policy_sha256": sha256(policy_path), "component_manifest_sha256": sha256(manifest_path),
        "source_release": manifest["release"], "magma_version": spec["magma_version"],
        "magma_binary_sha256": sha256(binary), "reference_member_sha256": reference_hashes,
        "annotation_sha256": sha256(annotation), "annotation_provenance_sha256": sha256(annotation_provenance_path),
        "gwas_path": str(gwas), "gwas_sha256": sha256(gwas), "pvalue_materialization": pvalue_details,
        "command": command, "source_lines_reported_by_magma": source_lines,
        "matched_reference_snps": matched_snps, "reference_overlap_fraction": overlap_fraction,
        "gene_rows": gene_rows, "genes_raw_path": str(genes_raw), "genes_raw_sha256": sha256(genes_raw),
        "genes_out_path": str(genes_out), "genes_out_sha256": sha256(genes_out),
        "log_path": str(log), "log_sha256": sha256(log),
    })
    print(
        f"MAGMA_GENE_RESULTS_OK trait={args.trait_id} genes={gene_rows} "
        f"matched_snps={matched_snps} overlap={overlap_fraction:.6g}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
