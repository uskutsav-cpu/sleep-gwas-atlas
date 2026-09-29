#!/usr/bin/env python3
"""Publish checksummed QC summaries from the verified Brain6 prep receipts."""
from __future__ import annotations

import argparse
import csv
import gzip
import io
import json
from pathlib import Path

from brain6.artifacts import verify_artifact
from brain6.io import atomic_text, sha256

ROOT = Path(__file__).resolve().parents[2]
TRAITS = ("adhd", "bipolar", "insomnia", "longsleep", "mdd", "parkinson", "scz")
PAIRS = ("insomnia__adhd", "insomnia__mdd", "longsleep__bipolar",
         "longsleep__parkinson", "longsleep__scz")
AVAILABILITY_FIELDS = (
    "source_archive_path", "source_archive_container_bytes", "source_archive_container_SHA256",
    "source_archive_content_transform", "source_archive_content_SHA256",
    "harmonized_input_path", "harmonized_input_SHA256", "source_archive_status",
)


def load_verified(root: Path, name: str, stage: str) -> tuple[dict, dict]:
    path = root / name
    receipt = verify_artifact(path)
    if receipt.get("stage") != stage or receipt.get("status") != "COMPLETE" or receipt.get("synthetic"):
        raise ValueError(f"Expected completed empirical {stage} receipt: {path}")
    qc_path = path / "qc.json"
    qc = json.loads(qc_path.read_text(encoding="utf-8"))
    output = {item["path"]: item["sha256"] for item in receipt["outputs"]}
    if output.get("qc.json") != sha256(qc_path):
        raise ValueError(f"QC file checksum differs from receipt: {qc_path}")
    return receipt, qc


def write_markdown(path: Path, row: dict) -> None:
    lines = [f"# Dense GWAS QC: {row['trait_id']}", "",
             "This report summarizes source normalization only; it is not an association result.", "",
             f"- Input: `{row['source_path']}`",
             f"- Input SHA256: `{row['source_sha256']}`",
             f"- Receipt SHA256: `{row['receipt_sha256']}`",
             f"- Genome build: {row['genome_build']}; ancestry: {row['ancestry']}",
             f"- Rows read / retained: {row['input_rows']} / {row['retained_rows']}",
             f"- Duplicate IDs: {row['duplicate_ids']}; ambiguous positions: {row['ambiguous_positions']}",
             f"- Source columns: `{row['source_header']}`",
             f"- EAF mapping: {row['frequency_field_status']}; INFO: {row['info_field_status']}.",
             f"- Brain6 normalization did not map an EAF field ({row['missing_eaf_seen']} rows) or INFO ({row['missing_info_seen']} rows), so its local MAF/INFO filters were not applied.",
             "- Source-level allele-frequency filtering is inherited only where an upstream, source-bound QC record documents it; do not infer it from the unused FRQ column.",
             "- This dataset-level QC does not establish participant-level cohort overlap.", ""]
    write_immutable_text(path, "\n".join(lines))


def write_immutable_text(path: Path, text: str) -> None:
    if path.exists():
        if path.read_text(encoding="utf-8") != text:
            raise ValueError(f"Existing QC report differs from verified inputs: {path}")
        return
    with atomic_text(path, immutable=True) as stream:
        stream.write(text)


def write_immutable_tsv(path: Path, fields: list[str], rows: list[dict]) -> None:
    buffer = io.StringIO()
    writer = csv.DictWriter(buffer, fieldnames=fields, delimiter="\t", lineterminator="\n")
    writer.writeheader(); writer.writerows(rows)
    write_immutable_text(path, buffer.getvalue())


def write_atomic_tsv(path: Path, fields: list[str], rows: list[dict]) -> None:
    with atomic_text(path, immutable=False) as stream:
        writer = csv.DictWriter(stream, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader(); writer.writerows(rows)


def compose_availability_rows(manifest_rows: list[dict[str, str]],
                              dense_by_trait: dict[str, dict[str, str]],
                              archive_by_trait: dict[str, dict[str, str]]) -> list[dict[str, str]]:
    """Build an external provenance overlay without modifying the run-bound ledger."""
    output = []
    for source_row in manifest_rows:
        row = dict(source_row)
        for field in AVAILABILITY_FIELDS:
            row[field] = "NA"
        trait = row["trait"]
        dense = dense_by_trait.get(trait)
        archive = archive_by_trait.get(trait)
        dense_status = "DENSE_INPUT_HASH_VERIFIED" if dense is not None else "DENSE_INPUT_UNAVAILABLE"
        if dense is not None:
            row.update({
                "variant_count": dense["retained_rows"],
                "beta": "YES_IN_DENSE_SOURCE", "SE": "YES_IN_DENSE_SOURCE", "P": "YES_IN_DENSE_SOURCE",
                "EAF": ("FRQ_PRESENT; EAF_SEMANTICS_NOT_MAPPED_IN_BRAIN6_SOURCE_CARD"
                        if "FRQ" in dense["source_header"].split(",")
                        else "NOT_SUPPLIED_IN_DENSE_SOURCE"),
                "per_variant_N": "YES_IN_DENSE_SOURCE",
                "INFO": "NOT_SUPPLIED_IN_DENSE_SOURCE",
                "harmonized_input_path": dense["source_path"],
                "harmonized_input_SHA256": dense["source_sha256"],
            })
        if archive is not None:
            status = archive["status"]
            row.update({
                "source_archive_path": archive["archive_path"],
                "source_archive_container_bytes": archive["container_bytes"],
                "source_archive_container_SHA256": archive["container_SHA256"],
                "source_archive_content_transform": archive["content_transform"],
                "source_archive_content_SHA256": archive["observed_source_SHA256"],
                "source_archive_status": status,
            })
            if status in {"VERIFIED_EXACT_ARCHIVE", "VERIFIED_DECOMPRESSED_SOURCE"}:
                row["raw_file_path"] = archive["archive_path"]
                row["raw_file_SHA256"] = archive["container_SHA256"]
                row["current_readiness"] = f"SOURCE_CONTENT_HASH_VERIFIED; {dense_status}; ACQUISITION_DATE_UNKNOWN"
            elif status == "MISMATCH":
                raise ValueError(f"Raw source archive failed checksum validation: {trait}")
            else:
                row["current_readiness"] = f"SOURCE_ARCHIVE_MISSING; {dense_status}; ACQUISITION_DATE_UNKNOWN"
        output.append(row)
    return output


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--preparation-root", type=Path, required=True)
    args = ap.parse_args()
    prep = args.preparation_root.resolve()
    out = ROOT / "brain6/qc"
    out.mkdir(parents=True, exist_ok=True)
    audit_path = ROOT / "brain6/manifests/locked_dense_input_audit.tsv"
    with audit_path.open(encoding="utf-8", newline="") as stream:
        audit = {row["trait_id"]: row for row in csv.DictReader(stream, delimiter="\t")}
    source_rows = []
    for trait in TRAITS:
        receipt, qc = load_verified(prep, f"normalize_{trait}", "normalize")
        if trait not in audit or audit[trait]["status"] != "VERIFIED":
            raise ValueError(f"No checksum-verified source audit row for {trait}")
        source = receipt["source"]
        raw = next(item for item in receipt["inputs"] if item["path"].endswith(".tsv.gz"))
        if raw["sha256"] != audit[trait]["observed_sha256"]:
            raise ValueError(f"Prep receipt and source audit disagree for {trait}")
        with gzip.open(raw["path"], "rt", encoding="utf-8", newline="") as stream:
            source_header = next(csv.reader(stream, delimiter="\t"))
        source_columns = set(source_header)
        required_columns = {"SNP", "CHR", "BP", "A1", "A2", "BETA", "SE", "P", "N"}
        if not required_columns <= source_columns:
            raise ValueError(f"Dense source schema is incomplete for {trait}: {sorted(required_columns-source_columns)}")
        row = {
            "trait_id": trait, "study_id": source["study_id"], "genome_build": source["genome_build"],
            "ancestry": source["ancestry"], "input_rows": qc["input_rows"],
            "retained_rows": qc["retained_rows"], "duplicate_ids": qc["duplicate_ids"],
            "ambiguous_positions": qc["ambiguous_positions"],
            "ambiguous_position_rows": qc["ambiguous_position_rows"],
            "missing_eaf_seen": qc.get("missing_eaf_seen", 0),
            "missing_info_seen": qc.get("missing_info_seen", 0),
            "source_header": ",".join(source_header),
            "frequency_field_status": ("FRQ_PRESENT; EAF_SEMANTICS_NOT_MAPPED_IN_BRAIN6_SOURCE_CARD"
                                        if "FRQ" in source_header else "NO_FREQUENCY_FIELD_IN_SOURCE_HEADER"),
            "info_field_status": ("INFO_PRESENT_IN_SOURCE_HEADER" if "INFO" in source_header
                                  else "NO_INFO_FIELD_IN_SOURCE_HEADER"),
            "min_maf_policy": receipt["parameters"]["min_maf"],
            "min_info_policy": receipt["parameters"]["min_info"],
            "eaf_required": receipt["parameters"]["require_eaf"],
            "info_required": receipt["parameters"]["require_info"],
            "source_path": raw["path"], "source_sha256": raw["sha256"],
            "source_receipt": str(prep / f"normalize_{trait}/receipt.json"),
            "receipt_sha256": sha256(prep / f"normalize_{trait}/receipt.json"),
            "normalized_sumstats_sha256": next(x["sha256"] for x in receipt["outputs"] if x["path"] == "sumstats.tsv.gz"),
            "status": "PASS",
        }
        source_rows.append(row)
        write_immutable_tsv(out / f"{trait}_dense_qc.tsv", list(row), [row])
        write_markdown(out / f"{trait}_dense_qc.md", row)

    pair_rows = []
    for pair in PAIRS:
        receipt, qc = load_verified(prep, f"pair_{pair}", "join_pair")
        pair_rows.append({"pair_id": pair, **qc,
                          "pair_receipt": str(prep / f"pair_{pair}/receipt.json"),
                          "receipt_sha256": sha256(prep / f"pair_{pair}/receipt.json"),
                          "pair_file_sha256": next(x["sha256"] for x in receipt["outputs"] if x["path"] == "pair.tsv.gz"),
                          "status": "PASS"})
    write_immutable_tsv(out / "dense_source_qc_master.tsv", list(source_rows[0]), source_rows)
    pair_fields = list(dict.fromkeys(key for row in pair_rows for key in row))
    write_immutable_tsv(out / "dense_pair_qc.tsv", pair_fields, pair_rows)

    # Preserve the source ledger SHA bound to the active LAVA materialization.
    # Publish current external availability as a separate derived overlay.
    manifest_path = ROOT / "brain6/manifests/gwas_master.tsv"
    rows = list(csv.DictReader(manifest_path.open(encoding="utf-8", newline=""), delimiter="\t"))
    source_by_trait = {r["trait_id"]: r for r in source_rows}
    for row in rows:
        qc_row = source_by_trait.get(row["trait"])
        if qc_row is None:
            continue
        row.update({
            "variant_count": qc_row["retained_rows"],
            "beta": "YES_IN_DENSE_SOURCE", "SE": "YES_IN_DENSE_SOURCE", "P": "YES_IN_DENSE_SOURCE",
            "EAF": ("FRQ_PRESENT; EAF_SEMANTICS_NOT_MAPPED_IN_BRAIN6_SOURCE_CARD"
                    if "FRQ" in source_by_trait[row["trait"]]["source_header"].split(",")
                    else "NOT_SUPPLIED_IN_DENSE_SOURCE"),
            "per_variant_N": "YES_IN_DENSE_SOURCE",
            "INFO": "NOT_SUPPLIED_IN_DENSE_SOURCE",
            "raw_file_path": qc_row["source_path"], "raw_file_SHA256": qc_row["source_sha256"],
            "current_readiness": "AVAILABLE_EXTERNAL; CHECKSUM_VERIFIED; PREPARATION_RECEIPTS_VERIFIED",
        })
    availability_path = ROOT / "brain6/manifests/gwas_external_availability.tsv"
    archive_path = ROOT / "brain6/qc/source_archive_audit.tsv"
    with archive_path.open(encoding="utf-8", newline="") as stream:
        archive_by_trait = {row["trait"]: row for row in csv.DictReader(stream, delimiter="\t")}
    availability_rows = compose_availability_rows(rows, source_by_trait, archive_by_trait)
    availability_fields = list(rows[0]) + [field for field in AVAILABILITY_FIELDS if field not in rows[0]]
    write_atomic_tsv(availability_path, availability_fields, availability_rows)

    table_path = ROOT / "brain6/results/supplement/table_S1_gwas_metadata.tsv"
    write_atomic_tsv(table_path, availability_fields, availability_rows)
    provenance = {
        "status": "PASS",
        "data_rows": len(availability_rows),
        "locked_gwas_master_sha256": sha256(manifest_path),
        "locked_dense_input_audit_sha256": sha256(ROOT / "brain6/manifests/locked_dense_input_audit.tsv"),
        "source_archive_audit_sha256": sha256(archive_path),
        "dense_source_qc_master_sha256": sha256(out / "dense_source_qc_master.tsv"),
        "external_availability_sha256": sha256(availability_path),
        "output_table_sha256": sha256(table_path),
        "interpretation": "External availability overlay; the run-bound GWAS master was not modified.",
    }
    write_immutable_text(
        ROOT / "brain6/results/supplement/table_S1_gwas_metadata.provenance.json",
        json.dumps(provenance, indent=2, sort_keys=True) + "\n",
    )

    print(json.dumps({"traits": len(source_rows), "pairs": len(pair_rows),
                      "qc_directory": str(out), "locked_gwas_ledger": str(manifest_path),
                      "external_availability_overlay": str(availability_path)}, indent=2))


if __name__ == "__main__":
    main()
