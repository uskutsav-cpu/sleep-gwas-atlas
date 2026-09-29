#!/usr/bin/env python3
"""Create isolated source-verified long-sleep pilot inputs; never edit canonical shards."""
from __future__ import annotations

import csv
import gzip
import hashlib
import io
import json
import math
import shutil
import zipfile
from collections import Counter
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
ARCHIVE = Path("/Volumes/Extreme SSD/brain6-work/atlas-raw/.archives/longsumstats.txt.zip")
CANONICAL = Path("/Volumes/Extreme SSD/brain6-work/preparation-v1/lava-inputs-v1")
PROVENANCE = CANONICAL / "provenance.json"
PILOT = ROOT / "brain6/results/lava_confirmatory_pilots_v1"
SELECTION = PILOT / "longsleep_linear_88_loci.tsv"
SELECTION_RECEIPT = PILOT / "longsleep_linear_88_loci.selection.json"
PROTOCOL = ROOT / "brain6/results/lava_confirmatory_protocol_v1/longsleep_linear_model_pilot_addendum_20260927.md"
OUTPUT = Path("/Volumes/Extreme SSD/brain6-work/lava-confirmatory-pilots-v1/longsleep_linear_88_v1")
EXPECTED = {
    ARCHIVE: "0afd9a3ccf459d283f1b6eb6a8b8b0780f54d9d634a05f816be31171a42ec885",
    PROVENANCE: "c50dae720237b548a282d931118f61348821a139f198c64d6d6876a015ce3d02",
    SELECTION: "9182ebca331a7362d96359cc41e130571ce274f7869477864ce7ca0a944a01e7",
    SELECTION_RECEIPT: "3d32855d7eaad85a05e798fc36620f7141b0ce8103051a06d74e4a4885eb7e07",
    PROTOCOL: "dc65468af30f17b3debc9223f94bb44c99a695f3bd2c79797f27fdf8fedd288f",
}
CASES, CONTROLS = 34_184, 305_742
TOTAL_N = CASES + CONTROLS
MIN_INFO, MIN_MAF = 0.95, 0.01
COMPLEMENT = str.maketrans("ACGT", "TGCA")


def sha(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(8 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as stream:
        return list(csv.DictReader(stream, delimiter="\t"))


def source_sign(a1: str, a2: str, raw1: str, raw0: str) -> int | None:
    a1, a2, raw1, raw0 = (x.upper() for x in (a1, a2, raw1, raw0))
    if any(len(x) != 1 or x not in "ACGT" for x in (a1, a2, raw1, raw0)):
        return None
    if a1 == raw1 and a2 == raw0:
        return 1
    if a1 == raw0 and a2 == raw1:
        return -1
    if a1 == raw1.translate(COMPLEMENT) and a2 == raw0.translate(COMPLEMENT):
        return 1
    if a1 == raw0.translate(COMPLEMENT) and a2 == raw1.translate(COMPLEMENT):
        return -1
    return None


def write_sumstats(path: Path, rows: list[dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("xb") as file, gzip.GzipFile(filename="", mode="wb", fileobj=file, mtime=0) as compressed:
        with io.TextIOWrapper(compressed, encoding="utf-8", newline="") as text:
            writer = csv.DictWriter(text, fieldnames=("SNP", "A1", "A2", "Z", "N"),
                                    delimiter="\t", lineterminator="\n")
            writer.writeheader()
            for row in rows:
                writer.writerow({"SNP": row["SNP"], "A1": row["A1"], "A2": row["A2"],
                                 "Z": row["Z"], "N": TOTAL_N})


def main() -> None:
    for path, expected in EXPECTED.items():
        if sha(path) != expected:
            raise ValueError(f"Frozen source hash mismatch: {path}")
    if OUTPUT.exists():
        raise FileExistsError(f"Pilot output already exists: {OUTPUT}")
    selected = tsv(SELECTION)
    selection_receipt = json.loads(SELECTION_RECEIPT.read_text())
    if len(selected) != 88 or selection_receipt["table_sha256"] != sha(SELECTION):
        raise ValueError("Prospective locus selection changed")
    canonical_prov = json.loads(PROVENANCE.read_text())
    old_records = {(str(row["locus_id"]), row.get("trait")): row for row in canonical_prov["records"]
                   if row.get("kind") == "sumstats" and row.get("trait") == "longsleep"}
    old_rows: dict[str, list[dict[str, str]]] = {}
    target: set[str] = set()
    for locus in selected:
        locus_id = locus["LOC"]
        receipt = old_records.get((locus_id, "longsleep"))
        expected_path = CANONICAL / f"locus_{locus_id}/longsleep.sumstats.tsv.gz"
        if not receipt or receipt["path"] != str(expected_path) or sha(expected_path) != receipt["sha256"]:
            raise ValueError(f"Original long-sleep shard receipt mismatch: {locus_id}")
        with gzip.open(expected_path, "rt", newline="", encoding="utf-8") as stream:
            reader = csv.DictReader(stream, delimiter="\t")
            if reader.fieldnames != ["SNP", "A1", "A2", "Z", "N"]:
                raise ValueError(f"Original shard schema changed: {locus_id}")
            rows = list(reader)
        if len(rows) != receipt["rows"] or len({row["SNP"] for row in rows}) != len(rows):
            raise ValueError(f"Original shard row identity changed: {locus_id}")
        old_rows[locus_id] = rows
        target.update(row["SNP"] for row in rows)

    raw_hits: dict[str, list[str] | None] = {}
    raw_rows = 0
    with zipfile.ZipFile(ARCHIVE) as archive:
        names = [item.filename for item in archive.infolist() if not item.is_dir()]
        if names != ["longsumstats.txt"]:
            raise ValueError("Original long-sleep archive members changed")
        with archive.open(names[0]) as raw:
            reader = csv.reader((line.decode("utf-8") for line in raw), delimiter="\t")
            expected_header = ["SNP", "CHR", "BP", "ALLELE1", "ALLELE0", "A1FREQ", "INFO",
                               "BETA_LONGSLEEP", "SE_LONGSLEEP", "P_LONGSLEEP"]
            if next(reader) != expected_header:
                raise ValueError("Original long-sleep source header changed")
            for row in reader:
                raw_rows += 1
                if len(row) != len(expected_header):
                    raise ValueError(f"Malformed original source row {raw_rows}")
                snp = row[0]
                if snp not in target:
                    continue
                if snp in raw_hits:
                    raw_hits[snp] = None  # Duplicate source SNP IDs are quarantined.
                else:
                    raw_hits[snp] = row

    outputs: list[dict] = []
    all_reasons: Counter[str] = Counter()
    prepared: dict[str, list[dict[str, str]]] = {}
    for locus in selected:
        locus_id = locus["LOC"]
        kept: list[dict[str, str]] = []
        reasons: Counter[str] = Counter()
        for row in old_rows[locus_id]:
            source = raw_hits.get(row["SNP"], "missing")
            if source == "missing":
                reasons["source_id_absent"] += 1
                continue
            if source is None:
                reasons["source_id_duplicate"] += 1
                continue
            try:
                chromosome, position = int(source[1]), int(source[2])
                frequency, info = float(source[5]), float(source[6])
                beta, se = float(source[7]), float(source[8])
                z_old = float(row["Z"])
            except (ValueError, OverflowError):
                reasons["source_non_numeric"] += 1
                continue
            if chromosome != int(locus["CHR"]) or not int(locus["START"]) <= position <= int(locus["STOP"]):
                reasons["source_position_mismatch"] += 1
                continue
            if not all(math.isfinite(x) for x in (frequency, info, beta, se, z_old)) or se <= 0:
                reasons["source_non_finite"] += 1
                continue
            if not (0 <= frequency <= 1 and MIN_MAF <= min(frequency, 1 - frequency)
                    and MIN_INFO <= info <= 1):
                reasons["fixed_quality_filter"] += 1
                continue
            sign = source_sign(row["A1"], row["A2"], source[3], source[4])
            if sign is None:
                reasons["allele_mismatch"] += 1
                continue
            if not math.isclose(z_old, sign * beta / se, rel_tol=1e-9, abs_tol=1e-9):
                reasons["source_z_mismatch"] += 1
                continue
            kept.append(row)
        prepared[locus_id] = kept
        all_reasons.update(reasons)
        outputs.append({"locus_id": locus_id, "chromosome": int(locus["CHR"]),
                        "canonical_status": locus["canonical_status"],
                        "original_shard_sha256": old_records[(locus_id, "longsleep")]["sha256"],
                        "original_rows": len(old_rows[locus_id]), "retained_rows": len(kept),
                        "excluded_reasons": dict(reasons)})
    if raw_rows != 14_661_601:
        raise ValueError("Original source row count changed")

    OUTPUT.mkdir(parents=True, exist_ok=False)
    loci_path = OUTPUT / "selected.loci"
    with loci_path.open("x", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=("LOC", "CHR", "START", "STOP"),
                                delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for row in selected:
            writer.writerow({key: row[key] for key in ("LOC", "CHR", "START", "STOP")})
    for output in outputs:
        locus_id = output["locus_id"]
        for arm in ("linear", "binary"):
            directory = OUTPUT / "inputs" / arm / f"locus_{locus_id}"
            shard = directory / "longsleep.sumstats.tsv.gz"
            if arm == "linear":
                write_sumstats(shard, prepared[locus_id])
            else:
                directory.mkdir(parents=True, exist_ok=False)
                shutil.copyfile(OUTPUT / "inputs/linear" / f"locus_{locus_id}" / shard.name, shard)
            info_path = directory / "input_info.tsv"
            cases, controls = ("NA", "NA") if arm == "linear" else (CASES, CONTROLS)
            info_path.write_text(f"phenotype\tcases\tcontrols\tfilename\n"
                                 f"longsleep\t{cases}\t{controls}\t{shard.name}\n")
            output[f"{arm}_shard_sha256"] = sha(shard)
            output[f"{arm}_info_sha256"] = sha(info_path)
        if output["linear_shard_sha256"] != output["binary_shard_sha256"]:
            raise ValueError("Primary and secondary arm SNP/Z/N files differ")
    receipt = {
        "analysis_id": "brain6_longsleep_linear_88_pilot_v1",
        "scope": "PREDECLARED_TRAIT_ONLY_DIAGNOSTIC_NO_PROMOTION",
        "source_sha256": {str(path): digest for path, digest in EXPECTED.items()},
        "materializer_sha256": sha(Path(__file__)),
        "loci_file_sha256": sha(loci_path),
        "original_source_rows_scanned": raw_rows,
        "target_source_ids": len(target),
        "matched_source_ids": sum(row is not None for row in raw_hits.values()),
        "duplicated_source_ids": sum(row is None for row in raw_hits.values()),
        "fixed_filters": {"minimum_info": MIN_INFO, "minimum_maf": MIN_MAF,
                          "N_proxy": TOTAL_N, "primary_model": "linear_BOLT_LMM",
                          "secondary_model": "binary_reconstruction"},
        "excluded_reasons_all_selected_loci": dict(all_reasons),
        "selected_loci": outputs,
        "canonical_modified": False,
        "lava_run_started": False,
    }
    receipt_path = OUTPUT / "materialization.receipt.json"
    receipt_path.write_text(json.dumps(receipt, indent=2, sort_keys=True) + "\n")
    print(json.dumps({"status": "MATERIALIZED", "loci": len(outputs),
                      "retained_rows": sum(row["retained_rows"] for row in outputs),
                      "receipt_sha256": sha(receipt_path)}, sort_keys=True))


if __name__ == "__main__":
    main()
