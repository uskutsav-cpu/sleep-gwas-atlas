#!/usr/bin/env python3
"""Create partial PLACO candidate loci using the frozen UKB/LAVA LD rule."""
from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
import math
import os
import subprocess
import tempfile
from collections import defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
RULE_PATH = ROOT / "brain6/config/shared_locus_rule_v1.json"
FAMILY_PATH = ROOT / "brain6/config/placo_family_v3/family_lock.json"
RUN_PATH = ROOT / "brain6/manifests/placo_family_v3_run.json"
MASTER_PATH = ROOT / "brain6/results/placo/placo_master.tsv"
VARIANT_COLUMNS = ["pair_id", "SNP", "CHR", "BP", "P_PLACO", "reference_status",
                   "candidate_status", "lead_SNP", "r2_to_lead", "locus_id"]
LOCUS_COLUMNS = ["locus_id", "pair_id", "sleep_trait", "brain_disorder", "CHR", "START", "STOP",
                 "lead_variant", "lead_P_PLACO", "n_lead_signals", "lead_variants",
                 "n_candidate_variants", "candidate_variants", "r2_threshold", "window_kb",
                 "reference", "locus_status", "interpretation"]


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def greedy_clumps(rows: list[dict], r2: dict[tuple[str, str], float], window_bp: int,
                  r2_threshold: float) -> list[dict]:
    """Return deterministic PLINK-style index clumps among exact-ref-matched hits."""
    by_chr: dict[int, list[dict]] = defaultdict(list)
    for row in rows:
        by_chr[int(row["CHR"])].append(row)
    clumps = []
    for chr_id in sorted(by_chr):
        pending = sorted(by_chr[chr_id], key=lambda x: (float(x["P_PLACO"]), int(x["BP"]), x["SNP"]))
        while pending:
            lead = pending[0]
            members = []
            remaining = []
            for candidate in pending:
                distance = abs(int(candidate["BP"]) - int(lead["BP"]))
                key = tuple(sorted((lead["SNP"], candidate["SNP"])))
                corr = float(r2.get(key, 1.0 if lead["SNP"] == candidate["SNP"] else 0.0))
                if distance <= window_bp and corr * corr >= r2_threshold:
                    members.append((candidate, corr * corr))
                else:
                    remaining.append(candidate)
            clumps.append({"chr": chr_id, "lead": lead, "members": members})
            pending = remaining
    return clumps


def merge_overlapping_regions(clumps: list[dict], flank_bp: int) -> list[dict]:
    """Merge overlapping same-pair intervals while retaining every LD-clump lead."""
    by_chr: dict[int, list[dict]] = defaultdict(list)
    for clump in clumps:
        bp = int(clump["lead"]["BP"])
        by_chr[clump["chr"]].append({"start": max(1, bp - flank_bp), "stop": bp + flank_bp,
                                     "clumps": [clump]})
    merged = []
    for chr_id in sorted(by_chr):
        intervals = sorted(by_chr[chr_id], key=lambda x: (x["start"], x["stop"]))
        for interval in intervals:
            if merged and merged[-1]["chr"] == chr_id and interval["start"] <= merged[-1]["stop"]:
                merged[-1]["stop"] = max(merged[-1]["stop"], interval["stop"])
                merged[-1]["clumps"].extend(interval["clumps"])
            else:
                merged.append({"chr": chr_id, **interval})
    return merged


def read_candidates(master: Path, threshold: float) -> tuple[dict[str, list[dict]], dict[str, str]]:
    candidates: dict[str, list[dict]] = {}
    input_hashes: dict[str, str] = {}
    with master.open(newline="", encoding="utf-8") as stream:
        rows = list(csv.DictReader(stream, delimiter="\t"))
    if not rows:
        raise ValueError("No published PLACO pair tables are available")
    for record in rows:
        if not record["stage_status"].startswith("COMPLETE_QC_PASS"):
            raise ValueError(f"Pair is not pair-QC passed: {record['pair_id']}")
        path = ROOT / record["output_path"]
        if sha256(path) != record["output_sha256"]:
            raise ValueError(f"Published variant hash differs from master: {path}")
        pair_rows = []
        seen: set[str] = set()
        with gzip.open(path, "rt", newline="", encoding="utf-8") as source:
            reader = csv.DictReader(source, delimiter="\t")
            required = {"SNP", "CHR", "BP", "P_PLACO", "status", "headline"}
            if not required.issubset(reader.fieldnames or []):
                raise ValueError(f"Unexpected published variant schema: {path}")
            for row in reader:
                if row["status"] != "TESTED":
                    continue
                p = float(row["P_PLACO"])
                if not math.isfinite(p) or not 0 <= p <= 1:
                    raise ValueError(f"Invalid PLACO P value in {path}")
                if row["headline"] != str(p < threshold):
                    raise ValueError(f"Headline flag disagrees with the frozen threshold in {path}: {row['SNP']}")
                if p >= threshold:
                    continue
                snp = row["SNP"].lower()
                chr_id, bp = int(row["CHR"]), int(row["BP"])
                if snp in seen:
                    raise ValueError(f"Duplicate candidate SNP in {path}: {snp}")
                seen.add(snp)
                pair_rows.append({"pair_id": record["pair_id"], "sleep_trait": record["sleep_trait"],
                                  "brain_disorder": record["brain_disorder"], "SNP": snp,
                                  "CHR": chr_id, "BP": bp, "P_PLACO": p})
        candidates[record["pair_id"]] = pair_rows
        input_hashes[record["pair_id"]] = record["output_sha256"]
    return candidates, input_hashes


def extract_reference_ld(candidates: list[dict], reference_prefix: Path, rscript: Path) -> tuple[dict, dict, int]:
    with tempfile.TemporaryDirectory(prefix="brain6-placo-ld-") as temp:
        tmp = Path(temp)
        candidate_path = tmp / "candidate_snps.tsv"
        with candidate_path.open("w", newline="", encoding="utf-8") as stream:
            writer = csv.DictWriter(stream, fieldnames=["SNP", "CHR", "BP"], delimiter="\t", lineterminator="\n")
            writer.writeheader()
            unique = {}
            for row in candidates:
                key = (row["CHR"], row["SNP"])
                prior = unique.get(key)
                if prior and prior["BP"] != row["BP"]:
                    raise ValueError(f"Same SNP has conflicting positions: {key}")
                unique[key] = row
            writer.writerows({k: row[k] for k in ("SNP", "CHR", "BP")} for row in unique.values())
        output = tmp / "ld"
        helper = ROOT / "brain6/scripts/extract_placo_ld.R"
        proc = subprocess.run([str(rscript), str(helper), str(reference_prefix), str(candidate_path), str(output)],
                              capture_output=True, text=True, check=False)
        if proc.returncode:
            raise RuntimeError("Pinned LAVA LD extraction failed: " + (proc.stdout + proc.stderr)[-4000:])
        ref_info = {}
        with (output / "reference_matches.tsv").open(newline="", encoding="utf-8") as stream:
            for row in csv.DictReader(stream, delimiter="\t"):
                ref_info[(int(row["CHR"]), row["SNP"].lower())] = {
                    "BP": int(row["POS"]), "A1": row["A1"], "A2": row["A2"]}
        ld = {}
        with (output / "ld_pairs.tsv").open(newline="", encoding="utf-8") as stream:
            for row in csv.DictReader(stream, delimiter="\t"):
                a, b = row["SNP1"].lower(), row["SNP2"].lower()
                ld[tuple(sorted((a, b)))] = float(row["R"])
        ref_n = int((output / "reference_sample_size.txt").read_text().strip())
        return ref_info, ld, ref_n


def build_outputs(candidates: dict[str, list[dict]], ref_info: dict,
                  ld: dict[tuple[str, str], float], rule: dict) -> tuple[list[dict], list[dict]]:
    clump_cfg = rule["clumping"]
    window = int(clump_cfg["window_kb"]) * 1000
    threshold = float(clump_cfg["r2_threshold"])
    flank = int(rule["region_definition"]["lead_flank_kb"]) * 1000
    variant_out, locus_out = [], []
    item_map = {}
    for pair in sorted(candidates):
        rows = candidates[pair]
        mapped, missing = [], []
        for row in rows:
            reference = ref_info.get((row["CHR"], row["SNP"]))
            item = {**row, "reference_status": "EXACT_MATCH", "candidate_status": "NOT_ASSIGNED",
                    "lead_SNP": "NA", "r2_to_lead": "NA", "locus_id": "NA"}
            if reference is None:
                item["reference_status"] = "NO_REFERENCE_VARIANT"
                item["candidate_status"] = "NOT_CLUMPED_NO_EXACT_REFERENCE_MATCH"
                missing.append(item)
            elif reference["BP"] != row["BP"]:
                item["reference_status"] = "POSITION_MISMATCH"
                item["candidate_status"] = "NOT_CLUMPED_NO_EXACT_REFERENCE_MATCH"
                missing.append(item)
            else:
                item["A1"], item["A2"] = reference["A1"], reference["A2"]
                mapped.append(row)
                item["_row"] = row
                item_map[(pair, row["SNP"])] = item
        clumps = greedy_clumps(mapped, ld, window, threshold)
        regions = merge_overlapping_regions(clumps, flank)
        for region in regions:
            signals = sorted((c["lead"] for c in region["clumps"]),
                             key=lambda x: (float(x["P_PLACO"]), int(x["BP"]), x["SNP"]))
            all_members = {}
            for clump in region["clumps"]:
                for member, r2 in clump["members"]:
                    all_members[member["SNP"]] = (member, r2, clump["lead"]["SNP"])
            lead = signals[0]
            locus_id = f"placo_{pair}_chr{region['chr']}_{region['start']}_{region['stop']}"
            for snp, (member, r2, lead_snp) in all_members.items():
                item = item_map[(pair, snp)]
                item["candidate_status"] = "LEAD" if snp == lead_snp else "LD_CLUMPED"
                item["lead_SNP"] = lead_snp
                item["r2_to_lead"] = f"{r2:.8g}"
                item["locus_id"] = locus_id
            locus_out.append({
                "locus_id": locus_id, "pair_id": pair, "sleep_trait": lead["sleep_trait"],
                "brain_disorder": lead["brain_disorder"], "CHR": region["chr"],
                "START": region["start"], "STOP": region["stop"], "lead_variant": lead["SNP"],
                "lead_P_PLACO": lead["P_PLACO"], "n_lead_signals": len(signals),
                "lead_variants": ";".join(x["SNP"] for x in signals),
                "n_candidate_variants": len(all_members),
                "candidate_variants": ";".join(sorted(all_members)),
                "r2_threshold": threshold, "window_kb": clump_cfg["window_kb"],
                "reference": clump_cfg["reference"], "locus_status": "PLACO_ONLY_CANDIDATE_LOCUS",
                "interpretation": "No local-rg, fine-mapping, replication, or functional evidence included",
            })
        variant_out.extend({k: v for k, v in item.items() if k in VARIANT_COLUMNS}
                           for item in item_map.values() if item.get("pair_id") == pair)
        variant_out.extend({k: v for k, v in item.items() if k in VARIANT_COLUMNS} for item in missing)
    return sorted(variant_out, key=lambda x: (x["pair_id"], int(x["CHR"]), int(x["BP"]), x["SNP"])), \
        sorted(locus_out, key=lambda x: (x["pair_id"], int(x["CHR"]), int(x["START"])))


def atomic_tsv(path: Path, fields: list[str], rows: list[dict]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fd, name = tempfile.mkstemp(prefix=f".{path.name}.", dir=path.parent)
    os.close(fd)
    temp = Path(name)
    try:
        with temp.open("w", newline="", encoding="utf-8") as stream:
            writer = csv.DictWriter(stream, fieldnames=fields, delimiter="\t", lineterminator="\n", extrasaction="ignore")
            writer.writeheader()
            writer.writerows(rows)
        os.replace(temp, path)
    finally:
        temp.unlink(missing_ok=True)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--reference-prefix", type=Path,
                        default=Path("/Volumes/Extreme SSD/brain6-work/lava-ukb-v1.1/lava-ukb-v1.1"))
    parser.add_argument("--rscript", type=Path,
                        default=Path("/Volumes/Extreme SSD/brain6-work/r-env-4.3.3-short/bin/Rscript"))
    args = parser.parse_args()
    rule = json.loads(RULE_PATH.read_text(encoding="utf-8"))
    family = json.loads(FAMILY_PATH.read_text(encoding="utf-8"))
    run = json.loads(RUN_PATH.read_text(encoding="utf-8"))
    candidates, input_hashes = read_candidates(MASTER_PATH, float(family["headline_threshold"]))
    all_candidates = [row for pair in sorted(candidates) for row in candidates[pair]]
    if not all_candidates:
        raise SystemExit("No PLACO headline variants in published QC-passing pair outputs")
    ref_info, ld, ref_n = extract_reference_ld(all_candidates, args.reference_prefix, args.rscript)
    if ref_n != int(rule["clumping"]["reference_sample_size"]):
        raise ValueError(f"LAVA reference N mismatch: {ref_n}")
    variants, loci = build_outputs(candidates, ref_info, ld, rule)
    output = ROOT / "brain6/results/loci"
    variant_path, locus_path = output / "placo_candidate_variants_partial.tsv", output / "placo_candidate_loci_partial.tsv"
    atomic_tsv(variant_path, VARIANT_COLUMNS, variants)
    atomic_tsv(locus_path, LOCUS_COLUMNS, loci)
    live_pairs = sorted(candidates)
    family_complete = (all(run["pairs"][p]["status"] == "COMPLETE" for p in run["pairs"])
                       and "insomnia__adhd" in input_hashes)
    reference_provenance = args.reference_prefix.parent / "reference.provenance.json"
    provenance = {
        "status": "PASS_PARTIAL_FAMILY" if not family_complete else "PASS",
        "rule_id": rule["rule_id"], "rule_sha256": sha256(RULE_PATH),
        "placo_family_lock_sha256": sha256(FAMILY_PATH),
        "placo_pair_outputs": input_hashes, "pairs_included": live_pairs,
        "required_v3_pairs": sorted(run["pairs"]), "protected_track_b_present": "insomnia__adhd" in input_hashes,
        "family_complete": family_complete, "candidate_threshold": family["headline_threshold"],
        "reference_name": rule["clumping"]["reference"], "reference_sample_size": ref_n,
        "reference_provenance_path": str(reference_provenance),
        "reference_provenance_sha256": sha256(reference_provenance),
        "reference_prefix": str(args.reference_prefix), "runtime_version": "LAVA 0.1.5",
        "runtime_binary_sha256": sha256(args.rscript), "r_extractor_sha256": sha256(ROOT / "brain6/scripts/extract_placo_ld.R"),
        "clumper_sha256": sha256(Path(__file__)),
        "n_candidate_variants": len(variants), "n_candidate_loci": len(loci),
        "n_exact_reference_matched_variants": sum(x["reference_status"] == "EXACT_MATCH" for x in variants),
        "n_unmatched_or_position_mismatched": sum(x["reference_status"] != "EXACT_MATCH" for x in variants),
        "outputs": [{"path": str(path.relative_to(ROOT)), "bytes": path.stat().st_size, "sha256": sha256(path)}
                    for path in (variant_path, locus_path)],
        "interpretation": "Partial PLACO-only candidate regions; not final independent/shared causal loci",
    }
    prov_path = output / "placo_candidate_loci_partial.provenance.json"
    prov_path.write_text(json.dumps(provenance, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(provenance, indent=2))


if __name__ == "__main__":
    main()
