#!/usr/bin/env python3
"""Read reference metadata only; no numerical imports, genotype decode or fits."""
from __future__ import annotations

import csv
import gzip
import hashlib
import itertools
import json
import math
import resource
import struct
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
V4 = ROOT / "sleep_unified_research_v4"
OUT = V4 / "statistical_validation/canonical_realistic_ld_feasibility_probe_receipt_v4.json"
GENO = Path("/Volumes/Extreme SSD/brain6-work/brain6-alternative-local-validation-v1/reference")
LAVA = Path("/Volumes/Extreme SSD/Utsav-Research-Archive/Sleep-GWAS/LAVA-UKB-v1.1-repair/extracted")
SOURCES = V4 / "statistical_validation/canonical_ld_feasibility_sources_v4"


def sha(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for b in iter(lambda: f.read(1024 * 1024), b""):
            h.update(b)
    return h.hexdigest()


def ident(path: Path, hash_content=True) -> dict:
    s = path.stat()
    return {"path": str(path), "bytes": s.st_size, "mtime_ns": s.st_mtime_ns,
            "inode": s.st_ino, "device": s.st_dev,
            "sha256": sha(path) if hash_content else None,
            "content_hash_verified": hash_content}


def main() -> None:
    if OUT.exists():
        raise SystemExit("preserve existing receipt; choose a distinct version")
    preparation = V4 / "manifests/canonical_200_interval_preparation_v4.json"
    j = json.loads(preparation.read_text())
    interval = Path(j["interval_tsv_path"])
    coordinate = Path(j["coordinate_map_path"])
    assert sha(interval) == j["interval_tsv_sha256"]
    assert sha(coordinate) == j["coordinate_map_sha256"]
    watched = [preparation, interval, coordinate, GENO / "1000G_EUR.bim",
               GENO / "1000G_EUR.fam", GENO / "prepared_v1/reference_prep_provenance.json",
               GENO / "per_chr_v1/split_reference_provenance.json", LAVA / "VERSION_INFO",
               LAVA / "lava-ukb-v1.1_chr22.info",
               SOURCES / "integrated_call_samples_v3.20130502.ALL.panel"]
    watched.extend(Path(p) for p in j["reference_sha256"])
    before = {str(p): ident(p) for p in watched}
    for p, expected in j["reference_sha256"].items():
        assert before[p]["sha256"] == expected
    prep = json.loads((GENO / "prepared_v1/reference_prep_provenance.json").read_text())
    for name in ["1000G_EUR.bim", "1000G_EUR.fam"]:
        assert before[str(GENO / name)]["sha256"] == prep["source_sha256"][name]

    panel = {r["sample"]: r for r in csv.DictReader(
        (SOURCES / "integrated_call_samples_v3.20130502.ALL.panel").open(), delimiter="\t")}
    fam = [s.split() for s in (GENO / "1000G_EUR.fam").read_text().splitlines()]
    samples = [r[1] for r in fam]
    assert len(samples) == len(set(samples)) == prep["n_samples"]
    ancestry = Counter(panel[s]["pop"] for s in samples if s in panel)
    ancestry_groups = Counter(panel[s]["super_pop"] for s in samples if s in panel)
    allele_rows = csv.DictReader(gzip.open(coordinate, "rt"), delimiter="\t")
    allele_groups = itertools.groupby(allele_rows, key=lambda r: int(r["CHR"]))
    allele_iterator = iter(allele_groups)
    current_allele_chr, current_allele_group = next(allele_iterator)
    with interval.open() as f:
        intervals = list(csv.DictReader(f, delimiter="\t"))
    # Actual builder schema is read, rather than inferred from positional fields.
    assert set(intervals[0]) == {"block_id", "chromosome", "start_inclusive", "end_exclusive", "construction_reference_snps"}
    block_counts = [0] * 200
    raw_seen = set()
    per_chr = []
    previous = (0, 0)
    counters = Counter()

    def raw_rows():
        nonlocal previous
        with (GENO / "1000G_EUR.bim").open() as f:
            for line in f:
                r = line.split()
                assert len(r) == 6
                c, bp = int(r[0]), int(r[3])
                assert 1 <= c <= 22 and (c, bp) >= previous
                previous = (c, bp)
                assert r[1] not in raw_seen
                raw_seen.add(r[1])
                counters["raw_variants"] += 1
                counters["zero_cm_rows"] += float(r[2]) == 0
                counters["non_acgt_or_equal_alleles"] += not (
                    r[4] in "ACGT" and len(r[4]) == 1 and r[5] in "ACGT" and len(r[5]) == 1 and r[4] != r[5])
                yield c, r[1], bp, r[4], r[5]

    for c, group in itertools.groupby(raw_rows(), key=lambda r: r[0]):
        raw = {r[1]: r[2:] for r in group}
        assert current_allele_chr == c
        alleles = {r["SNP"]: (int(r["BP"]), r["A1"], r["A2"]) for r in current_allele_group}
        try:
            current_allele_chr, current_allele_group = next(allele_iterator)
        except StopIteration:
            current_allele_chr = 23
        reference_path = next(Path(p) for p in j["reference_sha256"] if Path(p).name == f"{c}.l2.ldscore.gz")
        with gzip.open(reference_path, "rt") as f:
            reference = list(csv.DictReader(f, delimiter="\t"))
        count = Counter()
        for r in reference:
            count["consumed_reference_rows"] += 1
            snp, bp = r["SNP"], int(r["BP"])
            if snp not in raw:
                count["absent_from_raw_genotype"] += 1
                continue
            count["ID_present"] += 1
            rbp, a1, a2 = raw[snp]
            if rbp != bp:
                count["coordinate_conflict"] += 1
                continue
            count["ID_coordinate_match"] += 1
            if snp not in alleles:
                count["allele_map_absent"] += 1
                continue
            mapbp, ma1, ma2 = alleles[snp]
            if mapbp != bp or {a1, a2} != {ma1, ma2}:
                count["exact_allele_set_conflict"] += 1
                continue
            count["exact_coordinate_allele_set_match"] += 1
            count["palindromic_exact_match"] += {a1, a2} in ({"A", "T"}, {"C", "G"})
            for b in intervals:
                if int(b["chromosome"]) == c and int(b["start_inclusive"]) <= bp < int(b["end_exclusive"]):
                    block_counts[int(b["block_id"])] += 1
                    break
            else:
                raise AssertionError("matched reference SNP outside fixed intervals")
        per_chr.append({"CHR": c, "raw_variants": len(raw), **dict(count)})

    assert counters["raw_variants"] == prep["n_input_variants"]
    bed = GENO / "1000G_EUR.bed"
    bed_meta = ident(bed, False)
    with bed.open("rb") as f:
        bed_header = f.read(3)
    bed_meta["header_hex"] = bed_header.hex()
    bed_meta["expected_variant_major_bytes"] = 3 + counters["raw_variants"] * math.ceil(len(samples) / 4)
    bed_meta["size_and_header_match"] = bed_header == b"\x6c\x1b\x01" and bed_meta["bytes"] == bed_meta["expected_variant_major_bytes"]
    bed_meta["historical_expected_sha256_not_freshly_verified"] = prep["source_sha256"]["1000G_EUR.bed"]

    lava_family = []
    for c in range(1, 24):
        bcor = LAVA / f"lava-ukb-v1.1_chr{c}.bcor"
        info = LAVA / f"lava-ukb-v1.1_chr{c}.info"
        with bcor.open("rb") as f:
            header = f.read(48)
        values = struct.unpack("<12I", header)
        assert values[0] == 21775 and values[1] >= 12 and values[3:6] == (3, 4, 14)
        lava_family.append({"CHR": c, "bcor": ident(bcor, False), "info": ident(info, False),
                            "header_first_48_bytes_sha256": hashlib.sha256(header).hexdigest(),
                            "header_little_endian_uint32": values,
                            "format_interpreted_from_pinned_source": "EigenDecomposed / float"})
    lava22 = {}
    with (LAVA / "lava-ukb-v1.1_chr22.info").open() as f:
        for r in csv.DictReader(f, delimiter="\t"):
            assert r["SNP"] not in lava22
            lava22[r["SNP"]] = (int(r["POS"]), r["A1"], r["A2"])
    assert len(lava22) == lava_family[21]["header_little_endian_uint32"][2]
    with gzip.open(next(Path(p) for p in j["reference_sha256"] if Path(p).name == "22.l2.ldscore.gz"), "rt") as f:
        lava22_ref = Counter()
        for r in csv.DictReader(f, delimiter="\t"):
            lava22_ref["consumed_reference_rows"] += 1
            if r["SNP"] not in lava22:
                lava22_ref["absent_from_lava22"] += 1
            elif int(r["BP"]) != lava22[r["SNP"]][0]:
                lava22_ref["coordinate_conflict"] += 1
            else:
                lava22_ref["ID_coordinate_match"] += 1

    after = {str(p): ident(p) for p in watched}
    assert before == after
    result = {"recorded_utc": datetime.now(timezone.utc).isoformat(),
              "scope": "read-only metadata/source comparison, no BED decode/bcor payload read/numerical imports/fits/calibration/outcomes",
              "admitted_for_realistic_calibration": False,
              "consumed_reference_and_coordinate_identity_unchanged": True,
              "before_after_identities": before, "raw_genotype_counts": dict(counters),
              "fam_samples": len(samples), "fam_missing_from_official_panel": len(samples) - sum(ancestry.values()),
              "fam_super_pop_counts": dict(ancestry_groups), "fam_population_counts": dict(ancestry),
              "official_panel_matching_EUR_sample_set": set(samples) == {s for s, r in panel.items() if r["super_pop"] == "EUR"},
              "raw_genotype_per_chr_compatibility": per_chr,
              "raw_genotype_matched_canonical200_block_counts": block_counts,
              "raw_bed_metadata_only": bed_meta,
              "lava_family_metadata_only": lava_family, "lava_version_info": (LAVA / "VERSION_INFO").read_text(),
              "lava22_consumed_reference_compatibility": dict(lava22_ref),
              "lava22_info_rows": len(lava22),
              "lava_fresh_full_payload_hashes_verified": False,
              "peak_rss_bytes_macos": resource.getrusage(resource.RUSAGE_SELF).ru_maxrss}
    OUT.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    print(json.dumps({k: result[k] for k in ["admitted_for_realistic_calibration", "raw_genotype_counts", "fam_population_counts", "official_panel_matching_EUR_sample_set", "lava22_consumed_reference_compatibility", "peak_rss_bytes_macos"]}, indent=2))


if __name__ == "__main__":
    main()
