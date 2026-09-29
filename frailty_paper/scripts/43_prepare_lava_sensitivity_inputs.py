#!/usr/bin/env python3
"""Prepare and verify immutable chromosome-scoped inputs for locked FI×sleep LAVA sensitivity."""
from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
import math
import os
import re
import shutil
import sys
import uuid
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
LOCK = ROOT / "frailty_paper/config/lava_frailty_sensitivity_v1.yaml"
PANEL = ROOT / "config/analysis_panel.tsv"
RG = ROOT / "frailty_paper/results/frailty_v1/global_rg_sleep_frailty.tsv"
H2_SLEEP = ROOT / "frailty_paper/results/frailty_v1/h2_sleep_panel.tsv"
H2_FI = ROOT / "frailty_paper/results/frailty_v1/h2_frailty_verification.tsv"
EXPECTED = {
    "frailty": (6291635, 106450123, "241ed294ccc367b22787a0c051f8d377c020e16089e82920a555286f3148e6de"),
    "accel_sleep_duration": (6482514, 144197232, "2df99672307d103ac00d6f4ee280ee4ed788c8c777abd39f041177b48ae43bc4"),
    "chronotype": (6482405, 145654493, "c561f7951762cc8e98b7969f220150fafbd11d00b0dd470f18e0a245adab2286"),
    "insomnia": (6077635, 125531591, "81c71d838b19f0297e09ec027d00382732fee1bfd112f44cfd2feca44a4f8bf5"),
    "longsleep": (6549769, 147859672, "150b4136404adb77a8484c6326a6f233290762a71ca4c851af9462de820d0f4c"),
    "napping": (7042618, 155409529, "d9346be15b3f84b2afef1a9a1ce754dc776b90cc3aa7e78dad18d63607714c43"),
    "shortsleep": (6549791, 146959493, "47948b1e63eb15107545cde531e06f82227afce891f24b813039b930e574497d"),
    "sleep_apnea": (6841149, 167251776, "a9c37682a1da7687277781518a2d4b94b19a7985270941e7209d594f5f353aa5"),
    "sleep_efficiency": (6482540, 144204520, "18c766173bf0acca325b927ea4d7118ad51dd9e3310924132707dd12c0d34fb7"),
    "sleep_timing": (6482540, 144192547, "69f3f10605aa07017d108ad8c1b5f9c02a95743ab507823e7bf4c9f50a77d2d4"),
    "sleepdur": (6549809, 145984615, "f2c0b9eaa4623e9757de86d6542fb3267874a2def626ca548fccb3ab97fced41"),
    "sleepiness": (6549823, 144902667, "1a7002ce3679ab15cd3a6af2ba7bd5d54b937e19693ee6b5b1ec1e7c75cd3afe"),
    "snoring": (7168629, 163291287, "b2a2e6871941e56c5eee3799b1af63c6933329b68a6d8def16d4b503c08c7683"),
}
SOURCE_ROOT = Path("/Volumes/Extreme SSD/sleep-gwas-atlas-frailty-v1/analysis-workspace/frailty_v1/harmonized")
OUTPUT_ROOT = Path("/Volumes/Extreme SSD/sleep-gwas-atlas-frailty-v1/analysis-workspace/frailty_v1/lava_sensitivity_v1")
SOURCE_FIELDS = ["SNP", "CHR", "BP", "A1", "A2", "FRQ", "BETA", "SE", "P", "N"]
OUTPUT_FIELDS = ["SNP", "A1", "A2", "Z", "N"]
CHRS = range(1, 23)


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for block in iter(lambda: f.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def atomic_json(path: Path, value: object) -> None:
    tmp = path.with_name(path.name + f".{uuid.uuid4().hex}.partial")
    tmp.write_text(json.dumps(value, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    os.replace(tmp, path)


def read_lock_scalars() -> dict:
    """Read the two top-level/nested scalar lock fields used here without PyYAML."""
    text = LOCK.read_text(encoding="utf-8")
    analysis = re.search(r"^analysis_id:\s*([^\s#]+)\s*$", text, re.MULTILINE)
    outputs = re.search(r"^outputs:\s*\n((?:^[ \t]+.*\n?)*)", text, re.MULTILINE)
    root = re.search(r"^[ \t]+external_run_root:\s*(.+?)\s*$", outputs.group(1), re.MULTILINE) if outputs else None
    if not analysis or not root:
        raise RuntimeError("could not read analysis_id and outputs.external_run_root from lock")
    return {"analysis_id": analysis.group(1), "outputs": {"external_run_root": root.group(1)}}


def load() -> tuple[dict, list[dict[str, str]], dict[str, dict[str, str]], dict[str, dict[str, str]]]:
    cfg = read_lock_scalars()
    panel = list(csv.DictReader(PANEL.open(encoding="utf-8", newline=""), delimiter="\t"))
    sleep = [r for r in panel if r["domain"] == "sleep"]
    if len(panel) != 45 or len(sleep) != 12 or len({r["trait_id"] for r in sleep}) != 12:
        raise RuntimeError("frozen 45-trait panel / 12-trait sleep subset drifted")
    rg = {r["sleep_trait"]: r for r in csv.DictReader(RG.open(encoding="utf-8", newline=""), delimiter="\t")}
    h2s = {r["trait"]: r for r in csv.DictReader(H2_SLEEP.open(encoding="utf-8", newline=""), delimiter="\t")}
    h2fi_rows = list(csv.DictReader(H2_FI.open(encoding="utf-8", newline=""), delimiter="\t"))
    if len(rg) != 12 or set(rg) != {r["trait_id"] for r in sleep} or set(h2s) != set(rg):
        raise RuntimeError("locked global rg and h2 evidence do not cover all sleep traits")
    if len(h2fi_rows) != 1 or h2fi_rows[0]["trait"] != "frailty" or h2fi_rows[0]["verdict"] != "PASS":
        raise RuntimeError("Frailty Index h2 gate is not a unique PASS")
    if any(r["verdict"] != "PASS" for r in h2s.values()):
        raise RuntimeError("a locked sleep-trait h2 gate is not PASS")
    return cfg, sleep, rg, {"frailty": h2fi_rows[0], **h2s}


def expected_rho(trait: str, rg: dict[str, dict[str, str]], h2: dict[str, dict[str, str]]) -> float:
    cross = float(rg[trait]["gcov_int"])
    di = float(h2["frailty"]["intercept"])
    dj = float(h2[trait]["intercept"])
    if not all(map(math.isfinite, (cross, di, dj))) or di <= 0 or dj <= 0:
        raise RuntimeError(f"invalid LDSC intercepts for frailty × {trait}")
    rho = round(cross / math.sqrt(di * dj), 5)
    if not math.isfinite(rho) or abs(rho) >= 1:
        raise RuntimeError(f"invalid estimated overlap correlation for {trait}: {rho}")
    return rho


def write_gzip(path: Path, lines) -> None:
    with path.open("wb") as raw:
        with gzip.GzipFile(filename="", mode="wb", fileobj=raw, mtime=0, compresslevel=6) as gz:
            for line in lines:
                gz.write(line.encode("utf-8"))


def make_trait_shards(trait: str, rows_expected: int, byte_expected: int, hash_expected: str) -> dict:
    source = SOURCE_ROOT / f"{trait}.harmonized.tsv.gz"
    if not source.is_file() or source.stat().st_size != byte_expected or sha256(source) != hash_expected:
        raise RuntimeError(f"locked input identity mismatch: {source}")
    trait_root = OUTPUT_ROOT / "inputs" / trait
    manifest_path = trait_root / "trait_manifest.json"
    if trait_root.exists():
        if not manifest_path.is_file():
            raise RuntimeError(f"existing trait directory has no immutable manifest: {trait_root}")
        manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
        if manifest.get("source_sha256") != hash_expected or manifest.get("source_rows") != rows_expected:
            raise RuntimeError(f"existing shard manifest identity mismatch: {trait}")
        for item in manifest["shards"]:
            path = trait_root / item["file"]
            if not path.is_file() or path.stat().st_size != item["bytes"] or sha256(path) != item["sha256"]:
                raise RuntimeError(f"existing shard failed identity check: {path}")
        return manifest

    parent = trait_root.parent
    parent.mkdir(parents=True, exist_ok=True)
    staging = parent / f".{trait}.building-{uuid.uuid4().hex}"
    staging.mkdir()
    paths = {chrom: staging / f"chr{chrom:02d}.sumstats.tsv.gz" for chrom in CHRS}
    counts = {chrom: 0 for chrom in CHRS}
    writers = {}
    try:
        # Open all chromosome outputs in a trait-local staging directory. The directory is
        # renamed into place only after every source row and compressed member is verified.
        raw_files = {chrom: path.open("wb") for chrom, path in paths.items()}
        gzip_files = {chrom: gzip.GzipFile(filename="", mode="wb", fileobj=raw_files[chrom], mtime=0, compresslevel=6) for chrom in CHRS}
        writers = {chrom: gzip_files[chrom] for chrom in CHRS}
        for chrom in CHRS:
            writers[chrom].write(("\t".join(OUTPUT_FIELDS) + "\n").encode())
        total = 0
        with gzip.open(source, "rt", encoding="utf-8", newline="") as handle:
            reader = csv.reader(handle, delimiter="\t")
            header = next(reader, None)
            if header != SOURCE_FIELDS:
                raise RuntimeError(f"source schema drift for {trait}: {header}")
            for line_number, row in enumerate(reader, 2):
                if len(row) != len(SOURCE_FIELDS):
                    raise RuntimeError(f"malformed row {line_number} in {source}")
                rec = dict(zip(SOURCE_FIELDS, row))
                try:
                    chrom = int(rec["CHR"])
                    beta, se, p, n = map(float, (rec["BETA"], rec["SE"], rec["P"], rec["N"]))
                except ValueError as error:
                    raise RuntimeError(f"nonnumeric association row {trait}:{line_number}") from error
                if chrom not in counts or not all(map(math.isfinite, (beta, se, p, n))) or se <= 0 or not 0 <= p <= 1 or n <= 0:
                    raise RuntimeError(f"invalid association row {trait}:{line_number}")
                z = beta / se
                if not math.isfinite(z):
                    raise RuntimeError(f"non-finite Z at {trait}:{line_number}")
                values = (rec["SNP"], rec["A1"], rec["A2"], format(z, ".17g"), rec["N"])
                writers[chrom].write(("\t".join(values) + "\n").encode("utf-8"))
                counts[chrom] += 1
                total += 1
        for chrom in CHRS:
            writers[chrom].close()
            raw_files[chrom].close()
        writers.clear()
        if total != rows_expected:
            raise RuntimeError(f"row count mismatch for {trait}: observed={total}, expected={rows_expected}")
        shards = []
        for chrom, path in paths.items():
            # Empty chromosomes are permitted only if the frozen source contains no rows there.
            if counts[chrom] == 0:
                continue
            shards.append({"chromosome": chrom, "file": path.name, "rows": counts[chrom],
                           "bytes": path.stat().st_size, "sha256": sha256(path)})
        manifest = {"trait": trait, "source_path": str(source), "source_rows": total,
                    "source_bytes": byte_expected, "source_sha256": hash_expected,
                    "output_schema": OUTPUT_FIELDS, "shards": shards}
        atomic_json(staging / "trait_manifest.json", manifest)
        os.rename(staging, trait_root)
        return manifest
    except BaseException:
        for item in writers.values():
            try: item.close()
            except Exception: pass
        raise


def prepare_pair_files(cfg: dict, sleep: list[dict[str, str]], rg: dict[str, dict[str, str]], h2: dict[str, dict[str, str]], *, materialize: bool) -> dict:
    pair_root = OUTPUT_ROOT / "pairs"
    if materialize:
        pair_root.mkdir(parents=True, exist_ok=True)
    elif not pair_root.is_dir():
        raise RuntimeError("pair input directory is absent")
    records = {}
    for row in sleep:
        trait = row["trait_id"]
        target = pair_root / trait
        if materialize:
            target.mkdir(exist_ok=True)
        elif not target.is_dir():
            raise RuntimeError(f"pair input directory is absent: {target}")
        rho = expected_rho(trait, rg, h2)
        overlap_path = target / "sample_overlap.txt"
        overlap = f" frailty {trait}\nfrailty 1.00000 {rho:.5f}\n{trait} {rho:.5f} 1.00000\n"
        if overlap_path.exists():
            if overlap_path.read_text(encoding="utf-8") != overlap:
                raise RuntimeError(f"pair overlap file changed: {overlap_path}")
        else:
            if not materialize:
                raise RuntimeError(f"pair overlap file is absent: {overlap_path}")
            overlap_path.write_text(overlap, encoding="utf-8")
        for chrom in CHRS:
            info_path = target / f"input_info_chr{chrom:02d}.tsv"
            records_out = []
            for phenotype in ("frailty", trait):
                if phenotype == "frailty":
                    cases = controls = prevalence = "NA"
                else:
                    source_meta = next(r for r in sleep if r["trait_id"] == phenotype)
                    if source_meta["type"] == "binary":
                        cases, controls, prevalence = source_meta["ncase"], source_meta["ncontrol"], source_meta["pop_prev"]
                        if not cases.isdigit() or not controls.isdigit() or not 0 < float(prevalence) < 1:
                            raise RuntimeError(f"binary input metadata failed for {phenotype}")
                    else:
                        cases = controls = prevalence = "NA"
                shard = OUTPUT_ROOT / "inputs" / phenotype / f"chr{chrom:02d}.sumstats.tsv.gz"
                if not shard.is_file():
                    # Preserve a zero-row shard as a header-only gzip for input.info completeness.
                    raise RuntimeError(f"missing chromosome shard {shard}")
                records_out.append({"phenotype": phenotype, "cases": cases, "controls": controls,
                                    "prevalence": prevalence, "filename": str(shard)})
            buf = __import__("io").StringIO()
            writer = csv.DictWriter(buf, fieldnames=["phenotype", "cases", "controls", "prevalence", "filename"],
                                    delimiter="\t", lineterminator="\n")
            writer.writeheader(); writer.writerows(records_out)
            content = buf.getvalue()
            if info_path.exists():
                if info_path.read_text(encoding="utf-8") != content:
                    raise RuntimeError(f"pair input.info changed: {info_path}")
            else:
                if not materialize:
                    raise RuntimeError(f"pair input.info is absent: {info_path}")
                info_path.write_text(content, encoding="utf-8")
        records[trait] = {"rho": rho, "overlap_file": str(overlap_path),
                          "input_info_files": [str(target / f"input_info_chr{c:02d}.tsv") for c in CHRS]}
    return records


def verify_prepared(cfg: dict, sleep: list[dict[str, str]], rg: dict[str, dict[str, str]], h2: dict[str, dict[str, str]]) -> dict:
    manifest_path = OUTPUT_ROOT / "input_manifest.json"
    if not manifest_path.is_file():
        raise RuntimeError("prepared-input manifest is absent")
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    if manifest.get("analysis_id") != cfg["analysis_id"] or manifest.get("lock_sha256") != sha256(LOCK):
        raise RuntimeError("prepared-input manifest does not match current lock")
    if set(manifest.get("traits", {})) != set(EXPECTED):
        raise RuntimeError("prepared-input manifest trait set drifted")
    for trait, (rows, size, digest) in EXPECTED.items():
        source = SOURCE_ROOT / f"{trait}.harmonized.tsv.gz"
        if source.stat().st_size != size or sha256(source) != digest:
            raise RuntimeError(f"source input identity mismatch: {trait}")
        record = manifest["traits"][trait]
        if record.get("source_rows") != rows or record.get("source_sha256") != digest:
            raise RuntimeError(f"manifest source identity mismatch: {trait}")
        trait_root = OUTPUT_ROOT / "inputs" / trait
        local = json.loads((trait_root / "trait_manifest.json").read_text(encoding="utf-8"))
        if local != record:
            raise RuntimeError(f"trait manifest differs from run manifest: {trait}")
        total = 0
        for item in local["shards"]:
            path = trait_root / item["file"]
            if path.stat().st_size != item["bytes"] or sha256(path) != item["sha256"]:
                raise RuntimeError(f"shard identity mismatch: {path}")
            with gzip.open(path, "rt", encoding="utf-8") as handle:
                header = next(handle, "").rstrip("\n").split("\t")
                if header != OUTPUT_FIELDS:
                    raise RuntimeError(f"shard schema drift: {path}")
                count = sum(1 for _ in handle)
            if count != item["rows"]:
                raise RuntimeError(f"shard row count mismatch: {path}")
            total += count
        if total != rows:
            raise RuntimeError(f"total shard row count mismatch: {trait}")
    expected_pairs = {r["trait_id"] for r in sleep}
    if set(manifest.get("pairs", {})) != expected_pairs:
        raise RuntimeError("pair manifest set drifted")
    current_pairs = prepare_pair_files(cfg, sleep, rg, h2, materialize=False)
    if current_pairs != manifest["pairs"]:
        raise RuntimeError("pair-level overlap/input metadata differs from lock sources")
    for trait, pair in manifest["pairs"].items():
        overlap = Path(pair["overlap_file"])
        if not overlap.is_file():
            raise RuntimeError(f"overlap file missing: {overlap}")
        for path_text in pair["input_info_files"]:
            path = Path(path_text)
            if not path.is_file():
                raise RuntimeError(f"pair input.info missing: {path}")
    return manifest


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--verify", action="store_true", help="verify an existing immutable preparation")
    parser.add_argument("--minimum-free-gib", type=float, default=20.0)
    args = parser.parse_args()
    cfg, sleep, rg, h2 = load()
    if cfg["outputs"]["external_run_root"] != str(OUTPUT_ROOT):
        raise RuntimeError("output root differs from frozen analysis lock")
    if not args.verify:
        free_gib = shutil.disk_usage(OUTPUT_ROOT.parent).free / (1024**3) if OUTPUT_ROOT.parent.exists() else shutil.disk_usage(OUTPUT_ROOT.parent.parent).free / (1024**3)
        if free_gib < args.minimum_free_gib:
            raise RuntimeError(f"external output volume has only {free_gib:.1f} GiB free")
        OUTPUT_ROOT.mkdir(parents=True, exist_ok=True)
    if not SOURCE_ROOT.is_dir():
        raise RuntimeError(f"external harmonized input root missing: {SOURCE_ROOT}")
    if args.verify:
        manifest = verify_prepared(cfg, sleep, rg, h2)
        path = OUTPUT_ROOT / "input_manifest.json"
        print(f"PREPARED_INPUTS_VERIFIED traits={len(manifest['traits'])} pairs={manifest['pair_count']} loci={manifest['locus_count']} slots={manifest['candidate_pair_locus_slots']} manifest_sha256={sha256(path)}")
        return 0
    shard_manifests = {}
    for trait, (rows, size, digest) in EXPECTED.items():
        shard_manifests[trait] = make_trait_shards(trait, rows, size, digest)
    pairs = prepare_pair_files(cfg, sleep, rg, h2, materialize=True)
    manifest = {"analysis_id": cfg["analysis_id"], "lock_sha256": sha256(LOCK),
                "source_root": str(SOURCE_ROOT), "input_root": str(OUTPUT_ROOT / "inputs"),
                "pair_count": 12, "locus_count": 2495, "candidate_pair_locus_slots": 29940,
                "traits": shard_manifests, "pairs": pairs,
                "h2_source_sha256": {str(p): sha256(ROOT / p) for p in
                    ["frailty_paper/results/frailty_v1/h2_sleep_panel.tsv",
                     "frailty_paper/results/frailty_v1/h2_frailty_verification.tsv",
                     "frailty_paper/results/frailty_v1/global_rg_sleep_frailty.tsv"]}}
    manifest_path = OUTPUT_ROOT / "input_manifest.json"
    encoded = json.dumps(manifest, indent=2, sort_keys=True) + "\n"
    if manifest_path.exists():
        if manifest_path.read_text(encoding="utf-8") != encoded:
            raise RuntimeError("existing prepared-input manifest differs; refusing overwrite")
    else:
        atomic_json(manifest_path, manifest)
    print(f"PREPARED_INPUTS_PASS traits=13 shards={sum(len(v['shards']) for v in shard_manifests.values())} pairs=12 loci=2495 slots=29940 manifest_sha256={sha256(manifest_path)}")
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except Exception as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        raise SystemExit(1)
