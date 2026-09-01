#!/usr/bin/env python3
"""Run and verify the frozen Pair B FinnGen h2 gate and insomnia-ADHD LDSC."""

from __future__ import annotations

import argparse
import csv
import hashlib
import io
import json
import math
import re
import subprocess
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
PAIR_MANIFEST = ROOT / "results/track_b/pair_manifest.tsv"
SOURCE_LOCK = ROOT / "results/track_b/replication_source.lock.json"
INGEST_RECEIPT = ROOT / "results/track_b/replication/finngen_r13_F5_ADHD_ingest_receipt.json"
SLEEP_INPUT = ROOT / "data/munged/insomnia.sumstats.gz"
REPLICATION_INPUT = ROOT / "data/munged/track_b_finngen_r13_F5_ADHD.sumstats.gz"
LDSC_PYTHON = ROOT / ".ldsc-env/bin/python"
LDSC = ROOT / "ldsc/ldsc.py"
LD_PREFIX = ROOT / "ref/eur_w_ld_chr"
LOG_DIR = ROOT / "results/track_b/replication/logs"
H2_PREFIX = LOG_DIR / "h2_finngen_r13_F5_ADHD"
RG_PREFIX = LOG_DIR / "rg_insomnia__finngen_r13_F5_ADHD"
OUT = ROOT / "results/track_b/replication/pair_b_ldsc.tsv"
RECEIPT = ROOT / "results/track_b/replication/pair_b_ldsc.provenance.json"

H2_Z_MIN = 4.0
H2_INTERCEPT_MAX = 1.2
REPLICATION_ALPHA = 0.05

H2 = re.compile(r"Total (Observed|Liability) scale h2:\s*(-?[\d.eE+-]+)\s*\(([\d.eE+-]+)\)")
INTERCEPT = re.compile(r"^Intercept:\s*(-?[\d.eE+-]+)\s*\(([\d.eE+-]+)\)", re.M)
INPUT = re.compile(r"Read summary statistics for ([0-9]+) SNPs\.")
REGRESSION = re.compile(r"After merging with regression SNP LD, ([0-9]+) SNPs remain\.")
PAIR_QC = re.compile(
    r"Computing rg for phenotype [0-9]+/[0-9]+\s*\n"
    r"Reading summary statistics from (.+?) \.\.\.\s*\n"
    r"Read summary statistics for ([0-9]+) SNPs\.\s*\n"
    r"After merging with summary statistics, ([0-9]+) SNPs remain\.\s*\n"
    r"([0-9]+) SNPs with valid alleles\.", re.M,
)

FIELDS = [
    "pair_id", "sleep_trait", "external_trait", "replication_source_id",
    "h2_scale", "h2", "h2_se", "h2_z", "h2_intercept", "h2_intercept_se",
    "h2_input_snps", "h2_regression_snps", "h2_gate_status", "rg", "rg_se",
    "rg_z", "rg_p", "cross_trait_intercept", "cross_trait_intercept_se",
    "rg_input_snps", "rg_overlap_after_merge", "rg_valid_alleles",
    "direction_vs_discovery", "replication_class", "claim_limit",
]


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def finite(value: str, field: str) -> float:
    try:
        number = float(value)
    except ValueError as error:
        raise ValueError(f"non-numeric {field}: {value}") from error
    if not math.isfinite(number):
        raise ValueError(f"non-finite {field}: {value}")
    return number


def parse_h2_log(text: str) -> dict[str, object]:
    h2_match, intercept_match = H2.search(text), INTERCEPT.search(text)
    input_match, regression_match = INPUT.search(text), REGRESSION.search(text)
    if not all((h2_match, intercept_match, input_match, regression_match)):
        raise ValueError("required h2 diagnostics are missing")
    assert h2_match and intercept_match and input_match and regression_match
    h2, se = finite(h2_match.group(2), "h2"), finite(h2_match.group(3), "h2_se")
    intercept = finite(intercept_match.group(1), "h2_intercept")
    intercept_se = finite(intercept_match.group(2), "h2_intercept_se")
    if se <= 0:
        raise ValueError("h2_se must be positive")
    z = h2 / se
    passes = z >= H2_Z_MIN and intercept <= H2_INTERCEPT_MAX
    return {
        "h2_scale": h2_match.group(1).upper(), "h2": h2, "h2_se": se, "h2_z": z,
        "h2_intercept": intercept, "h2_intercept_se": intercept_se,
        "h2_input_snps": int(input_match.group(1)),
        "h2_regression_snps": int(regression_match.group(1)),
        "h2_gate_status": "PASS" if passes else "FAIL",
    }


def parse_rg_log(text: str) -> dict[str, object]:
    diagnostics = list(PAIR_QC.finditer(text))
    if len(diagnostics) != 1:
        raise ValueError(f"expected one cross-trait diagnostic block, found {len(diagnostics)}")
    lines = text.splitlines()
    try:
        start = next(i for i, line in enumerate(lines) if "Summary of Genetic Correlation Results" in line)
    except StopIteration as error:
        raise ValueError("rg result table is missing") from error
    header = lines[start + 1].split()
    if not {"p1", "p2", "rg", "se", "z", "p", "gcov_int", "gcov_int_se"}.issubset(header):
        raise ValueError(f"unexpected rg header: {header}")
    values = lines[start + 2].split()
    if len(values) != len(header):
        raise ValueError("rg result row width differs from header")
    parsed = dict(zip(header, values))
    if Path(parsed["p1"]).name != SLEEP_INPUT.name or Path(parsed["p2"]).name != REPLICATION_INPUT.name:
        raise ValueError("rg result does not match the frozen input order")
    rg = finite(parsed["rg"], "rg")
    se = finite(parsed["se"], "rg_se")
    z = finite(parsed["z"], "rg_z")
    if se <= 0:
        raise ValueError("rg_se must be positive")
    # LDSC's displayed p may be rounded; recompute the two-sided normal tail from z.
    p = math.erfc(abs(z) / math.sqrt(2))
    qc = diagnostics[0]
    return {
        "rg": rg, "rg_se": se, "rg_z": z, "rg_p": p,
        "cross_trait_intercept": finite(parsed["gcov_int"], "cross_trait_intercept"),
        "cross_trait_intercept_se": finite(parsed["gcov_int_se"], "cross_trait_intercept_se"),
        "rg_input_snps": int(qc.group(2)), "rg_overlap_after_merge": int(qc.group(3)),
        "rg_valid_alleles": int(qc.group(4)),
    }


def frozen_pair() -> dict[str, str]:
    rows = [r for r in read_tsv(PAIR_MANIFEST) if r["pair_id"] == "B"]
    if len(rows) != 1 or (rows[0]["sleep_trait"], rows[0]["external_trait"]) != ("insomnia", "adhd"):
        raise SystemExit("ERROR: frozen Pair B identity drifted")
    return rows[0]


def preflight() -> tuple[dict[str, str], dict[str, object]]:
    pair = frozen_pair()
    required = [SOURCE_LOCK, INGEST_RECEIPT, SLEEP_INPUT, REPLICATION_INPUT, LDSC_PYTHON, LDSC, LD_PREFIX / "1.l2.ldscore.gz"]
    missing = [str(path.relative_to(ROOT)) for path in required if not path.is_file()]
    if missing:
        raise SystemExit("ERROR: BLOCKED_BY_DATA missing=" + ",".join(missing))
    ingest = json.loads(INGEST_RECEIPT.read_text(encoding="utf-8"))
    if (
        ingest.get("pipeline_status") != "STREAM_HARMONIZE_DIRECT_MUNGE_PASS"
        or ingest.get("munged_output_sha256") != sha256(REPLICATION_INPUT)
        or ingest.get("source_lock_sha256") != sha256(SOURCE_LOCK)
    ):
        raise SystemExit("ERROR: frozen Pair B ingest receipt is invalid")
    return pair, ingest


def run_ldsc(arguments: list[str], prefix: Path) -> None:
    prefix.parent.mkdir(parents=True, exist_ok=True)
    command = [str(LDSC_PYTHON), str(LDSC), *arguments,
               "--ref-ld-chr", str(LD_PREFIX) + "/", "--w-ld-chr", str(LD_PREFIX) + "/",
               "--out", str(prefix)]
    result = subprocess.run(command, cwd=ROOT, check=False, text=True, capture_output=True)
    if result.returncode:
        raise SystemExit(f"ERROR: LDSC failed ({result.returncode}): {(result.stdout + result.stderr)[-2000:]}")


def classify(pair: dict[str, str], rg: dict[str, object]) -> tuple[str, str]:
    discovery = float(pair["discovery_rg"])
    observed = float(rg["rg"])
    concordant = observed == 0 or discovery == 0 or math.copysign(1, observed) == math.copysign(1, discovery)
    direction = "CONCORDANT" if concordant else "OPPOSITE"
    if float(rg["rg_p"]) < REPLICATION_ALPHA and concordant:
        status = "CONCORDANT_NOMINAL_REPLICATION"
    elif float(rg["rg_p"]) < REPLICATION_ALPHA:
        status = "SIGNIFICANT_OPPOSITE_DIRECTION_NO_GO"
    else:
        status = "NOT_SIGNIFICANT_EXTERNAL_SAMPLE"
    return direction, status


def output_text(record: dict[str, object]) -> str:
    buffer = io.StringIO(newline="")
    writer = csv.DictWriter(buffer, fieldnames=FIELDS, delimiter="\t", lineterminator="\n")
    writer.writeheader(); writer.writerow(record)
    return buffer.getvalue()


def atomic_text(path: Path, value: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(value, encoding="utf-8")
    temporary.replace(path)


def publish(record: dict[str, object], logs: list[Path], analysis_status: str) -> None:
    table = output_text(record)
    atomic_text(OUT, table)
    receipt = {
        "schema_version": 1, "analysis_status": analysis_status,
        "selection_timing": "SOURCE_FROZEN_BEFORE_RESULTS", "pair_id": "B",
        "thresholds": {"h2_z_min": H2_Z_MIN, "h2_intercept_max": H2_INTERCEPT_MAX, "replication_alpha": REPLICATION_ALPHA},
        "multiplicity": "ONE_PRESELECTED_PAIR_NO_POST_HOC_REPLACEMENT",
        "input_sha256": {
            str(SLEEP_INPUT.relative_to(ROOT)): sha256(SLEEP_INPUT),
            str(REPLICATION_INPUT.relative_to(ROOT)): sha256(REPLICATION_INPUT),
            str(SOURCE_LOCK.relative_to(ROOT)): sha256(SOURCE_LOCK),
            str(INGEST_RECEIPT.relative_to(ROOT)): sha256(INGEST_RECEIPT),
        },
        "log_sha256": {str(path.relative_to(ROOT)): sha256(path) for path in logs},
        "output_sha256": hashlib.sha256(table.encode()).hexdigest(),
        "script_sha256": sha256(Path(__file__)),
    }
    atomic_text(RECEIPT, json.dumps(receipt, indent=2, sort_keys=True) + "\n")


def verify_result() -> None:
    if not OUT.is_file() or not RECEIPT.is_file():
        raise SystemExit("ERROR: Pair B replication LDSC result is absent")
    receipt = json.loads(RECEIPT.read_text(encoding="utf-8"))
    expected_inputs = {
        str(SLEEP_INPUT.relative_to(ROOT)): sha256(SLEEP_INPUT),
        str(REPLICATION_INPUT.relative_to(ROOT)): sha256(REPLICATION_INPUT),
        str(SOURCE_LOCK.relative_to(ROOT)): sha256(SOURCE_LOCK),
        str(INGEST_RECEIPT.relative_to(ROOT)): sha256(INGEST_RECEIPT),
    }
    if (
        receipt.get("analysis_status") not in {
            "PAIR_B_REPLICATION_LDSC_COMPLETE", "PAIR_B_REPLICATION_H2_QC_FAIL",
        }
        or receipt.get("input_sha256") != expected_inputs
        or receipt.get("output_sha256") != sha256(OUT)
    ):
        raise SystemExit("ERROR: Pair B replication LDSC result or provenance drifted")
    rows = read_tsv(OUT)
    if len(rows) != 1 or rows[0]["pair_id"] != "B":
        raise SystemExit("ERROR: Pair B replication LDSC table is invalid")
    expected_h2 = "PASS" if receipt["analysis_status"] == "PAIR_B_REPLICATION_LDSC_COMPLETE" else "FAIL"
    if rows[0]["h2_gate_status"] != expected_h2:
        raise SystemExit("ERROR: Pair B replication h2 status differs from provenance")
    if expected_h2 == "FAIL" and any(rows[0][field] != "NA" for field in ("rg", "rg_se", "rg_z", "rg_p")):
        raise SystemExit("ERROR: h2-failed Pair B result contains rg statistics")
    print(f"verified Pair B replication LDSC: {receipt['output_sha256']}")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument("--preflight", action="store_true")
    mode.add_argument("--run", action="store_true")
    mode.add_argument("--verify", action="store_true")
    args = parser.parse_args()
    if args.verify:
        preflight(); verify_result(); return
    pair, _ = preflight()
    if args.preflight:
        print("TRACK_B_PAIR_B_LDSC_PREFLIGHT_PASS")
        return
    run_ldsc(["--h2", str(REPLICATION_INPUT)], H2_PREFIX)
    h2_log = H2_PREFIX.with_suffix(".log")
    h2 = parse_h2_log(h2_log.read_text(encoding="utf-8"))
    if h2["h2_gate_status"] != "PASS":
        record = {
            "pair_id": "B", "sleep_trait": "insomnia", "external_trait": "adhd",
            "replication_source_id": "finngen_r13_F5_ADHD", **h2,
            "rg": "NA", "rg_se": "NA", "rg_z": "NA", "rg_p": "NA",
            "cross_trait_intercept": "NA", "cross_trait_intercept_se": "NA",
            "rg_input_snps": "NA", "rg_overlap_after_merge": "NA", "rg_valid_alleles": "NA",
            "direction_vs_discovery": "NOT_APPLICABLE",
            "replication_class": "REPLICATION_H2_QC_FAIL",
            "claim_limit": "Replication source did not pass the predeclared h2 gate; no cross-trait LDSC result was run or inferred",
        }
        publish(record, [h2_log], "PAIR_B_REPLICATION_H2_QC_FAIL")
        print(f"TRACK_B_PAIR_B_REPLICATION_H2_QC_FAIL h2_z={h2['h2_z']} intercept={h2['h2_intercept']}")
        return
    run_ldsc(["--rg", f"{SLEEP_INPUT},{REPLICATION_INPUT}"], RG_PREFIX)
    rg_log = RG_PREFIX.with_suffix(".log")
    rg = parse_rg_log(rg_log.read_text(encoding="utf-8"))
    direction, replication_class = classify(pair, rg)
    record = {
        "pair_id": "B", "sleep_trait": "insomnia", "external_trait": "adhd",
        "replication_source_id": "finngen_r13_F5_ADHD", **h2, **rg,
        "direction_vs_discovery": direction, "replication_class": replication_class,
        "claim_limit": "External-cohort global genetic correlation replication does not establish a local shared locus, mechanism, or causal direction",
    }
    publish(record, [h2_log, rg_log], "PAIR_B_REPLICATION_LDSC_COMPLETE")
    print(f"TRACK_B_PAIR_B_LDSC_COMPLETE class={replication_class} rg={rg['rg']} p={rg['rg_p']}")


if __name__ == "__main__":
    main()
