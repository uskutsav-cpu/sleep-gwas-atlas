#!/usr/bin/env python3
"""Materialize the aggregate-only Brain6 record from the external LDSC receipt."""
from __future__ import annotations

import csv
import hashlib
import json
import re
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
RUN = Path("/Volumes/Extreme SSD/brain6-work/replication-pgc-mdd2025/insomnia_rg_sensitivity_v1")
TABLE = ROOT / "brain6/results/replication/pgc_mdd2025_insomnia_rg_sensitivity_v1.tsv"
PROVENANCE = TABLE.with_suffix(".provenance.json")
FIELDS = [
    "pair_id", "comparison_type", "sleep_trait", "disorder", "discovery_rg", "discovery_se",
    "external_source", "method", "genome_build", "ancestry", "mdd_cases", "mdd_controls",
    "mdd_effective_n", "mdd_cohorts", "mdd_source_variants", "insomnia_source_variants",
    "hm3_valid_allele_snps", "rg", "se", "p", "cross_trait_intercept", "cross_trait_intercept_se",
    "overlap_status", "status",
]


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def main() -> None:
    run_receipt = RUN / "provenance.json"
    run = json.loads(run_receipt.read_text(encoding="utf-8"))
    if run.get("status") != "PASS_LDSC_EXTERNAL_SENSITIVITY":
        raise ValueError("External LDSC receipt is not a successful completed sensitivity")
    result = run["result"]
    if (result.get("status") != "PASS_LDSC_EXTERNAL_SENSITIVITY" or
            result.get("rg") != 0.4771 or result.get("rg_se") != 0.0228 or
            result.get("p") != 4.6535e-97 or
            result.get("cross_trait_intercept") != 0.0087 or
            result.get("cross_trait_intercept_se") != 0.0069):
        raise ValueError("External receipt's result fields are missing or unexpected")
    source = run["sources"]["mdd"]
    insomnia = run["sources"]["insomnia"]
    mdd_rows = result["source_rows"]["mdd2025"]
    insomnia_rows = result["source_rows"]["insomnia"]
    rg_log = Path(result["log_path"])
    if sha256(rg_log) != result.get("rg_log_sha256"):
        raise ValueError("The external LDSC log does not match its run receipt")
    log_text = rg_log.read_text(encoding="utf-8", errors="replace")
    matched = re.search(r"(?m)^(\d+) SNPs with valid alleles\.?$", log_text)
    if matched is None:
        raise ValueError("Could not find the valid-allele SNP count in the LDSC rg log")
    row = {
        "pair_id": "insomnia__mdd_pgc2025_noUKBB",
        "comparison_type": "GLOBAL_RG_SOURCE_SENSITIVITY",
        "sleep_trait": "insomnia",
        "disorder": "mdd",
        "discovery_rg": "0.4425",
        "discovery_se": "0.0252",
        "external_source": source["study_id"],
        "method": "LDSC genome-wide rg; HapMap3 allele-merged, 1000G EUR LD reference",
        "genome_build": source["genome_build"],
        "ancestry": source["ancestry"],
        "mdd_cases": source["cases"],
        "mdd_controls": source["controls"],
        "mdd_effective_n": source["effective_n"],
        "mdd_cohorts": source["cohorts"],
        "mdd_source_variants": mdd_rows["input_rows_written"],
        "insomnia_source_variants": insomnia_rows["input_rows_written"],
        "hm3_valid_allele_snps": int(matched.group(1)),
        "rg": f"{result['rg']:.4f}",
        "se": f"{result['rg_se']:.4f}",
        "p": f"{result['p']:.6g}",
        "cross_trait_intercept": f"{result['cross_trait_intercept']:.4f}",
        "cross_trait_intercept_se": f"{result['cross_trait_intercept_se']:.4f}",
        "overlap_status": "PGC source excludes UK Biobank; exact participant overlap via other recruitment sources is unverified; LDSC cross-trait intercept 0.0087 (SE 0.0069)",
        "status": "SENSITIVITY_ONLY_NOT_INDEPENDENT_REPLICATION",
    }
    TABLE.parent.mkdir(parents=True, exist_ok=True)
    with TABLE.open("w", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=FIELDS, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerow(row)
    code_path = Path(__file__).resolve()
    data = {
        "schema_version": 1,
        "status": "PASS_AGGREGATE_ONLY_EXTERNAL_SENSITIVITY_RECORD",
        "scope": "One aggregate genome-wide LDSC sensitivity; not independent replication, new FDR family, local LAVA result, or a locus-level output.",
        "source_data_policy": "PGC raw and transformed summary statistics remain only on the external research drive and are not included in the repository.",
        "external_run_provenance": str(run_receipt),
        "external_run_provenance_sha256": sha256(run_receipt),
        "external_rg_log_sha256": result["rg_log_sha256"],
        "source_checksums": {
            "pgc_summary_stats_md5": source["source_md5"],
            "pgc_summary_stats_sha256": source["source_sha256"],
            "pgc_readme_md5": source["readme_md5"],
            "insomnia_source_sha256": insomnia["source_sha256"],
            "ldsc_hm3_snplist_sha256": run["reference"]["hm3_snplist_sha256"],
            "ldsc_reference_archive_md5": run["reference"]["archive_md5"],
        },
        "builder": str(code_path.relative_to(ROOT)),
        "builder_sha256": sha256(code_path),
        "output": str(TABLE.relative_to(ROOT)),
        "output_sha256": sha256(TABLE),
        "summary": {key: row[key] for key in ("rg", "se", "p", "cross_trait_intercept", "cross_trait_intercept_se", "status")},
    }
    PROVENANCE.write_text(json.dumps(data, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps({"status": data["status"], "output": str(TABLE), "output_sha256": data["output_sha256"]}, indent=2))


if __name__ == "__main__":
    main()
