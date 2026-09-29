#!/usr/bin/env python3
"""Record published rounded long-sleep LDSC estimates as non-independent context."""
from __future__ import annotations

import csv
import hashlib
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
GLOBAL = ROOT / "brain6/results/global/brain6_72_locked.tsv"
MASTER = ROOT / "brain6/manifests/gwas_master.tsv"
OUTPUT = ROOT / "brain6/results/novelty/sleepchart_2026_long_sleep_rg_context.tsv"
PROVENANCE = OUTPUT.with_suffix(".provenance.json")
FIELDS = [
    "sleep_trait", "brain_disorder", "sleepchart_rg_rounded", "sleepchart_se_rounded",
    "sleepchart_p_main_text", "brain6_rg", "brain6_se", "brain6_original_396_q",
    "direction_concordant", "sleepchart_phenotype", "sleepchart_sleep_n_cases",
    "sleepchart_sleep_n_controls", "sleepchart_sleep_source", "brain6_sleep_definition",
    "brain6_sleep_n_cases", "brain6_sleep_n_controls", "sleep_sample_overlap",
    "outcome_release_identity", "independent_replication", "power_replacement_eligible",
    "interpretation",
]

# Main-text rounded LDSC estimates from the article's long-vs-normal analysis.
# Exact coefficients and standard errors are provided in the linked supplement,
# which is intentionally not transcribed here.
PUBLISHED = {
    "adhd": ("0.28", "0.04", "2.24e-12"),
    "mdd": ("0.29", "0.04", "2.57e-11"),
    "scz": ("0.28", "0.03", "3.47e-16"),
    "bipolar": ("0.21", "0.03", "1.09e-07"),
}
DISORDER_NAMES = {"adhd": "ADHD", "mdd": "MDD", "scz": "SCZ", "bipolar": "bipolar disorder"}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as stream:
        return list(csv.DictReader(stream, delimiter="\t"))


def build() -> tuple[list[dict[str, str]], dict[str, object]]:
    global_rows = read_tsv(GLOBAL)
    master_rows = read_tsv(MASTER)
    global_index = {(row["sleep_trait"], row["brain_disorder"]): row for row in global_rows}
    long_sleep = next((row for row in master_rows if row["trait"] == "longsleep"), None)
    if long_sleep is None or long_sleep.get("cases") != "34184" or long_sleep.get("controls") != "305742":
        raise ValueError("Locked Brain6 long-sleep source counts changed")
    rows = []
    for disorder, (rg, se, p) in PUBLISHED.items():
        source = global_index.get(("longsleep", disorder))
        if source is None:
            raise ValueError(f"Locked Brain6 global map lacks longsleep-{disorder}")
        concordant = (float(source["rg"]) > 0) == (float(rg) > 0)
        rows.append({
            "sleep_trait": "longsleep",
            "brain_disorder": disorder,
            "sleepchart_rg_rounded": rg,
            "sleepchart_se_rounded": se,
            "sleepchart_p_main_text": p,
            "brain6_rg": source["rg"],
            "brain6_se": source["se"],
            "brain6_original_396_q": source["original_BH_FDR_q"],
            "direction_concordant": str(concordant),
            "sleepchart_phenotype": "Long sleep >8 h versus normal sleep 6–8 h; self-reported UK Biobank field 1160",
            "sleepchart_sleep_n_cases": "25049",
            "sleepchart_sleep_n_controls": "300420",
            "sleepchart_sleep_source": "UK Biobank; European-ancestry GWAS",
            "brain6_sleep_definition": "Dashti 2019 long sleep >=9 h versus 7–8 h; UK Biobank",
            "brain6_sleep_n_cases": long_sleep["cases"],
            "brain6_sleep_n_controls": long_sleep["controls"],
            "sleep_sample_overlap": "Same UK Biobank source population and field 1160; participant overlap is expected, exact overlap not reported",
            "outcome_release_identity": "PGC outcome sources; exact release identity versus Brain6 outcome GWAS unresolved",
            "independent_replication": "False",
            "power_replacement_eligible": "False",
            "interpretation": "Rounded article-text prior context only; different sleep threshold and shared UK Biobank source preclude independent replication. The smaller long-sleep case count is not a power upgrade. Do not change Brain6 estimates, novelty, or FDR labels.",
        })
    if len(rows) != 4 or any(row["direction_concordant"] != "True" for row in rows):
        raise ValueError("SleepChart source-context pair set or direction check failed")
    provenance: dict[str, object] = {
        "schema_version": 1,
        "status": "PASS_PRIOR_CONTEXT_ROUNDED_ESTIMATES_NOT_INDEPENDENT_REPLICATION",
        "analysis_id": "brain6_sleepchart_2026_longsleep_rg_context",
        "scope": {
            "rows": 4,
            "source_estimates": "Rounded values transcribed from article main text; supplementary workbook not reviewed",
            "independent_replication": False,
            "power_replacement_eligible": False,
            "novelty_or_fdr_changed": False,
            "pair_specific_loci_compared": False,
        },
        "builder": str(Path(__file__).resolve().relative_to(ROOT)),
        "builder_sha256": sha256(Path(__file__).resolve()),
        "input_hashes": {
            str(GLOBAL.relative_to(ROOT)): sha256(GLOBAL),
            str(MASTER.relative_to(ROOT)): sha256(MASTER),
        },
        "source": {
            "citation": "The MULTI Consortium, O'Toole CK, Song Z, et al. Sleep chart of biological ageing clocks in middle and late life. Nature. 2026. doi:10.1038/s41586-026-10524-5.",
            "doi": "10.1038/s41586-026-10524-5",
            "publisher_url": "https://www.nature.com/articles/s41586-026-10524-5",
            "correction_doi": "10.1038/s41586-026-10920-x",
            "correction_scope": "Adds two prior references; the correction does not report changes to the estimates transcribed here.",
            "accessed_date": "2026-09-26",
            "reviewed_sections": ["main-text long-sleep LDSC results", "UK Biobank sleep phenotype and sample sizes", "methods: UK Biobank and LDSC", "data availability"],
            "sleep_gwas": {
                "source": "UK Biobank field 1160",
                "ancestry": "European",
                "definition": "Long sleep >8 h versus normal sleep 6–8 h",
                "cases": 25049,
                "controls": 300420,
                "variant_count": 6477810,
                "summary_statistics_claimed_public": True,
                "summary_statistics_downloaded": False,
                "reasons_not_used_as_replacement": ["same UK Biobank sleep source population as Brain6", "different binary thresholds", "fewer long-sleep cases than Brain6 source"],
            },
            "estimate_precision": "Rounded rg and SE as stated in main text; the linked Supplementary Data 5a workbook was not downloaded or reviewed.",
            "outcome_overlap": "The paper states PGC outcome statistics exclude UK Biobank; outcome-release identity relative to each Brain6 disorder source is not established here. Shared UK Biobank sleep source prevents calling this independent two-trait replication.",
        },
        "output": {"path": str(OUTPUT.relative_to(ROOT)), "sha256": "PENDING", "rows": len(rows)},
    }
    return rows, provenance


def main() -> None:
    rows, provenance = build()
    OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    with OUTPUT.open("w", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=FIELDS, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    provenance["output"]["sha256"] = sha256(OUTPUT)  # type: ignore[index]
    PROVENANCE.write_text(json.dumps(provenance, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(json.dumps({"output": str(OUTPUT), "sha256": sha256(OUTPUT), "rows": len(rows), "status": provenance["status"]}, indent=2))


if __name__ == "__main__":
    main()
