#!/usr/bin/env python3
"""Check selected manuscript numbers against source tables and retain provenance."""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path


SOURCES = {
    "manuscript": "frailty_paper/paper/manuscript.md",
    "fi_rg": "frailty_paper/analysis/frozen_atlas_frailty_global_rg.tsv",
    "fi_correction": "frailty_paper/analysis/frozen_fi_correction_sensitivity.tsv",
    "latent_rg": "frailty_paper/results/frailty_v1/rg_sleep_latent_factors.tsv",
    "latent_correction": "frailty_paper/analysis/latent_factor_correction_sensitivity.tsv",
    "latent_h2": "frailty_paper/results/frailty_v1/h2_latent_factors.tsv",
    "aging_rg": "frailty_paper/analysis/frozen_atlas_aging_context_global_rg.tsv",
    "record_counts": "frailty_paper/review/source_record_counts.tsv",
    "screening_queue": "frailty_paper/review/screening/title_abstract_queue.tsv",
}
SCRIPT = Path(__file__).resolve()


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(4 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as stream:
        return list(csv.DictReader(stream, delimiter="\t"))


def normalized(path: Path) -> str:
    return " ".join(path.read_text(encoding="utf-8").replace("−", "-").split())


def audit(repo: Path) -> dict:
    paths = {name: repo / rel for name, rel in SOURCES.items()}
    missing = [str(path) for path in paths.values() if not path.is_file()]
    if missing:
        raise FileNotFoundError(f"missing manuscript-claim source(s): {missing}")
    text = normalized(paths["manuscript"])
    fi = read_tsv(paths["fi_rg"])
    fi_correction = read_tsv(paths["fi_correction"])
    latent = read_tsv(paths["latent_rg"])
    latent_correction = read_tsv(paths["latent_correction"])
    latent_h2 = read_tsv(paths["latent_h2"])
    aging = read_tsv(paths["aging_rg"])
    record_counts = read_tsv(paths["record_counts"])
    queue = read_tsv(paths["screening_queue"])

    claims: list[dict] = []

    def check(claim_id: str, manuscript_text: str, expected: object, observed: object, source: str) -> None:
        present = " ".join(manuscript_text.replace("−", "-").split()) in text
        def matches(left: object, right: object) -> bool:
            if isinstance(left, tuple) and isinstance(right, tuple):
                return len(left) == len(right) and all(matches(a, b) for a, b in zip(left, right))
            if isinstance(left, (int, float)) and isinstance(right, (int, float)):
                if isinstance(left, int) and isinstance(right, int):
                    return left == right
                tolerance = max(abs(float(left)) * 0.01, 1e-300) if abs(float(left)) < 1e-8 else 5e-5
                return math.isclose(float(left), float(right), rel_tol=0.0, abs_tol=tolerance)
            return left == right
        claims.append({
            "claim_id": claim_id,
            "expected_from_source": expected,
            "observed_from_source": observed,
            "manuscript_text_verified": present,
            "source": SOURCES[source],
            "status": "PASS" if matches(expected, observed) and present else "FAIL",
        })

    fi_n = len(fi)
    fi_q = sum(float(row["global_rg_fdr_all_396"]) <= 0.05 for row in fi)
    fi_bonf = sum(row["bonferroni_all_396_significant_at_0.05"] == "TRUE" for row in fi_correction)
    check("fi_bh_count", "9 of the 12 FI comparisons had Benjamini–Hochberg q≤0.05", (9, 12), (fi_q, fi_n), "fi_rg")
    check("fi_bonferroni_count", "All nine comparisons meeting the inherited BH q≤0.05 criterion also met the same-family Bonferroni sensitivity criterion", (9, 12), (fi_bonf, len(fi_correction)), "fi_correction")
    for trait, expected, snippet in [
        ("insomnia", (0.6405, 0.0235), "insomnia (r_g=0.6405, SE=0.0235)"),
        ("sleep_apnea", (0.4687, 0.0321), "sleep apnea (r_g=0.4687, SE=0.0321)"),
        ("shortsleep", (0.4543, 0.0272), "short sleep (r_g=0.4543, SE=0.0272)"),
        ("sleepdur", (-0.2220, 0.0275), "Sleep duration (r_g=-0.2220, SE=0.0275)"),
        ("sleep_efficiency", (-0.1520, 0.0365), "sleep efficiency (r_g=-0.1520, SE=0.0365)"),
    ]:
        row = next(row for row in fi if row["sleep_trait"] == trait)
        observed = (float(row["global_rg"]), float(row["global_rg_se"]))
        check(f"fi_{trait}_estimate", snippet, expected, observed, "fi_rg")

    latent_n = len(latent)
    latent_bh = sum(float(row["q_value"]) <= 0.05 for row in latent)
    latent_bonf = sum(row["bonferroni_significant_at_0.05"] == "TRUE" for row in latent_correction)
    check("latent_bh_count", "48 had Benjamini–Hochberg q≤0.05", (48, 84), (latent_bh, latent_n), "latent_rg")
    check("latent_bonferroni_count", "37 also met the same-family Bonferroni sensitivity", (37, 84), (latent_bonf, len(latent_correction)), "latent_correction")

    latent_examples = [
        ("insomnia", "frailty_general", (0.6658, 0.0256, 1.73e-147), "1.73×10^-147"),
        ("shortsleep", "frailty_general", (0.5278, 0.0319, 5.76e-60), "5.76×10^-60"),
        ("sleep_apnea", "frailty_general", (0.4272, 0.0375, 8.53e-29), "8.53×10^-29"),
        ("sleepdur", "frailty_general", (-0.2587, 0.0284, 6.66e-19), "6.66×10^-19"),
        ("sleep_apnea", "frailty_factor_4", (0.4276, 0.0419, 2.19e-23), "2.19×10^-23"),
    ]
    for trait, factor, expected, q_text in latent_examples:
        row = next(row for row in latent if row["sleep_trait"] == trait and row["disease_trait"] == factor)
        observed = (float(row["rg"]), float(row["se"]), float(row["q_value"]))
        snippet = f"r_g={expected[0]:.4f}, SE={expected[1]:.4f}, q={q_text}"
        check(f"latent_{trait}_{factor}", snippet, expected, observed, "latent_rg")

    passing_h2 = [row for row in latent_h2 if row["verdict"] == "PASS" and row["pass_z"] == "True" and row["pass_intercept"] == "True"]
    max_intercept = max(float(row["intercept"]) for row in latent_h2)
    factor3 = next(row for row in latent_h2 if row["trait"] == "frailty_factor_3")
    check("latent_h2_gates", "All seven factor GWAS passed the inherited h² input gate", (7, 7), (len(passing_h2), len(latent_h2)), "latent_h2")
    check("latent_max_intercept", "Factor 3 had the largest intercept (1.195)", ("frailty_factor_3", 1.195), (factor3["trait"], float(factor3["intercept"])), "latent_h2")
    if max_intercept != float(factor3["intercept"]):
        claims[-1]["status"] = "FAIL"

    aging_n = len(aging)
    aging_sig = [row for row in aging if float(row["global_rg_fdr_all_396"]) <= 0.05]
    check("aging_bh_count", "26 of 72 comparisons had q≤0.05", (26, 72), (len(aging_sig), aging_n), "aging_rg")
    aging_counts = Counter(row["non_sleep_trait"] for row in aging_sig)
    for trait, count, label in [
        ("parental_lifespan", 9, "9/12 for parental lifespan"),
        ("healthspan", 8, "8/12 for healthspan"),
        ("grip_strength", 7, "7/12 for hand-grip strength"),
        ("longevity", 1, "1/12 each for longevity and Parkinson disease"),
        ("parkinson", 1, "1/12 each for longevity and Parkinson disease"),
        ("alz", 0, "0/12 for Alzheimer disease"),
    ]:
        check(f"aging_{trait}_count", label, (count, 12), (aging_counts[trait], sum(row["non_sleep_trait"] == trait for row in aging)), "aging_rg")
    for trait, outcome, expected, q_text in [
        ("insomnia", "healthspan", (0.3432, 3.70e-15), "3.70×10^-15"),
        ("insomnia", "parental_lifespan", (-0.2740, 7.98e-13), "7.98×10^-13"),
    ]:
        row = next(row for row in aging if row["sleep_trait"] == trait and row["non_sleep_trait"] == outcome)
        observed = (float(row["global_rg"]), float(row["global_rg_fdr_all_396"]))
        snippet = f"r_g={expected[0]:.4f}, q={q_text}"
        check(f"aging_{trait}_{outcome}", snippet, expected, observed, "aging_rg")

    retrieved = next(int(row["records_imported_or_retrieved"]) for row in record_counts if row["database"] == "PubMed")
    retained = next(int(row["unique_records_retained_by_priority"]) for row in record_counts if row["database"] == "PubMed")
    queue_n = len(queue)
    check("pubmed_retrieved_count", "retrieval contains 61,009 source-record occurrences", (61009,), (retrieved,), "record_counts")
    check("pubmed_deduplicated_count", "after deduplication, 56,117 remain unscreened", (56117,), (retained,), "record_counts")
    check("screening_queue_count", "56,117 remain unscreened", (56117,), (queue_n,), "screening_queue")

    return {
        "schema_version": 1,
        "audited_at_utc": datetime.now(timezone.utc).isoformat(),
        "status": "PASS" if all(claim["status"] == "PASS" for claim in claims) else "FAIL",
        "claims_checked": len(claims),
        "claims_passed": sum(claim["status"] == "PASS" for claim in claims),
        "sources": {name: {"path": rel, "sha256": sha256(paths[name])} for name, rel in SOURCES.items()},
        "script_sha256": sha256(SCRIPT),
        "claims": claims,
        "interpretation": "This audit checks selected quantitative manuscript statements against current source tables. It does not validate scientific eligibility, screening decisions, causal interpretation, or complete manuscript accuracy.",
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--repo", type=Path, default=Path(__file__).resolve().parents[2])
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    repo = args.repo.resolve()
    output = args.output or repo / "frailty_paper/analysis/manuscript_quantitative_claims_audit_2026-09-24.json"
    if not output.is_absolute():
        output = repo / output
    result = audit(repo)
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(result, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(f"{result['status']} manuscript_claims={result['claims_checked']} passed={result['claims_passed']}")
    return 0 if result["status"] == "PASS" else 1


if __name__ == "__main__":
    raise SystemExit(main())
