import csv
import json
from pathlib import Path

from validate_power_optimized_sleep_audit import validate


ROOT = Path(__file__).resolve().parents[2]
AUDIT = ROOT / "results/power_optimized_sensitivity_v1"


def test_power_screen_audit_is_outcome_blinded_and_rejects_unavailable_candidates():
    result = validate()
    assert result["status"] == "PASS_BRAIN6_POWER_AUDIT"
    assert result["sleep_traits"] == 2
    assert result["candidates_screened"] == 6

    with (AUDIT / "candidate_replacement_screen_complete_v2.tsv").open(newline="", encoding="utf-8") as stream:
        complete_rows = {row["candidate"]: row for row in csv.DictReader(stream, delimiter="\t")}
    assert len(complete_rows) == 6
    assert complete_rows["Dashti 2019 continuous self-reported sleep duration"]["decision"].startswith(
        "NO FULL-SENSITIVITY PROMOTION"
    )
    assert complete_rows["Portas et al. 2026 device-measured sleep GWAS"]["decision"].startswith(
        "EXCLUDE from the long-sleep power screen"
    )

    with (AUDIT / "candidate_replacement_screen.tsv").open(newline="", encoding="utf-8") as stream:
        rows = {row["candidate"]: row for row in csv.DictReader(stream, delimiter="\t")}
    assert "Watanabe et al. 2022 2.365M insomnia meta-analysis" in rows
    assert "Austin-Zimmerman et al. 2025 multi-ancestry continuous sleep duration (preprint)" in rows
    assert "full meta unavailable" in rows[
        "Watanabe et al. 2022 2.365M insomnia meta-analysis"
    ]["decision"]
    assert rows[
        "Austin-Zimmerman et al. 2025 multi-ancestry continuous sleep duration (preprint)"
    ]["decision"].startswith("DO NOT USE NOW")

    provenance = json.loads((AUDIT / "provenance.json").read_text(encoding="utf-8"))
    assert provenance["decision_rule"]["downstream_association_results_consulted"] is False
    assert provenance["decision_rule"]["minimum_additional_loci"] == 250
    assert provenance["canonical_decision"] == "FAILED_QC_NOT_PROMOTED"
