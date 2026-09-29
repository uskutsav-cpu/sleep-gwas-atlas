import json
import sys
from pathlib import Path
from tempfile import TemporaryDirectory

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT / "brain6/scripts"))

import build_project_state_inventory as inventory_builder


def test_inventory_separates_technical_advancement_from_family_qc(monkeypatch, capsys):
    with TemporaryDirectory(dir=ROOT) as temporary_directory:
        output = Path(temporary_directory) / "inventory.json"
        monkeypatch.setattr(sys, "argv", ["build_project_state_inventory.py", "--output", str(output)])
        inventory_builder.main()
        snapshot = json.loads(output.read_text(encoding="utf-8"))

    screen = next(stage for stage in snapshot["stages"]
                  if stage["stage"] == "sleep_power_optimized_screen")
    metrics = screen["metrics"]
    assert metrics["technical_advancement_trigger_met"] is True
    assert metrics["technical_improvement_criteria_met"] == {
        "finite_local_h2_processability_gain": True,
        "low_local_h2_not_run_reduction": True,
    }
    assert metrics["candidate_specific_overlap_matrices"] == 3
    assert metrics["substitution_family_not_run"] == 3305
    assert metrics["substitution_family_cells"] == 17465
    assert metrics["substitution_family_maximum_not_run_cells"] == 873
    assert metrics["substitution_family_passes_frozen_5_percent_qc"] is False
    assert metrics["candidate_bivariate_lava_started"] is False
    assert "family fails the frozen 5% untested-cell ceiling" in screen["summary"]
    candidate_loci = next(stage for stage in snapshot["stages"]
                          if stage["stage"] == "partial_candidate_gene_context")
    candidate_metrics = candidate_loci["metrics"]
    assert candidate_metrics["regions_not_tiered"] == 19
    assert candidate_metrics["lead_signals_stable_across_ld_methods"] == 21
    assert candidate_metrics["range_filtered_method_specific_leads"] == 4
    assert candidate_metrics["method_specific_lead_ids"] == [
        "rs199534", "rs199535", "rs4938023", "rs6430538",
    ]
    assert candidate_metrics["candidate_placo_pvalues_recomputed"] == 2297
    assert candidate_metrics["candidate_placo_pvalues_matching"] == 2297
    assert candidate_metrics["candidate_placo_pvalue_max_relative_error"] < 1e-12
    assert "not genotype-level validated" in candidate_loci["summary"]
    assert json.loads(capsys.readouterr().out)["stages"] == 12
