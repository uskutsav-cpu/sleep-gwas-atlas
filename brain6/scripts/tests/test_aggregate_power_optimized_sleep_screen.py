import csv
import importlib.util
import json
import tempfile
import unittest
from pathlib import Path

SCRIPT = Path(__file__).resolve().parents[1] / "aggregate_power_optimized_sleep_screen_v1.py"
SPEC = importlib.util.spec_from_file_location("aggregate_power_screen", SCRIPT)
AGGREGATOR = importlib.util.module_from_spec(SPEC)
assert SPEC.loader is not None
SPEC.loader.exec_module(AGGREGATOR)


class AggregatePowerScreenTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.root = Path(self.temp.name)
        self.run_dir = self.root / AGGREGATOR.REL_DIR
        self.run_dir.mkdir(parents=True)
        self.base_ids = [f"L{i:04d}" for i in range(2495)]
        self.base = {
            "analysis_id": "brain6_power_optimized_sensitivity_v1_longsleep_trait_only_screen",
            "trait_id": "sleep_duration_continuous_dashti_2019",
            "strict_gate_p": 0.05 / 17465,
            "reference_prefix": "/reference/prefix",
            "input_root": "/candidate/inputs",
            "loci_file": "/frozen/loci.locfile",
            "random_seed": 20260922,
            "scope": "Brain6-only outcome-blinded test fixture",
            "execution_policy": {"worker_count": 4},
            "locus_ids": self.base_ids,
            "output_tsv": str(AGGREGATOR.REL_DIR / "combined.tsv"),
            "summary_json": str(AGGREGATOR.REL_DIR / "combined.summary.json"),
        }
        self._write(self.run_dir / "config.json", self.base)
        prov_dir = self.root / "brain6/results/power_optimized_sensitivity_v1"
        prov_dir.mkdir(parents=True, exist_ok=True)
        source_path = self.root / "brain6/config/source.json"
        source_path.parent.mkdir(parents=True)
        source_path.write_text("candidate-source-fixture\n", encoding="utf-8")
        inherited_output = prov_dir / "audit-fixture.tsv"
        inherited_output.write_text("audit-output-fixture\n", encoding="utf-8")
        self._write(prov_dir / "provenance.json", {
            "decision_rule": {
                "minimum_additional_loci": 250,
                "material_improvement_strict_gate_eligible_loci_percentage_points": 10,
                "downstream_association_results_consulted": False,
            },
            "inputs": {"brain6/config/source.json": AGGREGATOR.sha256(source_path)},
            "outputs": {"audit-fixture.tsv": AGGREGATOR.sha256(inherited_output)},
        })
        with (prov_dir / "current_sleep_trait_audit.tsv").open("w", newline="", encoding="utf-8") as f:
            writer = csv.DictWriter(
                f, fieldnames=["trait", "strict_univariate_gate_pass"], delimiter="\t"
            )
            writer.writeheader()
            writer.writerow({"trait": "longsleep", "strict_univariate_gate_pass": 1})
            writer.writerow({"trait": "insomnia", "strict_univariate_gate_pass": 0})
        self._make_workers()

    def tearDown(self):
        self.temp.cleanup()

    @staticmethod
    def _write(path, value):
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(json.dumps(value), encoding="utf-8")

    def _make_workers(self, gate_passes=0):
        for i, start in enumerate((0, 624, 1248, 1872), start=1):
            end = (624, 1248, 1872, 2495)[i - 1]
            ids = self.base_ids[start:end]
            cfg = {k: v for k, v in self.base.items() if k not in ("locus_ids", "output_tsv", "summary_json")}
            cfg.update({
                "locus_ids": ids,
                "expected_loci": len(ids),
                "output_tsv": str(AGGREGATOR.REL_DIR / f"worker_{i}.tsv"),
                "summary_json": str(AGGREGATOR.REL_DIR / f"worker_{i}.summary.json"),
            })
            config_path = self.run_dir / f"worker_{i}.config.json"
            self._write(config_path, cfg)
            result_path = self.root / cfg["output_tsv"]
            summary_path = self.root / cfg["summary_json"]
            worker_gate_passes = 0
            with result_path.open("w", newline="", encoding="utf-8") as f:
                writer = csv.DictWriter(f, fieldnames=AGGREGATOR.ROW_FIELDS, delimiter="\t")
                writer.writeheader()
                for global_index, locus_id in enumerate(ids, start=start):
                    passed = global_index < gate_passes
                    p_value = self.base["strict_gate_p"] / 2 if passed else 0.05
                    worker_gate_passes += int(passed)
                    writer.writerow({
                        "trait_id": self.base["trait_id"], "locus_id": locus_id,
                        "status": "TESTED", "n_snps": 100, "n_components": 2,
                        "h2_obs": 0.1, "p": p_value,
                        "strict_gate_pass": "TRUE" if passed else "FALSE", "reason": "",
                    })
            self._write(summary_path, {
                "analysis_id": self.base["analysis_id"], "scope": self.base["scope"],
                "rows": len(ids), "tested": len(ids), "not_run": 0, "failed": 0,
                "strict_gate_pass": worker_gate_passes,
                "threshold": self.base["strict_gate_p"],
            })

    def test_aggregates_complete_disjoint_batches_and_refuses_overwrite(self):
        result = AGGREGATOR.aggregate(self.root, write=False)
        self.assertEqual(result["status"], "COMPLETE_SCREEN_NO_FAILED_ROWS")
        self.assertEqual(result["rows"], 2495)
        self.assertEqual(result["tested"], 2495)
        self.assertEqual(result["strict_gate_pass"], 0)
        self.assertFalse(result["eligible_for_full_sensitivity_screen_gate"])
        written = AGGREGATOR.aggregate(self.root, write=True)
        self.assertTrue(written["write_requested"])
        output_path = self.root / self.base["output_tsv"]
        self.assertTrue(output_path.exists())
        self.assertEqual([row["locus_id"] for row in AGGREGATOR.read_tsv(output_path)], self.base_ids)
        with self.assertRaises(FileExistsError):
            AGGREGATOR.aggregate(self.root, write=True)

    def test_rejects_duplicate_or_wrong_worker_assignment(self):
        cfg_path = self.run_dir / "worker_2.config.json"
        cfg = json.loads(cfg_path.read_text(encoding="utf-8"))
        cfg["locus_ids"][0] = cfg["locus_ids"][1]
        self._write(cfg_path, cfg)
        with self.assertRaisesRegex(ValueError, "invalid locus assignment"):
            AGGREGATOR.aggregate(self.root)

    def test_rejects_failed_result_rows(self):
        cfg = json.loads((self.run_dir / "worker_1.config.json").read_text(encoding="utf-8"))
        result_path = self.root / cfg["output_tsv"]
        rows = AGGREGATOR.read_tsv(result_path)
        rows[0]["status"] = "FAILED"
        with result_path.open("w", newline="", encoding="utf-8") as f:
            writer = csv.DictWriter(f, fieldnames=AGGREGATOR.ROW_FIELDS, delimiter="\t")
            writer.writeheader()
            writer.writerows(rows)
        with self.assertRaisesRegex(ValueError, "Unacceptable worker status"):
            AGGREGATOR.aggregate(self.root)

    def test_rejects_unrecognized_not_run_reason(self):
        cfg = json.loads((self.run_dir / "worker_1.config.json").read_text(encoding="utf-8"))
        result_path = self.root / cfg["output_tsv"]
        rows = AGGREGATOR.read_tsv(result_path)
        rows[0].update({
            "status": "NOT_RUN", "n_snps": "NA", "n_components": "NA",
            "h2_obs": "NA", "p": "NA", "strict_gate_pass": "FALSE",
            "reason": "UNKNOWN_PROCESS_NULL_REASON",
        })
        with result_path.open("w", newline="", encoding="utf-8") as f:
            writer = csv.DictWriter(f, fieldnames=AGGREGATOR.ROW_FIELDS, delimiter="\t")
            writer.writeheader()
            writer.writerows(rows)
        summary_path = self.root / cfg["summary_json"]
        summary = json.loads(summary_path.read_text(encoding="utf-8"))
        summary.update({"tested": summary["tested"] - 1, "not_run": 1})
        summary_path.write_text(json.dumps(summary), encoding="utf-8")
        with self.assertRaisesRegex(ValueError, "NOT_RUN row has an invalid reason"):
            AGGREGATOR.aggregate(self.root)

    def test_rejects_changed_upstream_audit_input(self):
        source_path = self.root / "brain6/config/source.json"
        source_path.write_text("changed-source-fixture\n", encoding="utf-8")
        with self.assertRaisesRegex(ValueError, "Inherited audit input hash mismatch"):
            AGGREGATOR.aggregate(self.root)

    def test_accepts_worker_json_rounding_but_rejects_changed_threshold(self):
        summary_path = self.run_dir / "worker_1.summary.json"
        summary = json.loads(summary_path.read_text(encoding="utf-8"))
        summary["threshold"] = float(f"{self.base['strict_gate_p']:.5g}")
        summary_path.write_text(json.dumps(summary), encoding="utf-8")
        self.assertEqual(AGGREGATOR.aggregate(self.root)["rows"], 2495)

        summary["threshold"] *= 1.01
        summary_path.write_text(json.dumps(summary), encoding="utf-8")
        with self.assertRaisesRegex(ValueError, "Worker 1 threshold mismatch"):
            AGGREGATOR.aggregate(self.root)

    def test_requires_250_additional_strict_gate_loci_above_baseline(self):
        self._make_workers(gate_passes=250)
        below = AGGREGATOR.aggregate(self.root)
        self.assertEqual(below["additional_strict_gate_pass"], 249)
        self.assertFalse(below["eligible_for_full_sensitivity_screen_gate"])
        self._make_workers(gate_passes=251)
        at_threshold = AGGREGATOR.aggregate(self.root)
        self.assertEqual(at_threshold["additional_strict_gate_pass"], 250)
        self.assertTrue(at_threshold["eligible_for_full_sensitivity_screen_gate"])


if __name__ == "__main__":
    unittest.main()
