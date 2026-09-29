from __future__ import annotations

import hashlib
import importlib.util
import json
import tempfile
import unittest
from pathlib import Path


SCRIPT = Path(__file__).resolve().parents[1] / "scripts/13_validate_analysis_plan.py"
SPEC = importlib.util.spec_from_file_location("validate_analysis_plan", SCRIPT)
VALIDATOR = importlib.util.module_from_spec(SPEC)
assert SPEC and SPEC.loader
SPEC.loader.exec_module(VALIDATOR)


def digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


class AnalysisPlanLockTests(unittest.TestCase):
    def make_fixture(self, root: Path) -> Path:
        paper_config = root / "frailty_paper/config"
        config = root / "config"
        paper_config.mkdir(parents=True)
        config.mkdir()
        panel = "trait_id\tdomain\n" + "".join(
            f"t{i}\t{'sleep' if i < 12 else 'other'}\n" for i in range(45)
        )
        (config / "analysis_panel.tsv").write_text(panel, encoding="utf-8")
        basis = config / "basis.tsv"
        basis.write_text("locked\n", encoding="utf-8")
        plan = {
            "plan_id": "frailty_sleep_genetics_v1",
            "version": 1,
            "status": "FROZEN_PRE_NEW_FRAILTY_ANALYSES",
            "estimand_and_hierarchy": {"primary_endpoint": {"trait_id": "frailty", "source_id": "GCST90020053"}},
            "locked_inputs": {"sleep_trait_count": 12},
            "multiplicity": {"original_atlas_scope": 396},
        }
        plan_path = paper_config / "analysis_plan_v1.yaml"
        import yaml

        plan_path.write_text(yaml.safe_dump(plan, sort_keys=False), encoding="utf-8")
        lock = {
            "schema_version": "frailty-analysis-plan-lock.1",
            "plan_id": plan["plan_id"],
            "plan_path": "frailty_paper/config/analysis_plan_v1.yaml",
            "plan_sha256": digest(plan_path),
            "basis_file_sha256": {"config/basis.tsv": digest(basis)},
        }
        (paper_config / "analysis_plan_v1.lock.json").write_text(json.dumps(lock), encoding="utf-8")
        return root

    def test_valid_frozen_plan_and_basis_pass(self) -> None:
        with tempfile.TemporaryDirectory() as temp:
            self.assertEqual(VALIDATOR.validate(self.make_fixture(Path(temp))), [])

    def test_plan_tampering_is_reported(self) -> None:
        with tempfile.TemporaryDirectory() as temp:
            root = self.make_fixture(Path(temp))
            plan = root / "frailty_paper/config/analysis_plan_v1.yaml"
            plan.write_text(plan.read_text(encoding="utf-8") + "# post-freeze edit\n", encoding="utf-8")
            self.assertIn("frozen plan SHA-256 mismatch", VALIDATOR.validate(root))

    def test_basis_file_tampering_is_reported(self) -> None:
        with tempfile.TemporaryDirectory() as temp:
            root = self.make_fixture(Path(temp))
            (root / "config/basis.tsv").write_text("changed\n", encoding="utf-8")
            self.assertIn("locked basis file SHA-256 mismatch: config/basis.tsv", VALIDATOR.validate(root))


if __name__ == "__main__":
    unittest.main()
