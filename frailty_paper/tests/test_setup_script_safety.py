from __future__ import annotations

import os
import subprocess
import tempfile
import unittest
from pathlib import Path


REPO = Path(__file__).resolve().parents[2]
SETUP = REPO / "frailty_paper/scripts/00_setup.sh"


class SetupScriptSafetyTests(unittest.TestCase):
    def test_refuses_to_overlay_existing_environment_with_wrong_python(self):
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp) / "repo"
            (root / "config").mkdir(parents=True)
            (root / "config/traits.tsv").write_text("trait\n", encoding="utf-8")
            (root / "config/public_gwas_sources.tsv").write_text("id\n", encoding="utf-8")

            requested_python = Path(temp) / "python311"
            requested_python.write_text("#!/bin/sh\nprintf 'Python 3.11.11\\n'\n", encoding="utf-8")
            requested_python.chmod(0o755)

            existing = Path(temp) / "existing-venv"
            (existing / "bin").mkdir(parents=True)
            existing_python = existing / "bin/python"
            existing_python.write_text("#!/bin/sh\nprintf 'Python 3.13.5\\n'\n", encoding="utf-8")
            existing_python.chmod(0o755)
            sentinel = existing / "keep.txt"
            sentinel.write_text("preserve me\n", encoding="utf-8")

            env = os.environ.copy()
            env.update({"FRAILTY_PYTHON": str(requested_python), "FRAILTY_VENV_DIR": str(existing)})
            result = subprocess.run(
                ["bash", str(SETUP), root.name], capture_output=True, text=True, env=env,
                cwd=temp, check=False
            )

            self.assertEqual(result.returncode, 1)
            self.assertIn("refusing to layer Python 3.11.11", result.stderr)
            self.assertEqual(sentinel.read_text(encoding="utf-8"), "preserve me\n")
            self.assertIn("3.13.5", existing_python.read_text(encoding="utf-8"))

    def test_refuses_existing_non_venv_target(self):
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp) / "repo"
            (root / "config").mkdir(parents=True)
            (root / "config/traits.tsv").write_text("trait\n", encoding="utf-8")
            (root / "config/public_gwas_sources.tsv").write_text("id\n", encoding="utf-8")
            requested_python = Path(temp) / "python311"
            requested_python.write_text("#!/bin/sh\nprintf 'Python 3.11.11\\n'\n", encoding="utf-8")
            requested_python.chmod(0o755)
            target = Path(temp) / "occupied"
            target.mkdir()
            sentinel = target / "keep.txt"
            sentinel.write_text("preserve me\n", encoding="utf-8")

            env = os.environ.copy()
            env.update({"FRAILTY_PYTHON": str(requested_python), "FRAILTY_VENV_DIR": str(target)})
            result = subprocess.run(
                ["bash", str(SETUP), str(root)], capture_output=True, text=True, env=env, check=False
            )

            self.assertEqual(result.returncode, 1)
            self.assertIn("not a Python venv", result.stderr)
            self.assertEqual(sentinel.read_text(encoding="utf-8"), "preserve me\n")


if __name__ == "__main__":
    unittest.main()
