from __future__ import annotations

import importlib.util
import tempfile
import unittest
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import patch


SCRIPT = Path(__file__).resolve().parents[1] / "scripts/21_configure_external_storage.py"
SPEC = importlib.util.spec_from_file_location("configure_external_storage", SCRIPT)
STORAGE = importlib.util.module_from_spec(SPEC)
assert SPEC and SPEC.loader
SPEC.loader.exec_module(STORAGE)


class DerivedOutputStorageTests(unittest.TestCase):
    def make_layout(self, root: Path, output_root: Path) -> None:
        output_root.mkdir(parents=True)
        for relative in ("data/harmonized", "data/munged"):
            directory = root / relative
            directory.mkdir(parents=True)
            (directory / ".gitkeep").write_text("")

    def test_only_empty_placeholder_outputs_can_be_redirected(self) -> None:
        with tempfile.TemporaryDirectory() as temp:
            base = Path(temp)
            root, output_root = base / "repo", base / "external"
            self.make_layout(root, output_root)
            (root / "data/munged/existing.tsv.gz").write_bytes(b"preserve")
            with patch.object(STORAGE.shutil, "disk_usage", return_value=SimpleNamespace(free=21 * 1024**3)):
                with self.assertRaisesRegex(RuntimeError, "non-empty derived output"):
                    STORAGE.prepare_output_links(root, output_root)
            self.assertTrue((root / "data/munged/existing.tsv.gz").is_file())

    def test_redirect_preserves_empty_source_and_creates_external_links(self) -> None:
        with tempfile.TemporaryDirectory() as temp:
            base = Path(temp)
            root, output_root = base / "repo", base / "external"
            self.make_layout(root, output_root)
            with patch.object(STORAGE.shutil, "disk_usage", return_value=SimpleNamespace(free=21 * 1024**3)):
                specs = STORAGE.prepare_output_links(root, output_root)
            STORAGE.apply_output_links(specs)
            with patch.object(STORAGE.shutil, "disk_usage", return_value=SimpleNamespace(free=21 * 1024**3)):
                STORAGE.apply_output_links(STORAGE.prepare_output_links(root, output_root))
            for relative in ("data/harmonized", "data/munged"):
                linked = root / relative
                self.assertTrue(linked.is_symlink())
                self.assertEqual(linked.resolve(), (output_root / relative).resolve())
                self.assertTrue((linked.with_name(linked.name + ".local-preserved") / ".gitkeep").is_file())
                self.assertTrue((linked / ".gitkeep").is_file())


if __name__ == "__main__":
    unittest.main()
