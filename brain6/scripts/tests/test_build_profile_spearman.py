import sys
from pathlib import Path
import hashlib

import pytest

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT / "brain6" / "scripts"))

from build_profile_spearman import read_locked_profiles, render_table, write_exclusive


def make_fixture(tmp_path):
    config = tmp_path / "family.yaml"
    config.write_text(
        "family_name: original_locked_atlas_396\n"
        "recompute_brain6_fdr: false\n"
        "sleep_traits: [a, b, c, d]\n"
        "brain_disorders: [x, y, z]\n"
    )
    (tmp_path / "family.yaml.sha256").write_text(
        hashlib.sha256(config.read_bytes()).hexdigest() + "  family.yaml\n"
    )
    global_map = tmp_path / "global.tsv"
    global_map.write_text(
        "sleep_trait\tbrain_disorder\trg\n"
        "a\tx\t0.1\na\ty\t0.2\na\tz\t0.4\n"
        "b\tx\t0.2\nb\ty\t0.3\nb\tz\t0.3\n"
        "c\tx\t0.3\nc\ty\t0.4\nc\tz\t0.2\n"
        "d\tx\t0.4\nd\ty\t0.1\nd\tz\t0.1\n"
    )
    return config, global_map


def test_spearman_profiles_are_derived_from_config_order(tmp_path):
    config, global_map = make_fixture(tmp_path)
    sleeps, disorders, profiles = read_locked_profiles(config, global_map)
    assert sleeps == ["a", "b", "c", "d"]
    assert disorders == ["x", "y", "z"]
    assert list(profiles["x"]) == sleeps
    table = render_table(config, global_map)
    assert table.splitlines()[0] == "brain_disorder\tx\ty\tz"
    assert "x\t1\t-0.2\t-1" in table


def test_spearman_rejects_incomplete_or_duplicate_cells(tmp_path):
    config, global_map = make_fixture(tmp_path)
    global_map.write_text(global_map.read_text().replace("d\tz\t0.1\n", ""))
    with pytest.raises(ValueError, match="dimensions"):
        render_table(config, global_map)


def test_spearman_rejects_modified_locked_config(tmp_path):
    config, global_map = make_fixture(tmp_path)
    config.write_text(config.read_text() + "# modified\n")
    with pytest.raises(ValueError, match="checksum mismatch"):
        render_table(config, global_map)


def test_exclusive_write_preserves_existing_output(tmp_path):
    path = tmp_path / "out.tsv"
    write_exclusive(path, "original\n")
    with pytest.raises(FileExistsError):
        write_exclusive(path, "replacement\n")
    assert path.read_text() == "original\n"
