from __future__ import annotations

import csv
import hashlib
import json
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT / "brain6/scripts"))

from build_profile_leave_one_out import build, read_plan, render_table


def fixture_files(tmp_path: Path) -> tuple[Path, Path, Path]:
    config = tmp_path / "family.yaml"
    sleeps = [f"sleep_{index:02}" for index in range(12)]
    config.write_text(
        "family_name: original_locked_atlas_396\n"
        "recompute_brain6_fdr: false\n"
        f"sleep_traits: [{', '.join(sleeps)}]\n"
        "brain_disorders: [x, y, z]\n"
    )
    config.with_suffix(".yaml.sha256").write_text(
        hashlib.sha256(config.read_bytes()).hexdigest() + "  family.yaml\n"
    )
    global_map = tmp_path / "global.tsv"
    x = [index / 12 for index in range(12)]
    y = [((index * index + 3 * index) % 17) / 17 for index in range(12)]
    z = [((5 * index + index % 3) % 13) / 13 for index in range(12)]
    lines = ["sleep_trait\tbrain_disorder\trg"]
    for index, sleep in enumerate(sleeps):
        for disorder, values in (("x", x), ("y", y), ("z", z)):
            lines.append(f"{sleep}\t{disorder}\t{values[index]}\n")
    global_map.write_text("".join(line if line.endswith("\n") else line + "\n" for line in lines))

    plan = tmp_path / "plan.json"
    plan_data = {
        "analysis_id": "profile-loo-test",
        "status": "FROZEN_BEFORE_COMPUTATION",
        "method": {"expected_sleep_trait_count": 12, "traits_per_leave_one_out_profile": 11,
                   "inferential_p_values": False,
                   "statistic": "Pearson descriptive correlation"},
        "profiles": [
            {"disorder_a": "x", "disorder_b": "y", "claim_id": "xy"},
            {"disorder_a": "y", "disorder_b": "z", "claim_id": "yz"},
        ],
        "reporting": {"interpretation": "Descriptive only"},
    }
    plan.write_text(json.dumps(plan_data, sort_keys=True, indent=2) + "\n")
    plan.with_suffix(".json.sha256").write_text(
        hashlib.sha256(plan.read_bytes()).hexdigest() + "  plan.json\n"
    )
    return config, global_map, plan


def test_leave_one_out_covers_every_sleep_trait_and_emits_descriptive_ranges(tmp_path):
    config, global_map, plan_path = fixture_files(tmp_path)
    plan = read_plan(plan_path)
    table, summary = render_table(plan, config, global_map)
    rows = list(csv.DictReader(table.splitlines(), delimiter="\t"))
    assert len(rows) == 24
    assert len({(row["claim_id"], row["omitted_sleep_trait"]) for row in rows}) == 24
    assert {row["n_sleep_traits"] for row in rows} == {"11"}
    assert "p_value" not in rows[0] and "confidence_interval" not in rows[0]
    assert [row["claim_id"] for row in summary["profiles"]] == ["xy", "yz"]
    for group in summary["profiles"]:
        assert group["leave_one_out_min_r"] <= group["leave_one_out_max_r"]
        assert group["leave_one_out_range_width"] == pytest.approx(
            group["leave_one_out_max_r"] - group["leave_one_out_min_r"]
        )


def test_plan_checksum_and_locked_config_are_enforced(tmp_path):
    config, global_map, plan_path = fixture_files(tmp_path)
    plan_path.write_text(plan_path.read_text() + " ")
    with pytest.raises(ValueError, match="plan checksum mismatch"):
        read_plan(plan_path)
    fresh = tmp_path / "fresh"
    fresh.mkdir()
    config, global_map, plan_path = fixture_files(fresh)
    config.write_text(config.read_text() + "# changed\n")
    with pytest.raises(ValueError, match="configuration checksum mismatch"):
        render_table(read_plan(plan_path), config, global_map)


def test_build_is_reproducible_and_refuses_output_replacement(tmp_path):
    config, global_map, plan = fixture_files(tmp_path)
    output = tmp_path / "profile.tsv"
    first = build(plan, config, global_map, output)
    second = build(plan, config, global_map, output)
    assert first == second
    provenance = json.loads(output.with_suffix(".provenance.json").read_text())
    assert provenance["output_sha256"] == first["output_sha256"]
    assert provenance["plan_sha256"] == hashlib.sha256(plan.read_bytes()).hexdigest()
    output.write_text(output.read_text() + "changed\n")
    with pytest.raises(ValueError, match="refusing to replace different frozen output"):
        build(plan, config, global_map, output)
