"""Build a frozen, descriptive leave-one-sleep-trait-out profile sensitivity."""
from __future__ import annotations

import argparse
import hashlib
import json
import math
import os
import tempfile
from pathlib import Path
from typing import Any

from build_profile_spearman import read_locked_profiles, sha256
from scipy.stats import pearsonr

ROOT = Path(__file__).resolve().parents[2]
DEFAULT_PLAN = ROOT / "brain6/config/brain6_profile_leave_one_out_v1.json"
DEFAULT_CONFIG = ROOT / "brain6/config/brain6_locked_family.yaml"
DEFAULT_MAP = ROOT / "brain6/results/global/brain6_72_locked.tsv"
DEFAULT_OUTPUT = ROOT / "brain6/results/global/disorder_profile_leave_one_sleep_trait_out.tsv"
FIELDS = (
    "analysis_id", "claim_id", "disorder_a", "disorder_b", "omitted_sleep_trait",
    "n_sleep_traits", "full_profile_r", "leave_one_out_r", "delta_from_full",
    "leave_one_out_min_r", "leave_one_out_max_r", "leave_one_out_range_width",
)


def read_plan(path: Path) -> dict[str, Any]:
    checksum_path = path.with_suffix(path.suffix + ".sha256")
    if not checksum_path.is_file():
        raise ValueError(f"frozen sensitivity-plan checksum is missing: {checksum_path}")
    expected = checksum_path.read_text(encoding="utf-8").split()[0]
    if sha256(path) != expected:
        raise ValueError("frozen sensitivity-plan checksum mismatch")
    plan = json.loads(path.read_text(encoding="utf-8"))
    if plan.get("status") != "FROZEN_BEFORE_COMPUTATION":
        raise ValueError("profile sensitivity plan is not frozen")
    if plan.get("method", {}).get("inferential_p_values") is not False:
        raise ValueError("profile sensitivity plan must prohibit inferential p-values")
    return plan


def correlation(left: list[float], right: list[float]) -> float:
    if len(left) != len(right) or len(left) < 3:
        raise ValueError("Pearson profile correlation requires at least three paired traits")
    if len(set(left)) < 2 or len(set(right)) < 2:
        raise ValueError("degenerate profile: Pearson correlation is undefined")
    value = float(pearsonr(left, right).statistic)
    if not math.isfinite(value):
        raise ValueError("degenerate profile: Pearson correlation is not finite")
    return value


def render_table(plan: dict[str, Any], config_path: Path, map_path: Path) -> tuple[str, dict[str, Any]]:
    sleeps, disorders, profiles = read_locked_profiles(config_path, map_path)
    method = plan["method"]
    if len(sleeps) != int(method["expected_sleep_trait_count"]):
        raise ValueError("locked sleep-trait count differs from the frozen sensitivity plan")
    expected_remaining = len(sleeps) - 1
    if expected_remaining != int(method["traits_per_leave_one_out_profile"]):
        raise ValueError("leave-one-out trait count differs from the frozen sensitivity plan")

    pairs = plan.get("profiles")
    if not isinstance(pairs, list) or not pairs:
        raise ValueError("frozen sensitivity plan has no profile pairs")
    pair_ids: set[tuple[str, str]] = set()
    output_rows: list[dict[str, str]] = []
    summaries: list[dict[str, Any]] = []
    for pair in pairs:
        a, b, claim_id = pair["disorder_a"], pair["disorder_b"], pair["claim_id"]
        if a == b or a not in disorders or b not in disorders:
            raise ValueError(f"invalid disorder pair in frozen plan: {a}, {b}")
        if (a, b) in pair_ids or (b, a) in pair_ids:
            raise ValueError(f"duplicate profile pair in frozen plan: {a}, {b}")
        pair_ids.add((a, b))
        left = [profiles[a][trait] for trait in sleeps]
        right = [profiles[b][trait] for trait in sleeps]
        full = correlation(left, right)
        omitted_values: list[tuple[str, float]] = []
        for index, omitted_trait in enumerate(sleeps):
            loo_left = left[:index] + left[index + 1:]
            loo_right = right[:index] + right[index + 1:]
            omitted_values.append((omitted_trait, correlation(loo_left, loo_right)))
        loo_coefficients = [value for _, value in omitted_values]
        minimum, maximum = min(loo_coefficients), max(loo_coefficients)
        width = maximum - minimum
        summaries.append({
            "claim_id": claim_id,
            "disorder_a": a,
            "disorder_b": b,
            "full_profile_r": full,
            "leave_one_out_min_r": minimum,
            "leave_one_out_max_r": maximum,
            "leave_one_out_range_width": width,
        })
        for trait, value in omitted_values:
            output_rows.append({
                "analysis_id": plan["analysis_id"],
                "claim_id": claim_id,
                "disorder_a": a,
                "disorder_b": b,
                "omitted_sleep_trait": trait,
                "n_sleep_traits": str(expected_remaining),
                "full_profile_r": f"{full:.12g}",
                "leave_one_out_r": f"{value:.12g}",
                "delta_from_full": f"{value - full:.12g}",
                "leave_one_out_min_r": f"{minimum:.12g}",
                "leave_one_out_max_r": f"{maximum:.12g}",
                "leave_one_out_range_width": f"{width:.12g}",
            })

    lines = ["\t".join(FIELDS)]
    lines.extend("\t".join(row[field] for field in FIELDS) for row in output_rows)
    return "\n".join(lines) + "\n", {"profiles": summaries, "row_count": len(output_rows)}


def write_immutable(path: Path, payload: bytes) -> str:
    """Create a hash-bound file atomically; identical reruns are idempotent."""
    path.parent.mkdir(parents=True, exist_ok=True)
    if path.exists():
        if path.read_bytes() != payload:
            raise ValueError(f"refusing to replace different frozen output: {path}")
        return sha256(path)
    fd, name = tempfile.mkstemp(prefix=f".{path.name}.", suffix=".partial", dir=path.parent)
    temp = Path(name)
    try:
        with os.fdopen(fd, "wb") as stream:
            stream.write(payload)
            stream.flush()
            os.fsync(stream.fileno())
        try:
            os.link(temp, path)
        except FileExistsError:
            if path.read_bytes() != payload:
                raise ValueError(f"refusing to replace different frozen output: {path}")
        directory_fd = os.open(path.parent, os.O_RDONLY)
        try:
            os.fsync(directory_fd)
        finally:
            os.close(directory_fd)
        return sha256(path)
    finally:
        temp.unlink(missing_ok=True)


def build(plan_path: Path, config_path: Path, map_path: Path, output_path: Path) -> dict[str, Any]:
    plan = read_plan(plan_path)
    table, summary = render_table(plan, config_path, map_path)
    table_digest = hashlib.sha256(table.encode("utf-8")).hexdigest()
    provenance_path = output_path.with_suffix(".provenance.json")
    provenance = {
        "schema_version": 1,
        "status": "PASS",
        "analysis_id": plan["analysis_id"],
        "method": plan["method"]["statistic"],
        "interpretation": plan["reporting"]["interpretation"],
        "plan_sha256": sha256(plan_path),
        "locked_family_config_sha256": sha256(config_path),
        "global_map_sha256": sha256(map_path),
        "script_sha256": sha256(Path(__file__)),
        "output_sha256": table_digest,
        **summary,
    }
    output_digest = write_immutable(output_path, table.encode("utf-8"))
    provenance_digest = write_immutable(
        provenance_path,
        (json.dumps(provenance, indent=2, sort_keys=True) + "\n").encode("utf-8"),
    )
    return {"output": str(output_path), "output_sha256": output_digest,
            "provenance": str(provenance_path), "provenance_sha256": provenance_digest,
            **summary}


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--plan", type=Path, default=DEFAULT_PLAN)
    parser.add_argument("--config", type=Path, default=DEFAULT_CONFIG)
    parser.add_argument("--map", dest="map_path", type=Path, default=DEFAULT_MAP)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()
    print(json.dumps(build(args.plan, args.config, args.map_path, args.output), indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
