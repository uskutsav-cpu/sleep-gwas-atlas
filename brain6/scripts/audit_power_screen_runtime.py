"""Record a read-only, receipt-aware snapshot of the four trait-only workers."""
from __future__ import annotations

import argparse
import hashlib
import json
import subprocess
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parents[2]
SCREEN = ROOT / "brain6/results/power_optimized_sensitivity_v1/lava_trait_screen_v1"
SCRIPT = ROOT / "brain6/scripts/run_power_optimized_sleep_univariate_v1.R"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def parse_processes(ps_text: str) -> dict[int, dict[str, Any]]:
    found: dict[int, dict[str, Any]] = {}
    for line in ps_text.splitlines():
        parts = line.strip().split(maxsplit=4)
        if len(parts) != 5:
            continue
        pid_text, cpu, rss, elapsed, command = parts
        if "run_power_optimized_sleep_univariate_v1.R" not in command:
            continue
        for worker_id in range(1, 5):
            config = f"worker_{worker_id}.config.json"
            if config in command:
                found[worker_id] = {
                    "pid": int(pid_text), "cpu_percent": float(cpu),
                    "rss_kib": int(rss), "elapsed": elapsed,
                    "command": command,
                }
    return found


def snapshot() -> dict[str, Any]:
    process_inventory_error: dict[str, str] | None = None
    try:
        ps = subprocess.run(
            ["ps", "-axo", "pid=,pcpu=,rss=,etime=,command="],
            check=True, capture_output=True, text=True,
        )
    except (OSError, subprocess.SubprocessError) as error:
        # Sandboxed hosts may deny process-table access. Preserve the output
        # inventory, but never turn unavailable process evidence into a claim
        # that workers stopped or are active.
        processes: dict[int, dict[str, Any]] = {}
        process_inventory_error = {
            "type": type(error).__name__,
            "message": str(error),
        }
    else:
        processes = parse_processes(ps.stdout)
    worker_rows = []
    for worker_id in range(1, 5):
        config = SCREEN / f"worker_{worker_id}.config.json"
        data = json.loads(config.read_text(encoding="utf-8"))
        result = ROOT / data["output_tsv"]
        summary = ROOT / data["summary_json"]
        worker_rows.append({
            "worker_id": worker_id,
            "config_path": str(config.relative_to(ROOT)),
            "config_sha256": sha256(config),
            "expected_loci": int(data["expected_loci"]),
            "pid": processes.get(worker_id, {}).get("pid"),
            "process": processes.get(worker_id),
            "output_path": str(result.relative_to(ROOT)),
            "output_exists": result.is_file(),
            "output_sha256": sha256(result) if result.is_file() else None,
            "summary_path": str(summary.relative_to(ROOT)),
            "summary_exists": summary.is_file(),
            "summary_sha256": sha256(summary) if summary.is_file() else None,
        })
    active = (
        None if process_inventory_error is not None
        else sum(row["pid"] is not None for row in worker_rows)
    )
    durable = sum(row["output_exists"] and row["summary_exists"] for row in worker_rows)
    status = (
        "PROCESS_INVENTORY_UNAVAILABLE" if process_inventory_error is not None
        else "RUNNING" if active == 4
        else "PARTIAL" if active or durable
        else "NOT_RUNNING_OR_UNVERIFIED"
    )
    return {
        "schema_version": 1,
        "generated_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "status": status,
        "analysis_id": "brain6_power_optimized_sensitivity_v1_longsleep_trait_only_screen",
        "script_path": str(SCRIPT.relative_to(ROOT)),
        "script_sha256": sha256(SCRIPT),
        "configured_workers": 4,
        "active_workers": active,
        "process_inventory_error": process_inventory_error,
        "workers_with_durable_tsv_and_summary": durable,
        "outputs_are_final_write_only": True,
        "workers": worker_rows,
        "interpretation": "Runtime observation only; this screen is diagnostic and does not alter canonical LAVA, family gates, or promotion decisions.",
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    output = args.output if args.output.is_absolute() else ROOT / args.output
    if output.exists():
        raise FileExistsError(f"Refusing to overwrite runtime audit: {output}")
    data = snapshot()
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(data, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps({"output": str(output), "status": data["status"],
                      "active_workers": data["active_workers"],
                      "durable_worker_outputs": data["workers_with_durable_tsv_and_summary"],
                      "sha256": sha256(output)}, indent=2))


if __name__ == "__main__":
    main()
