#!/usr/bin/env python3
"""Resume the locked frailty LAVA family with atomic per-locus claims."""
from __future__ import annotations

import argparse
import importlib.util
import json
import os
import re
import subprocess
import sys
import tempfile
import time
from datetime import datetime, timezone
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
spec = importlib.util.spec_from_file_location("frailty_lava_serial", ROOT / "frailty_paper/scripts/45_run_lava_sensitivity.py")
serial = importlib.util.module_from_spec(spec)
assert spec and spec.loader
spec.loader.exec_module(serial)
OUTPUT_ROOT = serial.OUTPUT_ROOT
EXPECTED = 29940
WORKER_COUNT = 6
SUPPORTED_WORKER_COUNTS = {2, 3, 4, 6}
# Keep the six-worker trial requested by the intervention, then fall back to
# four before the previously observed ~84.7% swap failure point. Four workers
# remain the target fallback unless the host reaches the deeper emergency rung.
SWAP_FALLBACK_FRACTION = 0.82
DEEP_SWAP_FALLBACK_FRACTION = 0.95
# macOS grows the swapfile pool dynamically; keep absolute caps as a secondary
# guard, but set them above the current baseline so six can be measured safely.
SWAP_FALLBACK_USED_MIB = 5000.0
DEEP_SWAP_FALLBACK_USED_MIB = 5800.0
SEVERE_SWAP_FRACTION = 0.95
# Historical measurements show four-worker periods at 508–622 receipts/hour,
# while the initial six-worker interval produced about 392/hour. Treat a
# sustained 20% shortfall from the conservative four-worker floor as material;
# use a 30-minute window to avoid reacting to variable per-locus runtimes.
FOUR_WORKER_BASELINE_RECEIPTS_PER_HOUR = 508.0
THROUGHPUT_FALLBACK_FRACTION = 0.80
THROUGHPUT_WINDOW_SECONDS = 30 * 60


def now() -> str:
    return datetime.now(timezone.utc).isoformat()


def atomic_json(path: Path, value: dict) -> None:
    temp = path.with_name(f"{path.name}.{os.getpid()}.partial")
    temp.write_text(json.dumps(value, sort_keys=True, indent=2) + "\n", encoding="utf-8")
    os.replace(temp, path)


def swap_usage_mib() -> tuple[float, float] | None:
    """Return (used, total) swap MiB on macOS; unknown elsewhere is non-fatal."""
    try:
        result = subprocess.run(["sysctl", "vm.swapusage"], check=True, capture_output=True,
                                text=True, timeout=5)
    except (OSError, subprocess.SubprocessError):
        return None
    match = re.search(r"total\s*=\s*([\d.]+)M\s+used\s*=\s*([\d.]+)M", result.stdout)
    if not match:
        return None
    total, used = map(float, match.groups())
    return (used, total) if total > 0 else None


def maybe_reduce_workers(control_path: Path, control: dict, desired_workers: int,
                         failures: int, swap_usage: tuple[float, float] | None,
                         observed_receipts_per_hour: float | None = None) -> tuple[int, str | None]:
    """Step concurrency down only at safe locus boundaries and record every rung."""
    if desired_workers <= 2:
        return desired_workers, None

    reason = None
    if failures >= 2:
        reason = f"repeated worker failures ({failures})"
    elif desired_workers == 6 and observed_receipts_per_hour is not None:
        minimum_rate = FOUR_WORKER_BASELINE_RECEIPTS_PER_HOUR * THROUGHPUT_FALLBACK_FRACTION
        if observed_receipts_per_hour < minimum_rate:
            reason = (f"measured throughput {observed_receipts_per_hour:.0f} receipts/hour over the "
                      f"30-minute window fell below {minimum_rate:.0f}/hour (80% of the conservative "
                      f"four-worker baseline {FOUR_WORKER_BASELINE_RECEIPTS_PER_HOUR:.0f}/hour); "
                      "reducing six workers to four")
    if reason is None and swap_usage and swap_usage[1] > 0:
        fraction = swap_usage[0] / swap_usage[1]
        if desired_workers == 6 and (
            fraction >= SWAP_FALLBACK_FRACTION or swap_usage[0] >= SWAP_FALLBACK_USED_MIB
        ):
            trigger = (f"≥{SWAP_FALLBACK_FRACTION:.0%} of allocated swap" if fraction >= SWAP_FALLBACK_FRACTION
                       else f"≥{SWAP_FALLBACK_USED_MIB:.0f} MiB used")
            reason = (f"swap use reached {swap_usage[0]:.0f}/{swap_usage[1]:.0f} MiB ({trigger}); "
                      "reducing six workers to four")
        elif desired_workers == 4 and (
            fraction >= DEEP_SWAP_FALLBACK_FRACTION or swap_usage[0] >= DEEP_SWAP_FALLBACK_USED_MIB
        ):
            trigger = (f"≥{DEEP_SWAP_FALLBACK_FRACTION:.0%} of allocated swap" if fraction >= DEEP_SWAP_FALLBACK_FRACTION
                       else f"≥{DEEP_SWAP_FALLBACK_USED_MIB:.0f} MiB used")
            reason = (f"swap use reached {swap_usage[0]:.0f}/{swap_usage[1]:.0f} MiB ({trigger}); "
                      "reducing four workers to three")
        elif desired_workers == 3 and (
            fraction >= DEEP_SWAP_FALLBACK_FRACTION or swap_usage[0] >= DEEP_SWAP_FALLBACK_USED_MIB
        ):
            trigger = (f"≥{DEEP_SWAP_FALLBACK_FRACTION:.0%} of allocated swap" if fraction >= DEEP_SWAP_FALLBACK_FRACTION
                       else f"≥{DEEP_SWAP_FALLBACK_USED_MIB:.0f} MiB used")
            reason = (f"swap use reached {swap_usage[0]:.0f}/{swap_usage[1]:.0f} MiB ({trigger}); "
                      "reducing three workers to two")
    if not reason:
        return desired_workers, None

    next_workers = 4 if desired_workers == 6 else desired_workers - 1
    recorded_utc = now()
    try:
        latest_control = json.loads(control_path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError):
        latest_control = dict(control)
    if not isinstance(latest_control, dict):
        latest_control = dict(control)
    history = list(latest_control.get("resource_fallback_events", control.get("resource_fallback_events", [])))
    history.append({"from_workers": desired_workers, "to_workers": next_workers,
                    "reason": reason, "utc": recorded_utc})
    latest_control["desired_workers"] = next_workers
    latest_control["resource_fallback_reason"] = reason
    latest_control["resource_fallback_utc"] = recorded_utc
    latest_control["resource_fallback_events"] = history
    control.clear()
    control.update(latest_control)
    atomic_json(control_path, control)
    return next_workers, reason


def launch_resource_block(swap_usage: tuple[float, float] | None) -> str | None:
    """Do not release pause holds when the host is already critically swap-loaded."""
    if swap_usage and swap_usage[1] > 0 and swap_usage[0] / swap_usage[1] >= SEVERE_SWAP_FRACTION:
        return (f"pre-existing swap use {swap_usage[0]:.0f}/{swap_usage[1]:.0f} MiB "
                f"(≥{SEVERE_SWAP_FRACTION:.0%}); pause holds retained")
    return None


def process_group_alive(pid: int) -> bool:
    try:
        os.killpg(pid, 0)
        return True
    except ProcessLookupError:
        return False
    except PermissionError:
        return True


def reap_worker_processes(processes: dict[int, subprocess.Popen]) -> None:
    """Poll children started by this coordinator so exited workers are reaped."""
    for process in processes.values():
        process.poll()


def recover_stale_claims(claim_root: Path, live_worker_groups: set[int] | None = None) -> int:
    recovered = 0
    for claim in sorted(claim_root.glob("*/*.claim")):
        try:
            owner = json.loads(claim.read_text(encoding="utf-8"))
            pgid = int(owner["process_group"])
        except (OSError, ValueError, KeyError, json.JSONDecodeError):
            # O_EXCL publication is supported on the campaign's external volume,
            # but unlike hard-link publication it exposes the file while the
            # owner payload is being written. Keep malformed claims while any
            # registered worker can still be the writer; a later clean restart
            # can remove them after every worker group is confirmed dead.
            if live_worker_groups:
                continue
            claim.unlink(missing_ok=True)
            recovered += 1
            continue
        if owner.get("pause_hold"):
            # Pause holds are scheduling barriers, not in-flight work. They are
            # released only by explicit resume after the coordinator has stopped.
            continue
        if not process_group_alive(pgid):
            claim.unlink(missing_ok=True)
            recovered += 1
    return recovered


def release_pause_claims(output_root: Path, control_path: Path, control: dict) -> int:
    """Release only explicit pause holds after proving the runner is stopped."""
    if (output_root / "runner.lock").exists():
        raise RuntimeError("cannot resume while a coordinator lock exists")
    claim_root = output_root / "claims"
    non_pause = []
    pause_claims = []
    for claim in claim_root.glob("*/*.claim"):
        if claim.name.startswith("._"):
            continue
        try:
            value = json.loads(claim.read_text(encoding="utf-8"))
        except (OSError, json.JSONDecodeError) as error:
            raise RuntimeError(f"cannot verify claim before resume: {claim}") from error
        if value.get("pause_hold"):
            pause_claims.append(claim)
        else:
            non_pause.append(claim)
    if non_pause:
        raise RuntimeError(f"cannot resume with {len(non_pause)} non-pause claims present")
    for claim in pause_claims:
        claim.unlink()
    control["pause_requested"] = False
    control["desired_workers"] = 4
    control["resumed_utc"] = now()
    control["pause_holds_released"] = len(pause_claims)
    atomic_json(control_path, control)
    return len(pause_claims)


def jobs() -> list[tuple[int, str, int]]:
    return [(order, trait, locus) for locus in range(1, 2496)
            for order, trait in enumerate(serial.TRAITS, 1)]


def claim_job(claim_root: Path, pair_order: int, trait: str, locus: int) -> Path | None:
    path = claim_root / trait / f"locus_{locus:04d}.claim"
    path.parent.mkdir(parents=True, exist_ok=True)
    payload = {"trait": trait, "pair_order": pair_order, "locus_index": locus,
               "worker_pid": os.getpid(), "process_group": os.getpgrp(), "claimed_utc": now()}
    try:
        fd = os.open(path, os.O_CREAT | os.O_EXCL | os.O_WRONLY, 0o600)
    except FileExistsError:
        return None
    data = (json.dumps(payload, sort_keys=True) + "\n").encode("utf-8")
    try:
        offset = 0
        while offset < len(data):
            offset += os.write(fd, data[offset:])
        os.fsync(fd)
    finally:
        os.close(fd)
    return path


def wait_for_worker_registration(worker_index: int, timeout: float = 60.0) -> bool:
    """Do not let a worker claim work before the coordinator records its PGID."""
    record = OUTPUT_ROOT / "workers" / f"worker_{worker_index:02d}.json"
    deadline = time.monotonic() + timeout
    while time.monotonic() < deadline:
        try:
            value = json.loads(record.read_text(encoding="utf-8"))
            if (int(value["pid"]) == os.getpid()
                    and int(value["process_group"]) == os.getpgrp()):
                return True
        except (OSError, ValueError, KeyError, json.JSONDecodeError):
            pass
        time.sleep(0.1)
    return False


def prepare_worker(trait: str, info_root: Path) -> tuple[Path, Path]:
    pair_root = OUTPUT_ROOT / "pairs" / trait
    pair_root.mkdir(parents=True, exist_ok=True)
    info_root.mkdir(parents=True, exist_ok=True)
    return pair_root, info_root


def worker(worker_index: int, event_path: Path) -> int:
    if not wait_for_worker_registration(worker_index):
        raise RuntimeError(f"worker {worker_index} was not registered by its coordinator")
    lock = serial.read_lock()
    manifest = json.loads((OUTPUT_ROOT / "input_manifest.json").read_text(encoding="utf-8"))
    locus_rows = []
    with serial.LOCUS_FILE.open(encoding="utf-8") as handle:
        header = handle.readline().split()
        locus_rows = [dict(zip(header, line.split())) for line in handle if line.strip()]
    scratch = Path(tempfile.mkdtemp(prefix=f"frailtylava_worker{worker_index:02d}_"))
    (scratch / "inputs").symlink_to(OUTPUT_ROOT / "inputs", target_is_directory=True)
    claim_root = OUTPUT_ROOT / "claims"
    queue = jobs()
    # Rotate each worker's deterministic scan start to spread first claims across traits/loci.
    start = (worker_index * len(queue)) // WORKER_COUNT
    queue = queue[start:] + queue[:start]
    failures = 0
    try:
        for pair_order, trait, locus_index in queue:
            # Pause and scale-down requests take effect only between loci.
            try:
                control = json.loads((OUTPUT_ROOT / "concurrency.json").read_text(encoding="utf-8"))
            except (OSError, json.JSONDecodeError):
                control = {}
            if control.get("pause_requested") or worker_index > int(control.get("desired_workers", WORKER_COUNT)):
                break
            receipt = OUTPUT_ROOT / "pairs" / trait / "loci" / f"locus_{locus_index:04d}.json"
            if receipt.exists():
                serial.verify_receipt(receipt, trait, pair_order, locus_index)
                continue
            claim = claim_job(claim_root, pair_order, trait, locus_index)
            if claim is None:
                continue
            try:
                if receipt.exists():
                    serial.verify_receipt(receipt, trait, pair_order, locus_index)
                    continue
                pair_root, info_root = prepare_worker(trait, scratch / "input_info" / trait)
                chromosome = int(locus_rows[locus_index - 1]["CHR"])
                original_info = pair_root / f"input_info_chr{chromosome:02d}.tsv"
                info = info_root / f"input_info_chr{chromosome:02d}.tsv"
                import csv
                with original_info.open(encoding="utf-8", newline="") as source:
                    rows = list(csv.DictReader(source, delimiter="\t"))
                if [row["phenotype"] for row in rows] != ["frailty", trait]:
                    raise RuntimeError(f"unexpected pair input.info rows: {original_info}")
                for row in rows:
                    row["filename"] = f"{row['phenotype']}/chr{chromosome:02d}.sumstats.tsv.gz"
                with info.open("w", encoding="utf-8", newline="") as output:
                    writer = csv.DictWriter(output, fieldnames=["phenotype", "cases", "controls", "prevalence", "filename"], delimiter="\t", lineterminator="\n")
                    writer.writeheader()
                    writer.writerows(rows)
                log = pair_root / "logs" / f"locus_{locus_index:04d}.log"
                log.parent.mkdir(parents=True, exist_ok=True)
                command = [str(serial.R_SCRIPT), str(serial.WORKER), trait, str(pair_order), str(locus_index),
                           str(info), manifest["pairs"][trait]["overlap_file"], str(serial.LOCUS_FILE),
                           lock["reference_root"], str(scratch / "inputs"), str(receipt)]
                error = "worker did not produce a receipt"
                for attempt in (1, 2):
                    with log.open("a", encoding="utf-8") as stream:
                        stream.write(f"\nPARALLEL_WORKER={worker_index} ATTEMPT={attempt} {now()}\n")
                        result = subprocess.run(command, stdout=stream, stderr=subprocess.STDOUT, text=True)
                    if receipt.is_file():
                        serial.verify_receipt(receipt, trait, pair_order, locus_index)
                        error = ""
                        break
                    error = f"worker exit {result.returncode}; see {log}"
                if error:
                    failures += 1
                    record = {"worker": worker_index, "trait": trait, "pair_order": pair_order,
                              "locus_index": locus_index, "error": error, "utc": now()}
                    with (pair_root / "launch_failures.jsonl").open("a", encoding="utf-8") as handle:
                        handle.write(json.dumps(record, sort_keys=True) + "\n")
                    with event_path.open("a", encoding="utf-8") as handle:
                        handle.write(json.dumps({"type": "failure", **record}, sort_keys=True) + "\n")
                else:
                    with event_path.open("a", encoding="utf-8") as handle:
                        handle.write(json.dumps({"type": "receipt", "worker": worker_index,
                                                 "trait": trait, "locus_index": locus_index, "utc": now()}) + "\n")
            finally:
                claim.unlink(missing_ok=True)
    finally:
        with event_path.open("a", encoding="utf-8") as handle:
            handle.write(json.dumps({"type": "worker_finished", "worker": worker_index,
                                     "pid": os.getpid(), "utc": now()}) + "\n")
        import shutil
        shutil.rmtree(scratch, ignore_errors=True)
    return 2 if failures else 0


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--worker", type=int)
    parser.add_argument("--workers", type=int)
    parser.add_argument("--resume-paused", action="store_true",
                        help="explicitly resume a safely paused family after verifying no active runner")
    args = parser.parse_args()
    if args.worker is not None:
        return worker(args.worker, OUTPUT_ROOT / "parallel_events.jsonl")
    control_path = OUTPUT_ROOT / "concurrency.json"
    control = json.loads(control_path.read_text(encoding="utf-8")) if control_path.is_file() else {}
    if control.get("pause_requested") and not args.resume_paused:
        print("PARALLEL_PAUSED: explicit --resume-paused is required; no workers launched", flush=True)
        return 0
    desired_workers = args.workers if args.workers is not None else int(control.get("desired_workers", WORKER_COUNT))
    if desired_workers not in SUPPORTED_WORKER_COUNTS:
        parser.error("supported worker counts are six, four, three, or two")
    initial_requested_workers = desired_workers
    swap_snapshot = swap_usage_mib()
    desired_workers, initial_fallback = maybe_reduce_workers(
        control_path, control, desired_workers, 0, swap_snapshot)
    launch_block = launch_resource_block(swap_snapshot)
    if launch_block:
        print(f"PARALLEL_RESOURCE_PAUSED: {launch_block}; no workers launched", flush=True)
        return 0
    if args.resume_paused:
        if not control.get("pause_requested"):
            parser.error("--resume-paused requires a recorded pause request")
        release_pause_claims(OUTPUT_ROOT, control_path, control)
    control["desired_workers"] = desired_workers
    # A previous invocation's downshift reason remains in the event history;
    # it must not be reported as the active fallback after a clean restart.
    control["resource_fallback_reason"] = initial_fallback
    if initial_fallback is None:
        control.pop("resource_fallback_utc", None)
    atomic_json(control_path, control)
    lock = serial.read_lock()
    manifest = json.loads((OUTPUT_ROOT / "input_manifest.json").read_text(encoding="utf-8"))
    if manifest.get("analysis_id") != lock["analysis_id"] or manifest.get("lock_sha256") != serial.sha256(serial.LOCK):
        raise RuntimeError("prepared inputs do not match the locked analysis")
    if manifest.get("pair_count") != 12 or manifest.get("locus_count") != 2495 or manifest.get("candidate_pair_locus_slots") != EXPECTED:
        raise RuntimeError("prepared family dimensions drifted")
    subprocess.run([sys.executable, str(serial.PREPARE), "--verify"], check=True)
    validator = subprocess.run([str(serial.R_SCRIPT), str(serial.RUNTIME_VALIDATOR)], check=True, text=True,
                               stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    if "status=PASS" not in validator.stdout or "lava_commit=e729a245f7b6923967a96804fbf5246eadf2d6c6" not in validator.stdout:
        raise RuntimeError("pinned LAVA runtime validation failed")
    with serial.LOCUS_FILE.open(encoding="utf-8") as handle:
        header = handle.readline().split()
        locus_rows = [line for line in handle if line.strip()]
    if header != ["LOC", "CHR", "START", "STOP"] or len(locus_rows) != 2495:
        raise RuntimeError("frozen locus definition failed validation")
    OUTPUT_ROOT.mkdir(parents=True, exist_ok=True)
    global_lock = OUTPUT_ROOT / "runner.lock"
    workers_dir = OUTPUT_ROOT / "workers"
    workers_dir.mkdir(exist_ok=True)
    registered: dict[int, dict] = {}
    for record in workers_dir.glob("worker_*.json"):
        try:
            value = json.loads(record.read_text(encoding="utf-8"))
            wid = int(value["worker"])
            pgid = int(value["process_group"])
        except (OSError, ValueError, KeyError, json.JSONDecodeError):
            continue
        if process_group_alive(pgid):
            registered[wid] = value
    if len(registered) > desired_workers and not control.get("resource_fallback_reason"):
        raise RuntimeError(f"found {len(registered)} active worker groups, above requested {desired_workers}; retire excess workers at safe locus boundaries first")
    if global_lock.exists():
        try:
            owner = int(global_lock.read_text(encoding="utf-8").split("pid=", 1)[1].split()[0])
            os.kill(owner, 0)
        except ProcessLookupError:
            global_lock.unlink()
        except PermissionError:
            raise RuntimeError("cannot prove existing global runner is stale")
        else:
            raise RuntimeError(f"runner PID {owner} is still alive")
    try:
        fd = os.open(global_lock, os.O_CREAT | os.O_EXCL | os.O_WRONLY, 0o600)
    except FileExistsError as error:
        raise RuntimeError("another coordinator claimed runner.lock") from error
    os.write(fd, f"pid={os.getpid()} started={now()}\n".encode())
    os.close(fd)
    event_path = OUTPUT_ROOT / "parallel_events.jsonl"
    event_path.touch(exist_ok=True)
    claim_root = OUTPUT_ROOT / "claims"
    claim_root.mkdir(exist_ok=True)
    live_worker_groups = {int(value["process_group"]) for value in registered.values()}
    recovered = recover_stale_claims(claim_root, live_worker_groups)
    verified = 0
    for pair_order, trait in enumerate(serial.TRAITS, 1):
        for locus_index in range(1, 2496):
            receipt = OUTPUT_ROOT / "pairs" / trait / "loci" / f"locus_{locus_index:04d}.json"
            if receipt.exists():
                serial.verify_receipt(receipt, trait, pair_order, locus_index)
                verified += 1
    # Establish the receipt-event offset and rate window before any workers
    # start, so fast first jobs are counted and the denominator includes their
    # actual runtime rather than beginning after they have already completed.
    last_event_size = event_path.stat().st_size
    throughput_window_started = time.monotonic()
    throughput_window_start_count = verified
    worker_pids: dict[str, int] = {str(wid): int(value["pid"]) for wid, value in registered.items()}
    launched_processes: dict[int, subprocess.Popen] = {}
    for worker_id in range(1, desired_workers + 1):
        if worker_id in registered:
            continue
        log = (workers_dir / f"worker_{worker_id:02d}.log").open("a", encoding="utf-8")
        process = subprocess.Popen([sys.executable, str(Path(__file__).resolve()), "--worker", str(worker_id)],
                                   stdout=log, stderr=subprocess.STDOUT, start_new_session=True, close_fds=True)
        launched_processes[process.pid] = process
        worker_pids[str(worker_id)] = process.pid
        atomic_json(workers_dir / f"worker_{worker_id:02d}.json", {
            "worker": worker_id, "pid": process.pid, "process_group": process.pid, "started_utc": now()})
        log.close()
    atomic_json(OUTPUT_ROOT / "runner_state.json", {
        "analysis_id": lock["analysis_id"], "lock_sha256": serial.sha256(serial.LOCK),
        "expected_receipts": EXPECTED, "verified_receipts": verified,
        "launch_failures_this_invocation": 0, "parallel_workers": len(worker_pids),
        "requested_parallel_workers": initial_requested_workers,
        "effective_parallel_workers": desired_workers,
        "resource_fallback_reason": initial_fallback,
        "worker_pids": worker_pids,
        "stale_claims_recovered": recovered,
        "measured_receipts_per_hour": None, "throughput_window_seconds": 0.0,
        "started_utc": now(), "updated_utc": now()})
    print(f"PARALLEL_LAUNCHED workers={len(worker_pids)} pids={worker_pids} verified={verified}/{EXPECTED} stale_claims_recovered={recovered}", flush=True)
    try:
        completed = 0
        failures = 0
        measured_receipts_per_hour = None
        throughput_window_seconds = 0.0
        while True:
            time.sleep(10)
            reap_worker_processes(launched_processes)
            live_group_ids = set()
            for record in workers_dir.glob("worker_*.json"):
                try:
                    value = json.loads(record.read_text(encoding="utf-8"))
                    pgid = int(value["process_group"])
                    if process_group_alive(pgid):
                        live_group_ids.add(pgid)
                except (OSError, ValueError, KeyError, json.JSONDecodeError):
                    continue
            recovered += recover_stale_claims(claim_root, live_group_ids)
            current_size = event_path.stat().st_size
            if current_size > last_event_size:
                with event_path.open(encoding="utf-8") as stream:
                    stream.seek(last_event_size)
                    for line in stream:
                        try:
                            event = json.loads(line)
                        except json.JSONDecodeError:
                            continue
                        if event.get("type") == "receipt":
                            completed += 1
                        elif event.get("type") == "failure":
                            failures += 1
                last_event_size = current_size
            verified += completed
            completed = 0
            elapsed_window = time.monotonic() - throughput_window_started
            rate_for_window = None
            if elapsed_window >= THROUGHPUT_WINDOW_SECONDS:
                measured_receipts_per_hour = (
                    (verified - throughput_window_start_count) * 3600.0 / elapsed_window)
                throughput_window_seconds = elapsed_window
                rate_for_window = measured_receipts_per_hour
                throughput_window_started = time.monotonic()
                throughput_window_start_count = verified
            desired_workers, _ = maybe_reduce_workers(
                control_path, control, desired_workers, failures,
                swap_usage_mib() if desired_workers > 2 else None,
                observed_receipts_per_hour=rate_for_window)
            active_jobs = []
            for claim in sorted(claim_root.glob("*/*.claim")):
                if claim.name.startswith("._"):
                    continue
                try:
                    active_jobs.append(json.loads(claim.read_text(encoding="utf-8")))
                except (OSError, json.JSONDecodeError):
                    continue
            live = {}
            for record in workers_dir.glob("worker_*.json"):
                try:
                    value = json.loads(record.read_text(encoding="utf-8"))
                    wid = int(value["worker"])
                    if process_group_alive(int(value["process_group"])):
                        live[str(wid)] = int(value["pid"])
                except (OSError, ValueError, KeyError, json.JSONDecodeError):
                    continue
            atomic_json(OUTPUT_ROOT / "runner_state.json", {
                "analysis_id": lock["analysis_id"], "lock_sha256": serial.sha256(serial.LOCK),
                "expected_receipts": EXPECTED, "verified_receipts": min(verified, EXPECTED),
                "launch_failures_this_invocation": failures, "parallel_workers": len(live),
                "requested_parallel_workers": initial_requested_workers,
                "effective_parallel_workers": desired_workers,
                "resource_fallback_reason": control.get("resource_fallback_reason", initial_fallback),
                "worker_pids": worker_pids, "active_jobs": active_jobs,
                "stale_claims_recovered": recovered,
                "measured_receipts_per_hour": measured_receipts_per_hour,
                "throughput_window_seconds": throughput_window_seconds,
                "updated_utc": now()})
            if not live:
                print(f"PARALLEL_FINISHED verified={min(verified, EXPECTED)}/{EXPECTED} failures={failures}", flush=True)
                break
        return 0
    finally:
        global_lock.unlink(missing_ok=True)


if __name__ == "__main__":
    raise SystemExit(main())
