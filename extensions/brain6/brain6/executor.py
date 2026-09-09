"""Synchronous native-tool jobs with receipts, no shell strings or silent skips."""
from __future__ import annotations
import os
import re
import shutil
import subprocess
from pathlib import Path
from .artifacts import fingerprint, transaction, verify_artifact
from .io import (ContractError, check_hash, ensure_free, exclusive_lock, read_json,
                 require, safe_id, safe_write_path, write_json)


def validate_job(job: dict) -> None:
    require(job.get("schema_version") == 1, "Unsupported job schema")
    safe_id(job["job_id"])
    require(isinstance(job.get("argv"), list) and len(job["argv"]) > 0 and
            all(isinstance(x, str) and x for x in job["argv"]), "argv must be a nonempty string list")
    require(Path(job["argv"][0]).name not in {"bash", "sh", "zsh", "cmd", "powershell"},
            "Shell execution is not supported; use direct executable argv")
    require(job.get("reviewed") is True, "Native job must be explicitly reviewed")
    require(job.get("outputs"), "Expected outputs are required")
    for group in ["inputs", "outputs"]:
        require(isinstance(job[group], dict), f"{group} must be keyed by stable labels")
        for name in job[group]:
            safe_id(name)
    for path in job["outputs"].values():
        require(not Path(path).is_absolute() and ".." not in Path(path).parts, "Unsafe output")
    require(len(set(job["outputs"].values())) == len(job["outputs"]), "Output path collision")
    require("semantic_status" in job["outputs"], "Native adapter must provide semantic status JSON")
    require(int(job.get("threads", 1)) >= 1, "Invalid threads")


def run_job(job_path: str | Path, out_root: str | Path, *, timeout: int | None = None,
            allow_synthetic=False) -> dict:
    job = read_json(job_path)
    validate_job(job)
    synthetic = bool(job.get("synthetic", False))
    require(not synthetic or allow_synthetic, "Synthetic job requires explicit --synthetic")
    for record in job["inputs"].values():
        check_hash(record["path"], record["sha256"])
    for dep in job.get("dependency_receipts", []):
        receipt = verify_artifact(dep["path"], dep.get("fingerprint"))
        require(receipt["synthetic"] == synthetic, "Synthetic dependency contamination")
        require(receipt.get("scientific_status", "PASS") in {"PASS", "NO_SIGNAL"}, "Rejected scientific dependency")
    binary = shutil.which(job["argv"][0])
    require(binary is not None, f"BLOCKED_BY_SOFTWARE: {job['argv'][0]} is not installed")
    if job.get("executable_sha256"):
        check_hash(binary, job["executable_sha256"])
    root = Path(out_root).resolve()
    root.mkdir(parents=True, exist_ok=True)
    ensure_free(root, int(job.get("minimum_free_bytes", 2**30)))
    paths = [job_path, binary] + [r["path"] for r in job["inputs"].values()]
    key, _ = fingerprint(paths, job, "native_job")
    target = safe_write_path(root, job["job_id"])
    if target.exists():
        return verify_artifact(target, key)
    with exclusive_lock(root / ".native_worker.lock"):
        with transaction(root, job["job_id"], stage="native_job", inputs=paths,
                         parameters=job, synthetic=synthetic) as (work, meta):
            replacements = {"{work}": str(work)}
            for name, record in job["inputs"].items():
                replacements[f"{{input:{name}}}"] = str(Path(record["path"]).resolve())
            for name, relative in job["outputs"].items():
                output = safe_write_path(work, relative)
                output.parent.mkdir(parents=True, exist_ok=True)
                replacements[f"{{output:{name}}}"] = str(output)
            argv = []
            for argument in job["argv"]:
                for token, value in replacements.items():
                    argument = argument.replace(token, value)
                require(re.search(r"\{(?:input|output):", argument) is None, "Unresolved command placeholder")
                argv.append(argument)
            env = os.environ.copy()
            # Restrict numeric BLAS parallelism to avoid oversubscribing an 8 GB host.
            for name in ["OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "VECLIB_MAXIMUM_THREADS"]:
                env[name] = str(job.get("threads", 1))
            env["BRAIN6_SYNTHETIC"] = "1" if synthetic else "0"
            env["BRAIN6_JOB_DIR"] = str(work)
            with (work / "native.log").open("w", encoding="utf-8") as log:
                result = subprocess.run(argv, cwd=work, env=env, stdout=log, stderr=subprocess.STDOUT,
                                        timeout=timeout or job.get("timeout_seconds", 86400), check=False)
            if result.returncode:
                tail = (work / "native.log").read_text(errors="replace")[-5000:]
                raise ContractError(f"Native job failed ({result.returncode}); no outputs promoted:\n{tail}")
            for name, relative in job["outputs"].items():
                require((work / relative).is_file(), f"Native output missing: {name}")
            status = read_json(work / job["outputs"]["semantic_status"])
            require(status.get("status") in {"PASS", "NO_SIGNAL", "INSUFFICIENT_EVIDENCE"},
                    f"Native scientific QC rejected: {status}")
            # INSUFFICIENT_EVIDENCE is recorded, not relabelled as a successful signal.
            meta["scientific_status"] = status["status"]
            meta["native_status"] = status
            write_json(work / "invocation.json", {"argv": argv, "executable": binary})
    return verify_artifact(target)
