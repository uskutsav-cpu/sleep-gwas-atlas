"""Transactional artifact directories and content-verified resumability."""
from __future__ import annotations
import contextlib
import os
import shutil
import tempfile
import uuid
from pathlib import Path
from typing import Iterator
from .io import (ContractError, exclusive_lock, file_record, json_hash, read_json,
                 require, safe_write_path, sha256, utc_now, write_json)
from . import __version__


def verify_artifact(path: str | Path, fingerprint: str | None = None) -> dict:
    path = Path(path)
    receipt = read_json(path / "receipt.json")
    require(receipt.get("status") == "COMPLETE", f"Unfinished artifact: {path}")
    if fingerprint is not None:
        require(receipt.get("fingerprint") == fingerprint, f"Stale input/parameter fingerprint: {path}")
    require(isinstance(receipt.get("outputs"), list) and bool(receipt["outputs"]),
            f"Missing artifact output inventory: {path}")
    names = [r["path"] for r in receipt["outputs"]]
    require(len(names) == len(set(names)), "Duplicate artifact output records")
    for record in receipt["outputs"]:
        output = safe_write_path(path, record["path"])
        require(output.is_file() and output.stat().st_size == record["bytes"] and sha256(output) == record["sha256"],
                f"Corrupted output: {output}")
    return receipt


def fingerprint(inputs: list[str | Path], parameters: dict, stage: str) -> tuple[str, list[dict]]:
    records = [file_record(p) for p in inputs]
    # Paths not hashed: moving identical inputs is harmless; stage/code/bytes matter.
    content = [{"sha256": r["sha256"], "bytes": r["bytes"]} for r in records]
    code_root = Path(__file__).resolve().parent
    code = [(str(p.relative_to(code_root)), sha256(p))
            for p in sorted(code_root.rglob("*"))
            if p.is_file() and p.suffix in {".py", ".R", ".m"}]
    # Include checkout adapters too: editable installs use scripts/, wheels use resources/.
    scripts = code_root.parent / "scripts"
    if scripts.is_dir():
        code += [("scripts/" + p.name, sha256(p)) for p in sorted(scripts.iterdir())
                 if p.is_file() and p.suffix in {".py", ".R", ".sh"}]
    return json_hash({"inputs": content, "parameters": parameters, "stage": stage,
                      "code": code, "version": __version__}), records


@contextlib.contextmanager
def transaction(root: str | Path, relative: str, *, stage: str,
                inputs: list[str | Path], parameters: dict,
                synthetic: bool = False) -> Iterator[tuple[Path, dict]]:
    """Publish all outputs together or none. No clobbering existing stages.

    Explicitly re-run into a new run directory if code/parameters change.
    Call verify_artifact for read-only resumption of completed stages.
    """
    target = safe_write_path(root, relative)
    target.parent.mkdir(parents=True, exist_ok=True)
    key, records = fingerprint(inputs, parameters, stage)
    with exclusive_lock(target.with_name(target.name + ".lock")):
        require(not target.exists(), f"Artifact exists; verify/reuse it or choose a new run: {target}")
        work = Path(tempfile.mkdtemp(prefix=f".{target.name}.", suffix=".partial", dir=target.parent))
        meta = {"stage": stage, "fingerprint": key, "inputs": records,
                "parameters": parameters, "synthetic": synthetic, "created_at": utc_now(),
                "version": __version__}
        try:
            yield work, meta
            # Recheck bytes AND code after execution: do not seal a result whose
            # inputs were changed by another writer while this task ran.
            end_key, _ = fingerprint(inputs, parameters, stage)
            require(end_key == key, "Inputs or implementation changed during execution; artifact not published")
            outputs = []
            for p in sorted(work.rglob("*")):
                require(not p.is_symlink(), f"Output symlink forbidden: {p}")
                if p.is_file() and p.name != "receipt.json":
                    outputs.append({"path": str(p.relative_to(work)), "bytes": p.stat().st_size,
                                    "sha256": sha256(p)})
            require(bool(outputs), "Stage produced no outputs")
            meta.update(status="COMPLETE", completed_at=utc_now(), outputs=outputs)
            write_json(work / "receipt.json", meta)
            os.rename(work, target)
        except BaseException as exc:
            # Preserve failed computation for audit, but outside the consumable
            # artifact path. The caller must never consume these as null results.
            failure_root = safe_write_path(root, "_failed")
            failure_root.mkdir(parents=True, exist_ok=True)
            failed = failure_root / (target.name + "-" + uuid.uuid4().hex[:12])
            try:
                meta.update(status="FAILED", failed_at=utc_now(), error_type=type(exc).__name__,
                            error=str(exc), outputs=[])
                write_json(work / "receipt.json", meta, immutable=False)
                os.rename(work, failed)
            except Exception:
                # Never mask the original exception (including disk exhaustion).
                # Keep the .partial directory for manual inspection in this case.
                pass
            raise


def verify_current_artifact(path):
    """For resumed computations: verify inputs, parameters AND currently installed code."""
    previous=verify_artifact(path)
    key,_=fingerprint([r['path'] for r in previous['inputs']],previous['parameters'],previous['stage'])
    return verify_artifact(path,key)
