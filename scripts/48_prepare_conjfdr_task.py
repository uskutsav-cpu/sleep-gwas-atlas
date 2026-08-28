#!/usr/bin/env python3
"""Freeze one official pleioFDR conjunction-FDR task and its exact config."""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import re
import subprocess
from pathlib import Path


PAIR_ID = re.compile(r"[a-z0-9_]+__[a-z0-9_]+$")


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def atomic_text(path: Path, payload: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(payload, encoding="utf-8")
    temporary.replace(path)


def atomic_tsv(path: Path, row: dict[str, object]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".tmp")
    with temporary.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=list(row), delimiter="\t", lineterminator="\n"
        )
        writer.writeheader()
        writer.writerow(row)
    temporary.replace(path)


def relative(root: Path, path: Path) -> str:
    try:
        return str(path.resolve().relative_to(root))
    except ValueError as exc:
        raise SystemExit(f"ERROR: task path is outside repository root: {path}") from exc


def validate_trait_mat(root: Path, trait: str, policy_hash: str, template_hash: str) -> tuple[Path, str]:
    path = root / "data/pleiofdr" / f"{trait}.mat"
    provenance_path = root / "data/pleiofdr" / f"{trait}.provenance.json"
    if not path.is_file() or not provenance_path.is_file():
        raise SystemExit(f"ERROR: pleioFDR MAT input/provenance is absent for {trait}")
    provenance = json.loads(provenance_path.read_text(encoding="utf-8"))
    observed = sha256(path)
    checks = (
        provenance.get("trait_id") == trait,
        provenance.get("policy_sha256") == policy_hash,
        provenance.get("template_sha256") == template_hash,
        provenance.get("output_sha256") == observed,
        provenance.get("output") == relative(root, path),
        provenance.get("alignment_counts", {}).get("matched", 0) > 0,
    )
    if not all(checks):
        raise SystemExit(f"ERROR: pleioFDR MAT provenance failed for {trait}")
    return path, observed


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("pair_id")
    parser.add_argument("--root", default=".")
    parser.add_argument("--manifest", default="results/tables/pleiotropy_pair_manifest.tsv")
    parser.add_argument("--manifest-lock", default="results/tables/pleiotropy_pair_manifest.lock.json")
    parser.add_argument("--task-dir", default="results/pleiotropy/conjfdr_tasks")
    parser.add_argument("--result-dir", default="results/pleiotropy/conjfdr")
    args = parser.parse_args()
    if not PAIR_ID.fullmatch(args.pair_id):
        raise SystemExit("ERROR: invalid pair_id")
    root = Path(args.root).resolve()
    policy_path = root / "config/pleiotropy_analysis_policy.json"
    policy = json.loads(policy_path.read_text(encoding="utf-8"))
    policy_hash = sha256(policy_path)
    manifest_path = root / args.manifest
    manifest_lock_path = root / args.manifest_lock
    manifest_lock = json.loads(manifest_lock_path.read_text(encoding="utf-8"))
    if manifest_lock.get("manifest_sha256") != sha256(manifest_path):
        raise SystemExit("ERROR: pair manifest differs from its lock")
    if manifest_lock.get("policy_sha256") != policy_hash:
        raise SystemExit("ERROR: policy differs from the pair-manifest lock")
    rows = [row for row in read_tsv(manifest_path) if row["pair_id"] == args.pair_id]
    if len(rows) != 1:
        raise SystemExit("ERROR: pair_id must identify exactly one locked pair")
    pair = rows[0]
    if pair["input_status"] != "READY_FULL_SUMSTATS":
        raise SystemExit(f"ERROR: pair is blocked: {pair['blocker']}")

    template = root / "ref/pleiofdr/9545380.ref"
    reference = root / "ref/pleiofdr/ref9545380_1kgPhase3eur_LDr2p1.mat"
    code = root / "work/pleiofdr"
    patch = root / policy["pleiofdr_overlap_patch"]
    if (
        not template.is_file()
        or template.stat().st_size != policy["pleiofdr_variant_template_bytes"]
        or sha256(template) != policy["pleiofdr_variant_template_sha256"]
    ):
        raise SystemExit("ERROR: pinned pleioFDR variant template is absent or invalid")
    if not reference.is_file() or reference.stat().st_size != policy["pleiofdr_reference_bytes"]:
        raise SystemExit("ERROR: official pleioFDR LD reference is absent or byte-invalid")
    if not (code / ".git").is_dir():
        raise SystemExit("ERROR: pinned pleioFDR checkout is absent")
    commit = subprocess.run(
        ["git", "-C", str(code), "rev-parse", "HEAD"], capture_output=True, text=True, check=True
    ).stdout.strip()
    dirty = subprocess.run(
        ["git", "-C", str(code), "status", "--porcelain"], capture_output=True, text=True, check=True
    ).stdout.strip()
    if commit != policy["pleiofdr_commit"] or dirty:
        raise SystemExit("ERROR: pleioFDR checkout is not the clean pinned commit")
    if not patch.is_file() or sha256(patch) != policy["pleiofdr_overlap_patch_sha256"]:
        raise SystemExit("ERROR: sample-overlap patch is absent or checksum-invalid")
    sleep_mat, sleep_hash = validate_trait_mat(
        root, pair["sleep_trait"], policy_hash, policy["pleiofdr_variant_template_sha256"]
    )
    disease_mat, disease_hash = validate_trait_mat(
        root, pair["non_sleep_trait"], policy_hash, policy["pleiofdr_variant_template_sha256"]
    )

    task_dir = root / args.task_dir
    result_dir = root / args.result_dir / args.pair_id
    config_path = task_dir / f"{args.pair_id}.config.txt"
    for value in (root, template, reference, code, patch, sleep_mat, disease_mat, result_dir):
        if "'" in str(value) or "\n" in str(value):
            raise SystemExit("ERROR: paths containing quotes/newlines are unsupported by MATLAB task config")
    config = "\n".join([
        f"reffile={reference}",
        f"refinfo={template}",
        f"traitfolder={sleep_mat.parent}",
        f"traitfile1={sleep_mat.name}",
        f"traitname1={pair['sleep_trait']}",
        f"traitfiles={{'{disease_mat.name}'}}",
        f"traitnames={{'{pair['non_sleep_trait']}'}}",
        f"mlibrary={code}",
        f"outputdir={result_dir}",
        "stattype=conjfdr",
        f"fdrthresh={policy['pleiofdr_conjfdr_threshold']}",
        "randprune=true",
        f"randprune_n={policy['pleiofdr_random_prune_iterations']}",
        "reset_pruneidx=true",
        "randprune_repeats=default",
        "exclude_chr_pos=[6 25119106 33854733]",
        "exclude_from_discovery=false",
        "mafthresh=0.005",
        "exclude_ambiguous_snps=true",
        "perform_gc=true",
        "use_standard_gc=false",
        "randprune_gc=true",
        "onscreen=false",
        "manh_plot=false",
        "pthresh=1",
        "dummy_zscore=false",
        "exit_matlab_upon_completion=false",
        "",
    ])
    atomic_text(config_path, config)
    prefix = f"{pair['sleep_trait']}_{pair['non_sleep_trait']}_conjfdr_{policy['pleiofdr_conjfdr_threshold']}"
    task_path = task_dir / f"{args.pair_id}.tsv"
    lock_path = task_dir / f"{args.pair_id}.lock.tsv"
    task = {
        "analysis_id": policy["analysis_id"], "pair_id": args.pair_id,
        "sleep_trait": pair["sleep_trait"], "non_sleep_trait": pair["non_sleep_trait"],
        "analysis_tier": pair["analysis_tier"],
        "policy_sha256": policy_hash, "manifest_sha256": sha256(manifest_path),
        "config": relative(root, config_path), "config_sha256": sha256(config_path),
        "sleep_mat": relative(root, sleep_mat), "sleep_mat_sha256": sleep_hash,
        "non_sleep_mat": relative(root, disease_mat), "non_sleep_mat_sha256": disease_hash,
        "template": relative(root, template), "template_sha256": sha256(template),
        "reference": relative(root, reference), "reference_sha256": sha256(reference),
        "pleiofdr_code": relative(root, code), "pleiofdr_commit": commit,
        "overlap_patch": relative(root, patch), "overlap_patch_sha256": sha256(patch),
        "result_dir": relative(root, result_dir),
        "all_results": relative(root, result_dir / f"{prefix}_all.csv"),
        "locus_results": relative(root, result_dir / f"{prefix}_loci.csv"),
        "result_mat": relative(root, result_dir / "result.mat"),
        "completion": relative(root, result_dir / "atlas_completion.tsv"),
        "conjfdr_threshold": policy["pleiofdr_conjfdr_threshold"],
        "random_prune_iterations": policy["pleiofdr_random_prune_iterations"],
        "correct_sample_overlap": str(policy["pleiofdr_correct_sample_overlap"]).upper(),
    }
    atomic_tsv(task_path, task)
    atomic_tsv(lock_path, {
        "analysis_id": policy["analysis_id"], "pair_id": args.pair_id,
        "task": relative(root, task_path), "task_sha256": sha256(task_path),
    })
    print(f"ConjFDR task locked: {args.pair_id}")
    print(f"Run explicitly: python3 scripts/49_run_conjfdr_pair.py {relative(root, task_path)} {relative(root, lock_path)} --execute")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
