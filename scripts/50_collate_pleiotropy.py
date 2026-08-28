#!/usr/bin/env python3
"""Validate all 396 pair scans and collate cross-method shared LD blocks."""
from __future__ import annotations

import argparse
import bisect
import csv
import hashlib
import io
import json
import math
import re
from collections import defaultdict
from pathlib import Path

import pleiotropy_contract


PLACO_FIELDS = [
    "pair_id", "sleep_trait", "non_sleep_trait", "analysis_tier", "locus_id",
    "CHR", "START", "STOP", "lead_snp", "lead_p_placo_plus", "lead_z1", "lead_z2",
    "effect_direction",
    "family_significant_variant_count", "family_significant_variants", "method",
]
CONJFDR_FIELDS = [
    "pair_id", "sleep_trait", "non_sleep_trait", "analysis_tier", "locus_id",
    "CHR", "START", "STOP", "lead_snp", "lead_conjfdr",
    "significant_variant_count", "significant_variants", "method",
]
SHARED_FIELDS = [
    "shared_locus_id", "pair_id", "sleep_trait", "non_sleep_trait", "analysis_tier",
    "locus_id", "CHR", "START", "STOP", "placo_lead_snp", "placo_lead_p",
    "conjfdr_lead_snp", "conjfdr_lead_fdr", "same_lead_snp",
    "effect_direction",
    "placo_significant_variant_count", "conjfdr_significant_variant_count",
    "evidence_status", "claim_limit",
]


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def is_sha256(value: object) -> bool:
    return isinstance(value, str) and re.fullmatch(r"[0-9a-f]{64}", value) is not None


def validate_placo_input_identity(
    pair_id: str,
    pair_relative: str,
    pair_manifest: dict[str, str],
    pair_provenance: dict[str, object],
    placo_task: dict[str, str],
    placo_summary: dict[str, str],
    manifest_sha256: str,
    policy_sha256: str,
    live_pair_sha256: str | None,
) -> tuple[str, int]:
    """Validate the immutable PLACO input chain, with or without its temp payload."""
    bound_sha256 = placo_task.get("pair_input_sha256")
    rows_written = pair_provenance.get("alignment_counts", {}).get("written")
    expected_sources = [
        pair_manifest["sleep_sumstats"], pair_manifest["non_sleep_sumstats"],
    ]
    expected_source_hashes = [
        pair_manifest["sleep_sumstats_sha256"],
        pair_manifest["non_sleep_sumstats_sha256"],
    ]
    if (
        not is_sha256(bound_sha256)
        or (live_pair_sha256 is not None and live_pair_sha256 != bound_sha256)
        or placo_task.get("pair_input") != pair_relative
        or pair_provenance.get("pair_id") != pair_id
        or pair_provenance.get("output") != pair_relative
        or pair_provenance.get("output_sha256") != bound_sha256
        or pair_provenance.get("source_files") != expected_sources
        or pair_provenance.get("source_sha256") != expected_source_hashes
        or pair_provenance.get("manifest_sha256") != manifest_sha256
        or pair_provenance.get("policy_sha256") != policy_sha256
        or not isinstance(rows_written, int) or rows_written <= 0
        or placo_task.get("pair_input_rows") != str(rows_written)
        or placo_summary.get("input_sha256") != bound_sha256
    ):
        raise SystemExit(f"ERROR: PLACO+ input identity chain drifted for {pair_id}")
    return bound_sha256, rows_written


def validate_conjfdr_result_identity(
    root: Path,
    pair_id: str,
    task_path: Path,
    task_lock_path: Path,
    completion_path: Path,
    result_mat_path: Path,
    receipt_path: Path,
    task: dict[str, str],
    task_lock: dict[str, str],
    completion: dict[str, str],
) -> tuple[str, int, bool]:
    """Validate a sealed MAT result whether or not its temporary payload remains."""
    try:
        receipt = json.loads(receipt_path.read_text(encoding="utf-8"))
    except (OSError, UnicodeError, json.JSONDecodeError) as exc:
        raise SystemExit(
            f"ERROR: conjunction-FDR MAT receipt is unreadable for {pair_id}: {exc}"
        ) from exc
    if not isinstance(receipt, dict):
        raise SystemExit(
            f"ERROR: conjunction-FDR MAT receipt must be an object for {pair_id}"
        )
    result_sha256 = receipt.get("result_mat_sha256")
    result_bytes = receipt.get("result_mat_bytes")
    if result_mat_path.exists() and not result_mat_path.is_file():
        raise SystemExit(f"ERROR: conjunction-FDR MAT path is not a file for {pair_id}")
    retained = result_mat_path.is_file()
    if (
        receipt.get("schema_version") != "sleep-atlas-conjfdr-result-mat.1"
        or receipt.get("analysis_id") != task.get("analysis_id")
        or receipt.get("pair_id") != pair_id
        or receipt.get("task") != str(task_path.relative_to(root))
        or receipt.get("task_sha256") != digest(task_path)
        or receipt.get("task_lock") != str(task_lock_path.relative_to(root))
        or receipt.get("task_lock_sha256") != digest(task_lock_path)
        or receipt.get("completion") != str(completion_path.relative_to(root))
        or receipt.get("completion_sha256") != digest(completion_path)
        or receipt.get("result_mat") != str(result_mat_path.relative_to(root))
        or not is_sha256(result_sha256)
        or completion.get("result_mat_sha256") != result_sha256
        or type(result_bytes) is not int or result_bytes <= 0
        or receipt.get("sealer_sha256")
        != digest(root / "scripts/49_seal_conjfdr_result.py")
        or receipt.get("retention")
        != "TEMPORARY_RESULT_MAT_MAY_BE_EVICTED_AFTER_THIS_RECEIPT"
        or task_lock.get("task_sha256") != receipt.get("task_sha256")
    ):
        raise SystemExit(f"ERROR: conjunction-FDR MAT receipt chain drifted for {pair_id}")
    if retained and (
        result_mat_path.stat().st_size != result_bytes
        or digest(result_mat_path) != result_sha256
    ):
        raise SystemExit(f"ERROR: retained conjunction-FDR MAT drifted for {pair_id}")
    return result_sha256, result_bytes, retained


def rows(path: Path, delimiter: str = "\t") -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter=delimiter))


def rows_and_fields(path: Path, delimiter: str = "\t") -> tuple[list[dict[str, str]], set[str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter=delimiter)
        values = list(reader)
        return values, set(reader.fieldnames or [])


def one_row(path: Path, delimiter: str = "\t") -> dict[str, str]:
    values = rows(path, delimiter)
    if len(values) != 1:
        raise SystemExit(f"ERROR: expected one row in {path}")
    return values[0]


def probability(value: str, field: str, identity: str) -> float:
    try:
        number = float(value)
    except ValueError as exc:
        raise SystemExit(f"ERROR: invalid {field} for {identity}") from exc
    if not math.isfinite(number) or not 0 <= number <= 1:
        raise SystemExit(f"ERROR: invalid {field} for {identity}")
    return number


def table_text(fields: list[str], values: list[dict[str, object]]) -> str:
    output = io.StringIO(newline="")
    writer = csv.DictWriter(output, fieldnames=fields, delimiter="\t", lineterminator="\n")
    writer.writeheader()
    writer.writerows(values)
    return output.getvalue()


def effect_direction(z1: float, z2: float) -> str:
    if z1 == 0 or z2 == 0:
        raise SystemExit("ERROR: a PLACO+ locus lead has a zero signed effect")
    return ("+" if z1 > 0 else "-") + "/" + ("+" if z2 > 0 else "-")


def atomic_text(path: Path, payload: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(payload, encoding="utf-8")
    temporary.replace(path)


def load_blocks(path: Path, expected_hash: str, expected_count: int):
    if sha256(path) != expected_hash:
        raise SystemExit("ERROR: pleiotropy locus definition differs from its pin")
    blocks = rows(path, delimiter=" ")
    if len(blocks) != expected_count or set(blocks[0] if blocks else []) != {"LOC", "CHR", "START", "STOP"}:
        raise SystemExit("ERROR: pleiotropy locus definition has unexpected shape")
    by_chromosome: dict[int, list[tuple[int, int, str]]] = defaultdict(list)
    seen = set()
    for row in blocks:
        locus_id = row["LOC"]
        chromosome, start, stop = int(row["CHR"]), int(row["START"]), int(row["STOP"])
        if locus_id in seen or not (1 <= chromosome <= 22 and 0 < start <= stop):
            raise SystemExit("ERROR: invalid or duplicate locked locus")
        seen.add(locus_id)
        current = by_chromosome[chromosome]
        if current and start <= current[-1][1]:
            raise SystemExit("ERROR: locked locus definitions overlap or are unsorted")
        current.append((start, stop, locus_id))
    starts = {chromosome: [item[0] for item in values] for chromosome, values in by_chromosome.items()}
    return by_chromosome, starts


def assign_block(chromosome: int, position: int, blocks, starts) -> tuple[str, int, int]:
    candidates = blocks.get(chromosome, [])
    index = bisect.bisect_right(starts.get(chromosome, []), position) - 1
    if index < 0 or position > candidates[index][1]:
        raise SystemExit(f"ERROR: variant chr{chromosome}:{position} is outside locked loci")
    start, stop, locus_id = candidates[index]
    return locus_id, start, stop


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", default=".")
    parser.add_argument("--manifest", default="results/tables/pleiotropy_pair_manifest.tsv")
    parser.add_argument("--manifest-lock", default="results/tables/pleiotropy_pair_manifest.lock.json")
    parser.add_argument("--pair-dir", default="results/pleiotropy/inputs")
    parser.add_argument("--placo-task-dir", default="results/pleiotropy/tasks")
    parser.add_argument("--placo-dir", default="results/pleiotropy/placo")
    parser.add_argument("--conjfdr-task-dir", default="results/pleiotropy/conjfdr_tasks")
    parser.add_argument("--placo-out", default="results/tables/placo_loci.tsv")
    parser.add_argument("--conjfdr-out", default="results/tables/conjfdr_loci.tsv")
    parser.add_argument("--shared-out", default="results/atlas/shared_loci.tsv")
    parser.add_argument("--provenance-out", default="results/atlas/shared_loci.provenance.json")
    parser.add_argument("--report-only", action="store_true")
    parser.add_argument("--validate-only", action="store_true")
    parser.add_argument("--quiet", action="store_true")
    args = parser.parse_args()
    if args.report_only and args.validate_only:
        raise SystemExit("ERROR: choose report-only or validate-only")
    root = Path(args.root).resolve()
    policy_path = root / "config/pleiotropy_analysis_policy.json"
    policy = json.loads(policy_path.read_text(encoding="utf-8"))
    manifest_path = root / args.manifest
    manifest_lock_path = root / args.manifest_lock
    lock = json.loads(manifest_lock_path.read_text(encoding="utf-8"))
    if (
        lock.get("schema_version") != "sleep-atlas-pleiotropy-pairs.1"
        or lock.get("ready_pair_count") != policy["expected_sleep_non_sleep_pairs"]
        or lock.get("blocked_pair_count") != 0
        or lock.get("script_sha256") != sha256(root / "scripts/43_prepare_pleiotropy_pairs.py")
    ):
        raise SystemExit("ERROR: complete immutable pleiotropy pair family is required")
    if lock.get("manifest_sha256") != sha256(manifest_path):
        raise SystemExit("ERROR: pleiotropy manifest differs from its lock")
    if lock.get("policy_sha256") != sha256(policy_path):
        raise SystemExit("ERROR: pleiotropy policy differs from its manifest lock")
    manifest = rows(manifest_path)
    pair_ids = [row["pair_id"] for row in manifest]
    if (
        len(manifest) != policy["expected_sleep_non_sleep_pairs"]
        or pair_ids != lock.get("pair_ids_in_locked_order")
        or len(set(pair_ids)) != len(pair_ids)
    ):
        raise SystemExit("ERROR: pleiotropy pair family/order differs from lock")
    expected_paths: list[Path] = []
    for pair_id in pair_ids:
        expected_paths.extend([
            root / args.pair_dir / f"{pair_id}.provenance.json",
            root / args.placo_task_dir / f"{pair_id}.tsv",
            root / args.placo_task_dir / f"{pair_id}.lock.tsv",
            root / args.placo_dir / f"{pair_id}.summary.tsv",
            root / args.placo_dir / f"{pair_id}.hits.tsv",
            root / args.conjfdr_task_dir / f"{pair_id}.config.txt",
            root / args.conjfdr_task_dir / f"{pair_id}.tsv",
            root / args.conjfdr_task_dir / f"{pair_id}.lock.tsv",
        ])
    missing = [str(path.relative_to(root)) for path in expected_paths if not path.is_file()]
    if missing:
        message = f"Pleiotropy collation blocked: {len(missing)} required pair artifacts are absent"
        if not args.quiet:
            print(message)
            print("First missing artifacts: " + ", ".join(missing[:12]))
        if args.report_only:
            return 0
        raise SystemExit("ERROR: " + message)

    runtime_provenance = pleiotropy_contract.validate_runtime(
        root, rehash_reference=True,
    )

    locus_path = root / policy["locus_definition"]
    blocks, starts = load_blocks(
        locus_path, policy["locus_definition_sha256"], policy["expected_loci"]
    )
    placo_grouped: dict[tuple[str, str], list[dict[str, object]]] = defaultdict(list)
    conj_grouped: dict[tuple[str, str], list[dict[str, object]]] = defaultdict(list)
    input_hashes: dict[str, dict[str, object]] = {}
    digest_cache: dict[Path, str] = {}

    def digest(path: Path) -> str:
        if path not in digest_cache:
            digest_cache[path] = sha256(path)
        return digest_cache[path]

    for pair in manifest:
        pair_id = pair["pair_id"]
        pair_path = root / args.pair_dir / f"{pair_id}.tsv.gz"
        pair_provenance_path = root / args.pair_dir / f"{pair_id}.provenance.json"
        placo_task_path = root / args.placo_task_dir / f"{pair_id}.tsv"
        placo_lock_path = root / args.placo_task_dir / f"{pair_id}.lock.tsv"
        summary_path = root / args.placo_dir / f"{pair_id}.summary.tsv"
        hits_path = root / args.placo_dir / f"{pair_id}.hits.tsv"
        pair_provenance = json.loads(pair_provenance_path.read_text(encoding="utf-8"))
        placo_task, placo_lock = one_row(placo_task_path), one_row(placo_lock_path)
        summary = one_row(summary_path)
        if pair_path.exists() and not pair_path.is_file():
            raise SystemExit(f"ERROR: retained PLACO+ pair input is not a file for {pair_id}")
        pair_input_sha256, _pair_input_rows = validate_placo_input_identity(
            pair_id=pair_id,
            pair_relative=str(pair_path.relative_to(root)),
            pair_manifest=pair,
            pair_provenance=pair_provenance,
            placo_task=placo_task,
            placo_summary=summary,
            manifest_sha256=digest(manifest_path),
            policy_sha256=digest(policy_path),
            live_pair_sha256=digest(pair_path) if pair_path.is_file() else None,
        )
        if (
            placo_lock.get("schema_version") != "sleep-atlas-placo-task.1"
            or placo_lock.get("task_sha256") != digest(placo_task_path)
            or placo_lock.get("pair_id") != pair_id
            or placo_task.get("pair_id") != pair_id
            or placo_task.get("task_builder_sha256")
            != digest(root / "scripts/45_prepare_placo_task.py")
            or placo_task.get("runner_sha256") != digest(root / "scripts/46_run_placo_pair.R")
            or placo_task.get("pair_provenance_sha256") != digest(pair_provenance_path)
            or placo_task.get("variant_hits_out") != str(hits_path.relative_to(root))
            or placo_task.get("summary_out") != str(summary_path.relative_to(root))
        ):
            raise SystemExit(f"ERROR: PLACO+ input/task provenance drifted for {pair_id}")
        if (
            summary.get("pair_id") != pair_id
            or summary.get("analysis_status") != "PLACO_PLUS_COMPLETE"
            or summary.get("analysis_tier") != pair["analysis_tier"]
            or float(summary.get("family_threshold", "nan")) != policy["placo_locked_pair_family_threshold"]
            or float(summary.get("numerical_failure_fraction", "nan"))
            > policy["placo_maximum_numerical_failure_fraction"]
            or summary.get("task_sha256") != digest(placo_task_path)
            or summary.get("input_sha256") != pair_input_sha256
            or summary.get("placo_source_sha256") != policy["placo_source_sha256"]
            or summary.get("task_builder_sha256") != placo_task["task_builder_sha256"]
            or summary.get("runner_sha256") != placo_task["runner_sha256"]
        ):
            raise SystemExit(f"ERROR: invalid PLACO+ summary for {pair_id}")
        placo_hits, placo_fields = rows_and_fields(hits_path)
        required_placo = {
            "pair_id", "SNP", "CHR", "BP", "P_PLACO_PLUS",
            "Z1", "Z2", "conventional_significant", "locked_family_significant",
        }
        if not required_placo.issubset(placo_fields):
            raise SystemExit(f"ERROR: PLACO+ hit schema is invalid for {pair_id}")
        if int(summary["conventional_hit_count"]) != len(placo_hits):
            raise SystemExit(f"ERROR: PLACO+ hit count differs from summary for {pair_id}")
        seen_hits = set()
        family_count = 0
        for hit in placo_hits:
            identity = (hit.get("SNP"), hit.get("CHR"), hit.get("BP"))
            if (
                hit.get("pair_id") != pair_id or identity in seen_hits
                or hit.get("conventional_significant") != "TRUE"
            ):
                raise SystemExit(f"ERROR: invalid/duplicate PLACO+ hit for {pair_id}")
            seen_hits.add(identity)
            p_value = probability(hit["P_PLACO_PLUS"], "P_PLACO_PLUS", f"{pair_id}/{hit['SNP']}")
            if p_value > policy["placo_conventional_variant_threshold"]:
                raise SystemExit(f"ERROR: non-significant variant retained in PLACO+ hits for {pair_id}")
            family = p_value <= policy["placo_locked_pair_family_threshold"]
            if (hit.get("locked_family_significant") == "TRUE") != family:
                raise SystemExit(f"ERROR: PLACO+ family indicator differs from policy for {pair_id}")
            if family:
                family_count += 1
                chromosome, position = int(hit["CHR"]), int(hit["BP"])
                locus_id, start, stop = assign_block(chromosome, position, blocks, starts)
                placo_grouped[(pair_id, locus_id)].append({
                    "snp": hit["SNP"], "p": p_value, "chr": chromosome,
                    "start": start, "stop": stop, "z1": float(hit["Z1"]),
                    "z2": float(hit["Z2"]),
                })
        if family_count != int(summary["locked_family_hit_count"]):
            raise SystemExit(f"ERROR: PLACO+ family hit count differs from summary for {pair_id}")

        task_path = root / args.conjfdr_task_dir / f"{pair_id}.tsv"
        task_lock_path = root / args.conjfdr_task_dir / f"{pair_id}.lock.tsv"
        task, task_lock = one_row(task_path), one_row(task_lock_path)
        if (
            task_lock.get("schema_version") != "sleep-atlas-conjfdr-task.1"
            or task_lock.get("task_sha256") != digest(task_path)
            or task_lock.get("pair_id") != pair_id
            or task.get("pair_id") != pair_id
            or task.get("task_builder_sha256")
            != digest(root / "scripts/48_prepare_conjfdr_task.py")
            or task.get("runner_sha256") != digest(root / "scripts/49_run_conjfdr_pair.py")
            or task.get("runtime_provenance") != policy["runtime_provenance"]
            or task.get("runtime_provenance_sha256")
            != digest(root / str(policy["runtime_provenance"]))
            or task.get("reference_sha256")
            != runtime_provenance["pleiofdr_reference"]["sha256"]
            or int(task.get("reference_bytes", "0"))
            != runtime_provenance["pleiofdr_reference"]["bytes"]
            or task.get("pleiofdr_commit") != runtime_provenance["pleiofdr_code"]["commit"]
        ):
            raise SystemExit(f"ERROR: conjunction-FDR task/lock drifted for {pair_id}")
        for field, hash_field in (
            ("config", "config_sha256"), ("sleep_mat", "sleep_mat_sha256"),
            ("non_sleep_mat", "non_sleep_mat_sha256"), ("template", "template_sha256"),
            ("overlap_patch", "overlap_patch_sha256"),
        ):
            if digest(root / task[field]) != task[hash_field]:
                raise SystemExit(f"ERROR: conjunction-FDR input drifted for {pair_id}: {field}")
        completion_path = root / task["completion"]
        all_results_path = root / task["all_results"]
        locus_results_path = root / task["locus_results"]
        result_mat_path = root / task["result_mat"]
        result_mat_receipt_path = result_mat_path.with_name("result_mat.receipt.json")
        matlab_log_path = root / task["result_dir"] / "matlab.log"
        if not all(path.is_file() and path.stat().st_size for path in (
            completion_path, all_results_path, locus_results_path,
            result_mat_receipt_path, matlab_log_path,
        )):
            raise SystemExit(f"ERROR: conjunction-FDR result is absent for {pair_id}")
        completion = one_row(completion_path)
        if (
            completion.get("analysis_status") != "CONJFDR_COMPLETE"
            or completion.get("pair_id") != pair_id
            or completion.get("task_sha256") != digest(task_path)
            or completion.get("all_results_sha256") != digest(all_results_path)
            or completion.get("locus_results_sha256") != digest(locus_results_path)
            or completion.get("matlab_log_sha256") != digest(matlab_log_path)
            or completion.get("runtime_provenance_sha256")
            != digest(root / str(policy["runtime_provenance"]))
            or completion.get("task_builder_sha256") != task["task_builder_sha256"]
            or completion.get("runner_sha256") != task["runner_sha256"]
            or completion.get("correct_sample_overlap") != "TRUE"
            or int(completion.get("random_prune_iterations", "0"))
            != policy["pleiofdr_random_prune_iterations"]
        ):
            raise SystemExit(f"ERROR: conjunction-FDR completion drifted for {pair_id}")
        result_mat_sha256, result_mat_bytes, result_mat_retained = (
            validate_conjfdr_result_identity(
                root, pair_id, task_path, task_lock_path, completion_path,
                result_mat_path, result_mat_receipt_path, task, task_lock, completion,
            )
        )
        conj_rows, conj_fields = rows_and_fields(all_results_path, delimiter=",")
        required_conj = {"snpid", "chrnum", "chrpos", "min_conjfdr"}
        if not required_conj.issubset(conj_fields):
            raise SystemExit(f"ERROR: conjunction-FDR result schema is invalid for {pair_id}")
        seen_conj = set()
        for hit in conj_rows:
            fdr = probability(hit["min_conjfdr"], "min_conjfdr", f"{pair_id}/{hit['snpid']}")
            if fdr > policy["pleiofdr_conjfdr_threshold"]:
                continue
            identity = (hit["snpid"], hit["chrnum"], hit["chrpos"])
            if not hit["snpid"].startswith("rs") or identity in seen_conj:
                raise SystemExit(f"ERROR: invalid/duplicate conjunction-FDR hit for {pair_id}")
            seen_conj.add(identity)
            chromosome, position = int(hit["chrnum"]), int(hit["chrpos"])
            locus_id, start, stop = assign_block(chromosome, position, blocks, starts)
            conj_grouped[(pair_id, locus_id)].append({
                "snp": hit["snpid"], "fdr": fdr, "chr": chromosome,
                "start": start, "stop": stop,
            })
        input_hashes[pair_id] = {
            "pair_input": pair_input_sha256,
            "pair_input_retained_at_collation": str(pair_path.is_file()).upper(),
            "pair_provenance": digest(pair_provenance_path),
            "placo_task": digest(placo_task_path), "placo_lock": digest(placo_lock_path),
            "placo_summary": digest(summary_path), "placo_hits": digest(hits_path),
            "conjfdr_task": digest(task_path), "conjfdr_lock": digest(task_lock_path),
            "conjfdr_completion": digest(completion_path),
            "conjfdr_all": digest(all_results_path),
            "conjfdr_loci": digest(locus_results_path),
            "conjfdr_result_mat": result_mat_sha256,
            "conjfdr_result_mat_bytes": result_mat_bytes,
            "conjfdr_result_mat_retained_at_collation": str(result_mat_retained).upper(),
            "conjfdr_result_mat_receipt": digest(result_mat_receipt_path),
            "conjfdr_matlab_log": digest(matlab_log_path),
        }

    manifest_by_pair = {row["pair_id"]: row for row in manifest}
    pair_order = {pair_id: index for index, pair_id in enumerate(pair_ids)}
    locus_order = {item[2]: index for values in blocks.values() for index, item in enumerate(values)}
    key_order = lambda key: (pair_order[key[0]], int(key[1]))
    placo_loci = []
    for key in sorted(placo_grouped, key=key_order):
        pair_id, locus_id = key
        values = placo_grouped[key]
        lead = min(values, key=lambda item: (item["p"], item["snp"]))
        pair = manifest_by_pair[pair_id]
        placo_loci.append({
            "pair_id": pair_id, "sleep_trait": pair["sleep_trait"],
            "non_sleep_trait": pair["non_sleep_trait"], "analysis_tier": pair["analysis_tier"],
            "locus_id": locus_id, "CHR": lead["chr"], "START": lead["start"], "STOP": lead["stop"],
            "lead_snp": lead["snp"], "lead_p_placo_plus": f"{lead['p']:.15g}",
            "lead_z1": f"{lead['z1']:.15g}", "lead_z2": f"{lead['z2']:.15g}",
            "effect_direction": effect_direction(lead["z1"], lead["z2"]),
            "family_significant_variant_count": len(values),
            "family_significant_variants": ";".join(sorted(item["snp"] for item in values)),
            "method": "PLACO_PLUS_FAMILY_CORRECTED",
        })
    conj_loci = []
    for key in sorted(conj_grouped, key=key_order):
        pair_id, locus_id = key
        values = conj_grouped[key]
        lead = min(values, key=lambda item: (item["fdr"], item["snp"]))
        pair = manifest_by_pair[pair_id]
        conj_loci.append({
            "pair_id": pair_id, "sleep_trait": pair["sleep_trait"],
            "non_sleep_trait": pair["non_sleep_trait"], "analysis_tier": pair["analysis_tier"],
            "locus_id": locus_id, "CHR": lead["chr"], "START": lead["start"], "STOP": lead["stop"],
            "lead_snp": lead["snp"], "lead_conjfdr": f"{lead['fdr']:.15g}",
            "significant_variant_count": len(values),
            "significant_variants": ";".join(sorted(item["snp"] for item in values)),
            "method": "CONJFDR_0.05",
        })
    placo_index = {(row["pair_id"], row["locus_id"]): row for row in placo_loci}
    conj_index = {(row["pair_id"], row["locus_id"]): row for row in conj_loci}
    shared = []
    for pair_id, locus_id in sorted(set(placo_index) & set(conj_index), key=key_order):
        placo, conj = placo_index[(pair_id, locus_id)], conj_index[(pair_id, locus_id)]
        pair = manifest_by_pair[pair_id]
        shared.append({
            "shared_locus_id": f"{pair_id}__loc{locus_id}", "pair_id": pair_id,
            "sleep_trait": pair["sleep_trait"], "non_sleep_trait": pair["non_sleep_trait"],
            "analysis_tier": pair["analysis_tier"], "locus_id": locus_id,
            "CHR": placo["CHR"], "START": placo["START"], "STOP": placo["STOP"],
            "placo_lead_snp": placo["lead_snp"], "placo_lead_p": placo["lead_p_placo_plus"],
            "conjfdr_lead_snp": conj["lead_snp"], "conjfdr_lead_fdr": conj["lead_conjfdr"],
            "same_lead_snp": str(placo["lead_snp"] == conj["lead_snp"]).upper(),
            "effect_direction": placo["effect_direction"],
            "placo_significant_variant_count": placo["family_significant_variant_count"],
            "conjfdr_significant_variant_count": conj["significant_variant_count"],
            "evidence_status": "PLACO_PLUS_AND_CONJFDR_SAME_LOCKED_LD_BLOCK",
            "claim_limit": "STATISTICAL_PLEIOTROPY_NOT_SHARED_CAUSAL_VARIANT_OR_CAUSAL_DIRECTION",
        })

    payloads = {
        root / args.placo_out: table_text(PLACO_FIELDS, placo_loci),
        root / args.conjfdr_out: table_text(CONJFDR_FIELDS, conj_loci),
        root / args.shared_out: table_text(SHARED_FIELDS, shared),
    }
    provenance = {
        "schema_version": "sleep-atlas-pleiotropic-loci.1",
        "analysis_id": policy["analysis_id"], "policy_sha256": sha256(policy_path),
        "manifest_sha256": sha256(manifest_path), "manifest_lock_sha256": sha256(manifest_lock_path),
        "runtime_provenance_sha256": digest(root / str(policy["runtime_provenance"])),
        "locus_definition_sha256": sha256(locus_path), "completed_pair_scans": len(pair_ids),
        "placo_locus_count": len(placo_loci), "conjfdr_locus_count": len(conj_loci),
        "canonical_shared_locus_count": len(shared), "pair_input_hashes": input_hashes,
        "canonical_rule": policy["canonical_rule"], "claim_limit": policy["claim_limit"],
        "canonical_output_sha256": {
            str(path.relative_to(root)): hashlib.sha256(payload.encode("utf-8")).hexdigest()
            for path, payload in payloads.items()
        },
        "script_sha256": {
            relative: digest(root / relative) for relative in (
                "scripts/43_prepare_pleiotropy_pairs.py",
                "scripts/44_materialize_pleiotropy_pair.py",
                "scripts/45_prepare_placo_task.py",
                "scripts/46_run_placo_pair.R",
                "scripts/47_prepare_pleiofdr_trait.py",
                "scripts/48_prepare_conjfdr_task.py",
                "scripts/49_run_conjfdr_pair.py",
                "scripts/49_seal_conjfdr_result.py",
                "scripts/50_collate_pleiotropy.py",
                "scripts/pleiotropy_contract.py",
            )
        },
    }
    provenance_text = json.dumps(provenance, indent=2, sort_keys=True) + "\n"
    provenance_path = root / args.provenance_out
    if args.report_only:
        if not args.quiet:
            print(
                f"Pleiotropy collation ready: {len(pair_ids)} pairs; {len(shared)} "
                "cross-method shared loci; no canonical files written"
            )
        return 0
    if args.validate_only:
        for path, payload in payloads.items():
            if not path.is_file() or path.read_text(encoding="utf-8") != payload:
                raise SystemExit(f"ERROR: canonical pleiotropy output drifted: {path}")
        if not provenance_path.is_file() or provenance_path.read_text(encoding="utf-8") != provenance_text:
            raise SystemExit("ERROR: canonical pleiotropy provenance drifted")
    else:
        if provenance_path.exists() or any(path.exists() for path in payloads):
            raise SystemExit("ERROR: immutable canonical pleiotropy publication already exists")
        for path, payload in payloads.items():
            atomic_text(path, payload)
        atomic_text(provenance_path, provenance_text)
    if not args.quiet:
        print(
            f"Pleiotropy complete: {len(pair_ids)} pairs; {len(placo_loci)} PLACO+ loci; "
            f"{len(conj_loci)} conjFDR loci; {len(shared)} cross-method shared loci"
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
