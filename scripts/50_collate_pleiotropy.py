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
from collections import defaultdict
from pathlib import Path


PLACO_FIELDS = [
    "pair_id", "sleep_trait", "non_sleep_trait", "analysis_tier", "locus_id",
    "CHR", "START", "STOP", "lead_snp", "lead_p_placo_plus",
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
    "placo_significant_variant_count", "conjfdr_significant_variant_count",
    "evidence_status", "claim_limit",
]


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


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
    root = Path(args.root).resolve()
    policy_path = root / "config/pleiotropy_analysis_policy.json"
    policy = json.loads(policy_path.read_text(encoding="utf-8"))
    manifest_path = root / args.manifest
    manifest_lock_path = root / args.manifest_lock
    lock = json.loads(manifest_lock_path.read_text(encoding="utf-8"))
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
            root / args.placo_dir / f"{pair_id}.summary.tsv",
            root / args.placo_dir / f"{pair_id}.hits.tsv",
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

    locus_path = root / policy["locus_definition"]
    blocks, starts = load_blocks(
        locus_path, policy["locus_definition_sha256"], policy["expected_loci"]
    )
    placo_grouped: dict[tuple[str, str], list[dict[str, object]]] = defaultdict(list)
    conj_grouped: dict[tuple[str, str], list[dict[str, object]]] = defaultdict(list)
    input_hashes: dict[str, dict[str, str]] = {}
    for pair in manifest:
        pair_id = pair["pair_id"]
        summary_path = root / args.placo_dir / f"{pair_id}.summary.tsv"
        hits_path = root / args.placo_dir / f"{pair_id}.hits.tsv"
        summary = one_row(summary_path)
        if (
            summary.get("pair_id") != pair_id
            or summary.get("analysis_status") != "PLACO_PLUS_COMPLETE"
            or summary.get("analysis_tier") != pair["analysis_tier"]
            or float(summary.get("family_threshold", "nan")) != policy["placo_locked_pair_family_threshold"]
            or float(summary.get("numerical_failure_fraction", "nan"))
            > policy["placo_maximum_numerical_failure_fraction"]
        ):
            raise SystemExit(f"ERROR: invalid PLACO+ summary for {pair_id}")
        placo_hits, placo_fields = rows_and_fields(hits_path)
        required_placo = {
            "pair_id", "SNP", "CHR", "BP", "P_PLACO_PLUS",
            "conventional_significant", "locked_family_significant",
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
                    "start": start, "stop": stop,
                })
        if family_count != int(summary["locked_family_hit_count"]):
            raise SystemExit(f"ERROR: PLACO+ family hit count differs from summary for {pair_id}")

        task_path = root / args.conjfdr_task_dir / f"{pair_id}.tsv"
        task_lock_path = root / args.conjfdr_task_dir / f"{pair_id}.lock.tsv"
        task, task_lock = one_row(task_path), one_row(task_lock_path)
        if task_lock.get("task_sha256") != sha256(task_path) or task.get("pair_id") != pair_id:
            raise SystemExit(f"ERROR: conjunction-FDR task/lock drifted for {pair_id}")
        completion_path = root / task["completion"]
        all_results_path = root / task["all_results"]
        if not completion_path.is_file() or not all_results_path.is_file():
            raise SystemExit(f"ERROR: conjunction-FDR result is absent for {pair_id}")
        completion = one_row(completion_path)
        if (
            completion.get("analysis_status") != "CONJFDR_COMPLETE"
            or completion.get("pair_id") != pair_id
            or completion.get("task_sha256") != sha256(task_path)
            or completion.get("all_results_sha256") != sha256(all_results_path)
            or completion.get("correct_sample_overlap") != "TRUE"
            or int(completion.get("random_prune_iterations", "0"))
            != policy["pleiofdr_random_prune_iterations"]
        ):
            raise SystemExit(f"ERROR: conjunction-FDR completion drifted for {pair_id}")
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
            "placo_summary": sha256(summary_path), "placo_hits": sha256(hits_path),
            "conjfdr_task": sha256(task_path), "conjfdr_completion": sha256(completion_path),
            "conjfdr_all": sha256(all_results_path),
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
        "analysis_id": policy["analysis_id"], "policy_sha256": sha256(policy_path),
        "manifest_sha256": sha256(manifest_path), "manifest_lock_sha256": sha256(manifest_lock_path),
        "locus_definition_sha256": sha256(locus_path), "completed_pair_scans": len(pair_ids),
        "placo_locus_count": len(placo_loci), "conjfdr_locus_count": len(conj_loci),
        "canonical_shared_locus_count": len(shared), "pair_input_hashes": input_hashes,
        "canonical_rule": policy["canonical_rule"], "claim_limit": policy["claim_limit"],
    }
    provenance_text = json.dumps(provenance, indent=2, sort_keys=True) + "\n"
    provenance_path = root / args.provenance_out
    if args.validate_only:
        for path, payload in payloads.items():
            if not path.is_file() or path.read_text(encoding="utf-8") != payload:
                raise SystemExit(f"ERROR: canonical pleiotropy output drifted: {path}")
        if not provenance_path.is_file() or provenance_path.read_text(encoding="utf-8") != provenance_text:
            raise SystemExit("ERROR: canonical pleiotropy provenance drifted")
    else:
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
