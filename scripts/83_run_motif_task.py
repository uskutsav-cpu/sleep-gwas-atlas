#!/usr/bin/env python3
"""Run one context-gated HOCOMOCO H14CORE allele-specific motif task."""
from __future__ import annotations

import argparse
from bisect import bisect_right
import csv
import hashlib
import io
import json
import math
from pathlib import Path
import tarfile
from urllib.error import HTTPError, URLError
from urllib.request import Request, urlopen

from liftover_chain import load_chain, reverse_complement


BASE_INDEX = {"A": 0, "C": 1, "G": 2, "T": 3}
MISSING = {"", "NA"}


def fail(message: str) -> None:
    raise SystemExit(f"ERROR: {message}")


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def read_tsv(path: Path) -> tuple[list[str], list[dict[str, str]]]:
    if not path.is_file() or path.stat().st_size == 0:
        fail(f"missing real non-empty artifact: {path}")
    with path.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        return reader.fieldnames or [], list(reader)


def table_text(fields: list[str], rows: list[dict[str, str]]) -> str:
    output = io.StringIO(newline="")
    writer = csv.DictWriter(output, fieldnames=fields, delimiter="\t", lineterminator="\n")
    writer.writeheader()
    writer.writerows(rows)
    return output.getvalue()


def atomic_text(path: Path, payload: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(payload, encoding="utf-8")
    temporary.replace(path)


def fmt(value: float) -> str:
    return format(value, ".12g")


def evidence_id(*values: str) -> str:
    return "AUTOMOTIF::" + hashlib.sha256("\x1f".join(values).encode()).hexdigest()[:24]


def eligible_variants(rows: list[dict[str, str]], locus_id: str) -> list[dict[str, str]]:
    return [
        row for row in rows
        if row["locus_id"] == locus_id and row["qc_status"] == "PASS"
        and any(row[field] not in MISSING for field in (
            "credible_set_sleep", "credible_set_non_sleep", "shared_signal_posterior",
        ))
    ]


def hocomoco_components(root: Path, policy: dict[str, object]) -> tuple[dict[str, Path], dict[str, object]]:
    spec = policy["hocomoco_v14"]
    manifest_path = root / spec["component_manifest"]
    if not manifest_path.is_file() or sha256(manifest_path) != spec["component_manifest_sha256"]:
        fail("HOCOMOCO component manifest is absent or differs from the policy pin")
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    components: dict[str, Path] = {}
    for component in manifest["components"]:
        path = root / component["path"]
        if (
            not path.is_file() or path.stat().st_size != component["bytes"]
            or sha256(path) != component["sha256"]
        ):
            fail(f"HOCOMOCO component differs from its exact pin: {component['component_id']}")
        components[component["component_id"]] = path
    return components, manifest


def load_motifs(
    pwm_archive: Path, threshold_archive: Path, annotation_path: Path,
) -> dict[str, dict[str, object]]:
    annotations = {
        row["name"]: row
        for row in (
            json.loads(line) for line in annotation_path.read_text(encoding="utf-8").splitlines() if line
        )
    }
    motifs: dict[str, dict[str, object]] = {}
    with tarfile.open(pwm_archive, "r:gz") as archive:
        for member in archive.getmembers():
            if not member.isfile() or not member.name.endswith(".pwm"):
                continue
            handle = archive.extractfile(member)
            if handle is None:
                fail(f"could not read HOCOMOCO PWM member: {member.name}")
            lines = handle.read().decode("ascii").strip().splitlines()
            name = lines[0].removeprefix(">").strip()
            matrix = [[float(value) for value in line.split("\t")] for line in lines[1:]]
            if not matrix or any(len(row) != 4 or any(not math.isfinite(value) for value in row) for row in matrix):
                fail(f"invalid HOCOMOCO PWM: {name}")
            motifs[name] = {"pwm": matrix}
    with tarfile.open(threshold_archive, "r:gz") as archive:
        for member in archive.getmembers():
            if not member.isfile() or not member.name.endswith(".thr"):
                continue
            name = Path(member.name).stem
            handle = archive.extractfile(member)
            if handle is None or name not in motifs:
                fail(f"threshold has no matching HOCOMOCO PWM: {name}")
            values = [tuple(map(float, line.split("\t"))) for line in handle.read().decode("ascii").splitlines() if line]
            scores, probabilities = zip(*values)
            if (
                list(scores) != sorted(scores) or any(not 0 <= value <= 1 for value in probabilities)
                or any(probabilities[index] < probabilities[index + 1] for index in range(len(probabilities) - 1))
            ):
                fail(f"invalid HOCOMOCO threshold map: {name}")
            motifs[name]["threshold_scores"] = list(scores)
            motifs[name]["threshold_p"] = list(probabilities)
    if set(motifs) != set(annotations) or any("threshold_scores" not in value for value in motifs.values()):
        fail("HOCOMOCO annotation/PWM/threshold motif identities differ")
    for name, motif in motifs.items():
        motif["tf"] = str(annotations[name].get("tf", "NA"))
    return motifs


def score(sequence: str, matrix: list[list[float]]) -> float | None:
    if len(sequence) != len(matrix):
        return None
    total = 0.0
    for base, weights in zip(sequence, matrix):
        index = BASE_INDEX.get(base)
        if index is None:
            return None
        total += weights[index]
    return total


def threshold_p(observed: float, scores: list[float], probabilities: list[float]) -> float:
    index = bisect_right(scores, observed + 1e-12) - 1
    return 1.0 if index < 0 else probabilities[index]


def best_allele_hit(
    reference: str, alternate: str, center: int, motif: dict[str, object], p_threshold: float,
) -> dict[str, object] | None:
    matrix = motif["pwm"]
    width = len(matrix)
    candidates: list[dict[str, object]] = []
    for start in range(center - width + 1, center + 1):
        if start < 0 or start + width > len(reference):
            continue
        ref_slice, alt_slice = reference[start:start + width], alternate[start:start + width]
        for strand, ref_sequence, alt_sequence in (
            ("+", ref_slice, alt_slice),
            ("-", reverse_complement(ref_slice), reverse_complement(alt_slice)),
        ):
            ref_score, alt_score = score(ref_sequence, matrix), score(alt_sequence, matrix)
            if ref_score is None or alt_score is None or ref_score == alt_score:
                continue
            ref_p = threshold_p(ref_score, motif["threshold_scores"], motif["threshold_p"])
            alt_p = threshold_p(alt_score, motif["threshold_scores"], motif["threshold_p"])
            if min(ref_p, alt_p) > p_threshold:
                continue
            candidates.append({
                "start": start, "strand": strand, "ref_score": ref_score, "alt_score": alt_score,
                "ref_p": ref_p, "alt_p": alt_p, "minimum_p": min(ref_p, alt_p),
                "delta": alt_score - ref_score,
            })
    if not candidates:
        return None
    return min(candidates, key=lambda row: (
        row["minimum_p"], -abs(row["delta"]), row["start"], row["strand"],
    ))


def screen_accessible_variants(
    root: Path, task: dict[str, str], policy: dict[str, object], tasks: list[dict[str, str]],
) -> set[str]:
    screen_ids = set(policy["screen_registry_v4"]["source_ids"])
    dependencies = [
        row for row in tasks
        if row["analysis_family"] == "regulatory" and row["locus_id"] == task["locus_id"]
        and row["domain"] == task["domain"] and row["source_id"] in screen_ids
    ]
    if len(dependencies) != len(screen_ids):
        fail("motif task does not have the exact same-domain SCREEN dependency family")
    accessible: set[str] = set()
    for dependency in dependencies:
        result_path, provenance_path = root / dependency["normalized_result_path"], root / dependency["provenance_path"]
        if not result_path.is_file() or not provenance_path.is_file():
            fail(f"motif task SCREEN dependency is incomplete: {dependency['task_id']}")
        provenance = json.loads(provenance_path.read_text(encoding="utf-8"))
        if (
            provenance.get("task_id") != dependency["task_id"]
            or provenance.get("result_sha256") != sha256(result_path)
            or provenance.get("terminal_status") not in {"COMPLETED", "NO_EVIDENCE_FOUND"}
        ):
            fail(f"motif task SCREEN dependency is invalid or blocked: {dependency['task_id']}")
        _, rows = read_tsv(result_path)
        accessible.update(row["variant_id"] for row in rows)
    return accessible


def fetch_sequence_snapshot(
    snapshot_path: Path, requests: list[dict[str, object]], url_template: str,
) -> tuple[list[dict[str, object]], bool, str]:
    if snapshot_path.is_file():
        snapshot = json.loads(snapshot_path.read_text(encoding="utf-8"))
        if snapshot.get("requests") != requests or not isinstance(snapshot.get("responses"), list):
            fail("existing motif sequence snapshot differs from the locked request family")
        return snapshot["responses"], True, ""
    responses: list[dict[str, object]] = []
    try:
        for item in requests:
            url = url_template.format(**item)
            request = Request(url, headers={"User-Agent": "sleep-gwas-atlas/atlas-v1.0"})
            with urlopen(request, timeout=60) as handle:
                raw = handle.read()
            response = json.loads(raw)
            sequence = str(response.get("dna", "")).upper()
            expected_length = int(item["end0"]) - int(item["start0"])
            if (
                response.get("genome") != "hg38" or response.get("chrom") != f"chr{item['chromosome']}"
                or response.get("start") != item["start0"] or response.get("end") != item["end0"]
                or len(sequence) != expected_length or set(sequence) - set("ACGTN")
            ):
                fail(f"invalid UCSC sequence response for chr{item['chromosome']}:{item['start0']}-{item['end0']}")
            responses.append({
                "request": item, "url": url, "raw_response": raw.decode("utf-8"),
                "raw_response_sha256": hashlib.sha256(raw).hexdigest(), "sequence": sequence,
            })
    except (HTTPError, URLError, TimeoutError, json.JSONDecodeError) as error:
        return [], False, f"UCSC hg38 sequence access failed: {type(error).__name__}: {error}"
    snapshot = {"schema_version": "atlas-v1.0-motif-sequences.1", "requests": requests, "responses": responses}
    atomic_text(snapshot_path, json.dumps(snapshot, indent=2, sort_keys=True) + "\n")
    return responses, True, ""


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("task_id")
    parser.add_argument("--root", default=".")
    parser.add_argument("--policy", default="config/interpretation_analysis_policy.json")
    parser.add_argument("--manifest", default="results/tables/interpretation_task_manifest.tsv")
    parser.add_argument("--out", required=True)
    parser.add_argument("--provenance-out", required=True)
    args = parser.parse_args()
    root = Path(args.root).resolve()
    policy_path, manifest_path = root / args.policy, root / args.manifest
    policy = json.loads(policy_path.read_text(encoding="utf-8"))
    _, tasks = read_tsv(manifest_path)
    selected = [row for row in tasks if row["task_id"] == args.task_id]
    if len(selected) != 1:
        fail(f"unknown or duplicate interpretation task: {args.task_id}")
    task = selected[0]
    if (
        task["analysis_family"] != "regulatory" or task["method"] != "motif_overlap"
        or task["source_id"] != policy["hocomoco_v14"]["source_id"]
    ):
        fail("task is not the locked HOCOMOCO regulatory task")
    components, source_manifest = hocomoco_components(root, policy)
    motifs = load_motifs(
        components["H14CORE_PWM"], components["H14CORE_THRESHOLDS"], components["H14CORE_ANNOTATION"],
    )
    if len(motifs) != int(policy["hocomoco_v14"]["expected_motifs"]):
        fail("HOCOMOCO motif family is incomplete")
    _, variants = read_tsv(root / "results/atlas/variants.tsv")
    eligible = eligible_variants(variants, task["locus_id"])
    accessible = screen_accessible_variants(root, task, policy, tasks)
    eligible = [row for row in eligible if row["variant_id"] in accessible]
    chain_spec = policy["regulatory_build_harmonization"]
    chain, chain_provenance = load_chain(
        root / chain_spec["chain_path"], chain_spec["chain_sha256"], chain_spec["chain_bytes"],
    )
    flank = int(policy["hocomoco_v14"]["maximum_motif_width"]) - 1
    mapped: list[dict[str, object]] = []
    exclusions: list[dict[str, str]] = []
    for variant in eligible:
        alleles = [variant["effect_allele"].upper(), variant["other_allele"].upper()]
        if any(len(allele) != 1 or allele not in BASE_INDEX for allele in alleles):
            exclusions.append({"variant_id": variant["variant_id"], "reason": "non_biallelic_snp"})
            continue
        status, target = chain.map_point(variant["chromosome"], variant["position_bp"])
        if status != "mapped" or target is None:
            exclusions.append({"variant_id": variant["variant_id"], "reason": status})
            continue
        if target[2] == "-":
            alleles = [reverse_complement(allele) for allele in alleles]
        position = int(target[1])
        if position <= flank:
            exclusions.append({"variant_id": variant["variant_id"], "reason": "insufficient_left_flank"})
            continue
        mapped.append({
            "variant": variant, "chromosome": int(target[0]), "position": position,
            "alleles_grch38": alleles, "target_strand": target[2],
        })
    request_by_coordinate: dict[tuple[int, int], dict[str, object]] = {}
    for item in mapped:
        key = (int(item["chromosome"]), int(item["position"]))
        request_by_coordinate[key] = {
            "chromosome": key[0], "position1": key[1],
            "start0": key[1] - 1 - flank, "end0": key[1] + flank,
        }
    sequence_requests = [request_by_coordinate[key] for key in sorted(request_by_coordinate)]
    snapshot_path = root / "results/sources/interpretation/motif_sequences" / f"{task['task_id']}.json"
    if sequence_requests:
        responses, access_ok, access_error = fetch_sequence_snapshot(
            snapshot_path, sequence_requests, policy["hocomoco_v14"]["sequence_api_template"],
        )
    else:
        responses, access_ok, access_error = [], True, ""
    sequence_by_coordinate = {
        (int(row["request"]["chromosome"]), int(row["request"]["position1"])): row["sequence"]
        for row in responses
    }
    rows: list[dict[str, str]] = []
    if access_ok:
        center = flank
        for item in mapped:
            variant = item["variant"]
            sequence = sequence_by_coordinate[(int(item["chromosome"]), int(item["position"]))]
            reference_base = sequence[center]
            alleles = item["alleles_grch38"]
            if alleles.count(reference_base) != 1:
                exclusions.append({"variant_id": variant["variant_id"], "reason": "alleles_do_not_identify_unique_hg38_reference"})
                continue
            alternate_base = alleles[1] if alleles[0] == reference_base else alleles[0]
            alternate = sequence[:center] + alternate_base + sequence[center + 1:]
            for motif_name, motif in motifs.items():
                hit = best_allele_hit(
                    sequence, alternate, center, motif,
                    float(policy["hocomoco_v14"]["primary_motif_p_threshold"]),
                )
                if hit is None:
                    continue
                start1 = int(item["position"]) - center + int(hit["start"])
                end1 = start1 + len(motif["pwm"]) - 1
                element = f"H14CORE::{motif_name}::chr{item['chromosome']}:{start1}-{end1}:{hit['strand']}"
                rows.append({
                    "regulatory_evidence_id": evidence_id(task["task_id"], variant["variant_id"], element),
                    "regulatory_element_id": element, "variant_id": variant["variant_id"],
                    "locus_id": task["locus_id"], "element_type": "TF_motif",
                    "annotation": (
                        f"TF::{motif['tf']}::REF::{reference_base}::ALT::{alternate_base}::"
                        f"REF_P::{fmt(hit['ref_p'])}::ALT_P::{fmt(hit['alt_p'])}"
                    ),
                    "effect": fmt(hit["delta"]), "p_value": fmt(hit["minimum_p"]),
                    "biosample": "SCREEN_CORE_COLLECTION_DOMAIN_AGGREGATE", "tissue": task["domain"],
                    "cell_type": "NA", "context_domain": task["domain"], "target_gene_id": "NA",
                    "link_method": "HOCOMOCO_H14CORE_ALLELE_PWM",
                    "source_dataset": task["source_id"], "source_version": task["source_release"],
                    "evidence_level": "PREDICTED_TF_MOTIF_ALLELE_EFFECT", "provenance_id": "ADAPTER_PENDING",
                })
    rows.sort(key=lambda row: (row["variant_id"], row["regulatory_element_id"], row["regulatory_evidence_id"]))
    fields = policy["regulatory_mapping"]["canonical_fields"]
    if any(list(row) != fields for row in rows):
        fail("automatic motif row differs from the locked regulatory schema")
    payload = table_text(fields, rows)
    out_path = root / args.out
    atomic_text(out_path, payload)
    if not access_ok:
        status, reason = "ACCESS_BLOCKED", access_error
    elif rows:
        status = "COMPLETED"
        reason = f"H14CORE allele scanning produced {len(rows)} context-gated motif evidence rows"
    else:
        status = "NO_EVIDENCE_FOUND"
        reason = "No SCREEN-accessible eligible variant produced an allele-dependent H14CORE hit at the locked P-value threshold"
    provenance = {
        "schema_version": policy["schema_version"], "analysis_id": policy["analysis_id"],
        "task_id": task["task_id"], "terminal_status": status, "terminal_reason": reason,
        "policy_sha256": sha256(policy_path), "task_manifest_sha256": sha256(manifest_path),
        "task_input_scope_sha256": task["input_scope_sha256"], "source_id": task["source_id"],
        "source_release": task["source_release"],
        "source_manifest_release": source_manifest["release"],
        "source_component_sha256": {identity: sha256(path) for identity, path in sorted(components.items())},
        "chain": chain_provenance, "eligible_variant_count": len(eligible),
        "screen_accessible_variant_count": len(accessible), "mapped_variant_count": len(mapped),
        "liftover_or_allele_exclusions": exclusions, "motif_count": len(motifs),
        "sequence_snapshot": str(snapshot_path.relative_to(root)) if snapshot_path.is_file() else "NA",
        "sequence_snapshot_sha256": sha256(snapshot_path) if snapshot_path.is_file() else "NA",
        "result_rows": len(rows), "normalized_result_path": str(out_path),
        "normalized_result_sha256": hashlib.sha256(payload.encode()).hexdigest(),
        "claim_limit": policy["hocomoco_v14"]["claim_limit"],
    }
    provenance_path = root / args.provenance_out
    atomic_text(provenance_path, json.dumps(provenance, indent=2, sort_keys=True) + "\n")
    print(f"MOTIF_ADAPTER_OK task={task['task_id']} status={status} rows={len(rows)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
