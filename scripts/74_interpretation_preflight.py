#!/usr/bin/env python3
"""Audit locked regulatory, cell, pathway, and causal interpretation inputs."""
from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
import re
import tarfile
from pathlib import Path


REGISTRY_FIELDS = [
    "source_id", "analysis_family", "method", "layer_or_resource",
    "coverage_domains", "exact_release", "build", "access_mode", "source_url",
    "local_path", "expected_bytes", "expected_sha256", "source_status", "notes",
]
READINESS_FIELDS = [
    "source_id", "analysis_family", "method", "layer_or_resource",
    "coverage_domains", "exact_release", "local_path", "observed_bytes",
    "observed_sha256", "readiness_status", "blocker",
]
SHA256 = re.compile(r"^[0-9a-f]{64}$")
FAMILIES = {"regulatory", "cell_type", "pathway", "causal"}
SCREEN_METADATA = re.compile(
    r'\{"name":"([^"]*)","collection":"Core","ontology":"([^"]*)",'
    r'"lifeStage":"([^"]*)","sampleType":"([^"]*)","displayName":"([^"]*)",'
    r'"assays":\[\{"id":"dnase-[^"]*","assay":"dnase","url":"[^"]*",'
    r'"experimentAccession":"(ENCSR[^"]*)"'
)


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


def atomic_tsv(path: Path, fields: list[str], rows: list[dict[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".tmp")
    with temporary.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    temporary.replace(path)


def atomic_json(path: Path, value: object) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(json.dumps(value, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    temporary.replace(path)


def split_domains(value: str) -> set[str]:
    return {item for item in value.split(";") if item and item != "NA"}


def validate_screen_bundle(root: Path, policy: dict[str, object]) -> tuple[bool, str, dict[str, object]]:
    spec = policy["screen_registry_v4"]
    manifest_path = root / spec["component_manifest"]
    if not manifest_path.is_file():
        return False, "SCREEN component manifest is absent", {}
    if sha256(manifest_path) != spec["component_manifest_sha256"]:
        return False, "SCREEN component manifest differs from the policy pin", {}
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    observed: dict[str, object] = {"manifest_sha256": sha256(manifest_path), "components": {}}
    components: dict[str, Path] = {}
    for component in manifest.get("components", []):
        identity = component.get("component_id", "UNKNOWN")
        relative = Path(str(component.get("path", "")))
        if relative.is_absolute() or ".." in relative.parts:
            return False, f"unsafe SCREEN component path for {identity}", observed
        path = root / relative
        if not path.is_file():
            return False, f"SCREEN component is absent: {identity}", observed
        actual_bytes, actual_hash = path.stat().st_size, sha256(path)
        observed["components"][identity] = {"bytes": actual_bytes, "sha256": actual_hash}
        if actual_bytes != component.get("bytes") or actual_hash != component.get("sha256"):
            return False, f"SCREEN component differs from its exact pin: {identity}", observed
        components[identity] = path
    required = {
        "GRCH38_CCRE_COORDINATES", "CORE_COLLECTION_CLASS_MATRIX",
        "DOWNLOADS_BIOSAMPLE_METADATA_SNAPSHOT",
    }
    if set(components) != required:
        return False, "SCREEN component manifest is not the exact required three-file bundle", observed
    bed_rows = sum(1 for _ in components["GRCH38_CCRE_COORDINATES"].open("rb"))
    if bed_rows != manifest.get("ccre_count"):
        return False, "SCREEN coordinate row count differs from the release manifest", observed
    with gzip.open(components["CORE_COLLECTION_CLASS_MATRIX"], "rb") as handle:
        header = handle.readline()
        header_fields = header.rstrip(b"\n").decode("ascii").split("\t")
        matrix_lines = 1 + sum(1 for _ in handle)
    if (
        len(header_fields) != manifest.get("matrix_column_count")
        or matrix_lines != manifest.get("matrix_line_count_including_header")
        or hashlib.sha256(header).hexdigest() != manifest.get("matrix_header_sha256")
        or len(header_fields[1:]) != manifest.get("core_biosample_count")
    ):
        return False, "SCREEN Core Collection matrix shape/header differs from its release pin", observed
    metadata_text = components["DOWNLOADS_BIOSAMPLE_METADATA_SNAPSHOT"].read_text(encoding="utf-8")
    metadata_rows = SCREEN_METADATA.findall(metadata_text)
    metadata_accessions = [row[5] for row in metadata_rows]
    if (
        len(metadata_rows) != manifest.get("core_biosample_count")
        or len(metadata_accessions) != len(set(metadata_accessions))
        or set(metadata_accessions) != set(header_fields[1:])
    ):
        return False, "SCREEN Core Collection metadata does not exactly match the matrix header", observed
    observed.update({"coordinate_rows": bed_rows, "matrix_lines": matrix_lines, "core_biosamples": len(metadata_rows)})
    return True, "", observed


def validate_hocomoco_bundle(root: Path, policy: dict[str, object]) -> tuple[bool, str, dict[str, object]]:
    spec = policy["hocomoco_v14"]
    manifest_path = root / spec["component_manifest"]
    if not manifest_path.is_file() or sha256(manifest_path) != spec["component_manifest_sha256"]:
        return False, "HOCOMOCO component manifest is absent or differs from the policy pin", {}
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    components: dict[str, Path] = {}
    observed: dict[str, object] = {"manifest_sha256": sha256(manifest_path), "components": {}}
    for component in manifest.get("components", []):
        identity = component.get("component_id", "UNKNOWN")
        relative = Path(str(component.get("path", "")))
        if relative.is_absolute() or ".." in relative.parts:
            return False, f"unsafe HOCOMOCO component path for {identity}", observed
        path = root / relative
        if not path.is_file():
            return False, f"HOCOMOCO component is absent: {identity}", observed
        actual_bytes, actual_hash = path.stat().st_size, sha256(path)
        observed["components"][identity] = {"bytes": actual_bytes, "sha256": actual_hash}
        if actual_bytes != component.get("bytes") or actual_hash != component.get("sha256"):
            return False, f"HOCOMOCO component differs from its exact pin: {identity}", observed
        components[identity] = path
    required = {"H14CORE_MEME", "H14CORE_ANNOTATION", "H14CORE_PWM", "H14CORE_THRESHOLDS"}
    if set(components) != required:
        return False, "HOCOMOCO manifest is not the exact required four-file bundle", observed
    motif_count = sum(
        line.startswith("MOTIF ")
        for line in components["H14CORE_MEME"].read_text(encoding="utf-8").splitlines()
    )
    annotation = [
        json.loads(line) for line in components["H14CORE_ANNOTATION"].read_text(encoding="utf-8").splitlines()
        if line
    ]
    annotation_names = [row.get("name") for row in annotation]
    with tarfile.open(components["H14CORE_PWM"], "r:gz") as archive:
        pwm_names = [Path(member.name).stem for member in archive.getmembers() if member.isfile() and member.name.endswith(".pwm")]
    with tarfile.open(components["H14CORE_THRESHOLDS"], "r:gz") as archive:
        threshold_names = [Path(member.name).stem for member in archive.getmembers() if member.isfile() and member.name.endswith(".thr")]
    expected = int(spec["expected_motifs"])
    if (
        motif_count != expected or len(annotation_names) != expected
        or len(set(annotation_names)) != expected or set(annotation_names) != set(pwm_names)
        or set(annotation_names) != set(threshold_names)
    ):
        return False, "HOCOMOCO motif, annotation, PWM, and threshold identities are not the exact H14CORE family", observed
    observed["motif_count"] = motif_count
    return True, "", observed


def validate_abc_bundle(root: Path, policy: dict[str, object]) -> tuple[bool, str, dict[str, object]]:
    spec = policy["abc_2021"]
    manifest_path = root / spec["component_manifest"]
    if not manifest_path.is_file() or sha256(manifest_path) != spec["component_manifest_sha256"]:
        return False, "ABC component manifest is absent or differs from the policy pin", {}
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    components = manifest.get("components", [])
    observed: dict[str, object] = {"manifest_sha256": sha256(manifest_path), "components": {}}
    if len(components) != 1 or components[0].get("component_id") != "ABC_ALL_PREDICTIONS":
        return False, "ABC manifest is not the exact required one-file bundle", observed
    component = components[0]
    relative = Path(str(component.get("path", "")))
    if relative.is_absolute() or ".." in relative.parts:
        return False, "unsafe ABC component path", observed
    path = root / relative
    if not path.is_file():
        return False, "ABC all-predictions source is absent", observed
    actual_bytes, actual_hash = path.stat().st_size, sha256(path)
    observed["components"]["ABC_ALL_PREDICTIONS"] = {
        "bytes": actual_bytes, "sha256": actual_hash,
    }
    if actual_bytes != component.get("bytes") or actual_hash != component.get("sha256"):
        return False, "ABC all-predictions source differs from its exact pin", observed
    expected_header = [
        "chr", "start", "end", "name", "class", "activity_base", "TargetGene",
        "TargetGeneTSS", "TargetGeneExpression", "TargetGenePromoterActivityQuantile",
        "TargetGeneIsExpressed", "distance", "isSelfPromoter", "hic_contact",
        "powerlaw_contact", "powerlaw_contact_reference", "hic_contact_pl_scaled",
        "hic_pseudocount", "hic_contact_pl_scaled_adj", "ABC.Score.Numerator",
        "ABC.Score", "powerlaw.Score.Numerator", "powerlaw.Score", "CellType",
    ]
    rows = 0
    biosamples: set[str] = set()
    genes: set[str] = set()
    minimum_score, maximum_score = float("inf"), float("-inf")
    try:
        with gzip.open(path, "rt", encoding="utf-8", newline="") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            if reader.fieldnames != expected_header:
                return False, "ABC header differs from the frozen 24-column schema", observed
            raw_header = ("\t".join(reader.fieldnames) + "\n").encode()
            if hashlib.sha256(raw_header).hexdigest() != manifest.get("header_sha256"):
                return False, "ABC header differs from its release pin", observed
            for row in reader:
                if len(row) != len(expected_header) or None in row:
                    return False, "ABC source contains a malformed row", observed
                score = float(row["ABC.Score"])
                if not 0 <= score <= 1:
                    return False, "ABC source contains an invalid score", observed
                rows += 1
                biosamples.add(row["CellType"])
                genes.add(row["TargetGene"])
                minimum_score = min(minimum_score, score)
                maximum_score = max(maximum_score, score)
    except (OSError, UnicodeError, ValueError, KeyError) as exc:
        return False, f"ABC source is not a valid compressed prediction table: {exc}", observed
    biosample_hash = hashlib.sha256(("\n".join(sorted(biosamples)) + "\n").encode()).hexdigest()
    domain_map = manifest.get("domain_biosamples", {})
    mapped = [value for domain in policy["regulatory_mapping"]["required_context_domains"] for value in domain_map.get(domain, [])]
    if (
        rows != manifest.get("row_count") or rows != spec["expected_rows"]
        or len(biosamples) != manifest.get("biosample_count")
        or len(biosamples) != spec["expected_biosamples"]
        or len(genes) != manifest.get("gene_symbol_count")
        or biosample_hash != manifest.get("sorted_biosample_list_sha256")
        or minimum_score != manifest.get("minimum_released_abc_score")
        or maximum_score != manifest.get("maximum_observed_abc_score")
        or set(domain_map) != set(policy["regulatory_mapping"]["required_context_domains"])
        or not mapped or len(mapped) != len(set(mapped)) or not set(mapped).issubset(biosamples)
    ):
        return False, "ABC rows, scores, genes, biosamples, or domain map differ from the frozen release", observed
    observed.update({
        "rows": rows, "biosamples": len(biosamples), "mapped_biosamples": len(mapped),
        "gene_symbols": len(genes), "minimum_score": minimum_score,
        "maximum_score": maximum_score, "sorted_biosample_list_sha256": biosample_hash,
    })
    return True, "", observed


def validate_policy_alignment(policy: dict[str, object], downstream: dict[str, object]) -> None:
    pairs = (
        (policy["regulatory_mapping"]["required_layers"], downstream["regulatory_mapping"]["required_layers"]),
        (policy["regulatory_mapping"]["required_context_domains"], downstream["regulatory_mapping"]["required_context_domains"]),
        (policy["cell_types"]["required_strategies"], downstream["cell_types"]["required_strategies"]),
        (policy["cell_types"]["required_domains"], downstream["cell_types"]["required_domains"]),
        (policy["pathways"]["required_resources"], downstream["pathways"]["required_resources"]),
        (policy["causal_inference"]["directions"], downstream["causal_inference"]["directions"]),
        (policy["causal_inference"]["methods"], downstream["causal_inference"]["methods"]),
        (policy["robustness"]["required_families"], downstream["robustness"]["required_families"]),
    )
    if any(left != right for left, right in pairs):
        fail("interpretation and downstream policies define different locked method families")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", default=".")
    parser.add_argument("--policy", default="config/interpretation_analysis_policy.json")
    parser.add_argument("--downstream-policy", default="config/downstream_analysis_policy.json")
    parser.add_argument("--readiness-out", default="results/tables/interpretation_source_readiness.tsv")
    parser.add_argument("--report-out", default="results/tables/interpretation_preflight.json")
    parser.add_argument("--report-only", action="store_true")
    parser.add_argument("--quiet", action="store_true")
    args = parser.parse_args()
    root = Path(args.root).resolve()
    policy_path = root / args.policy
    policy = json.loads(policy_path.read_text(encoding="utf-8"))
    downstream_path = root / args.downstream_policy
    downstream = json.loads(downstream_path.read_text(encoding="utf-8"))
    validate_policy_alignment(policy, downstream)

    panel_fields, panel = read_tsv(root / "config/analysis_panel.tsv")
    del panel_fields
    traits = [row["trait_id"] for row in panel]
    sleep_count = sum(row["domain"] == "sleep" for row in panel)
    if (
        len(panel) != policy["expected_traits"] or len(set(traits)) != len(traits)
        or sleep_count != policy["expected_sleep_traits"]
        or len(panel) - sleep_count != policy["expected_non_sleep_traits"]
    ):
        fail("interpretation workflow requires the exact locked 45/12/33 panel")

    registry_path = root / policy["source_registry"]
    references_path = root / policy["method_references"]
    reference_fields, references = read_tsv(references_path)
    if reference_fields != ["reference_id", "analysis_layer", "title", "doi", "pmid", "url", "role"] or not references:
        fail("interpretation method-reference registry has an unexpected schema")
    if any(row["analysis_layer"] not in {"regulatory", "cell_type", "pathway", "causal"} for row in references):
        fail("interpretation method-reference registry contains an unknown layer")
    fields, registry = read_tsv(registry_path)
    if fields != REGISTRY_FIELDS or not registry:
        fail("interpretation source registry has an unexpected schema")
    source_ids = [row["source_id"] for row in registry]
    if len(source_ids) != len(set(source_ids)):
        fail("interpretation source registry contains duplicate source IDs")
    if {row["analysis_family"] for row in registry} != FAMILIES:
        fail("interpretation source registry does not cover the exact four analysis families")

    regulatory = [row for row in registry if row["analysis_family"] == "regulatory"]
    cells = [row for row in registry if row["analysis_family"] == "cell_type"]
    pathways = [row for row in registry if row["analysis_family"] == "pathway"]
    causal = [row for row in registry if row["analysis_family"] == "causal"]
    required_domains = set(policy["regulatory_mapping"]["required_context_domains"])
    if {row["layer_or_resource"] for row in regulatory} != set(policy["regulatory_mapping"]["required_layers"]):
        fail("source registry does not cover each locked regulatory layer exactly")
    if any(not required_domains.issubset(split_domains(row["coverage_domains"])) for row in regulatory):
        fail("a regulatory source omits a locked context domain")
    if {row["method"] for row in cells} != set(policy["cell_types"]["required_strategies"]):
        fail("source registry does not cover each locked cell-type strategy exactly")
    if any(not set(policy["cell_types"]["required_domains"]).issubset(split_domains(row["coverage_domains"])) for row in cells):
        fail("a cell-type source omits a locked biological domain")
    if {row["layer_or_resource"] for row in pathways} != set(policy["pathways"]["required_resources"]):
        fail("source registry does not cover each locked pathway resource exactly")
    if {row["method"] for row in causal} != set(policy["causal_inference"]["methods"]):
        fail("source registry does not cover each locked causal estimator family exactly")

    readiness: list[dict[str, object]] = []
    screen_ready, screen_blocker, screen_validation = validate_screen_bundle(root, policy)
    screen_source_ids = set(policy["screen_registry_v4"]["source_ids"])
    hocomoco_ready, hocomoco_blocker, hocomoco_validation = validate_hocomoco_bundle(root, policy)
    hocomoco_source_id = policy["hocomoco_v14"]["source_id"]
    abc_ready, abc_blocker, abc_validation = validate_abc_bundle(root, policy)
    abc_source_id = policy["abc_2021"]["source_id"]
    allowed_source_status = {
        "SOURCE_VERIFIED", "CURATION_REQUIRED", "DERIVED_UPSTREAM",
        "DERIVED_WITHIN_WORKFLOW",
    }
    for row in registry:
        identity = row["source_id"]
        if row["source_status"] not in allowed_source_status:
            fail(f"unknown source status for {identity}")
        relative = Path(row["local_path"])
        if relative.is_absolute() or ".." in relative.parts:
            fail(f"unsafe source path for {identity}")
        path = root / relative
        observed_bytes = path.stat().st_size if path.is_file() else 0
        observed_hash = sha256(path) if path.is_file() else "NA"
        ready = False
        blocker = ""
        if row["source_status"] == "SOURCE_VERIFIED":
            if not row["expected_bytes"].isdigit() or not SHA256.fullmatch(row["expected_sha256"]):
                fail(f"verified source lacks an exact byte/SHA256 pin: {identity}")
            ready = (
                path.is_file() and observed_bytes == int(row["expected_bytes"])
                and observed_hash == row["expected_sha256"]
            )
            if not ready:
                blocker = "local source is absent or differs from exact bytes/SHA256"
            elif identity in screen_source_ids and not screen_ready:
                ready = False
                blocker = screen_blocker
            elif identity == hocomoco_source_id and not hocomoco_ready:
                ready = False
                blocker = hocomoco_blocker
            elif identity == abc_source_id and not abc_ready:
                ready = False
                blocker = abc_blocker
        elif row["source_status"] == "DERIVED_UPSTREAM":
            ready = path.is_file() and observed_bytes > 0
            if not ready:
                blocker = "required upstream atlas artifact is not complete"
        elif row["source_status"] == "DERIVED_WITHIN_WORKFLOW":
            ready = True
        else:
            blocker = "source release, exact bytes, and SHA256 require pre-result curation"
        readiness.append({
            "source_id": identity, "analysis_family": row["analysis_family"],
            "method": row["method"], "layer_or_resource": row["layer_or_resource"],
            "coverage_domains": row["coverage_domains"], "exact_release": row["exact_release"],
            "local_path": row["local_path"], "observed_bytes": observed_bytes,
            "observed_sha256": observed_hash, "readiness_status": "READY" if ready else "BLOCKED",
            "blocker": blocker,
        })

    upstream_paths = [
        "results/atlas/traits.tsv", "results/atlas/trait_pairs.tsv", "results/atlas/loci.tsv",
        "results/atlas/variants.tsv", "results/atlas/genes.tsv",
        "results/tables/molecular_locus_coverage.tsv",
    ]
    upstream = {
        relative: {
            "ready": (root / relative).is_file() and (root / relative).stat().st_size > 0,
            "bytes": (root / relative).stat().st_size if (root / relative).is_file() else 0,
        }
        for relative in upstream_paths
    }
    sources_ready = all(row["readiness_status"] == "READY" for row in readiness)
    upstream_ready = all(value["ready"] for value in upstream.values())
    report = {
        "schema_version": policy["schema_version"], "analysis_id": policy["analysis_id"],
        "policy_sha256": sha256(policy_path), "downstream_policy_sha256": sha256(downstream_path),
        "source_registry_sha256": sha256(registry_path),
        "method_references_sha256": sha256(references_path), "code_ready": True,
        "screen_source_bundle": screen_validation,
        "hocomoco_source_bundle": hocomoco_validation,
        "abc_source_bundle": abc_validation,
        "sources_ready": sources_ready, "upstream_ready": upstream_ready,
        "production_ready": sources_ready and upstream_ready,
        "ready_source_count": sum(row["readiness_status"] == "READY" for row in readiness),
        "source_count": len(readiness), "upstream": upstream,
    }
    atomic_tsv(root / args.readiness_out, READINESS_FIELDS, readiness)
    atomic_json(root / args.report_out, report)
    if not args.quiet:
        print(
            "INTERPRETATION_PREFLIGHT_"
            + ("READY" if report["production_ready"] else "BLOCKED")
            + f" sources={report['ready_source_count']}/{report['source_count']} "
            + f"upstream={sum(value['ready'] for value in upstream.values())}/{len(upstream)}"
        )
        for row in readiness:
            if row["readiness_status"] != "READY":
                print(f"BLOCKED {row['source_id']}: {row['blocker']}")
    return 0 if report["production_ready"] or args.report_only else 1


if __name__ == "__main__":
    raise SystemExit(main())
