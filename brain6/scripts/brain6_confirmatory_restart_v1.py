#!/usr/bin/env python3
"""Fail-closed admission checkpoint for a future exact-phenotype LAVA rescue.

This tool freezes source and review evidence before a new pilot. It never runs
LAVA, writes a protected result, or treats a numerical gain as source admission.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import os
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
ANCHORS = {
    "canonical_decision": (
        ROOT / "work/lava-canonical-v3-production/d730debf45266d298401564f3260bdecb14739d1c1f1835a5aebd615c83fa60b/canonical_family_decision.json",
        "3569ce81d33479321f9b3adec7bd09c9b54a10a01f1fa4abd15888d9d2e66b22",
    ),
    "canonical_family_lock": (ROOT / "brain6/config/lava_family_canonical_v3.json",
                              "f30c46f1eaaaf50d4697cf8af97f3acd44bef2ad8e067c016f3b2568b424b417"),
    "exact_phenotype_contract": (
        ROOT / "brain6/results/lava_longsleep_source_rescue_v1/phenotype_contract_v1.freeze.json",
        "91aa1e55e3a52270aec7280de184395e807c48a03f7e62760bd227c46ae0d3bf",
    ),
    "prespecified_88_loci": (ROOT / "brain6/results/lava_confirmatory_pilots_v1/longsleep_linear_88_loci.tsv",
                             "9182ebca331a7362d96359cc41e130571ce274f7869477864ce7ca0a944a01e7"),
}
PROTECTED = (
    ROOT / "work/lava-canonical-v3-production",
    ROOT / "work/lava-family-v2",
    ROOT / "work/lava-roundoff-v1",
    ROOT / "brain6/results/lava",
    Path("/Volumes/Extreme SSD/brain6-work/lava-confirmatory-pilots-v1"),
)


def sha(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def verified_file(item: object, label: str, errors: list[str]) -> dict[str, str] | None:
    if not isinstance(item, dict) or not item.get("path") or not item.get("sha256"):
        errors.append(f"{label}: path and SHA-256 required")
        return None
    path = Path(str(item["path"])).expanduser()
    expected = str(item["sha256"])
    if not path.is_absolute() or not path.is_file() or len(expected) != 64:
        errors.append(f"{label}: absolute existing file and 64-character SHA-256 required")
        return None
    observed = sha(path)
    if observed != expected:
        errors.append(f"{label}: SHA-256 mismatch")
        return None
    return {"path": str(path.resolve()), "sha256": observed}


def preflight(manifest: dict, manifest_path: Path, *, allow_existing_output: bool = False) -> tuple[list[str], dict]:
    errors: list[str] = []
    anchors = {}
    for name, (path, expected) in ANCHORS.items():
        if not path.is_file() or sha(path) != expected:
            errors.append(f"frozen anchor mismatch: {name}")
        else:
            anchors[name] = {"path": str(path), "sha256": expected}
    decision_path = ANCHORS["canonical_decision"][0]
    if decision_path.is_file():
        decision = json.loads(decision_path.read_text())
        if decision.get("overall_status") != "FAILED_QC_NOT_PROMOTED" or decision.get("maximum_allowed_untested_cells") != 873:
            errors.append("frozen family decision/status/ceiling changed")
    if manifest.get("schema_version") != 1:
        errors.append("manifest schema_version must be 1")
    analysis_id = manifest.get("analysis_id")
    if not isinstance(analysis_id, str) or not analysis_id.startswith("brain6-confirmatory-rescue-") or "REPLACE" in analysis_id or analysis_id in (
        "brain6-confirmatory-rescue-v1", "brain6-lava-canonical-v3"):
        errors.append("new, specific brain6-confirmatory-rescue-* analysis_id required")
    if manifest.get("phenotype") != "Dashti_2019_long_sleep_ge9h_vs_7to8h":
        errors.append("exact frozen long-sleep phenotype identity required")
    if manifest.get("source_kind") not in ("ORIGINAL_RELEASE_WITH_VERIFIED_VARIANT_N", "SOURCE_VERIFIED_EXACT_PHENOTYPE_RERUN"):
        errors.append("original-plus-variant-N or exact-phenotype rerun source kind required")
    if manifest.get("ancestry") != "EUR" or manifest.get("genome_build") != "GRCh37":
        errors.append("EUR and GRCh37 source contract required before this workflow")
    if manifest.get("n_semantics") not in ("LITERAL_PER_VARIANT_ANALYZED_N", "SOURCE_VERIFIED_CONSTANT_ANALYZED_N"):
        errors.append("source-verified literal analyzed-N semantics required")
    if manifest.get("source_association_model") in (None, "", "UNKNOWN"):
        errors.append("source association model required")
    if manifest.get("worker_count") != 4:
        errors.append("frozen four-worker scheduler setting required")
    if manifest.get("lava_statistic_mapping") in (None, "", "UNKNOWN"):
        errors.append("reviewed source-to-LAVA statistic mapping required")
    files = []
    source_files = manifest.get("source_files")
    if not isinstance(source_files, list) or not source_files:
        errors.append("at least one source file and hash required")
    else:
        for index, item in enumerate(source_files):
            found = verified_file(item, f"source_files[{index}]", errors)
            if found:
                files.append(found)
    if len({item["sha256"] for item in files}) != len(files):
        errors.append("source file hashes must be unique")
    if manifest.get("source_kind") == "ORIGINAL_RELEASE_WITH_VERIFIED_VARIANT_N":
        if "0afd9a3ccf459d283f1b6eb6a8b8b0780f54d9d634a05f816be31171a42ec885" not in {
            item["sha256"] for item in files
        } or len(files) < 2:
            errors.append("original-source route requires the frozen Dashti archive plus a distinct N source")
    custodian = verified_file(manifest.get("custodian_or_primary_source_document"), "custodian_or_primary_source_document", errors)
    review = verified_file(manifest.get("signed_method_review"), "signed_method_review", errors)
    if review:
        try:
            review_json = json.loads(Path(review["path"]).read_text())
        except (OSError, UnicodeError, json.JSONDecodeError):
            errors.append("signed_method_review must be structured JSON")
        else:
            if review_json.get("status") != "ADMITTED_BEFORE_PILOT" or not review_json.get("reviewer") or not review_json.get(
                "reviewer_signature_or_attestation"):
                errors.append("method review requires named reviewer, attestation, and pre-pilot admission")
            if sorted(review_json.get("source_file_sha256", [])) != sorted(item["sha256"] for item in files):
                errors.append("method review does not bind all exact source file hashes")
            if review_json.get("n_semantics") != manifest.get("n_semantics") or review_json.get(
                "lava_statistic_mapping") != manifest.get("lava_statistic_mapping"):
                errors.append("method review N/mapping differs from manifest")
            if not review_json.get("case_fraction_or_continuous_model_rationale") or not review_json.get(
                "phenotype_equivalence_rationale") or not review_json.get("model_to_lava_rationale"):
                errors.append("method review needs phenotype, model, and case-fraction rationale")
    required_review = (
        "exact_phenotype_pass", "source_rights_pass", "variant_build_allele_pass",
        "effect_statistic_model_pass", "analyzed_n_pass", "case_fraction_or_continuous_model_pass",
        "lava_mapping_pass", "source_independent_of_candidate_results_pass",
    )
    if any(manifest.get(flag) is not True for flag in required_review):
        errors.append("all pre-output scientific review flags must be true; numerical gain cannot substitute")
    output_root = manifest.get("new_output_root")
    output = None
    if not isinstance(output_root, str) or not Path(output_root).is_absolute():
        errors.append("new absolute output root required")
    else:
        output = Path(output_root).expanduser().resolve()
        if "REPLACE" in output_root or (isinstance(analysis_id, str) and output.name != analysis_id):
            errors.append("output directory name must equal the final source-specific analysis_id")
        drive = Path("/Volumes/Extreme SSD")
        if (output == drive or drive in output.parents) and not os.path.ismount(drive):
            errors.append("Extreme SSD is not mounted; refuse to create an internal-disk lookalike path")
        if output.exists() and not allow_existing_output:
            errors.append("new output root already exists; never reuse or overwrite a prior run")
        for protected in PROTECTED:
            p = protected.resolve()
            if output == p or output.is_relative_to(p) or p.is_relative_to(output):
                errors.append("new output root overlaps protected outputs")
                break
    freeze = {
        "schema_version": 1, "analysis_id": analysis_id, "admission": "SOURCE_AND_METHOD_REVIEW_PASS" if not errors else "BLOCKED",
        "manifest_path": str(manifest_path.resolve()), "manifest_sha256": sha(manifest_path),
        "frozen_anchors": anchors, "verified_source_files": files,
        "custodian_or_primary_source_document": custodian,
        "signed_method_review": review,
        "new_output_root": str(output) if output else None,
        "next_stage": "PROSPECTIVE_88_LOCUS_PILOT_WITH_FROZEN_RULE" if not errors else "RESOLVE_ADMISSION_ERRORS",
        "note": "Admission is not a LAVA result or authorization for a full-family run.",
    }
    return errors, freeze


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", required=True, type=Path)
    parser.add_argument("--freeze-admission", action="store_true", help="write a new admission receipt only when all checks pass")
    args = parser.parse_args()
    manifest = json.loads(args.manifest.read_text())
    errors, freeze = preflight(manifest, args.manifest)
    result = {"status": "BLOCKED" if errors else "READY_FOR_PROSPECTIVE_PILOT", "errors": errors, "freeze": freeze}
    print(json.dumps(result, indent=2, sort_keys=True))
    if errors:
        return 2
    if args.freeze_admission:
        out = Path(freeze["new_output_root"])
        out.mkdir(parents=True, exist_ok=False)
        receipt = out / "source_admission.freeze.json"
        flags = os.O_CREAT | os.O_EXCL | os.O_WRONLY
        descriptor = os.open(receipt, flags, 0o444)
        with os.fdopen(descriptor, "w") as stream:
            stream.write(json.dumps(freeze, indent=2, sort_keys=True) + "\n")
            stream.flush()
            os.fsync(stream.fileno())
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
