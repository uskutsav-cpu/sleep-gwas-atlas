#!/usr/bin/env python3
"""Freeze or verify the result-blind P14/P16 broad-cell execution amendment.

Revision 4 is additive and scoped: revisions 1--3 remain byte-exact audit
records and revision 3 remains authoritative outside the two amended broad
cell-class analyses.  No future fine-mapping, PLACO, cell-enrichment, or other
scientific result is read by this program.
"""
from __future__ import annotations

import argparse
import hashlib
import importlib.util
import json
import os
import stat
import subprocess
import sys
import tempfile
from pathlib import Path
from typing import Callable, Mapping, Sequence

sys.dont_write_bytecode = True

ROOT_AT_IMPORT = Path(__file__).resolve().parents[1]
ENRICHMENT_REL = "scripts/164_run_track_b_cell_class_enrichment.py"
ENRICHMENT_SHA256 = "9d5bb622b5c27c6f38b1cf6aff4402bdcfa847e748465266207f07e6dee06619"
CLASSIFIER_REL = "scripts/165_classify_track_b_cell_replication.py"
CLASSIFIER_SHA256 = "a2cc8fd53f2053a2f9acf44b971959950b861e14675e297095d4568f03c093cc"
SELF_REL = "scripts/163_build_track_b_cell_class_contract.py"
V1_LOCK_REL = "results/track_b/mechanism_followup/PRE_RESULT_READINESS.lock.json"
V1_LOCK_SHA256 = "4ab3c004d109fa856936a5411f4495975b176333916bbe3f70b7974aac9de40f"
V2_LOCK_REL = "results/track_b/mechanism_followup/v2/PRE_RESULT_READINESS.lock.json"
V2_LOCK_SHA256 = "6643561657f69d13e6e5c8e0352d9c095a401943e4028585562f81e3de3ac530"
READINESS_FIELDS = [
    "analysis_id", "phase", "resource_status", "status", "status_reason",
    "required_upstream_gate", "ram_decomposition", "planned_output", "claim_limit",
]


def _bootstrap_hash(path: Path) -> str:
    before_path = os.lstat(path)
    if stat.S_ISLNK(before_path.st_mode) or not stat.S_ISREG(before_path.st_mode):
        raise RuntimeError(f"bootstrap dependency is not a real regular file: {path}")
    fd = os.open(path, os.O_RDONLY | getattr(os, "O_NOFOLLOW", 0))
    try:
        before = os.fstat(fd)
        digest = hashlib.sha256()
        for block in iter(lambda: os.read(fd, 1024 * 1024), b""):
            digest.update(block)
        after = os.fstat(fd)
    finally:
        os.close(fd)
    final = os.lstat(path)
    identity = lambda value: (
        value.st_dev, value.st_ino, value.st_mode, value.st_size,
        value.st_mtime_ns, value.st_ctime_ns,
    )
    if identity(before_path) != identity(before) or identity(before) != identity(after) or identity(after) != identity(final):
        raise RuntimeError("bootstrap dependency drifted while hashing")
    return digest.hexdigest()


for _relative, _digest in (
    (ENRICHMENT_REL, ENRICHMENT_SHA256), (CLASSIFIER_REL, CLASSIFIER_SHA256),
):
    if _bootstrap_hash(ROOT_AT_IMPORT / _relative) != _digest:
        raise RuntimeError(f"revision-4 dependency differs from its bootstrap pin: {_relative}")
_SPEC = importlib.util.spec_from_file_location("track_b_cell_class_enrichment_v4_builder", ROOT_AT_IMPORT / ENRICHMENT_REL)
if _SPEC is None or _SPEC.loader is None:
    raise RuntimeError("cannot load the revision-4 enrichment implementation")
E = importlib.util.module_from_spec(_SPEC)
sys.modules[_SPEC.name] = E
_SPEC.loader.exec_module(E)


def validate_prior_lock(root: Path, relative: str, digest: str, revision: int) -> dict[str, object]:
    record = E.stable_file(root, relative)
    if record.sha256 != digest:
        raise E.ContractError(f"frozen revision-{revision} lock differs from its byte-exact pin")
    try:
        lock = json.loads(E.read_stable_bytes(root, relative).decode("utf-8"))
    except (UnicodeError, json.JSONDecodeError) as exc:
        raise E.ContractError(f"frozen revision-{revision} lock is invalid JSON") from exc
    if not isinstance(lock, dict):
        raise E.ContractError(f"frozen revision-{revision} lock is not an object")
    if (
        lock.get("schema_version") != f"track-b-mechanism-pre-result-readiness.{revision}"
        or (revision == 1 and "contract_revision" in lock)
        or (revision > 1 and lock.get("contract_revision") != revision)
        or lock.get("contract_kind") != "PRE_RESULT_LOCAL_RESOURCE_READINESS_LOCK"
        or lock.get("result_blind") is not True
        or lock.get("future_results_accessed") is not False
    ):
        raise E.ContractError(f"frozen revision-{revision} lock semantics differ")
    outputs = lock.get("outputs")
    if not isinstance(outputs, dict) or not outputs:
        raise E.ContractError(f"frozen revision-{revision} output family is absent")
    for output, expected in outputs.items():
        if not isinstance(expected, dict) or set(expected) != {"bytes", "sha256"}:
            raise E.ContractError(f"frozen revision-{revision} output identity is malformed")
        observed = E.stable_file(root, str(output))
        if observed.bytes != expected["bytes"] or observed.sha256 != expected["sha256"]:
            raise E.ContractError(f"frozen revision-{revision} output differs: {output}")
    return lock


def validate_fuma_manifest(root: Path, policy: Mapping[str, object]) -> dict[str, object]:
    manifest_rel = str(E._static_by_id(policy)["FUMA_MANIFEST"]["path"])
    try:
        manifest = json.loads(E.read_stable_bytes(root, manifest_rel).decode("utf-8"))
    except (UnicodeError, json.JSONDecodeError) as exc:
        raise E.ContractError("FUMA manifest is invalid JSON") from exc
    if (
        not isinstance(manifest, dict)
        or manifest.get("schema_version") != "fuma-scrna-magma-celltyping-source-bundle.1"
        or manifest.get("commit") != "dd526163ea80af1a80a6cdc80db167144500694b"
        or manifest.get("gene_annotation") != "Ensembl_v92_GRCh37_20260_coding_genes"
        or manifest.get("magma_version") != "1.10"
    ):
        raise E.ContractError("FUMA manifest release/build/model semantics differ")
    source_rows = manifest.get("datasets")
    if not isinstance(source_rows, list):
        raise E.ContractError("FUMA manifest dataset family is absent")
    by_id = {row.get("dataset_id"): row for row in source_rows if isinstance(row, dict)}
    prefix = "FUMA_scRNA_data-dd526163ea80af1a80a6cdc80db167144500694b/"
    for dataset_id, spec in E._dataset_specs(policy).items():
        row = by_id.get(dataset_id)
        if not isinstance(row, dict) or (
            spec["archive_member"] != prefix + str(row.get("member"))
            or spec["species"] != row.get("species") or spec["tissue"] != row.get("tissue")
            or spec["gene_rows"] != row.get("gene_rows")
            or spec["cell_type_count"] != row.get("cell_type_count")
            or spec["compressed_bytes"] != row.get("compressed_bytes")
            or spec["compressed_sha256"] != row.get("compressed_sha256")
            or spec["text_bytes"] != row.get("text_bytes")
            or spec["text_sha256"] != row.get("text_sha256")
        ):
            raise E.ContractError(f"V4 source pin is not an exact FUMA manifest member: {dataset_id}")
    return manifest


def read_v3_scope_rows(root: Path) -> dict[str, dict[str, str]]:
    fields, rows, _record = E._read_tsv_stable(
        root, "results/track_b/mechanism_followup/v3/PRE_RESULT_READINESS.tsv",
    )
    if not {"analysis_id", "resource_status", "status", "status_reason"}.issubset(fields):
        raise E.ContractError("V3 readiness table lacks scoped amendment fields")
    wanted = {
        "P14_BRAIN_CELL_CLASS_DISCOVERY", "P14_HUMAN_BRAIN_PERICYTE_DISCOVERY",
        "P15_BRAIN_CELL_SUBTYPE_DISCOVERY", "P16_HUMAN_CELL_CLASS_REPLICATION",
        "P16_HUMAN_CELL_SUBTYPE_REPLICATION",
    }
    selected = {row["analysis_id"]: row for row in rows if row.get("analysis_id") in wanted}
    if set(selected) != wanted:
        raise E.ContractError("V3 P14--P16 readiness scope is incomplete")
    expected = {
        "P14_BRAIN_CELL_CLASS_DISCOVERY": ("BLOCKED_BY_SOFTWARE", "BLOCKED_BY_SOFTWARE"),
        "P14_HUMAN_BRAIN_PERICYTE_DISCOVERY": ("BLOCKED_BY_DATA", "BLOCKED_BY_DATA"),
        "P15_BRAIN_CELL_SUBTYPE_DISCOVERY": ("READY", "NOT_APPLICABLE_UNTIL_UPSTREAM"),
        "P16_HUMAN_CELL_CLASS_REPLICATION": ("BLOCKED_BY_SOFTWARE", "BLOCKED_BY_SOFTWARE"),
        "P16_HUMAN_CELL_SUBTYPE_REPLICATION": ("BLOCKED_BY_DATA", "BLOCKED_BY_DATA"),
    }
    for identity, statuses in expected.items():
        if (selected[identity]["resource_status"], selected[identity]["status"]) != statuses:
            raise E.ContractError(f"V3 scoped status drifted: {identity}")
    return selected


def readiness_rows() -> list[dict[str, object]]:
    common_gate = (
        "An exact canonical track-b-p14-p16-upstream-gate.1 must bind the complete candidate "
        "locus/gene family and its deterministic G4/G5 plus brain-PASS subset, with upstream "
        "failures preserved and full policy-bound genome-wide MAGMA gene results present."
    )
    common_ram = (
        "Safe sequentially by one eligible target trait and one complete atlas at a time; retain "
        "all seven hypotheses in each atlas, never split an atlas family, and apply BH over n=7 "
        "after failed rows are retained with p=1 for correction only."
    )
    return [
        {
            "analysis_id": "P14_BRAIN_CELL_CLASS_DISCOVERY", "phase": 14,
            "resource_status": "READY", "status": "NOT_APPLICABLE_UNTIL_UPSTREAM",
            "status_reason": (
                "The exact 75-column Allen Human MTG source, exhaustive 75-to-7 map, arithmetic-mean "
                "aggregator, MAGMA v1.10 model, and seven-class correction family are frozen; no "
                "eligible upstream family has been supplied to this result-blind contract."
            ),
            "required_upstream_gate": common_gate, "ram_decomposition": common_ram,
            "planned_output": "results/track_b/mechanism_followup/science/cell_class/{pair_id}/{trait_id}/enrichment/P14_ALLEN_BROAD_CLASS_RESULTS.tsv",
            "claim_limit": (
                "A corrected positive Allen MTG gene-property association is broad-class expression "
                "enrichment, not a causal-cell, subtype, locus-mediation, or pericyte claim."
            ),
        },
        {
            "analysis_id": "P16_HUMAN_CELL_CLASS_REPLICATION", "phase": 16,
            "resource_status": "READY", "status": "NOT_APPLICABLE_UNTIL_UPSTREAM",
            "status_reason": (
                "The independent human GSE67835 seven-column matrix and frozen alignment/classifier "
                "are available; hybrid is denominator-retained but unaligned, while excitatory and "
                "inhibitory Allen classes share one coarse neurons target without double credit."
            ),
            "required_upstream_gate": common_gate + " P14 and GSE67835 must both finish with their full seven-row families.",
            "ram_decomposition": common_ram,
            "planned_output": "results/track_b/mechanism_followup/science/cell_class/{pair_id}/{trait_id}/replication/P16_CLASS_REPLICATION.tsv",
            "claim_limit": (
                "Replication requires positive direction and BH FDR <=0.05 in both complete atlas "
                "families; shared neurons support is PARTIAL_COARSE_NEURON, hybrid is never aligned, "
                "and neither pericyte nor mouse evidence is implied."
            ),
        },
    ]


def render_report(policy_sha: str, mapping_sha: str, alignment_sha: str) -> bytes:
    text = f"""# Track B P14/P16 pre-result execution amendment (V4)

This immutable, result-blind amendment resolves only the software/resource blockers for `P14_BRAIN_CELL_CLASS_DISCOVERY` and `P16_HUMAN_CELL_CLASS_REPLICATION`. Revision 3 and all revision-1/2/3 files remain byte-exact; revision 3 remains authoritative for every other analysis.

## Frozen source and model

- Allen Human MTG level 2: human middle temporal gyrus, 29,155 gene rows and exactly 75 subtype columns. Every subtype maps exactly once to seven broad classes (2 astrocyte, 1 endothelial, 24 excitatory, 45 inhibitory, 1 microglial, 1 oligodendrocyte and 1 OPC column). Broad columns are arithmetic means computed in source order with `math.fsum`.
- The source `Average` column is preserved only as the MAGMA covariate. It is never an eighth hypothesis.
- GSE67835 Human Cortex without fetal samples: independent human cerebral-cortex matrix with exactly seven tested columns. `hybrid` remains tested and in the seven-test BH denominator but is explicitly unaligned.
- MAGMA v1.10 uses `condition-hide=Average direction=greater`; the recorded statistic is beta/SE and its one-sided p-value.

## Complete-family inference

BH is recomputed separately over all seven Allen hypotheses and all seven GSE67835 hypotheses for each target trait. Failed, missing, malformed, zero-SE and NA rows remain in the family; failed rows use p=1 only for BH. A literal p=0 remains valid. Replication requires positive effect direction and BH FDR <=0.05 in both atlases. GSE `neurons` is a single many-to-one target for Allen excitatory and inhibitory classes: either can reach `PARTIAL_COARSE_NEURON`, but neither receives independent replication credit and the shared target is counted at most once.

## Upstream and claim boundary

Execution is inert until a canonical upstream contract supplies the complete candidate locus/Ensembl-gene family and an exact deterministic G4/G5 plus brain-tissue-PASS subset, with all failures retained, for a source-verified European target trait. The gene-property model still consumes the complete genome-wide policy-bound MAGMA gene-result family. No future result was read to choose an atlas, class, threshold, direction or correction family.

This amendment creates no human-brain pericyte resource and does not unblock subtype replication. Allen endothelial cells cannot be relabelled pericytes, GSE `hybrid` cannot be silently mapped, and mouse matrices cannot substitute for independent human replication.

## RAM boundary

The executor is safely decomposable by target trait and then by complete atlas. It never splits an atlas family and always reconstructs all seven rows before BH. The largest frozen expression text is 30,505,558 bytes, so source transformation is streamed; actual scientific execution has not been run or benchmarked by this pre-result freeze.

Policy SHA-256: `{policy_sha}`  
Mapping SHA-256: `{mapping_sha}`  
Alignment SHA-256: `{alignment_sha}`
"""
    return text.encode("utf-8")


def construct_artifacts(root: Path) -> tuple[dict[str, bytes], dict[str, object]]:
    policy, policy_record = E.load_policy(root)
    records = E.validate_static_inputs(root, policy)
    validate_prior_lock(root, V1_LOCK_REL, V1_LOCK_SHA256, 1)
    validate_prior_lock(root, V2_LOCK_REL, V2_LOCK_SHA256, 2)
    E.validate_parent_v3_family(root)
    validate_fuma_manifest(root, policy)
    v3_scope = read_v3_scope_rows(root)
    if E.stable_file(root, ENRICHMENT_REL).sha256 != ENRICHMENT_SHA256:
        raise E.ContractError("V4 enrichment executor differs from its pin")
    if E.stable_file(root, CLASSIFIER_REL).sha256 != CLASSIFIER_SHA256:
        raise E.ContractError("V4 replication classifier differs from its pin")
    headers = E.archive_headers(root, policy)
    mapping = E.class_mapping_rows(headers["Allen_Human_MTG_level2"][1:-1])
    alignment = E.alignment_rows()
    mapping_payload = E.tsv_bytes(E.MAPPING_FIELDS, mapping)
    alignment_payload = E.tsv_bytes(E.ALIGNMENT_FIELDS, alignment)
    magma = root / records["MAGMA_BINARY"].path
    version = subprocess.run([magma, "--version"], cwd=root, capture_output=True, text=True, check=False, timeout=30)
    if version.returncode != 0 or (version.stdout + version.stderr).strip() != "MAGMA version: v1.10 (custom)":
        raise E.ContractError("MAGMA binary does not report the frozen v1.10 build")
    rows = readiness_rows()
    mapping_sha = hashlib.sha256(mapping_payload).hexdigest()
    alignment_sha = hashlib.sha256(alignment_payload).hexdigest()
    script_record = E.stable_file(root, SELF_REL)
    amendment = {
        "schema_version": E.LOCK_SCHEMA,
        "contract_revision": 4,
        "contract_kind": "PRE_RESULT_P14_P16_EXECUTION_AMENDMENT",
        "authoritative_scope": ["P14_BRAIN_CELL_CLASS_DISCOVERY", "P16_HUMAN_CELL_CLASS_REPLICATION"],
        "result_blind": True, "future_results_accessed": False,
        "parent_v3": {"path": E.PARENT_LOCK_REL, "sha256": E.PARENT_LOCK_SHA256, "retention": "BYTE_EXACT"},
        "retained_predecessors": [
            {"revision": 1, "path": V1_LOCK_REL, "sha256": V1_LOCK_SHA256},
            {"revision": 2, "path": V2_LOCK_REL, "sha256": V2_LOCK_SHA256},
            {"revision": 3, "path": E.PARENT_LOCK_REL, "sha256": E.PARENT_LOCK_SHA256},
        ],
        "policy": {"path": E.POLICY_REL, "bytes": policy_record.bytes, "sha256": policy_record.sha256},
        "builder": {"path": SELF_REL, "bytes": script_record.bytes, "sha256": script_record.sha256},
        "executors": policy["executors"],
        "source_evidence": [
            {"input_id": name, "path": record.path, "bytes": record.bytes, "sha256": record.sha256}
            for name, record in sorted(records.items())
        ],
        "mapping": {"path": E.MAPPING_REL, "rows": 75, "bytes": len(mapping_payload), "sha256": mapping_sha},
        "alignment": {"path": E.ALIGNMENT_REL, "rows": 7, "bytes": len(alignment_payload), "sha256": alignment_sha},
        "statistics": policy["statistics"], "upstream_gate": policy["upstream_gate"],
        "production_outputs": policy["production_outputs"], "readiness": rows,
        "previous_scoped_statuses": {
            name: {"resource_status": row["resource_status"], "status": row["status"]}
            for name, row in sorted(v3_scope.items())
        },
        "preserved_blockers": {
            "P14_HUMAN_BRAIN_PERICYTE_DISCOVERY": "BLOCKED_BY_DATA",
            "P16_HUMAN_CELL_SUBTYPE_REPLICATION": "BLOCKED_BY_DATA",
            "P15_BRAIN_CELL_SUBTYPE_DISCOVERY": "NOT_APPLICABLE_UNTIL_UPSTREAM",
        },
        "ram_decomposition": {
            "safe": True, "unit": "ONE_TARGET_TRAIT_THEN_ONE_COMPLETE_ATLAS",
            "maximum_source_text_bytes": 30505558,
            "never_split": "SEVEN_CLASS_ATLAS_HYPOTHESIS_FAMILY",
            "actual_science_benchmark_status": "NOT_RUN_PRE_RESULT",
        },
        "claim_limit": (
            "Resource readiness is not scientific-result readiness. No causal-cell, subtype, pericyte, "
            "locus-mediation, or mouse-as-human-replication claim is authorized."
        ),
    }
    amendment_payload = E.canonical_json(amendment)
    readiness_payload = E.tsv_bytes(READINESS_FIELDS, rows)
    report_payload = render_report(policy_record.sha256, mapping_sha, alignment_sha)
    receipt = {
        "schema_version": E.LOCK_SCHEMA,
        "record_kind": "SCOPED_IMMUTABLE_EXECUTION_AMENDMENT",
        "superseded_scope": ["P14_BRAIN_CELL_CLASS_DISCOVERY", "P16_HUMAN_CELL_CLASS_REPLICATION"],
        "parent_contract": {
            "lock_path": E.PARENT_LOCK_REL, "lock_sha256": E.PARENT_LOCK_SHA256,
            "retention": "PRESERVE_REVISION_3_BYTE_EXACT_AND_KEEP_AUTHORITATIVE_OUTSIDE_SUPERSEDED_SCOPE",
        },
        "replacement_contract": {
            "lock_path": E.CONTRACT_LOCK_REL, "policy_path": E.POLICY_REL,
            "policy_sha256": policy_record.sha256, "builder_path": SELF_REL,
            "builder_sha256": script_record.sha256,
        },
        "scientific_correction": (
            "Adds the predeclared Allen 75-to-7 aggregation and Allen-to-GSE67835 replication "
            "classifier; it does not change any V3 analysis outside the two broad-class rows."
        ),
        "result_blind": True, "future_results_accessed": False,
    }
    artifacts: dict[str, bytes] = {
        E.MAPPING_REL: mapping_payload,
        E.ALIGNMENT_REL: alignment_payload,
        E.AMENDMENT_REL: amendment_payload,
        E.READINESS_REL: readiness_payload,
        E.REPORT_REL: report_payload,
        E.SUPERSESSION_REL: E.canonical_json(receipt),
    }
    lock = {
        "schema_version": E.LOCK_SCHEMA, "contract_revision": 4,
        "contract_kind": "PRE_RESULT_P14_P16_EXECUTION_AMENDMENT_LOCK",
        "authoritative_for_execution": True, "result_blind": True,
        "future_results_accessed": False, "policy_path": E.POLICY_REL,
        "policy_sha256": policy_record.sha256, "builder_path": SELF_REL,
        "builder_sha256": script_record.sha256, "parent_v3_lock_path": E.PARENT_LOCK_REL,
        "parent_v3_lock_sha256": E.PARENT_LOCK_SHA256,
        "analysis_ids": ["P14_BRAIN_CELL_CLASS_DISCOVERY", "P16_HUMAN_CELL_CLASS_REPLICATION"],
        "mapping_row_count": 75, "alignment_row_count": 7,
        "outputs": {
            relative: {"bytes": len(payload), "sha256": hashlib.sha256(payload).hexdigest()}
            for relative, payload in sorted(artifacts.items())
        },
    }
    artifacts[E.CONTRACT_LOCK_REL] = E.canonical_json(lock)
    return artifacts, lock


def _same_payload(root: Path, relative: str, payload: bytes) -> bool:
    try:
        record = E.stable_file(root, relative)
    except FileNotFoundError:
        return False
    return record.bytes == len(payload) and record.sha256 == hashlib.sha256(payload).hexdigest()


def _validate_no_unknown_v4_files(root: Path, allowed: set[str]) -> None:
    directory = root / "results/track_b/mechanism_followup/v4"
    try:
        info = os.lstat(directory)
    except FileNotFoundError:
        return
    if stat.S_ISLNK(info.st_mode) or not stat.S_ISDIR(info.st_mode):
        raise E.ContractError("V4 output namespace is not a real directory")
    observed = {entry.name for entry in os.scandir(directory)}
    expected_names = {Path(relative).name for relative in allowed if "/v4/" in relative}
    if not observed.issubset(expected_names):
        raise E.ContractError(f"V4 output namespace contains unknown entries: {sorted(observed - expected_names)}")


def freeze_no_replace(
    root: Path,
    artifacts: Mapping[str, bytes],
    *,
    after_link: Callable[[str], None] | None = None,
) -> None:
    ordered = [
        E.MAPPING_REL, E.ALIGNMENT_REL, E.AMENDMENT_REL, E.READINESS_REL,
        E.REPORT_REL, E.SUPERSESSION_REL, E.CONTRACT_LOCK_REL,
    ]
    if set(artifacts) != set(ordered):
        raise E.ContractError("V4 artifact family differs from the exact canonical family")
    _validate_no_unknown_v4_files(root, set(ordered))
    for relative in ordered:
        try:
            os.lstat(root / relative)
        except FileNotFoundError:
            continue
        if not _same_payload(root, relative, artifacts[relative]):
            raise E.ContractError(f"no-replace V4 freeze refused a nonidentical target: {relative}")
    stage_parent = E._secure_mkdirs(root, "results/track_b/mechanism_followup")
    E._secure_mkdirs(root, "results/track_b/mechanism_followup/v4")
    stage = Path(tempfile.mkdtemp(prefix=".p14-p16-v4-freeze-", dir=stage_parent))
    staged: dict[str, Path] = {}
    try:
        for index, relative in enumerate(ordered):
            if _same_payload(root, relative, artifacts[relative]):
                continue
            path = stage / f"{index:02d}.stage"
            fd = os.open(path, os.O_WRONLY | os.O_CREAT | os.O_EXCL | getattr(os, "O_NOFOLLOW", 0), 0o600)
            try:
                view = memoryview(artifacts[relative])
                while view:
                    written = os.write(fd, view)
                    if written <= 0:
                        raise E.ContractError(f"short write while staging V4 artifact: {relative}")
                    view = view[written:]
                os.fsync(fd)
            finally:
                os.close(fd)
            os.chmod(path, 0o444)
            staged[relative] = path
        for relative in ordered:
            if not _same_payload(root, relative, artifacts[relative]):
                parent = (root / relative).parent
                E._secure_mkdirs(root, parent.relative_to(root).as_posix())
                try:
                    os.link(staged[relative], root / relative, follow_symlinks=False)
                except FileExistsError as exc:
                    if not _same_payload(root, relative, artifacts[relative]):
                        raise E.ContractError(f"V4 publication race produced a nonidentical target: {relative}") from exc
                if not _same_payload(root, relative, artifacts[relative]):
                    raise E.ContractError(f"published V4 artifact differs: {relative}")
                if after_link is not None:
                    after_link(relative)
        for parent in sorted({(root / relative).parent for relative in ordered}, key=str):
            fd = os.open(parent, os.O_RDONLY)
            try:
                os.fsync(fd)
            finally:
                os.close(fd)
    finally:
        for path in staged.values():
            try:
                os.unlink(path)
            except FileNotFoundError:
                pass
        try:
            os.rmdir(stage)
        except FileNotFoundError:
            pass


def verify_frozen(root: Path) -> None:
    lock_payload = E.read_stable_bytes(root, E.CONTRACT_LOCK_REL)
    lock = E.read_canonical_json(root, E.CONTRACT_LOCK_REL)
    if (
        set(lock) != {
            "schema_version", "contract_revision", "contract_kind", "authoritative_for_execution",
            "result_blind", "future_results_accessed", "policy_path", "policy_sha256",
            "builder_path", "builder_sha256", "parent_v3_lock_path", "parent_v3_lock_sha256",
            "analysis_ids", "mapping_row_count", "alignment_row_count", "outputs",
        }
        or lock.get("schema_version") != E.LOCK_SCHEMA or lock.get("contract_revision") != 4
        or lock.get("contract_kind") != "PRE_RESULT_P14_P16_EXECUTION_AMENDMENT_LOCK"
        or lock.get("authoritative_for_execution") is not True
        or lock.get("result_blind") is not True or lock.get("future_results_accessed") is not False
        or lock.get("parent_v3_lock_path") != E.PARENT_LOCK_REL
        or lock.get("parent_v3_lock_sha256") != E.PARENT_LOCK_SHA256
        or lock.get("analysis_ids") != ["P14_BRAIN_CELL_CLASS_DISCOVERY", "P16_HUMAN_CELL_CLASS_REPLICATION"]
        or lock.get("mapping_row_count") != 75 or lock.get("alignment_row_count") != 7
    ):
        raise E.ContractError("V4 lock schema/lineage/scope is invalid")
    outputs = lock.get("outputs")
    if not isinstance(outputs, dict) or set(outputs) != set(E.CONTRACT_OUTPUT_RELS):
        raise E.ContractError("V4 lock output family is incomplete")
    expected_paths = set(E.CONTRACT_OUTPUT_RELS) | {E.CONTRACT_LOCK_REL}
    _validate_no_unknown_v4_files(root, expected_paths)
    for relative in E.CONTRACT_OUTPUT_RELS:
        expected = outputs[relative]
        if not isinstance(expected, dict) or set(expected) != {"bytes", "sha256"}:
            raise E.ContractError("V4 lock output identity is malformed")
        record = E.stable_file(root, relative)
        if record.bytes != expected["bytes"] or record.sha256 != expected["sha256"]:
            raise E.ContractError(f"V4 frozen output differs: {relative}")
    recomputed, expected_lock = construct_artifacts(root)
    for relative in E.CONTRACT_OUTPUT_RELS:
        if E.read_stable_bytes(root, relative) != recomputed[relative]:
            raise E.ContractError(f"current immutable inputs no longer reproduce V4 output: {relative}")
    if lock != expected_lock or lock_payload != recomputed[E.CONTRACT_LOCK_REL]:
        raise E.ContractError("current immutable inputs no longer reproduce the V4 lock")


def parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument("--freeze", action="store_true", help="publish the V4 family without replacement")
    mode.add_argument("--verify", action="store_true", help="verify V4 and all retained predecessors without writes")
    parser.add_argument("--root", type=Path, default=ROOT_AT_IMPORT)
    return parser.parse_args(argv)


def main(argv: Sequence[str] | None = None) -> int:
    args = parse_args(argv)
    root = args.root.resolve()
    if not root.is_dir():
        raise E.ContractError("repository root is not a directory")
    if args.freeze:
        artifacts, _lock = construct_artifacts(root)
        freeze_no_replace(root, artifacts)
        print("FROZEN_P14_P16_EXECUTION_AMENDMENT_V4 analyses=2 mapping_rows=75 alignment_rows=7")
    else:
        verify_frozen(root)
        print("VERIFIED_P14_P16_EXECUTION_AMENDMENT_V4 analyses=2 mapping_rows=75 alignment_rows=7")
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except E.ContractError as exc:
        raise SystemExit(f"ERROR: {exc}")
