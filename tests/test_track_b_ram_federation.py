from __future__ import annotations

import hashlib
import importlib.util
import json
import os
import sys
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import patch

import pytest


ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "scripts/153_finalize_track_b_ram_aware_execution.py"
SPEC = importlib.util.spec_from_file_location("track_b_ram_federation_153_test", SCRIPT)
assert SPEC is not None and SPEC.loader is not None
RAM = importlib.util.module_from_spec(SPEC)
PREVIOUS = sys.dont_write_bytecode
sys.dont_write_bytecode = True
try:
    SPEC.loader.exec_module(RAM)
finally:
    sys.dont_write_bytecode = PREVIOUS


def identity(root: Path, path: Path) -> dict[str, object]:
    content = path.read_bytes()
    return {
        "path": str(path.relative_to(root)), "bytes": len(content),
        "sha256": hashlib.sha256(content).hexdigest(),
    }


def write_json(path: Path, value: object) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, indent=2, sort_keys=True) + "\n", encoding="utf-8")


def row(
    *, analysis: str = "TEST", pair: str = "A", locus: str = "L1",
    chromosome: str = "1", n_snps: str = "2", peak: str = "0.5",
    runtime: str = "1.25", status: str = "0", output_hash: str = "a" * 64,
) -> dict[str, str]:
    return {
        "analysis": analysis, "pair": pair, "locus": locus,
        "chromosome": chromosome, "n_snps": n_snps,
        "peak_ram_gb": peak, "runtime_sec": runtime,
        "exit_status": status, "output_hash": output_hash,
    }


def stale_audit(*, observed: int = 2, declared: int = 1) -> dict[str, object]:
    return {
        "fully_consistent": False, "row_count_match": False,
        "benchmark_bytes_match": False, "benchmark_sha256_match": False,
        "observed_benchmark_rows": observed, "declared_attempt_rows": declared,
        "handling": (
            "PRESERVED_AS_PROVISIONAL_HISTORY_ONLY;"
            "FINAL_FEDERATION_MUST_REVALIDATE_RAW_RECEIPTS"
        ),
    }


def make_v1_fixture(root: Path) -> tuple[Path, dict[str, object], dict[str, object], Path]:
    fingerprint = RAM.V1_FINGERPRINT
    attempt_rel = (
        Path("results/track_b/checkpoints/lava") / fingerprint
        / "discovery/.attempts/locus_0001.test"
    )
    attempt = root / attempt_rel
    attempt.mkdir(parents=True)
    result = attempt / "result.rds"
    semantic = attempt / "semantic_validation.txt"
    log = attempt / "worker.log"
    result.write_bytes(b"sealed-result")
    semantic.write_text("VALID\n", encoding="utf-8")
    log.write_text("worker\n", encoding="utf-8")
    receipt = {
        "schema_version": 3, "analysis": "LAVA_DISCOVERY", "phase": "discovery",
        "locus_index": 1, "execution_fingerprint": fingerprint, "exit_status": 0,
        "peak_rss_bytes": 1024**3, "runtime_sec": 1.25,
        "marker": {
            "index": "1", "locus": "L1", "chromosome": "1", "qc": "PROCESSED",
            "pair": "NONE", "n_snps": "NA",
        },
        "artifacts": [
            {"path": "result.rds", "bytes": result.stat().st_size,
             "sha256": hashlib.sha256(result.read_bytes()).hexdigest()},
            {"path": "semantic_validation.txt", "bytes": semantic.stat().st_size,
             "sha256": hashlib.sha256(semantic.read_bytes()).hexdigest()},
        ],
        "log_bytes": log.stat().st_size,
        "log_sha256": hashlib.sha256(log.read_bytes()).hexdigest(),
    }
    receipt_path = attempt / "receipt.json"
    write_json(receipt_path, receipt)
    ready_path = attempt / "READY"
    write_json(ready_path, {
        "schema_version": 1, "state": "READY",
        "receipt_sha256": hashlib.sha256(receipt_path.read_bytes()).hexdigest(),
    })
    bundle_rel = attempt_rel.parent.parent / "locus_0001"
    bundle = root / bundle_rel
    os.symlink(os.path.relpath(attempt, bundle.parent), bundle, target_is_directory=True)
    entry = {
        "locus_index": 1, "locus": "L1", "chromosome": "1", "qc": "PROCESSED",
        "bundle_path": str(bundle_rel),
        "attempt_path": str(attempt_rel),
        "ready": identity(root, ready_path), "receipt": identity(root, receipt_path),
        "result": identity(root, result), "log": identity(root, log),
        "semantic_validation": identity(root, semantic),
    }
    source_lock = {
        "source_discovery_fingerprint": fingerprint, "checkpoint_count": 1,
        "checkpoints": [entry], "checkpoint_family_sha256": RAM.digest_json([entry]),
    }
    failed_root = (
        root / "results/track_b/checkpoints/lava" / fingerprint / "failed_attempts"
    )
    failed_root.mkdir(parents=True)
    failed_log = failed_root / "aggregate-discovery_all_1.log"
    failed_log.write_text(
        f"command --phase aggregate-discovery --fingerprint {fingerprint}\nfailed\n",
        encoding="utf-8",
    )
    discovery = row(
        analysis="LAVA_DISCOVERY", pair="NONE", locus="L1", chromosome="1",
        n_snps="NA", peak="1.000000", runtime="1.250", status="0",
        output_hash=hashlib.sha256(result.read_bytes()).hexdigest(),
    )
    failed = row(
        analysis="LAVA_DISCOVERY_AGGREGATION", pair="UNKNOWN", locus="ALL",
        chromosome="ALL", n_snps="NA", peak="0.25", runtime="2",
        status="1", output_hash=hashlib.sha256(failed_log.read_bytes()).hexdigest(),
    )
    benchmark = root / "history/RAM_BENCHMARK.tsv"
    benchmark.parent.mkdir(parents=True)
    benchmark.write_bytes(RAM.tsv_bytes(RAM.BENCHMARK_FIELDS, [discovery, failed]))
    history_receipt = {"historical_provenance_audit": stale_audit()}
    return benchmark, history_receipt, source_lock, result


def make_publication_candidate(root: Path) -> dict[str, object]:
    for role, relative in RAM.TOP_LEVEL.items():
        path = root / relative
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_bytes(f"old-{role}\n".encode())
    contents = {
        "benchmark": RAM.tsv_bytes(RAM.BENCHMARK_FIELDS, [row()]),
        "namespace": b'{"new":"namespace"}\n',
        "report": b"# final report\n",
        "provenance": b'{"new":"commit-marker"}\n',
    }
    return {
        "fingerprint": "f" * 64, "rows": [row()], "contents": contents,
        "old_sources": {
            role: identity(root, root / relative) for role, relative in RAM.TOP_LEVEL.items()
        },
        "package_rel": RAM.PACKAGE_ROOT_REL / ("f" * 64),
    }


def file_snapshot(root: Path) -> dict[str, tuple[int, int, str]]:
    return {
        str(path.relative_to(root)): (
            path.stat().st_ino, path.stat().st_size,
            hashlib.sha256(path.read_bytes()).hexdigest(),
        )
        for path in sorted(root.rglob("*")) if path.is_file()
    }


def test_known_stale_v1_provenance_is_quarantined_not_used_as_row_truth() -> None:
    RAM.validate_stale_v1_audit(
        stale_audit(observed=2496, declared=1376), 2496,
        expected_discovery_count=2495,
    )
    RAM.validate_stale_v1_audit(stale_audit(), 2, expected_discovery_count=1)
    consistent = {**stale_audit(), "fully_consistent": True}
    with pytest.raises(RAM.FederationError, match="quarantined"):
        RAM.validate_stale_v1_audit(consistent, 2, expected_discovery_count=1)
    with pytest.raises(RAM.FederationError, match="discovery plus one"):
        RAM.validate_stale_v1_audit(stale_audit(observed=3), 3, expected_discovery_count=1)


def test_v1_rows_are_reconstructed_from_exact_raw_ready_receipt_and_result(tmp_path: Path) -> None:
    benchmark, history, lock, result = make_v1_fixture(tmp_path)
    verified = RAM.validate_v1_snapshot_rows(
        tmp_path, benchmark, history, lock, expected_count=1,
    )
    assert verified["actual_rows"] == 2
    assert verified["stale_declared_rows"] == 1
    assert verified["rows"][0]["n_snps"] == "NA"
    result.write_bytes(b"tamper-result")
    with pytest.raises(RAM.FederationError, match="immutable identity"):
        RAM.validate_v1_snapshot_rows(tmp_path, benchmark, history, lock, expected_count=1)


def test_v1_raw_receipt_mismatch_is_refused(tmp_path: Path) -> None:
    benchmark, history, lock, _ = make_v1_fixture(tmp_path)
    receipt = tmp_path / lock["checkpoints"][0]["receipt"]["path"]
    payload = json.loads(receipt.read_text(encoding="utf-8"))
    payload["runtime_sec"] = 99
    write_json(receipt, payload)
    with pytest.raises(RAM.FederationError, match="immutable identity"):
        RAM.validate_v1_snapshot_rows(tmp_path, benchmark, history, lock, expected_count=1)


def test_superseded_continuations_are_never_admitted() -> None:
    with pytest.raises(RAM.FederationError, match="superseded continuation"):
        RAM.merge_component_rows([{
            "name": "old-continuation", "execution_fingerprint": RAM.SUPERSEDED_CONTINUATIONS[0],
            "rows": [row()],
        }])


def test_incomplete_c139_family_is_refused() -> None:
    conditional = lambda index: {"receipt": {"locus_index": index}}
    with pytest.raises(RAM.FederationError, match="conditional=1/2"):
        RAM.require_complete_continuation_family(
            {}, [conditional(1)], {}, None, expected_count=2,
        )
    with pytest.raises(RAM.FederationError, match="exactly one terminal"):
        RAM.require_complete_continuation_family(
            {}, [conditional(1), conditional(2)], {}, {}, expected_count=2,
        )
    assert RAM.require_complete_continuation_family(
        {}, [conditional(1), conditional(2)], None, {}, expected_count=2,
    ) == "TERMINAL_FAILED_QC"


def test_placo_and_finemapping_incomplete_families_are_refused() -> None:
    with pytest.raises(RAM.FederationError, match="PLACO full-family cleanup is incomplete"):
        RAM.require_placo_terminal_family({"A": {}, "B": {}})
    with pytest.raises(RAM.FederationError, match="fine-mapping family is incomplete"):
        RAM.require_finemap_terminal_state(
            {"state": "MATERIALIZED_OR_PARTIAL_RUN_FAMILY"},
            "NOT_APPLICABLE_ZERO_ELIGIBLE_LOCUS_FAMILY",
        )
    assert RAM.require_finemap_terminal_state(
        {"state": "COMPLETE_CANONICAL_FAMILY"},
        "NOT_APPLICABLE_ZERO_ELIGIBLE_LOCUS_FAMILY",
    ) == "COMPLETE_CANONICAL_FAMILY"


def test_placo_sealed_contract_is_reconstructed_without_legacy_deep_gate(tmp_path: Path) -> None:
    reference = tmp_path / "reference.provenance.json"
    write_json(reference, {"sealed": True})
    reference_record = identity(tmp_path, reference)
    gate_path = Path("terminal_gate.lock.json")
    readiness_path = Path("readiness_gate.tsv")
    (tmp_path / readiness_path).write_text("state\nREADY\n", encoding="utf-8")
    write_json(tmp_path / gate_path, {
        "analysis_id": "track-b", "reserve": 1,
        "readiness_gate": str(readiness_path),
        "readiness_gate_sha256": hashlib.sha256((tmp_path / readiness_path).read_bytes()).hexdigest(),
    })
    guard = {
        "accepted_states": ["PRESENT", "EVICTED"],
        "guard_implementation": {"path": "guard.py", "bytes": 1, "sha256": "a" * 64},
        "reference_provenance": reference_record,
        "archive_family_sha256": "b" * 64, "extracted_payload_count": 44,
        "extracted_family_sha256": "c" * 64, "verification": "SEALED",
    }

    class Gate:
        TERMINAL_GATE_LOCK = gate_path
        READINESS_GATE = readiness_path

    class Engine:
        PAIR_IDS = RAM.EXPECTED_PLACO_PAIR_ORDER

    class Module:
        EXECUTION_CONTRACT = Path("sequential_execution.contract.json")
        GATE = Gate
        ENGINE = Engine
        PAIR_ORDER = RAM.EXPECTED_PLACO_PAIR_ORDER

        @staticmethod
        def read_json(path: Path):
            return json.loads(path.read_text(encoding="utf-8"))

        @staticmethod
        def execution_contract_payload(gate, reference_state):
            reconstructed_guard = {
                **guard,
                "reference_provenance": reference_state["reference_provenance_identity"],
                "archive_family_sha256": reference_state["archive_family_sha256"],
                "extracted_payload_count": reference_state["extracted_payload_count"],
                "extracted_family_sha256": reference_state["extracted_family_sha256"],
            }
            return {
                "gate": gate, "pair_order": list(Module.PAIR_ORDER),
                "lava_reference_archive_state_guard": reconstructed_guard,
            }

        @staticmethod
        def verify_reference_archive_state(_gate):
            return {
                "reference_provenance_identity": reference_record,
                "archive_family_sha256": guard["archive_family_sha256"],
                "extracted_payload_count": guard["extracted_payload_count"],
                "extracted_family_sha256": guard["extracted_family_sha256"],
            }

    contract = Module.execution_contract_payload(
        Module.read_json(tmp_path / gate_path), {
            "reference_provenance_identity": reference_record,
            "archive_family_sha256": guard["archive_family_sha256"],
            "extracted_payload_count": guard["extracted_payload_count"],
            "extracted_family_sha256": guard["extracted_family_sha256"],
        },
    )
    write_json(tmp_path / Module.EXECUTION_CONTRACT, contract)
    assert RAM.validate_placo_execution_contract(Module, tmp_path) == contract
    Module.PAIR_ORDER = ("A", "B")
    Engine.PAIR_IDS = ("A", "B")
    with pytest.raises(RAM.FederationError, match="module pair order"):
        RAM.validate_placo_execution_contract(Module, tmp_path)
    Module.PAIR_ORDER = RAM.EXPECTED_PLACO_PAIR_ORDER
    Engine.PAIR_IDS = RAM.EXPECTED_PLACO_PAIR_ORDER
    contract["pair_order"] = ["A", "B"]
    write_json(tmp_path / Module.EXECUTION_CONTRACT, contract)
    with pytest.raises(RAM.FederationError, match="sealed contract pair order"):
        RAM.validate_placo_execution_contract(Module, tmp_path)
    contract["pair_order"] = list(RAM.EXPECTED_PLACO_PAIR_ORDER)
    contract["gate"] = {"analysis_id": "tampered"}
    write_json(tmp_path / Module.EXECUTION_CONTRACT, contract)
    with pytest.raises(RAM.FederationError, match="missing or drifted"):
        RAM.validate_placo_execution_contract(Module, tmp_path)


def test_placo_pair_context_uses_verified_gate_without_recursive_gate_walk(tmp_path: Path) -> None:
    paths = {
        name: Path(name) for name in (
            "policy.json", "v1-contract.json", "v1-input.json", "terminal-gate.json",
            "engine.py", "bridge.py", "runner-v1.R", "runner-v2.R",
        )
    }
    for name, relative in paths.items():
        target = tmp_path / relative
        write_json(target, {"artifact": name})

    gate_lock = {
        "pair_scope": ["A", "B", "CONTROL"],
        "resource_envelope": {"ram_safety_reserve_bytes": 1024**3},
        "lava_scientific_evidence": {
            "status": "FAILED_QC_NOT_CONSUMED", "scientific_validation_passed": False,
            "passed_lava_results_used_by_placo": False,
            "consumption": "ORDERING_TERMINALITY_ATTESTATION_ONLY",
        },
    }

    class Engine:
        PAIR_IDS = ("A", "B", "CONTROL")
        POLICY = paths["policy.json"]
        CONTRACT_LOCK = paths["v1-contract.json"]
        INPUT_LOCK = Path("original-input.json")
        RUNNER = Path("original-runner.R")

        @staticmethod
        def assert_upstream_ready(_row):
            raise AssertionError("original readiness verifier must be temporarily replaced")

        @staticmethod
        def sha256(path: Path) -> str:
            return hashlib.sha256(Path(path).read_bytes()).hexdigest()

        @classmethod
        def production_context(cls, root: Path, pair: str):
            assert root == tmp_path
            assert pair == "A"
            assert cls.INPUT_LOCK == paths["v1-input.json"]
            assert cls.RUNNER == paths["runner-v2.R"]
            return {
                "policy": {
                    "analysis_id": "track-b", "placo_plus": {"source_sha256": "a" * 64},
                },
                "pair_id": pair, "trait1": "t1", "trait2": "t2",
                "family_role": "PRIMARY", "sources": [
                    {"path": "dense.tsv.gz", "sha256": "b" * 64, "absolute_path": "/tmp/dense"},
                ],
            }

    original_assert = Engine.assert_upstream_ready

    class Gate:
        V1_INPUT_LOCK = paths["v1-input.json"]
        TERMINAL_GATE_LOCK = paths["terminal-gate.json"]

        @staticmethod
        def verify_gate():
            raise AssertionError("legacy recursive gate verifier must not be called")

    class Bridge:
        __file__ = str(tmp_path / paths["bridge.py"])
        ENGINE_PATH = tmp_path / paths["engine.py"]
        V1_RUNNER = paths["runner-v1.R"]
        V2_RUNNER = paths["runner-v2.R"]
        inherited_v1_upstream = staticmethod(lambda _row: None)

        @staticmethod
        def validate_v2_paths(context):
            assert context["pair_id"] == "A"
            assert context["input_lock"] == gate_lock

    class Module:
        ENGINE = Engine
        GATE = Gate
        BRIDGE = Bridge

    context, observed_gate = RAM.derive_placo_context_from_verified_contract(
        Module, tmp_path, "A", gate_lock,
    )
    assert observed_gate == gate_lock
    assert context["fingerprint"] == RAM.digest_json(context["fingerprint_payload"])
    assert context["fingerprint_payload"]["sources"] == [
        {"path": "dense.tsv.gz", "sha256": "b" * 64},
    ]
    assert Engine.assert_upstream_ready == original_assert
    assert Engine.INPUT_LOCK == Gate.TERMINAL_GATE_LOCK
    assert Engine.RUNNER == Bridge.V2_RUNNER


def test_finemap_import_uses_only_public_terminal_and_component_verifiers(tmp_path: Path) -> None:
    verifier = tmp_path / RAM.FINEMAP_REL
    verifier.parent.mkdir(parents=True)
    verifier.write_bytes((ROOT / RAM.FINEMAP_REL).read_bytes())
    report_identity = {
        "path": "results/track_b/finemapping/RAM_BY_LOCUS." + "a" * 20 + ".tsv",
        "bytes": 1, "sha256": "a" * 64,
    }
    namespace_identity = {
        "path": "results/track_b/finemapping/RAM_BY_LOCUS." + "a" * 20 + ".namespace.json",
        "bytes": 1, "sha256": "b" * 64,
    }
    provenance_identity = {
        "path": "results/track_b/finemapping/RAM_BY_LOCUS." + "a" * 20 + ".provenance.json",
        "bytes": 1, "sha256": "c" * 64,
    }
    hint = {
        "report_identity": report_identity, "namespace_identity": namespace_identity,
        "provenance_identity": provenance_identity,
    }
    complete = row(
        analysis="TRACK_B_FINEMAPPING_TRAIT_COLOC_V1", status="COMPLETE_SEALED",
    )
    blocked = row(
        analysis="TRACK_B_FINEMAPPING_TRAIT_COLOC_V1", pair="B", locus="UNMAPPED_1",
        n_snps="0", peak="0", runtime="0",
        status="BLOCKED_BY_DATA:NO_PREDECLARED_SIGNED_LD_BLOCK", output_hash="d" * 64,
    )

    class PublicOnlyModule:
        @staticmethod
        def verify_current_state(_root):
            return {
                "state": "COMPLETE_CANONICAL_FAMILY",
                    "canonical": {
                        "state": "COMPLETE_CANONICAL_FAMILY", "component_report": hint,
                        "runs": [{
                            "row": {"locus_entry_id": "L1", "pair_id": "A", "CHR": "1"},
                            "resource_row": {"peak_rss_bytes": RAM.GIB // 2},
                            "materialization": {"receipt": {
                                "peak_rss_bytes": RAM.GIB // 4,
                                "larger_host_continuation": True,
                            }},
                            "provenance": {
                                "larger_host_continuation": True,
                                "engine_continuation_predecessors": [{"terminal_status": "FAILED_RESOURCE_OOM"}],
                            },
                        }],
                        "provenance_identity": {
                        "path": "results/track_b/finemapping/results.provenance.json",
                        "bytes": 1, "sha256": "e" * 64,
                    },
                },
            }

        @staticmethod
        def verify_component_ram_report(_root, report_rel):
            assert str(report_rel) == report_identity["path"]
            return {
                "state": "VERIFIED_COMPONENT_RAM_REPORT", "rows": [complete, blocked],
                "collision_key": RAM.COLLISION_KEY,
                "report_identity": report_identity, "namespace_identity": namespace_identity,
                "provenance_identity": provenance_identity,
                "paths": {
                    "report": report_identity["path"], "namespace": namespace_identity["path"],
                    "provenance": provenance_identity["path"],
                },
            }

    with patch.object(RAM, "_load_module", return_value=PublicOnlyModule):
        observed = RAM.validate_finemap_component(tmp_path)
    assert observed["rows"] == [complete, blocked]
    assert observed["peak_rss_bytes"] == [RAM.GIB // 2, 0]
    assert observed["evidence"]["blocked_by_data_rows"] == 1
    assert [
        stage["stage"] for stage in observed["evidence"]["genuine_larger_host_stages"]
    ] == ["MATERIALIZATION", "R_ENGINE"]
    assert all(
        stage["locus"] == "L1"
        for stage in observed["evidence"]["genuine_larger_host_stages"]
    )
    assert observed["evidence"]["component_verified_only_via_script_151_public_api"] is True


def test_finemap_zero_state_uses_public_terminal_api_without_component_placeholder(
    tmp_path: Path,
) -> None:
    verifier = tmp_path / RAM.FINEMAP_REL
    verifier.parent.mkdir(parents=True)
    verifier.write_bytes((ROOT / RAM.FINEMAP_REL).read_bytes())

    class ZeroModule:
        @staticmethod
        def verify_current_state(_root):
            return {"state": RAM.FINEMAP_ZERO_STATE, "family": {}, "canonical": {}}

        @staticmethod
        def verify_component_ram_report(*_args, **_kwargs):
            raise AssertionError("zero-state API must not require a component placeholder")

    with patch.object(RAM, "_load_module", return_value=ZeroModule):
        observed = RAM.validate_finemap_component(tmp_path)
    assert observed["rows"] == []
    assert observed["peak_rss_bytes"] == []
    assert observed["execution_fingerprint"] == RAM.FINEMAP_ZERO_STATE


def test_exact_schema_collision_key_and_failure_retention() -> None:
    failed = row(status="17:PROCESS_FAILED", peak="7.25", output_hash="b" * 64)
    merged = RAM.merge_component_rows([{
        "name": "component", "execution_fingerprint": "c" * 64, "rows": [failed],
    }])
    assert merged == [failed]
    rendered = RAM.tsv_bytes(RAM.BENCHMARK_FIELDS, merged)
    assert rendered.decode().splitlines()[0].split("\t") == RAM.BENCHMARK_FIELDS
    assert RAM.read_tsv_bytes(rendered, RAM.BENCHMARK_FIELDS, "candidate") == [failed]
    with pytest.raises(RAM.FederationError, match="exact required schema"):
        RAM.read_tsv_bytes(rendered.replace(b"output_hash", b"wrong_hash"), RAM.BENCHMARK_FIELDS, "bad")
    with pytest.raises(RAM.FederationError, match="duplicate federation collision key"):
        RAM.merge_component_rows([{
            "name": "duplicates", "execution_fingerprint": "c" * 64,
            "rows": [failed, {**failed, "output_hash": "d" * 64}],
        }])
    report = RAM.build_report("e" * 64, merged, [{
        "name": "component", "rows": merged, "execution_fingerprint": "c" * 64,
        "peak_rss_bytes": [int(7.25 * RAM.GIB)], "evidence": {},
    }]).decode()
    assert "17:PROCESS_FAILED" in report
    assert "PLACO global nuisance" in report
    assert "1-GiB reserve" in report
    assert "genuine larger-server" in report.lower()
    assert "| TEST |" in report and "| NO |" in report


def test_eight_gib_fit_uses_exact_bytes_not_rounded_display() -> None:
    display_boundary = row(peak="7.000000")
    exact = RAM._analysis_statistics([display_boundary], [7 * RAM.GIB])[0]
    one_byte_over = RAM._analysis_statistics([display_boundary], [7 * RAM.GIB + 1])[0]
    unknown = RAM._analysis_statistics([display_boundary], [None])[0]
    assert exact["fits"] is True
    assert one_byte_over["fits"] is False
    assert unknown["fits"] is None


def test_conditional_n_snps_is_inherited_from_exact_source_locus() -> None:
    bundle = {
        "receipt": {
            "phase": "conditional", "locus_index": 1, "analysis": "LAVA_CONDITIONAL_V2",
            "peak_rss_bytes": 1024, "runtime_sec": 0.5, "exit_status": 0,
            "marker": {
                "locus": "L1", "chromosome": "1", "pair": "NONE",
                "n_snps": "NA", "qc": "NOT_BIVARIATE_FDR_ELIGIBLE",
            },
        },
        "result_identity": {"sha256": "a" * 64},
    }
    observed = RAM._continuation_federation_row(
        bundle, {1: {"locus": "L1", "chromosome": "1", "pair": "NONE", "n_snps": "NA"}},
    )
    assert observed["n_snps"] == "NA"
    assert observed["exit_status"] == "0:NOT_BIVARIATE_FDR_ELIGIBLE"


def test_c139_unbenchmarkable_failed_attempt_is_explicitly_fail_closed(tmp_path: Path) -> None:
    failure_root = tmp_path / "results/failed_attempts"
    failure_root.mkdir(parents=True)
    RAM._require_no_unbenchmarkable_failed_attempts(tmp_path, failure_root)
    (failure_root / "conditional_1_123.log").write_text("failed\n", encoding="utf-8")
    with pytest.raises(RAM.FederationError, match="without benchmarkable immutable resource receipts"):
        RAM._require_no_unbenchmarkable_failed_attempts(tmp_path, failure_root)


def test_c139_invariant_hashes_seed_the_one_pass_active_input_cache(tmp_path: Path) -> None:
    paths = {
        "lineage": tmp_path / "lineage.json",
        "source": tmp_path / "source.lock.json",
        "semantic": tmp_path / "checkpoint_validator.py",
        "terminal": tmp_path / "terminal_validator.py",
    }
    for name, path in paths.items():
        path.write_bytes(f"{name}\n".encode())

    contract = SimpleNamespace(
        lineage_path=lambda _fingerprint: paths["lineage"],
        source_lock_path=lambda _fingerprint: paths["source"],
    )

    class Engine:
        CHECKPOINT_VALIDATOR = paths["semantic"]
        RESULT_VALIDATOR = paths["terminal"]

        @staticmethod
        def _safe_relative(value: str) -> Path:
            return tmp_path / value

        @staticmethod
        def stable_stat(path: Path) -> dict[str, int]:
            status = path.stat()
            return {"bytes": status.st_size, "mtime_ns": status.st_mtime_ns}

    cache: dict[str, tuple[int, str]] = {}
    with patch.object(RAM, "stable_bytes", wraps=RAM.stable_bytes) as reader:
        observed = RAM._precompute_continuation_invariants(tmp_path, Engine, contract, cache)
        assert reader.call_count == 4
        relative = str(paths["lineage"].relative_to(tmp_path))
        current = Engine.stable_stat(paths["lineage"])
        item = {
            "path": relative, "bytes": paths["lineage"].stat().st_size,
            "sha256": observed["lineage_sha256"], "role": "CONTINUATION_LINEAGE",
            "pre_stat": current, "post_stat": current,
        }
        RAM._active_content_identity(Engine, item, cache, "cached lineage")
        assert reader.call_count == 4


def test_c139_official_loci_are_loaded_once_and_factory_is_restored() -> None:
    calls = {"factory": 0, "load": 0}
    loci = list(range(RAM.EXPECTED_DISCOVERY_COUNT))

    class Legacy:
        sentinel = "delegated"

        @staticmethod
        def load_loci():
            calls["load"] += 1
            return loci

    def factory():
        calls["factory"] += 1
        return Legacy()

    contract = SimpleNamespace(legacy_supervisor=factory)
    engine = SimpleNamespace(CONTRACT=contract)
    with RAM._cache_continuation_legacy_loci(engine) as cached_loci:
        assert cached_loci is loci
        assert contract.legacy_supervisor().load_loci() is loci
        assert contract.legacy_supervisor().sentinel == "delegated"
    assert calls == {"factory": 1, "load": 1}
    assert contract.legacy_supervisor is factory


def test_preflight_and_verify_are_side_effect_free(tmp_path: Path) -> None:
    candidate = make_publication_candidate(tmp_path)
    lazy = tmp_path / "lazy_nested_verifier.py"
    lazy.write_text("VALUE = 1\n", encoding="utf-8")

    def collector(_root: Path):
        specification = importlib.util.spec_from_file_location(
            f"lazy_nested_verifier_{id(object())}", lazy,
        )
        assert specification is not None and specification.loader is not None
        module = importlib.util.module_from_spec(specification)
        specification.loader.exec_module(module)
        assert module.VALUE == 1
        return candidate

    before = file_snapshot(tmp_path)
    result = RAM.preflight(tmp_path, collector=collector)
    assert result["publication_state"] == "UNPUBLISHED_FROZEN_V1"
    assert result["writes_performed"] is False
    assert file_snapshot(tmp_path) == before
    assert not (tmp_path / "__pycache__").exists()

    RAM.publish_candidate(tmp_path, candidate)
    committed = file_snapshot(tmp_path)
    verified = RAM.verify(tmp_path, collector=collector)
    assert verified["publication_state"] == "COMMITTED"
    assert verified["writes_performed"] is False
    assert file_snapshot(tmp_path) == committed


@pytest.mark.parametrize("crash_after", [1, 2, 3, 4])
def test_partial_publication_recovers_at_every_promotion_boundary(
    tmp_path: Path, crash_after: int,
) -> None:
    candidate = make_publication_candidate(tmp_path)
    unrelated = tmp_path / "results/track_b/user-owned.keep"
    unrelated.write_bytes(b"preserve-me\n")
    unrelated_identity = identity(tmp_path, unrelated)
    old_contents = {
        role: (tmp_path / relative).read_bytes() for role, relative in RAM.TOP_LEVEL.items()
    }
    original = RAM._promote_file_preserving
    first_attempt: list[str] = []

    def crash_after_boundary(root, candidate_value, package, role) -> None:
        original(root, candidate_value, package, role)
        first_attempt.append(role)
        if len(first_attempt) == crash_after:
            raise RuntimeError("simulated SIGKILL boundary")

    with patch.object(RAM, "_promote_file_preserving", side_effect=crash_after_boundary):
        with pytest.raises(RuntimeError, match="simulated"):
            RAM.publish_candidate(tmp_path, candidate)
    assert first_attempt == list(RAM.PROMOTION_ORDER[:crash_after])
    state = RAM.publication_state(tmp_path, candidate)
    assert state["promoted_prefix"] == crash_after
    assert identity(tmp_path, unrelated) == unrelated_identity

    recovered: list[str] = []

    def record_promotion(root, candidate_value, package, role) -> None:
        recovered.append(role)
        original(root, candidate_value, package, role)

    with patch.object(RAM, "_promote_file_preserving", side_effect=record_promotion):
        final = RAM.publish_candidate(tmp_path, candidate)
    assert recovered == list(RAM.PROMOTION_ORDER[crash_after:])
    assert final["publication"]["state"] == "COMMITTED"
    assert identity(tmp_path, unrelated) == unrelated_identity
    for role, relative in RAM.TOP_LEVEL.items():
        assert (tmp_path / relative).read_bytes() == candidate["contents"][role]
        assert (tmp_path / relative).stat().st_mode & 0o777 == 0o444
        backup = RAM._transaction_backup(tmp_path, candidate, role)
        assert backup.read_bytes() == old_contents[role]


def test_withdrawal_crash_is_visible_and_resumable(tmp_path: Path) -> None:
    candidate = make_publication_candidate(tmp_path)
    original_link = os.link
    crashed = False

    def crash_before_link(source, destination, *args, **kwargs):
        nonlocal crashed
        if not crashed and Path(destination) == tmp_path / RAM.TOP_LEVEL["benchmark"]:
            crashed = True
            raise RuntimeError("crash after preservation-first withdrawal")
        return original_link(source, destination, *args, **kwargs)

    with patch.object(os, "link", side_effect=crash_before_link):
        with pytest.raises(RuntimeError, match="preservation-first"):
            RAM.publish_candidate(tmp_path, candidate)
    state = RAM.publication_state(tmp_path, candidate)
    assert state["withdrawn_role"] == "benchmark"
    assert RAM._transaction_backup(tmp_path, candidate, "benchmark").is_file()
    assert RAM.publish_candidate(tmp_path, candidate)["publication"]["state"] == "COMMITTED"


@pytest.mark.parametrize("raced_role", RAM.PROMOTION_ORDER)
def test_hard_crash_after_withdrawing_raced_bytes_restores_without_overwrite(
    tmp_path: Path, raced_role: str,
) -> None:
    candidate = make_publication_candidate(tmp_path)
    target = tmp_path / RAM.TOP_LEVEL[raced_role]
    unknown = f"external-before-hard-crash-{raced_role}\n".encode()
    target.write_bytes(unknown)
    backup = RAM._transaction_backup(tmp_path, candidate, raced_role)
    backup.parent.mkdir(parents=True, exist_ok=True)
    RAM._rename_no_replace(target, backup)
    assert not target.exists() and backup.read_bytes() == unknown

    with pytest.raises(RAM.FederationError, match="restored without overwriting"):
        RAM.publish_candidate(tmp_path, candidate)
    assert target.read_bytes() == unknown
    assert not backup.exists()
    assert not (tmp_path / candidate["package_rel"]).exists()


@pytest.mark.parametrize("raced_role", RAM.PROMOTION_ORDER)
def test_each_promotion_preserves_concurrent_unknown_bytes(
    tmp_path: Path, raced_role: str,
) -> None:
    candidate = make_publication_candidate(tmp_path)
    unknown = f"external-{raced_role}\n".encode()
    original = RAM._rename_no_replace
    injected = False

    def inject_before_withdrawal(source: Path, destination: Path) -> None:
        nonlocal injected
        if not injected and source == tmp_path / RAM.TOP_LEVEL[raced_role]:
            replacement = source.with_name(f".{source.name}.external")
            replacement.write_bytes(unknown)
            os.replace(replacement, source)
            injected = True
        original(source, destination)

    with patch.object(RAM, "_rename_no_replace", side_effect=inject_before_withdrawal):
        with pytest.raises(RAM.FederationError, match="unexpected bytes were preserved"):
            RAM.publish_candidate(tmp_path, candidate)
    assert injected is True
    assert (tmp_path / RAM.TOP_LEVEL[raced_role]).read_bytes() == unknown


def test_nonprefix_partial_publication_and_unrelated_drift_fail_closed(tmp_path: Path) -> None:
    candidate = make_publication_candidate(tmp_path)
    (tmp_path / RAM.TOP_LEVEL["namespace"]).write_bytes(candidate["contents"]["namespace"])
    (tmp_path / RAM.TOP_LEVEL["namespace"]).chmod(0o444)
    with pytest.raises(RAM.FederationError, match="not a valid recoverable prefix"):
        RAM.publication_state(tmp_path, candidate)
    (tmp_path / RAM.TOP_LEVEL["namespace"]).chmod(0o644)
    (tmp_path / RAM.TOP_LEVEL["namespace"]).write_bytes(b"unrelated replacement\n")
    with pytest.raises(RAM.FederationError, match="neither frozen V1 nor candidate"):
        RAM.publication_state(tmp_path, candidate)
    with pytest.raises(RAM.FederationError, match="neither frozen V1 nor candidate"):
        RAM.publish_candidate(tmp_path, candidate)
    assert not (tmp_path / candidate["package_rel"]).exists()
