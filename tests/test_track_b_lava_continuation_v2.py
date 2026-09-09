import argparse
import hashlib
import importlib.util
import json
import subprocess
import tempfile
import unittest
from pathlib import Path
from unittest.mock import Mock, patch


ROOT = Path(__file__).resolve().parents[1]


def load_module(name: str, relative: str):
    specification = importlib.util.spec_from_file_location(name, ROOT / relative)
    if specification is None or specification.loader is None:
        raise RuntimeError(f"could not load {relative}")
    module = importlib.util.module_from_spec(specification)
    specification.loader.exec_module(module)
    return module


CONTRACT = load_module(
    "track_b_lava_continuation_contract_test",
    "scripts/134_track_b_lava_continuation_contract.py",
)
SUPERVISOR = load_module(
    "track_b_lava_continuation_supervisor_test",
    "scripts/136_run_track_b_lava_continuation.py",
)
RESULTS = load_module(
    "track_b_lava_continuation_results_test",
    "scripts/138_validate_track_b_lava_results_v2.py",
)


class TrackBLavaContinuationContractTests(unittest.TestCase):
    def test_legacy_fingerprint_inputs_remain_byte_exact(self):
        for relative, expected in CONTRACT.LEGACY_CODE_SHA256.items():
            observed = hashlib.sha256((ROOT / relative).read_bytes()).hexdigest()
            self.assertEqual(observed, expected, relative)

    def test_v2_runner_is_post_discovery_only_and_fixes_integer_collision(self):
        text = (ROOT / "scripts/135_run_track_b_lava_postdiscovery_v2.R").read_text(
            encoding="utf-8"
        )
        self.assertIn(
            'allowed_phases <- c("aggregate-discovery", "conditional", "finalize")',
            text,
        )
        self.assertNotIn('allowed_phases <- c("discovery"', text)
        self.assertIn('policy_int <- function(key)', text)
        self.assertIn('locus_index = integer(0)', text)
        self.assertNotIn('integer <- function(key)', text)

    def test_no_replace_publisher_refuses_even_identical_replacement(self):
        with tempfile.TemporaryDirectory() as temporary:
            destination = Path(temporary) / "lock.json"
            CONTRACT.publish_no_replace(destination, b"first\n")
            with self.assertRaises(CONTRACT.ContinuationError):
                CONTRACT.publish_no_replace(destination, b"first\n")
            self.assertEqual(destination.read_bytes(), b"first\n")

    def test_successor_contract_binds_exact_predecessor_final_lock(self):
        source = CONTRACT.PINNED_SOURCE_FINGERPRINT
        record = CONTRACT.validate_predecessor_supersession_lock()
        self.assertEqual(
            record,
            {
                "path": (
                    "results/track_b/lava_continuations/supersessions/final_locks/"
                    f"{CONTRACT.PREDECESSOR_CONTINUATION_FINGERPRINT}.lock.json"
                ),
                "bytes": CONTRACT.PREDECESSOR_SUPERSESSION_LOCK.stat().st_size,
                "sha256": "547aa0aa82fe6abce986fdc13b3172374dd2621cb1efbf6a8afca821c6e97930",
            },
        )
        payload = CONTRACT.continuation_contract_payload(source)
        self.assertEqual(
            payload["predecessor_supersession"],
            {
                "relationship": "SUCCESSOR_CONTRACT_BINDS_FINAL_PREDECESSOR_LOCK",
                "superseded_continuation_fingerprint": (
                    CONTRACT.PREDECESSOR_CONTINUATION_FINGERPRINT
                ),
                "final_lock": record,
            },
        )
        self.assertIn(
            "scripts/144_freeze_track_b_lava_supersession.py",
            payload["continuation_code"],
        )
        self.assertNotEqual(
            CONTRACT.continuation_fingerprint(source),
            CONTRACT.PREDECESSOR_CONTINUATION_FINGERPRINT,
        )

    def test_successor_contract_rejects_predecessor_lock_scope_tampering(self):
        lock_module = CONTRACT.supersession_validator()
        payload = lock_module.verify(CONTRACT.PREDECESSOR_CONTINUATION_FINGERPRINT)
        tampered = json.loads(json.dumps(payload))
        tampered["verified_snapshot"]["conditional_partial_family"]["ready_count"] = 21
        validator = Mock()
        validator.SupersessionError = lock_module.SupersessionError
        validator.DIRECT_LOCK_SCHEMA = lock_module.DIRECT_LOCK_SCHEMA
        validator.verify.return_value = tampered
        validator.lock_path.return_value = CONTRACT.PREDECESSOR_SUPERSESSION_LOCK
        with patch.object(
            CONTRACT, "supersession_validator", return_value=validator
        ):
            with self.assertRaisesRegex(
                CONTRACT.ContinuationError, "scope drifted"
            ):
                CONTRACT.validate_predecessor_supersession_lock()

    def test_source_lock_rejects_decoy_path_even_with_recorded_identity(self):
        canonical = {
            "locus_index": 1,
            "bundle_path": "results/canonical-link",
            "attempt_path": "results/canonical-attempt",
            "result": {
                "path": "results/canonical-attempt/result.rds",
                "bytes": 12,
                "sha256": "a" * 64,
            },
        }
        decoy = json.loads(json.dumps(canonical))
        decoy["result"]["path"] = "results/decoy/result.rds"
        with patch.object(CONTRACT, "_source_entry", return_value=canonical):
            with self.assertRaises(CONTRACT.ContinuationError):
                CONTRACT._require_canonical_source_entry(
                    decoy, 1, object(), "a" * 64, revalidate_semantics=False
                )
            CONTRACT._require_canonical_source_entry(
                canonical, 1, object(), "a" * 64, revalidate_semantics=False
            )

    def test_supervisor_forbids_discovery_and_requires_dual_marker_identity(self):
        source = "a" * 64
        continuation = "b" * 64
        with self.assertRaises(SUPERVISOR.ExecutionError):
            SUPERVISOR.bundle_name("discovery", 1)
        with tempfile.TemporaryDirectory() as temporary:
            result = Path(temporary) / "result.rds"
            result.write_bytes(b"not-an-rds-but-a-real-file")
            marker_line = "\t".join(
                [
                    "TRACK_B_LAVA_CONTINUATION_WORKER",
                    "phase=aggregate-discovery",
                    "index=0",
                    "locus=ALL",
                    "chromosome=ALL",
                    "n_snps=123",
                    "pair=A;B;CONTROL",
                    "qc=FULL_FAMILY_BH_COMPLETE",
                    "construction_complete=TRUE",
                    f"output={result}",
                    f"source_fingerprint={source}",
                    f"continuation_fingerprint={continuation}",
                ]
            )
            marker = SUPERVISOR.parse_marker(marker_line + "\n")
            SUPERVISOR.validate_marker(
                marker, "aggregate-discovery", None, result, source, continuation
            )
            marker["continuation_fingerprint"] = source
            with self.assertRaises(SUPERVISOR.ExecutionError):
                SUPERVISOR.validate_marker(
                    marker, "aggregate-discovery", None, result, source, continuation
                )

    def test_symlinked_attempt_namespace_is_rejected(self):
        fingerprint = "b" * 64
        with tempfile.TemporaryDirectory() as temporary:
            temporary_root = Path(temporary)
            runs = temporary_root / "runs"
            run = runs / fingerprint
            run.mkdir(parents=True)
            outside = temporary_root / "outside"
            outside.mkdir()
            (run / ".attempts").symlink_to(outside, target_is_directory=True)
            target = outside / "conditional_locus_0001.1"
            target.mkdir()
            bundle = run / "conditional" / "locus_0001"
            bundle.parent.mkdir()
            bundle.symlink_to(target, target_is_directory=True)
            with (
                patch.object(SUPERVISOR.CONTRACT, "ROOT", temporary_root),
                patch.object(SUPERVISOR.CONTRACT, "RUN_ROOT", runs),
            ):
                with self.assertRaises(SUPERVISOR.CONTRACT.ContinuationError):
                    SUPERVISOR._safe_attempt_target(bundle, fingerprint)

    def test_conditional_binds_chromosome_inputs_only_for_sealed_candidate(self):
        entries = [
            {
                "path": "common.tsv", "bytes": 1, "sha256": "a" * 64,
                "role": "COMMON_INPUT", "chromosome": None,
            }
        ]
        for chromosome in (3, 4):
            entries.extend(
                {
                    "path": f"chr{chromosome}/ref-{offset}", "bytes": 1,
                    "sha256": f"{chromosome}" * 64,
                    "role": "REFERENCE_CHROMOSOME", "chromosome": chromosome,
                }
                for offset in range(2)
            )
            entries.extend(
                {
                    "path": f"chr{chromosome}/shard-{offset}", "bytes": 1,
                    "sha256": f"{chromosome}" * 64,
                    "role": "CHROMOSOME_SHARD", "chromosome": chromosome,
                }
                for offset in range(8)
            )
            entries.append({
                "path": f"chr{chromosome}/input.tsv", "bytes": 1,
                "sha256": f"{chromosome}" * 64,
                "role": "CHROMOSOME_INPUT_INFO", "chromosome": chromosome,
            })
        legacy = Mock()
        legacy.validate_input_baseline.return_value = {"entries": entries}
        legacy.load_loci.return_value = [{"CHR": 3}]
        with patch.object(SUPERVISOR.CONTRACT, "legacy_supervisor", return_value=legacy):
            noncandidate = SUPERVISOR._scientific_input_specs(
                "conditional", 1, "a" * 64,
                include_conditional_chromosome=False,
            )
            observed = SUPERVISOR._scientific_input_specs(
                "conditional", 1, "a" * 64,
                include_conditional_chromosome=True,
            )
        self.assertEqual([item["role"] for item in noncandidate], ["COMMON_INPUT"])
        self.assertEqual(sum(item["role"] == "COMMON_INPUT" for item in observed), 1)
        self.assertEqual(
            sum(item["role"] == "REFERENCE_CHROMOSOME" for item in observed), 2
        )
        self.assertEqual(sum(item["role"] == "CHROMOSOME_SHARD" for item in observed), 8)
        self.assertEqual(
            sum(item["role"] == "CHROMOSOME_INPUT_INFO" for item in observed), 1
        )
        self.assertTrue(all("chr4/" not in str(item["path"]) for item in observed))

    def test_active_specs_uses_exact_sealed_candidate_membership(self):
        scientific = Mock(side_effect=lambda phase, index, source, **options: [{
            "path": "chromosome" if options["include_conditional_chromosome"] else "common",
            "bytes": 1, "sha256": "1" * 64,
            "role": "REFERENCE_CHROMOSOME" if options["include_conditional_chromosome"] else "COMMON_INPUT",
        }])
        aggregate_specs = [
            {"path": "aggregate.rds", "bytes": 1, "sha256": "2" * 64,
             "role": "CONTINUATION_CHECKPOINT"},
            {"path": "candidates.tsv", "bytes": 1, "sha256": "3" * 64,
             "role": "CONDITIONAL_CANDIDATE_ATTESTATION"},
            {"path": "aggregate-receipt.json", "bytes": 1, "sha256": "4" * 64,
             "role": "CONTINUATION_CHECKPOINT_RECEIPT"},
        ]
        artifact = lambda path, role, expected=None: {
            "path": str(path), "bytes": 1, "sha256": "5" * 64, "role": role,
        }
        source_spec = {
            "path": "source.rds", "bytes": 1, "sha256": "6" * 64,
            "role": "SOURCE_DISCOVERY_CHECKPOINT",
        }
        with (
            patch.object(SUPERVISOR, "_artifact_spec", side_effect=artifact),
            patch.object(SUPERVISOR, "_execution_code_specs", return_value=[]),
            patch.object(
                SUPERVISOR, "_aggregate_dependency_specs",
                return_value=(aggregate_specs, frozenset({2})),
            ),
            patch.object(SUPERVISOR, "_scientific_input_specs", scientific),
            patch.object(SUPERVISOR, "_source_discovery_spec", return_value=source_spec),
            patch.object(
                SUPERVISOR.CONTRACT, "lineage_path", return_value=Path("lineage.json")
            ),
            patch.object(
                SUPERVISOR.CONTRACT, "source_lock_path", return_value=Path("source-lock.json")
            ),
        ):
            locus_one = SUPERVISOR.active_specs("conditional", 1, "a" * 64, "b" * 64)
            locus_two = SUPERVISOR.active_specs("conditional", 2, "a" * 64, "b" * 64)
        self.assertFalse(any(item["role"] == "REFERENCE_CHROMOSOME" for item in locus_one))
        self.assertTrue(any(item["role"] == "REFERENCE_CHROMOSOME" for item in locus_two))
        for observed in (locus_one, locus_two):
            self.assertTrue(any(
                item["role"] == "CONDITIONAL_CANDIDATE_ATTESTATION" for item in observed
            ))
            self.assertTrue(any(
                item["role"] == "CONTINUATION_CHECKPOINT_RECEIPT" for item in observed
            ))
            self.assertTrue(any(
                item["role"] == "SOURCE_DISCOVERY_CHECKPOINT" for item in observed
            ))

    def test_candidate_attestation_header_only_and_candidate_rows(self):
        header = "\t".join(SUPERVISOR.CONDITIONAL_CANDIDATE_FIELDS) + "\n"
        legacy = Mock()
        loci = [
            {"LOC": index, "CHR": 1, "START": index * 10, "STOP": index * 10 + 9}
            for index in range(1, 2496)
        ]
        loci[463] = {"LOC": 464, "CHR": 3, "START": 47588462, "STOP": 50387742}
        legacy.load_loci.return_value = loci
        with tempfile.TemporaryDirectory() as temporary:
            path = Path(temporary) / "conditional_candidates.tsv"
            path.write_text(header, encoding="utf-8")
            with patch.object(SUPERVISOR.CONTRACT, "legacy_supervisor", return_value=legacy):
                self.assertEqual(SUPERVISOR._candidate_loci(path), frozenset())
                path.write_text(
                    header + "464\t464\t3\t47588462\t50387742\t312\t1\tB\n",
                    encoding="utf-8",
                )
                self.assertEqual(SUPERVISOR._candidate_loci(path), frozenset({464}))
                path.write_text(
                    header + "464\t464\t4\t47588462\t50387742\t312\t1\tB\n",
                    encoding="utf-8",
                )
                with self.assertRaises(SUPERVISOR.ExecutionError):
                    SUPERVISOR._candidate_loci(path)

    def test_aggregate_candidate_dependency_rechecks_receipted_hash(self):
        source = "a" * 64
        continuation = "b" * 64
        header = "\t".join(SUPERVISOR.CONDITIONAL_CANDIDATE_FIELDS) + "\n"
        with tempfile.TemporaryDirectory() as temporary:
            temporary_root = Path(temporary)
            target = temporary_root / "attempt"
            target.mkdir()
            result = target / "result.rds"
            candidates = target / "conditional_candidates.tsv"
            receipt_path = target / "receipt.json"
            result.write_bytes(b"result")
            candidates.write_text(header, encoding="utf-8")
            receipt_path.write_text("{}\n", encoding="utf-8")
            candidate_size, candidate_hash = SUPERVISOR.file_identity(candidates)
            result_size, result_hash = SUPERVISOR.file_identity(result)
            receipt = {"artifacts": [
                {"path": "result.rds", "bytes": result_size, "sha256": result_hash},
                {"path": "conditional_candidates.tsv", "bytes": candidate_size,
                 "sha256": candidate_hash},
            ]}
            legacy = Mock()
            legacy.load_loci.return_value = [
                {"LOC": index, "CHR": 1, "START": index, "STOP": index + 1}
                for index in range(1, 2496)
            ]
            with (
                patch.object(SUPERVISOR, "ROOT", temporary_root),
                patch.object(SUPERVISOR, "validate_bundle", return_value=receipt),
                patch.object(SUPERVISOR, "_safe_attempt_target", return_value=target),
                patch.object(SUPERVISOR.CONTRACT, "legacy_supervisor", return_value=legacy),
            ):
                specs, loci = SUPERVISOR._aggregate_dependency_specs(source, continuation)
                self.assertEqual(loci, frozenset())
                self.assertEqual(
                    {item["role"] for item in specs},
                    {"CONTINUATION_CHECKPOINT", "CONDITIONAL_CANDIDATE_ATTESTATION",
                     "CONTINUATION_CHECKPOINT_RECEIPT"},
                )
                candidates.write_text(header + "tampered\n", encoding="utf-8")
                with self.assertRaises(SUPERVISOR.ExecutionError):
                    SUPERVISOR._aggregate_dependency_specs(source, continuation)

    def test_aggregate_result_receipt_and_changed_input_tampering_fail_closed(self):
        source = "a" * 64
        continuation = "b" * 64
        header = "\t".join(SUPERVISOR.CONDITIONAL_CANDIDATE_FIELDS) + "\n"
        with tempfile.TemporaryDirectory() as temporary:
            temporary_root = Path(temporary)
            target = temporary_root / "attempt"
            target.mkdir()
            result = target / "result.rds"
            candidates = target / "conditional_candidates.tsv"
            receipt_path = target / "receipt.json"
            result.write_bytes(b"sealed aggregate")
            candidates.write_text(header, encoding="utf-8")
            receipt_path.write_text("{}\n", encoding="utf-8")
            result_size, result_hash = SUPERVISOR.file_identity(result)
            candidate_size, candidate_hash = SUPERVISOR.file_identity(candidates)
            receipt = {"artifacts": [
                {"path": "result.rds", "bytes": result_size, "sha256": result_hash},
                {"path": "conditional_candidates.tsv", "bytes": candidate_size,
                 "sha256": candidate_hash},
            ]}
            legacy = Mock()
            legacy.load_loci.return_value = [
                {"LOC": index, "CHR": 1, "START": index, "STOP": index + 1}
                for index in range(1, 2496)
            ]
            with (
                patch.object(SUPERVISOR, "ROOT", temporary_root),
                patch.object(SUPERVISOR, "validate_bundle", return_value=receipt),
                patch.object(SUPERVISOR, "_safe_attempt_target", return_value=target),
                patch.object(SUPERVISOR.CONTRACT, "legacy_supervisor", return_value=legacy),
            ):
                specs, _ = SUPERVISOR._aggregate_dependency_specs(source, continuation)
                result.write_bytes(b"tampered aggregate")
                with self.assertRaises(SUPERVISOR.ExecutionError):
                    SUPERVISOR.capture_active_inputs(specs)

                result.write_bytes(b"sealed aggregate")
                specs, _ = SUPERVISOR._aggregate_dependency_specs(source, continuation)
                receipt_path.write_text('{"tampered":true}\n', encoding="utf-8")
                with self.assertRaises(SUPERVISOR.ExecutionError):
                    SUPERVISOR.capture_active_inputs(specs)

                live_input = temporary_root / "live-input.tsv"
                live_input.write_text("before\n", encoding="utf-8")
                captured = SUPERVISOR.capture_active_inputs([
                    SUPERVISOR._artifact_spec(live_input, "COMMON_INPUT")
                ])
                live_input.write_text("after!\n", encoding="utf-8")
                with self.assertRaises(SUPERVISOR.ExecutionError):
                    SUPERVISOR.complete_active_inputs(captured)

    def test_header_only_candidate_attestation_cannot_omit_an_eligible_locus(self):
        script = r'''
expressions <- parse(file = "scripts/135_run_track_b_lava_postdiscovery_v2.R")
environment <- new.env(parent = globalenv())
wanted <- c(
  "build_candidates", "canonical_candidates", "read_conditional_candidates",
  "candidate_fields"
)
for (expression in expressions) {
  if (is.call(expression) && identical(expression[[1L]], as.name("<-")) &&
      is.symbol(expression[[2L]]) && as.character(expression[[2L]]) %in% wanted) {
    eval(expression, environment)
  }
}
environment$loci <- data.frame(LOC = 1L, CHR = 3L, START = 10L, STOP = 20L)
environment$conditioners <- data.frame(
  pair_id = "B", conditional_model_id = "B_MDD_ONLY", covariates = "MDD"
)
environment$policy_num <- function(key) 0.05
family <- list(
  bivar = data.frame(
    LOC = 1L, pair_id = "B", fdr_significant = TRUE,
    stringsAsFactors = FALSE
  ),
  univ = data.frame(
    LOC = 1L, phen = "MDD", analysis_status = "TESTED", p = 0.01,
    stringsAsFactors = FALSE
  ),
  status = data.frame(LOC = 1L, n_snps = 12L)
)
path <- tempfile()
writeLines(paste(environment$candidate_fields, collapse = "\t"), path)
omission_failed <- tryCatch({
  environment$read_conditional_candidates(path, family)
  FALSE
}, error = function(error) TRUE)
stopifnot(omission_failed)
write.table(
  environment$build_candidates(family), path, sep = "\t", quote = FALSE,
  row.names = FALSE, na = "NA"
)
observed <- environment$read_conditional_candidates(path, family)
stopifnot(identical(as.integer(observed$locus_index), 1L))
'''
        result = subprocess.run(
            [str(ROOT / ".r-env/bin/Rscript"), "-e", script],
            cwd=ROOT, text=True, capture_output=True, check=False,
        )
        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)

    def test_conditional_prerequisite_uses_sealed_aggregate_without_recomputing_semantics(self):
        validator = Mock(return_value={"phase": "aggregate-discovery"})
        with patch.object(SUPERVISOR, "validate_bundle", validator):
            SUPERVISOR.validate_phase_prerequisites(
                "conditional", "a" * 64, "b" * 64
            )
        validator.assert_called_once_with(
            "aggregate-discovery", None, "a" * 64, "b" * 64,
            revalidate_semantics=False, reconcile_active=False,
        )

    def test_quick_contract_checks_live_fingerprint_once_then_exact_lineage(self):
        source = "a" * 64
        continuation = "b" * 64
        with (
            patch.object(
                SUPERVISOR.CONTRACT, "continuation_fingerprint",
                return_value=continuation,
            ) as fingerprint,
            patch.object(
                SUPERVISOR.CONTRACT, "validate_lineage_record",
                return_value={"continuation_execution_fingerprint": continuation},
            ) as lineage,
        ):
            SUPERVISOR.quick_contract(source, continuation)
        fingerprint.assert_called_once_with(source)
        lineage.assert_called_once_with(source, continuation)

        with (
            patch.object(
                SUPERVISOR.CONTRACT, "continuation_fingerprint",
                return_value="c" * 64,
            ),
            patch.object(SUPERVISOR.CONTRACT, "validate_lineage_record") as lineage,
        ):
            with self.assertRaisesRegex(
                SUPERVISOR.ExecutionError, "live continuation code/runtime"
            ):
                SUPERVISOR.quick_contract(source, continuation)
        lineage.assert_not_called()

    def test_internal_measure_rejects_stale_identity_before_launch_or_output(self):
        source = "a" * 64
        stale = "b" * 64
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            arguments = argparse.Namespace(
                source_fingerprint=source,
                continuation_fingerprint=stale,
                phase="conditional",
                locus_index=1,
                worker_output=str(root / "worker.rds"),
                internal_result=str(root / "metrics.json"),
            )
            with (
                patch.object(
                    SUPERVISOR, "quick_contract",
                    side_effect=SUPERVISOR.ExecutionError("stale fingerprint"),
                ) as contract_check,
                patch.object(SUPERVISOR.subprocess, "run") as worker,
            ):
                with self.assertRaisesRegex(
                    SUPERVISOR.ExecutionError, "stale fingerprint"
                ):
                    SUPERVISOR.internal_measure(arguments)
            contract_check.assert_called_once_with(source, stale)
            worker.assert_not_called()
            self.assertFalse((root / "worker.rds").exists())
            self.assertFalse((root / "metrics.json").exists())

    def test_source_lock_cache_is_identity_keyed_and_local_lookup_is_exact(self):
        source = "a" * 64
        entries = [
            {
                "locus_index": index,
                "result": {
                    "path": f"source/{index:04d}/result.rds",
                    "bytes": index,
                    "sha256": format(index, "064x"),
                },
            }
            for index in range(1, 2496)
        ]
        lock = {"checkpoints": entries}
        SUPERVISOR._SOURCE_LOCK_CACHE.clear()
        with (
            patch.object(
                SUPERVISOR.CONTRACT, "source_lock_path",
                return_value=Path("source-lock.json"),
            ),
            patch.object(
                SUPERVISOR, "file_identity", return_value=(123, "c" * 64)
            ),
            patch.object(
                SUPERVISOR.CONTRACT, "_validate_legacy_code", return_value={"v1": "exact"}
            ),
            patch.object(
                SUPERVISOR.CONTRACT, "validate_source_lock", return_value=lock
            ) as validate,
            patch.object(
                SUPERVISOR, "_artifact_spec",
                side_effect=lambda path, role, expected=None: {
                    "path": str(path), "bytes": expected[0],
                    "sha256": expected[1], "role": role,
                },
            ),
        ):
            first = SUPERVISOR._source_discovery_spec(source, 22)
            second = SUPERVISOR._source_discovery_spec(source, 22)
        self.assertEqual(first, second)
        self.assertEqual(first["bytes"], 22)
        self.assertEqual(first["sha256"], format(22, "064x"))
        validate.assert_called_once_with(
            source, validate_bundles=False, require_live_fingerprint=False,
        )

        malformed = {"checkpoints": list(entries)}
        malformed["checkpoints"][21] = dict(malformed["checkpoints"][21])
        malformed["checkpoints"][21]["locus_index"] = 23
        with patch.object(SUPERVISOR, "_source_lock", return_value=malformed):
            with self.assertRaisesRegex(
                SUPERVISOR.ExecutionError, "local checkpoint identity"
            ):
                SUPERVISOR._source_discovery_spec(source, 22)

    def test_source_lock_cache_rejects_change_during_initial_validation(self):
        source = "a" * 64
        SUPERVISOR._SOURCE_LOCK_CACHE.clear()
        with (
            patch.object(
                SUPERVISOR.CONTRACT, "source_lock_path",
                return_value=Path("source-lock.json"),
            ),
            patch.object(
                SUPERVISOR, "file_identity",
                side_effect=[(123, "c" * 64), (123, "d" * 64)],
            ),
            patch.object(
                SUPERVISOR.CONTRACT, "_validate_legacy_code", return_value={"v1": "exact"}
            ),
            patch.object(
                SUPERVISOR.CONTRACT, "validate_source_lock",
                return_value={"checkpoints": []},
            ),
        ):
            with self.assertRaisesRegex(
                SUPERVISOR.ExecutionError, "changed while entering"
            ):
                SUPERVISOR._source_lock(source)
        self.assertEqual(SUPERVISOR._SOURCE_LOCK_CACHE, {})

    def test_post_publication_check_reuses_proven_semantics_but_binds_capture(self):
        source = CONTRACT.PINNED_SOURCE_FINGERPRINT
        continuation = "b" * 64
        active = [{
            "path": "input.tsv", "bytes": 1, "sha256": "c" * 64,
            "role": "COMMON_INPUT", "pre_stat": {"bytes": 1},
            "post_stat": {"bytes": 1},
        }]
        marker = {"qc": "NOT_BIVARIATE_FDR_ELIGIBLE"}
        metrics = {"command": ["Rscript", "worker.R"], "runtime_sec": 2.0,
                   "peak_rss_bytes": 1024}

        with tempfile.TemporaryDirectory() as temporary:
            attempt = Path(temporary) / "attempt"
            attempt.mkdir()
            (attempt / "result.rds").write_bytes(b"sealed-result")
            semantic_artifacts = SUPERVISOR.capture_semantic_artifacts(
                attempt, "conditional"
            )
            semantic_binding = SUPERVISOR.bind_semantic_validation(
                semantic_artifacts, semantic_artifacts, "conditional"
            )
            receipted_artifacts = [
                {key: item[key] for key in ("path", "bytes", "sha256")}
                for item in semantic_artifacts
            ]
            validator = Mock(return_value={
                "active_inputs": active,
                "semantic_artifact_binding": semantic_binding,
                "artifacts": receipted_artifacts,
            })
            with (
                patch.object(SUPERVISOR, "sha256", return_value="d" * 64),
                patch.object(
                    SUPERVISOR, "bundle_path", return_value=Path(temporary) / "ready"
                ),
                patch.object(SUPERVISOR, "_relative", return_value="conditional/locus_0001"),
                patch.object(SUPERVISOR, "_publish_ready_bundle"),
                patch.object(SUPERVISOR, "validate_bundle", validator),
                patch.object(SUPERVISOR, "update_benchmark") as benchmark,
                patch.object(SUPERVISOR, "_semantic_validation") as semantic,
            ):
                observed = SUPERVISOR._publish_receipt(
                    attempt, "conditional", 1, source, continuation, marker,
                    metrics, active, semantic_binding,
                    "SEMANTIC_PASS\n", "worker log\n",
                )
            self.assertEqual(observed["active_inputs"], active)
            validator.assert_called_once_with(
                "conditional", 1, source, continuation,
                revalidate_semantics=False, reconcile_active=False,
            )
            semantic.assert_not_called()
            benchmark.assert_called_once_with(observed)

        with tempfile.TemporaryDirectory() as temporary:
            attempt = Path(temporary) / "attempt"
            attempt.mkdir()
            (attempt / "result.rds").write_bytes(b"sealed-result")
            semantic_artifacts = SUPERVISOR.capture_semantic_artifacts(
                attempt, "conditional"
            )
            semantic_binding = SUPERVISOR.bind_semantic_validation(
                semantic_artifacts, semantic_artifacts, "conditional"
            )
            receipted_artifacts = [
                {key: item[key] for key in ("path", "bytes", "sha256")}
                for item in semantic_artifacts
            ]
            with (
                patch.object(SUPERVISOR, "sha256", return_value="d" * 64),
                patch.object(
                    SUPERVISOR, "bundle_path", return_value=Path(temporary) / "ready"
                ),
                patch.object(SUPERVISOR, "_relative", return_value="conditional/locus_0001"),
                patch.object(SUPERVISOR, "_publish_ready_bundle"),
                patch.object(
                    SUPERVISOR, "validate_bundle",
                    return_value={
                        "active_inputs": [],
                        "semantic_artifact_binding": semantic_binding,
                        "artifacts": receipted_artifacts,
                    },
                ),
                patch.object(SUPERVISOR, "update_benchmark"),
            ):
                with self.assertRaisesRegex(
                    SUPERVISOR.ExecutionError, "active-input receipt differs"
                ):
                    SUPERVISOR._publish_receipt(
                        attempt, "conditional", 1, source, continuation, marker,
                        metrics, active, semantic_binding,
                        "SEMANTIC_PASS\n", "worker log\n",
                    )

    def test_semantic_artifact_binding_detects_content_and_inode_replacement(self):
        with tempfile.TemporaryDirectory() as temporary:
            attempt = Path(temporary) / "attempt"
            attempt.mkdir()
            result = attempt / "result.rds"
            result.write_bytes(b"sealed-result")
            sealed = SUPERVISOR.capture_semantic_artifacts(attempt, "conditional")

            result.write_bytes(b"forged-result")
            with self.assertRaisesRegex(
                SUPERVISOR.ExecutionError, "changed across semantic validation"
            ):
                SUPERVISOR.require_semantic_artifacts_unchanged(
                    attempt, "conditional", sealed
                )

            result.write_bytes(b"sealed-result")
            sealed = SUPERVISOR.capture_semantic_artifacts(attempt, "conditional")
            replacement = attempt / "replacement.rds"
            replacement.write_bytes(b"sealed-result")
            replacement.replace(result)
            with self.assertRaisesRegex(
                SUPERVISOR.ExecutionError, "changed across semantic validation"
            ):
                SUPERVISOR.require_semantic_artifacts_unchanged(
                    attempt, "conditional", sealed
                )

    def test_publication_rejects_mutation_after_semantic_validation(self):
        source = CONTRACT.PINNED_SOURCE_FINGERPRINT
        continuation = "b" * 64
        active = [{
            "path": "input.tsv", "bytes": 1, "sha256": "c" * 64,
            "role": "COMMON_INPUT", "pre_stat": {"bytes": 1},
            "post_stat": {"bytes": 1},
        }]
        marker = {"qc": "NOT_BIVARIATE_FDR_ELIGIBLE"}
        metrics = {
            "command": ["Rscript", "worker.R"], "runtime_sec": 2.0,
            "peak_rss_bytes": 1024,
        }
        with tempfile.TemporaryDirectory() as temporary:
            attempt = Path(temporary) / "attempt"
            attempt.mkdir()
            result = attempt / "result.rds"
            result.write_bytes(b"sealed-result")
            semantic_artifacts = SUPERVISOR.capture_semantic_artifacts(
                attempt, "conditional"
            )
            semantic_binding = SUPERVISOR.bind_semantic_validation(
                semantic_artifacts, semantic_artifacts, "conditional"
            )
            real_publish = SUPERVISOR.publish_no_replace

            def publish_then_mutate(path, value):
                real_publish(path, value)
                if path.name == "worker.log":
                    result.write_bytes(b"forged-result")

            with (
                patch.object(SUPERVISOR, "publish_no_replace", side_effect=publish_then_mutate),
                patch.object(
                    SUPERVISOR, "bundle_path", return_value=Path(temporary) / "ready"
                ),
                patch.object(SUPERVISOR, "_relative", return_value="conditional/locus_0001"),
            ):
                with self.assertRaisesRegex(
                    SUPERVISOR.ExecutionError,
                    "receipt artifacts differ from the semantically validated bytes",
                ):
                    SUPERVISOR._publish_receipt(
                        attempt, "conditional", 1, source, continuation, marker,
                        metrics, active, semantic_binding,
                        "SEMANTIC_PASS\n", "worker log\n",
                    )
            self.assertFalse((attempt / "receipt.json").exists())
            self.assertFalse((attempt / "READY").exists())

    def test_terminal_semantic_binding_covers_result_qc_and_all_staged_tables(self):
        with tempfile.TemporaryDirectory() as temporary:
            attempt = Path(temporary) / "attempt"
            attempt.mkdir()
            for name in SUPERVISOR.SEMANTIC_ARTIFACT_NAMES["terminal-qc"]:
                (attempt / name).write_bytes((name + "\n").encode())
            sealed = SUPERVISOR.capture_semantic_artifacts(attempt, "terminal-qc")
            self.assertEqual(
                {item["path"] for item in sealed},
                set(SUPERVISOR.SEMANTIC_ARTIFACT_NAMES["terminal-qc"]),
            )
            (attempt / "terminal_qc.json").write_bytes(b"{}\n")
            with self.assertRaisesRegex(
                SUPERVISOR.ExecutionError, "changed across semantic validation"
            ):
                SUPERVISOR.require_semantic_artifacts_unchanged(
                    attempt, "terminal-qc", sealed
                )

    def test_finalize_keeps_full_family_validation_and_identity_binding(self):
        source = "a" * 64
        continuation = "b" * 64
        validator = Mock(return_value={"phase": "ready"})
        with patch.object(SUPERVISOR, "validate_bundle", validator):
            SUPERVISOR.validate_phase_prerequisites(
                "finalize", source, continuation
            )
        calls = validator.call_args_list
        self.assertEqual(len(calls), 2496)
        self.assertEqual(
            calls[0].args,
            ("aggregate-discovery", None, source, continuation),
        )
        self.assertEqual(calls[0].kwargs, {})
        self.assertEqual(
            [call.args[1] for call in calls[1:]], list(range(1, 2496))
        )
        self.assertTrue(all(
            call.kwargs == {
                "revalidate_semantics": False, "reconcile_active": False,
            }
            for call in calls[1:]
        ))

        runner = (ROOT / "scripts/135_run_track_b_lava_postdiscovery_v2.R").read_text(
            encoding="utf-8"
        )
        validator_r = (ROOT / "scripts/137_validate_track_b_lava_checkpoint_v2.R").read_text(
            encoding="utf-8"
        )
        self.assertIn(
            "conditional_reads <- lapply(seq_along(conditional_paths)", runner
        )
        self.assertIn(
            "aggregate_read <- load_aggregate(full_source_validation = TRUE)",
            validator_r,
        )
        self.assertIn(
            "conditional_reads <- lapply(seq_along(conditional_paths)", validator_r
        )

    def test_rds_manifest_hashes_the_exact_bytes_deserialized(self):
        script = r'''
expressions <- parse(file = "scripts/135_run_track_b_lava_postdiscovery_v2.R")
environment <- new.env(parent = globalenv())
wanted <- c("sha256_files", "read_rds_with_sha256")
for (expression in expressions) {
  if (is.call(expression) && identical(expression[[1L]], as.name("<-")) &&
      is.symbol(expression[[2L]]) && as.character(expression[[2L]]) %in% wanted) {
    eval(expression, environment)
  }
}
source <- tempfile()
saveRDS(list(value = 1L), source, version = 3)
expected <- environment$sha256_files(source)[[1L]]
real_sha256 <- environment$sha256_files
environment$sha256_files <- function(paths) {
  saveRDS(list(value = 999L), source, version = 3)
  real_sha256(paths)
}
observed <- environment$read_rds_with_sha256(source, "mutation fixture")
stopifnot(identical(observed$value, list(value = 1L)))
stopifnot(identical(observed$sha256, expected))
'''
        for filename in (
            "scripts/135_run_track_b_lava_postdiscovery_v2.R",
            "scripts/137_validate_track_b_lava_checkpoint_v2.R",
        ):
            with self.subTest(filename=filename):
                current = script.replace(
                    "scripts/135_run_track_b_lava_postdiscovery_v2.R", filename
                )
                result = subprocess.run(
                    [str(ROOT / ".r-env/bin/Rscript"), "-e", current],
                    cwd=ROOT, text=True, capture_output=True, check=False,
                )
                self.assertEqual(result.returncode, 0, result.stdout + result.stderr)

    def test_conditional_numeric_anomalies_never_become_tested(self):
        script = r'''
expressions <- parse(file = "scripts/135_run_track_b_lava_postdiscovery_v2.R")
environment <- new.env(parent = globalenv())
wanted <- c("empty_conditional", "finite_or_na_scalar", "process_conditional_locus")
for (expression in expressions) {
  if (is.call(expression) && identical(expression[[1L]], as.name("<-")) &&
      is.symbol(expression[[2L]]) && as.character(expression[[2L]]) %in% wanted) {
    eval(expression, environment)
  }
}
environment$loci <- data.frame(LOC = "1", CHR = 1, START = 1, STOP = 2)
environment$conditioners <- data.frame(
  pair_id = "A", conditional_model_id = "A_BMI_ONLY", covariates = "BMI"
)
environment$values <- list(source_fingerprint = "source", continuation_fingerprint = "continuation")
environment$schema_version <- "fixture"
environment$aggregate_sha256 <- "aggregate"
environment$policy_num <- function(key) if (key == "univariate_p_threshold") 0.05 else 0.95
environment$policy_int <- function(key) 1L
environment$policy_bool <- function(key) TRUE
environment$scoped_input <- function(locus) list()
environment$process.locus <- function(...) list()
family <- list(
  bivar = data.frame(
    LOC = "1", CHR = 1, START = 1, STOP = 2, pair_id = "A",
    trait1 = "T1", trait2 = "T2", fdr_significant = TRUE
  ),
  univ = data.frame(LOC = "1", phen = "BMI", analysis_status = "TESTED", p = 0.01),
  status = data.frame(LOC = "1", n_snps = 10)
)
run_case <- function(partial, expected) {
  environment$run.pcor <- function(...) partial
  observed <- environment$process_conditional_locus(1L, family, TRUE)$conditional
  stopifnot(identical(as.character(observed$analysis_status), expected))
  observed
}
valid <- data.frame(
  pcor = 0.2, ci.lower = 0.1, ci.upper = 0.3, p = 0.01,
  r2.phen1_z = 0.1, r2.phen2_z = 0.2
)
run_case(valid, "TESTED")
bad_pcor <- valid; bad_pcor$pcor <- 1.2
stopifnot(is.na(run_case(bad_pcor, "CONDITIONAL_FAILED")$pcor))
bad_ci <- valid; bad_ci$ci.lower <- 0.25
stopifnot(is.na(run_case(bad_ci, "CONDITIONAL_FAILED")$pcor))
bad_p <- valid; bad_p$p <- 1.1
stopifnot(is.na(run_case(bad_p, "CONDITIONAL_FAILED")$p))
bad_r2 <- valid; bad_r2$r2.phen1_z <- -0.1
stopifnot(is.na(run_case(bad_r2, "CONDITIONAL_FAILED")$r2.trait1_z))
high_r2 <- valid; high_r2$r2.phen1_z <- 0.96
stopifnot(is.na(run_case(high_r2, "CONDITIONAL_UNSTABLE_MAX_R2")$p))
no_target <- family
no_target$bivar <- no_target$bivar[0, , drop = FALSE]
environment$scoped_input <- function(...) stop("unused chromosome input was evaluated")
noncandidate <- environment$process_conditional_locus(1L, no_target, FALSE)
stopifnot(identical(noncandidate$state, "NOT_BIVARIATE_FDR_ELIGIBLE"))
candidate_mismatch_failed <- tryCatch({
  environment$process_conditional_locus(1L, family, FALSE)
  FALSE
}, error = function(error) TRUE)
stopifnot(candidate_mismatch_failed)
'''
        result = subprocess.run(
            [str(ROOT / ".r-env/bin/Rscript"), "-e", script],
            cwd=ROOT, text=True, capture_output=True, check=False,
        )
        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)


class TrackBLavaTerminalQCTests(unittest.TestCase):
    def test_incomplete_family_is_error_not_terminal_qc(self):
        with tempfile.TemporaryDirectory() as temporary:
            directory = Path(temporary) / "finalize_unit_all.1"
            directory.mkdir()
            with (
                patch.object(RESULTS, "_quick_contract"),
                patch.object(RESULTS, "_attempt_directory", return_value=directory),
                patch.object(RESULTS, "_staged_paths", return_value=[]),
                patch.object(
                    RESULTS,
                    "_diagnostic_validate",
                    side_effect=RESULTS.ValidationError("incomplete family"),
                ),
            ):
                with self.assertRaises(RESULTS.ValidationError):
                    RESULTS.validate_staged(directory, "a" * 64, "b" * 64)
            self.assertFalse((directory / "terminal_qc.json").exists())

    def test_diagnostic_pass_collects_all_gate_failures_and_restores_policy(self):
        policy = {
            "maximum_locus_failure_fraction": 0.10,
            "maximum_univariate_untested_fraction": 0.25,
            "maximum_bivariate_failure_fraction": 0.05,
            "maximum_conditional_failure_fraction": 0.05,
        }
        original_loader = lambda: dict(policy)
        observed_diagnostic_thresholds = []

        def fake_validate_low_level(_paths):
            diagnostic = RESULTS.LEGACY.contract.load_policy()
            observed_diagnostic_thresholds.append(
                tuple(
                    diagnostic[key]
                    for key in (
                        "maximum_locus_failure_fraction",
                        "maximum_univariate_untested_fraction",
                        "maximum_bivariate_failure_fraction",
                        "maximum_conditional_failure_fraction",
                    )
                )
            )
            return {
                "status": [
                    {"status": "PROCESSED"},
                    {"status": "PROCESS_FAILED"},
                    {"status": "PROCESSED"},
                    {"status": "UNIVARIATE_FAILED"},
                ],
                "univ": [
                    {"analysis_status": "TESTED"},
                    {"analysis_status": "PHENOTYPE_DROPPED"},
                    {"analysis_status": "LOCUS_PROCESS_FAILED"},
                    {"analysis_status": "TESTED"},
                ],
                "bivar": [],
                "conditional": [
                    {"analysis_status": "CONDITIONER_LOCAL_H2_INELIGIBLE"}
                ],
            }

        with tempfile.TemporaryDirectory() as temporary:
            temporary_root = Path(temporary)
            paths = [temporary_root / f"result-{index}.tsv" for index in range(4)]
            for path in paths:
                path.write_text("header\n", encoding="utf-8")
            with (
                patch.object(RESULTS, "ROOT", temporary_root),
                patch.object(RESULTS.LEGACY.contract, "load_policy", original_loader),
                patch.object(
                    RESULTS.LEGACY, "validate_low_level", side_effect=fake_validate_low_level
                ),
            ):
                _, qc = RESULTS._diagnostic_validate(paths)
                self.assertIs(RESULTS.LEGACY.contract.load_policy, original_loader)

        self.assertEqual(observed_diagnostic_thresholds, [(1.0, 1.0, 1.0, 1.0)])
        self.assertEqual(
            [reason["code"] for reason in qc["reasons"]],
            [
                "LOCUS_FAILURE_FRACTION_EXCEEDED",
                "UNIVARIATE_UNTESTED_FRACTION_EXCEEDED",
            ],
        )
        self.assertEqual(qc["reasons"][0]["numerator"], 2)
        self.assertEqual(qc["reasons"][0]["denominator"], 4)

    def test_passing_diagnostic_is_rechecked_by_unmodified_validator(self):
        validator = Mock(return_value={})
        with tempfile.TemporaryDirectory() as temporary:
            directory = Path(temporary) / "finalize_unit_all.1"
            directory.mkdir()
            with (
                patch.object(RESULTS, "_quick_contract"),
                patch.object(RESULTS, "_attempt_directory", return_value=directory),
                patch.object(RESULTS, "_staged_paths", return_value=[]),
                patch.object(RESULTS, "_diagnostic_validate", return_value=({}, {"reasons": []})),
                patch.object(RESULTS.LEGACY, "validate", validator),
            ):
                RESULTS.validate_staged(directory, "a" * 64, "b" * 64)
        validator.assert_called_once_with([], build_summary_files=True)

    def test_terminal_qc_is_no_replace_and_never_becomes_a_pass(self):
        payload = {
            "schema_version": RESULTS.TERMINAL_QC_SCHEMA,
            "state": "TERMINAL_FAILED_QC",
            "full_family_complete": True,
            "scientific_validation_passed": False,
            "canonical_publication_allowed": False,
            "qc": {"reasons": [{"code": "LOCUS_FAILURE_FRACTION_EXCEEDED"}]},
        }
        diagnostic = ({}, {"reasons": payload["qc"]["reasons"]})
        with tempfile.TemporaryDirectory() as temporary:
            directory = Path(temporary) / "finalize_unit_all.1"
            directory.mkdir()
            patches = (
                patch.object(RESULTS, "_quick_contract"),
                patch.object(RESULTS, "_attempt_directory", return_value=directory),
                patch.object(RESULTS, "_staged_paths", return_value=[]),
                patch.object(RESULTS, "_diagnostic_validate", return_value=diagnostic),
                patch.object(RESULTS, "_terminal_payload", return_value=payload),
            )
            with patches[0], patches[1], patches[2], patches[3], patches[4]:
                with self.assertRaises(RESULTS.TerminalQCOutcome):
                    RESULTS.validate_staged(directory, "a" * 64, "b" * 64)
                attestation = directory / "terminal_qc.json"
                first = attestation.read_bytes()
                self.assertEqual(json.loads(first)["state"], "TERMINAL_FAILED_QC")
                with self.assertRaises(RESULTS.TerminalQCOutcome):
                    RESULTS.validate_staged(directory, "a" * 64, "b" * 64)
                self.assertEqual(attestation.read_bytes(), first)
                attestation.write_text("{}\n", encoding="utf-8")
                with self.assertRaises(RESULTS.ValidationError):
                    RESULTS.validate_staged(directory, "a" * 64, "b" * 64)


if __name__ == "__main__":
    unittest.main()
