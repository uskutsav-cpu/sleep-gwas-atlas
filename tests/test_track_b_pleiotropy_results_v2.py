#!/usr/bin/env python3

from __future__ import annotations

import csv
from contextlib import redirect_stderr
import gzip
import importlib.util
import io
import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest
from unittest import mock


ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "scripts/147_build_track_b_pleiotropy_results_v2.py"


def load_collator():
    specification = importlib.util.spec_from_file_location("track_b_results_v2_test", SCRIPT)
    if specification is None or specification.loader is None:
        raise RuntimeError("could not load post-PLACO collator")
    module = importlib.util.module_from_spec(specification)
    specification.loader.exec_module(module)
    return module


COLLATOR = load_collator()


def placo_row(
    pair: str, snp: str, position: int, p_value: str, q_value: str,
    *, status: str = "TESTED", family_n: int = 1,
) -> dict[str, str]:
    numerical = status == "NUMERICAL_FAILURE_P_SET_TO_ONE"
    return {
        "analysis_id": "track-b-v1.0-pleiotropy", "pair_id": pair,
        "SNP": snp, "CHR": "1", "BP": str(position), "A1": "A", "A2": "C",
        "Z1": "1.25", "Z2": "-0.75", "P1": "0.1", "P2": "0.2",
        "T_PLACO_PLUS": "NA" if numerical else "0.5",
        "P_PLACO_PLUS": p_value, "PLACO_BH_Q": q_value,
        "within_pair_family_n": str(family_n), "analysis_status": status,
        "numerical_error": "synthetic_underflow" if numerical else "NA",
    }


def write_pair(root: Path, pair_id: str, rows: list[dict[str, str]]):
    directory = root / "inputs"
    directory.mkdir(parents=True, exist_ok=True)
    for row in rows:
        row["within_pair_family_n"] = str(len(rows))
    ledger = directory / f"{pair_id}.full.tsv.gz"
    with gzip.open(ledger, "wt", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=COLLATOR.PLACO_LEDGER_FIELDS, delimiter="\t",
            lineterminator="\n",
        )
        writer.writeheader()
        writer.writerows(rows)
    identity = COLLATOR.portable_identity(root, ledger.relative_to(root))
    primary = sum(
        pair_id != "CONTROL" and float(row["P_PLACO_PLUS"]) <= 2.5e-8
        for row in rows
    )
    pairwise = sum(float(row["P_PLACO_PLUS"]) <= 5e-8 for row in rows)
    bh = sum(float(row["PLACO_BH_Q"]) <= 0.05 for row in rows)
    counts = {
        "rows": len(rows),
        "failures": sum(row["analysis_status"].startswith("NUMERICAL") for row in rows),
        "primary": primary, "pairwise": pairwise, "bh": bh,
        "minimum_p": min(float(row["P_PLACO_PLUS"]) for row in rows),
        "terminal_status": "COMPLETE_WITH_HITS" if primary or pairwise or bh else "TESTED_NO_HIT",
    }
    provenance_path = directory / f"{pair_id}.provenance.json"
    provenance = {
        "schema_version": "sleep-atlas-track-b-placo-result.1",
        "analysis_id": "track-b-v1.0-pleiotropy", "pair_id": pair_id,
        "terminal_result_state": counts["terminal_status"], "qc_status": "PASS",
        "output_rows": len(rows), "within_pair_bh_family_n": len(rows),
        "output_sha256": identity["sha256"], "output_bytes": identity["bytes"],
        "exact_schema": ",".join(COLLATOR.PLACO_LEDGER_FIELDS),
        "reference_sha256": "NOT_APPLICABLE_TO_PLACO_FULL_P_SCAN_LD_REQUIRED_FOR_LATER_CLUMPING",
        "complete_family_counts": counts,
        "nuisance_estimation": {
            "input_rows": len(rows),
            "scope": "ALL_VALID_ALIGNED_GENOME_WIDE_VARIANTS_BEFORE_ANY_SINGLE_VARIANT_SHARD",
        },
    }
    provenance_path.write_text(json.dumps(provenance) + "\n", encoding="utf-8")
    provenance_identity = COLLATOR.portable_identity(root, provenance_path.relative_to(root))
    return COLLATOR.PairInput(
        pair_id, ledger, provenance_path, counts["terminal_status"], len(rows),
        identity, provenance_identity,
    )


def insert_eligible(
    connection, pair: str, snp: str, position: int, p_value: str = "1e-9",
    q_value: str = "0.01",
) -> None:
    role = COLLATOR.PAIR_IDENTITIES[pair][2]
    connection.execute(
        "INSERT INTO eligible (pair_id,family_role,source_row_index,evidence_id,snp,chr,bp,"
        "a1,a2,z1,z2,p1,p2,p,q,eligibility_basis,evidence_labels,reference_status,"
        "reference_chr,reference_bp,reference_a1,reference_a2) "
        "VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?)",
        (
            pair, role, position, f"PLACO:{pair}:{snp}:1:{position}", snp, 1, position,
            "A", "C", "1.2", "-0.8", "0.1", "0.2", p_value, q_value,
            "PLACO_P_LE_5E_8;PLACO_COMPLETE_WITHIN_PAIR_BH_Q_LE_0_05",
            ";".join(COLLATOR.evidence_labels(pair, float(p_value), float(q_value))),
            "MATCHED", "1", str(position), "C", "A",
        ),
    )
    connection.commit()


def fake_zero_package(root: Path, fingerprint: str) -> tuple[Path, dict[Path, Path]]:
    package = root / COLLATOR.PACKAGE_ROOT / fingerprint
    package.mkdir(parents=True)
    mapping = COLLATOR.package_source_map("ZERO_ELIGIBLE_SIGNAL_FAMILY")
    for index, packaged in enumerate(mapping.values(), start=1):
        path = package / packaged
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(f"artifact-{index}\n", encoding="utf-8")
    files = COLLATOR._inventory_files(package)
    manifest = {
        "schema_version": "sleep-atlas-track-b-post-placo-package.2",
        "fingerprint": fingerprint, "science_state": "ZERO_ELIGIBLE_SIGNAL_FAMILY",
        "primary_has_loci": False, "control_has_loci": False,
        "canonical_mapping": {str(key): str(value) for key, value in mapping.items()},
        "files": files,
        "publication_order": [str(path) for path in COLLATOR.publication_order(mapping)],
    }
    (package / "package.manifest.json").write_text(
        json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8",
    )
    return package, mapping


def synthetic_contract(root: Path) -> dict:
    (root / COLLATOR.COLLATION_CONTRACT.parent).mkdir(parents=True, exist_ok=True)
    (root / COLLATOR.COLLATION_CONTRACT).write_text("synthetic contract\n", encoding="utf-8")
    identity = {"path": "unused", "bytes": 1, "sha256": "e" * 64}
    upstream = {
        str(path): {**identity, "path": str(path)}
        for path in (
            COLLATOR.POLICY, COLLATOR.PLEIOTROPY_CONTRACT_LOCK,
            COLLATOR.INPUT_GATE_LOCK, COLLATOR.SCRIPT, COLLATOR.SEQUENTIAL_SCRIPT,
        )
    }
    reference = {
        suffix: {"path": f"ref/g1000_eur.{suffix}", "bytes": 1, "sha256": suffix[0] * 64}
        for suffix in ("bed", "bim", "fam")
    }
    reference.update({
        "variant_count": 1, "autosomal_variant_count": 1,
        "chromosome_lengths": {"1": 1000}, "source_release": "synthetic",
        "ancestry": "EUR", "build": "GRCh37/hg19",
        "plink": {
            "path": "ref/plink", "bytes": 1, "sha256": "9" * 64,
            "version": "PLINK v1.9.0-b.7.11 synthetic",
        },
    })
    return {
        "upstream": upstream, "ld_reference": reference,
        "thresholds": {
            "eligible_union": "P_PLACO_PLUS<=5e-8 OR COMPLETE_WITHIN_PAIR_BH_Q<=0.05",
            "primary_headline_A_B_only": 2.5e-8, "clump_window_kb": 1000,
            "clump_r2": 0.1, "lead_independence": "STRICT_R2_LT_0.1_WITHIN_1_MB",
        },
        "official_lava_partition": {
            "path": str(COLLATOR.OFFICIAL_BLOCKS), "bytes": 1, "sha256": "8" * 64,
            "block_count": 2495, "role": "DOWNSTREAM_WORK_PARTITION_ONLY",
        },
        "lava_failed_qc_accounting": {"scientific_state": "TERMINAL_FAILED_QC"},
        "archive_reference_state_guard": {"synthetic": "invariant"},
    }


class TrackBPleiotropyResultsV2Tests(unittest.TestCase):
    def test_threshold_union_and_primary_vs_control_labels(self) -> None:
        self.assertEqual(
            COLLATOR.eligibility_basis("A", 6e-8, 0.04),
            "PLACO_COMPLETE_WITHIN_PAIR_BH_Q_LE_0_05",
        )
        self.assertEqual(
            COLLATOR.eligibility_basis("A", 5e-8, 0.05),
            "PLACO_P_LE_5E_8;PLACO_COMPLETE_WITHIN_PAIR_BH_Q_LE_0_05",
        )
        self.assertIsNone(COLLATOR.eligibility_basis("CONTROL", 5.01e-8, 0.0501))
        self.assertEqual(
            COLLATOR.evidence_labels("A", 2.5e-8, 0.05),
            ["PLACO_PRIMARY_HEADLINE", "PLACO_PAIRWISE_GWS", "PLACO_WITHIN_PAIR_FDR"],
        )
        self.assertEqual(
            COLLATOR.evidence_labels("CONTROL", 2e-8, 0.01),
            ["CONTROL_PLACO_GWS_RECOVERED", "CONTROL_PLACO_WITHIN_PAIR_FDR"],
        )

    def test_every_mutating_cli_route_requires_explicit_execute(self) -> None:
        for action in ("--materialize-ld", "--seal-contract", "--run"):
            with self.subTest(action=action), mock.patch.object(
                COLLATOR, "materialize_ld_reference"
            ) as materialize, mock.patch.object(
                COLLATOR, "seal_pre_result_contract"
            ) as seal, mock.patch.object(
                COLLATOR, "build_and_publish_results"
            ) as run:
                with redirect_stderr(io.StringIO()):
                    self.assertEqual(COLLATOR.main([action]), 1)
                materialize.assert_not_called()
                seal.assert_not_called()
                run.assert_not_called()

    def test_complete_07_preserves_all_rows_and_numerical_failures(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            pairs = [
                write_pair(root, "A", [
                    placo_row("A", "rs1", 100, "1e-9", "0.01"),
                    placo_row(
                        "A", "rs2", 200, "1", "1",
                        status="NUMERICAL_FAILURE_P_SET_TO_ONE",
                    ),
                ]),
                write_pair(root, "B", [placo_row("B", "rs3", 300, "1e-4", "0.04")]),
                write_pair(
                    root, "CONTROL", [placo_row("CONTROL", "rs4", 400, "0.5", "0.7")],
                ),
            ]
            database = root / "eligible.sqlite"
            connection = COLLATOR.create_index_database(database)
            counts, projected = COLLATOR.index_complete_placo_family(root, pairs, connection)
            output = root / "07.tsv"
            identity = COLLATOR.write_complete_07(root, pairs, output, projected)
            with output.open(newline="", encoding="utf-8") as handle:
                rows = list(csv.DictReader(handle, delimiter="\t"))
            self.assertEqual([row["pair_id"] for row in rows], ["A", "A", "B", "CONTROL"])
            self.assertEqual(rows[1]["analysis_status"], "NUMERICAL_FAILURE_P_SET_TO_ONE")
            self.assertEqual(rows[1]["P_PLACO_PLUS"], "1")
            self.assertEqual(rows[1]["PLACO_BH_Q"], "1")
            self.assertEqual(rows[1]["numerical_error"], "synthetic_underflow")
            self.assertEqual(rows[-1]["family_role"], "POSITIVE_CONTROL")
            self.assertEqual(identity["bytes"], projected)
            self.assertEqual(counts["A"]["rows"], 2)
            self.assertEqual(connection.execute("SELECT COUNT(*) FROM eligible").fetchone()[0], 2)
            connection.close()

    def test_complete_family_denominator_mismatch_fails(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            pair = write_pair(root, "A", [placo_row("A", "rs1", 100, "0.5", "0.5")])
            with gzip.open(pair.ledger, "rt", encoding="utf-8") as handle:
                lines = handle.read().replace("\t1\tTESTED", "\t2\tTESTED")
            with gzip.open(pair.ledger, "wt", encoding="utf-8") as handle:
                handle.write(lines)
            changed = pair._replace(ledger_identity=COLLATOR.portable_identity(root, pair.ledger.relative_to(root)))
            provenance = json.loads(pair.provenance.read_text())
            provenance["output_sha256"] = changed.ledger_identity["sha256"]
            provenance["output_bytes"] = changed.ledger_identity["bytes"]
            pair.provenance.write_text(json.dumps(provenance) + "\n")
            changed = changed._replace(
                provenance_identity=COLLATOR.portable_identity(root, pair.provenance.relative_to(root)),
            )
            with self.assertRaisesRegex(COLLATOR.CollationError, "variant identity"):
                COLLATOR._stream_pair_rows(root, changed)

    def test_exact_reference_match_retains_every_failure_class(self) -> None:
        self.assertEqual(
            COLLATOR.exact_reference_match(1, 10, "A", "C", [])[0], "REFERENCE_ABSENT",
        )
        self.assertEqual(
            COLLATOR.exact_reference_match(1, 10, "A", "C", [(1, 11, "A", "C")])[0],
            "REFERENCE_COORDINATE_MISMATCH",
        )
        self.assertEqual(
            COLLATOR.exact_reference_match(1, 10, "A", "C", [(1, 10, "A", "G")])[0],
            "REFERENCE_ALLELE_MISMATCH",
        )
        self.assertEqual(
            COLLATOR.exact_reference_match(1, 10, "A", "C", [(1, 10, "T", "A")])[0],
            "REFERENCE_ALLELES_INVALID_OR_PALINDROMIC",
        )
        self.assertEqual(
            COLLATOR.exact_reference_match(1, 10, "A", "C", [(1, 10, "C", "A")])[0],
            "MATCHED",
        )
        self.assertEqual(
            COLLATOR.exact_reference_match(
                1, 10, "A", "C", [(1, 10, "A", "C"), (1, 10, "C", "A")],
            )[0],
            "REFERENCE_RSID_AMBIGUOUS",
        )

    def test_bim_scan_preserves_absent_coordinate_and_allele_mismatch_rows(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            database = root / "eligible.sqlite"
            connection = COLLATOR.create_index_database(database)
            for snp, position in (("rsMatch", 10), ("rsAbsent", 20), ("rsCoord", 30), ("rsAllele", 40)):
                insert_eligible(connection, "A", snp, position)
            bim = root / "ref.bim"
            bim.write_text(
                "1 rsMatch 0 10 C A\n1 rsCoord 0 31 A C\n1 rsAllele 0 40 A G\n",
                encoding="utf-8",
            )
            observed = COLLATOR.match_eligible_to_bim(root, connection, "ref.bim")
            self.assertEqual(observed["eligible_rows"], 4)
            self.assertEqual(observed["MATCHED"], 1)
            self.assertEqual(observed["REFERENCE_ABSENT"], 1)
            self.assertEqual(observed["REFERENCE_COORDINATE_MISMATCH"], 1)
            self.assertEqual(observed["REFERENCE_ALLELE_MISMATCH"], 1)
            self.assertEqual(
                connection.execute(
                    "SELECT COUNT(*) FROM eligible WHERE failure_reason!='NONE'"
                ).fetchone()[0],
                3,
            )
            connection.close()

    def test_real_plink_output_parser_and_strict_pairwise_r2(self) -> None:
        clumped = COLLATOR.parse_plink_clumped_text(
            "CHR F SNP BP P TOTAL NSIG S05 S01 S001 S0001 SP2\n"
            "1 1 rs1 100 1e-9 2 0 0 0 1 1 rs2(1)\n"
        )
        self.assertEqual(clumped[0]["members"], ["rs1", "rs2"])
        passed = COLLATOR.parse_pairwise_r2_text(
            "CHR_A BP_A SNP_A CHR_B BP_B SNP_B R2\n1 100 rs1 1 200 rs2 0.099999\n",
            {"rs1": (1, 100), "rs2": (1, 200)},
        )
        self.assertEqual(passed["expected_pairs"], 1)
        with self.assertRaisesRegex(COLLATOR.CollationError, "strict r2<0.1"):
            COLLATOR.parse_pairwise_r2_text(
                "CHR_A BP_A SNP_A CHR_B BP_B SNP_B R2\n1 100 rs1 1 200 rs2 0.1\n",
                {"rs1": (1, 100), "rs2": (1, 200)},
            )
        with self.assertRaisesRegex(COLLATOR.CollationError, "omitted"):
            COLLATOR.parse_pairwise_r2_text(
                "CHR_A BP_A SNP_A CHR_B BP_B SNP_B R2\n",
                {"rs1": (1, 100), "rs2": (1, 200)},
            )

    def test_mocked_plink_commands_use_actual_p_and_pairwise_audit(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            connection = COLLATOR.create_index_database(root / "eligible.sqlite")
            insert_eligible(connection, "A", "rs1", 100, "1e-9")
            insert_eligible(connection, "A", "rs2", 200, "2e-9")

            def runner(argv, _cwd):
                out = Path(argv[argv.index("--out") + 1])
                if "--clump" in argv:
                    out.with_suffix(".clumped").write_text(
                        "CHR F SNP BP P TOTAL NSIG S05 S01 S001 S0001 SP2\n"
                        "1 1 rs1 100 1e-9 1 0 0 0 0 0 NONE\n"
                        "1 1 rs2 200 2e-9 1 0 0 0 0 0 NONE\n",
                        encoding="utf-8",
                    )
                else:
                    out.with_suffix(".ld").write_text(
                        "CHR_A BP_A SNP_A CHR_B BP_B SNP_B R2\n"
                        "1 100 rs1 1 200 rs2 0.09\n",
                        encoding="utf-8",
                    )
                return subprocess.CompletedProcess(argv, 0, stdout="ok", stderr="")

            result = COLLATOR.run_plink_clump_for_pair(
                connection, "A", root / "plink", root / "g1000_eur",
                root / "run", root=root, runner=runner,
            )
            self.assertEqual(result["leads"], ["rs1", "rs2"])
            self.assertEqual(result["pairwise"]["maximum_r2"], 0.09)
            command = result["commands"][0]["argv"]
            self.assertEqual(command[command.index("--clump-kb") + 1], "1000")
            self.assertEqual(command[command.index("--clump-r2") + 1], "0.1")
            self.assertEqual(command[command.index("--clump-p1") + 1], "1")
            statuses = connection.execute(
                "SELECT DISTINCT clump_status FROM eligible ORDER BY clump_status"
            ).fetchall()
            self.assertEqual(statuses, [("INDEPENDENT_LEAD_VERIFIED_R2_LT_0.1",)])
            connection.close()

    def test_mocked_plink_rejects_lead_not_ordered_by_actual_p(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            connection = COLLATOR.create_index_database(root / "eligible.sqlite")
            insert_eligible(connection, "A", "rs1", 100, "1e-9")
            insert_eligible(connection, "A", "rs2", 200, "2e-9")

            def runner(argv, _cwd):
                out = Path(argv[argv.index("--out") + 1])
                out.with_suffix(".clumped").write_text(
                    "CHR F SNP BP P TOTAL NSIG S05 S01 S001 S0001 SP2\n"
                    "1 1 rs2 200 2e-9 2 0 0 0 1 1 rs1(1)\n",
                    encoding="utf-8",
                )
                return subprocess.CompletedProcess(argv, 0, stdout="ok", stderr="")

            with self.assertRaisesRegex(COLLATOR.CollationError, "actual PLACO P"):
                COLLATOR.run_plink_clump_for_pair(
                    connection, "A", root / "plink", root / "g1000_eur",
                    root / "run", root=root, runner=runner,
                )
            connection.close()

    def test_multiple_independent_leads_in_one_block_and_no_cross_pair_collapse(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            connection = COLLATOR.create_index_database(root / "eligible.sqlite")
            insert_eligible(connection, "A", "rs1", 100)
            insert_eligible(connection, "A", "rs2", 200)
            insert_eligible(connection, "B", "rs1", 100)
            insert_eligible(connection, "CONTROL", "rsC", 300)
            _, block_index = COLLATOR.parse_blocks("LOC CHR START STOP\n1 1 1 1000\n")
            pairwise = {
                "expected_pairs": 1, "observed_pairs": 1, "maximum_r2": 0.02,
                "status": "PASS_STRICT_R2_LT_0.1_WITHIN_1_MB",
            }
            singleton = {
                "expected_pairs": 0, "observed_pairs": 0, "maximum_r2": None,
                "status": "PASS_NO_WITHIN_WINDOW_LEAD_PAIR",
            }
            results = {
                "A": {
                    "leads": ["rs1", "rs2"],
                    "clumps": [
                        {"lead": "rs1", "members": ["rs1"]},
                        {"lead": "rs2", "members": ["rs2"]},
                    ],
                    "pairwise": pairwise,
                },
                "B": {"leads": ["rs1"], "clumps": [{"lead": "rs1", "members": ["rs1"]}], "pairwise": singleton},
                "CONTROL": {"leads": ["rsC"], "clumps": [{"lead": "rsC", "members": ["rsC"]}], "pairwise": singleton},
            }
            lava = {
                (pair, "1"): {
                    "locus_status": "FAILED_QC", "bivariate_status": "NOT_TESTED",
                    "conditional_failures": "FAILED_QC",
                }
                for pair in COLLATOR.PAIR_ORDER
            }
            primary, control = COLLATOR.build_shared_locus_rows(
                connection, results, block_index, {"1": 1000}, lava,
            )
            self.assertEqual(len(primary), 3)
            self.assertEqual(len(control), 1)
            a_rows = [row for row in primary if row["pair_id"] == "A"]
            self.assertEqual(len(a_rows), 2)
            self.assertEqual({row["lava_block_id"] for row in a_rows}, {"LOC1"})
            self.assertTrue(all(len(row["pair_block_evidence_ids"].split(";")) == 2 for row in a_rows))
            self.assertEqual(len({row["clump_id"] for row in [*primary, *control]}), 4)
            self.assertTrue(all(row["family_role"] == "PRIMARY_DISCOVERY" for row in primary))
            self.assertEqual(control[0]["primary_headline"], "NOT_APPLICABLE_CONTROL")
            connection.close()

    def test_blocked_conjfdr_never_emits_only_or_both_labels(self) -> None:
        shared = [{field: "x" for field in COLLATOR.SHARED_LOCI_08_FIELDS}]
        shared[0].update({
            "analysis_id": "track-b-v1.0-pleiotropy", "pair_id": "A",
            "family_role": "PRIMARY_DISCOVERY", "clump_id": "PLACO:A:rs1",
            "chr": "1", "locus_start": "1", "locus_end": "2", "lead_variant": "rs1",
            "PLACO_P": "1e-9", "FDR": "0.01", "evidence_labels": "PLACO_PRIMARY_HEADLINE",
            "method_availability": f"PLACO:COMPLETE_WITH_HITS;CONJFDR:{COLLATOR.CONJFDR_BLOCKED}",
            "claim_limit": "NO_MEDIATION",
        })
        row = list(COLLATOR.comparison_rows(shared))[0]
        self.assertEqual(row["comparison_label"], "NOT_EMITTED_OTHER_METHOD_NOT_TESTED")
        self.assertEqual(row["conjfdr"], "NA")
        self.assertNotIn(row["comparison_label"], {"PLACO_ONLY", "CONJFDR_ONLY", "PLACO_AND_CONJFDR"})

    def test_zero_and_ld_blocked_mappings_have_no_header_only_science_tables(self) -> None:
        zero = COLLATOR.package_source_map("ZERO_ELIGIBLE_SIGNAL_FAMILY")
        blocked = COLLATOR.package_source_map("BLOCKED_BY_LD_REFERENCE_COVERAGE")
        for mapping in (zero, blocked):
            self.assertNotIn(COLLATOR.SHARED_LOCI_08_PATH, mapping)
            self.assertNotIn(COLLATOR.COMPARISON_09_PATH, mapping)
            self.assertNotIn(COLLATOR.CONTROL_08_PATH, mapping)
            self.assertNotIn(COLLATOR.CONTROL_09_PATH, mapping)
        self.assertIn(COLLATOR.ZERO_FAMILY_PATH, zero)
        self.assertIn(COLLATOR.BLOCKED_COVERAGE_PATH, blocked)

    def test_partial_publication_is_prefix_only_and_crash_resumes_without_replace(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            fingerprint = "a" * 64
            package, mapping = fake_zero_package(root, fingerprint)
            with self.assertRaisesRegex(COLLATOR.CollationError, "SIMULATED_CRASH"):
                COLLATOR.publish_package(root, fingerprint, crash_after_links=2)
            order = COLLATOR.publication_order(mapping)
            self.assertTrue(all((root / path).exists() for path in order[:2]))
            self.assertTrue(all(not (root / path).exists() for path in order[2:]))
            result = COLLATOR.publish_package(root, fingerprint)
            self.assertEqual(result["linked_this_call"], len(order) - 2)
            for canonical, packaged in mapping.items():
                self.assertTrue(COLLATOR._same_regular_inode(root / canonical, package / packaged))

    def test_no_replace_rejects_same_bytes_different_inode_and_symlink(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            fingerprint = "b" * 64
            package, mapping = fake_zero_package(root, fingerprint)
            first = COLLATOR.publication_order(mapping)[0]
            destination = root / first
            destination.parent.mkdir(parents=True, exist_ok=True)
            destination.write_bytes((package / mapping[first]).read_bytes())
            with self.assertRaisesRegex(COLLATOR.CollationError, "not the immutable package inode"):
                COLLATOR.publish_package(root, fingerprint)
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            fingerprint = "c" * 64
            package, mapping = fake_zero_package(root, fingerprint)
            first = COLLATOR.publication_order(mapping)[0]
            destination = root / first
            destination.parent.mkdir(parents=True, exist_ok=True)
            destination.symlink_to(package / mapping[first])
            with self.assertRaisesRegex(COLLATOR.CollationError, "not the immutable package inode"):
                COLLATOR.publish_package(root, fingerprint)

    def test_own_result_verifier_allows_later_authorized_finemapping_outputs(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            fingerprint = "d" * 64
            fake_zero_package(root, fingerprint)
            COLLATOR.publish_package(root, fingerprint)
            for relative in COLLATOR.COMPETING_PATHS:
                path = root / relative
                path.parent.mkdir(parents=True, exist_ok=True)
                path.write_text("authorized downstream artifact\n", encoding="utf-8")
            manifest = COLLATOR.verify_package(root, fingerprint)
            mapping = COLLATOR.verify_published_link_family(root, fingerprint, manifest)
            self.assertIn(COLLATOR.RESULT_PROVENANCE_PATH, mapping)
            with self.assertRaisesRegex(COLLATOR.CollationError, "competing fine-mapping"):
                COLLATOR.assert_no_competing_finemapping(root)

    def test_deep_zero_family_verifier_remains_valid_after_downstream_publication(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            pairs = [
                write_pair(root, pair, [placo_row(pair, f"rs{index}", index, "0.5", "0.5")])
                for index, pair in enumerate(COLLATOR.PAIR_ORDER, start=1)
            ]
            contract = synthetic_contract(root)
            archive = {"state": COLLATOR.ARCHIVES_PRESENT, "synthetic": "archive"}
            conjfdr = {
                "family_state": COLLATOR.CONJFDR_BLOCKED,
                "pair_states": {pair: COLLATOR.CONJFDR_BLOCKED for pair in COLLATOR.PAIR_ORDER},
                "evidence_rows": 0, "semantics": "blocked", "frozen_readiness": {},
            }
            runtime_audit = {
                "reference": {
                    suffix: contract["ld_reference"][suffix] for suffix in ("bed", "bim", "fam")
                },
                "plink": {
                    key: contract["ld_reference"]["plink"][key]
                    for key in ("path", "bytes", "sha256")
                },
                "snapshot_mode": (
                    "HARDLINKED_EXACT_INODES_BEFORE_BIM_SCAN_AND_PLINK;FULL_HASH_AND_"
                    "INODE_MTIME_CTIME_RECHECK_AFTER_ALL_COMMANDS"
                ),
            }
            frozen = {
                "reference": {"bim": {"snapshot": {"path": "unused", "bytes": 1, "sha256": "1" * 64}}},
                "plink_path": root / "unused-plink", "reference_prefix": root / "unused-ref",
            }
            fingerprint = "f" * 64
            coverage = {"eligible_rows": 0, "bim_rows_streamed": 1, "bloom_candidate_rows": 0}
            blocks = ([COLLATOR.Block(1, 1, 1, 1000)], {1: ([1], [COLLATOR.Block(1, 1, 1, 1000)])})
            patches = (
                mock.patch.object(COLLATOR, "archive_state_invariant", return_value={"synthetic": "invariant"}),
                mock.patch.object(COLLATOR, "freeze_runtime_family", return_value=frozen),
                mock.patch.object(COLLATOR, "verify_and_release_runtime_family", return_value=runtime_audit),
                mock.patch.object(COLLATOR, "match_eligible_to_bim", return_value=coverage),
                mock.patch.object(COLLATOR, "load_official_blocks", return_value=blocks),
                mock.patch.object(COLLATOR, "load_lava_diagnostic_index", return_value={}),
            )
            with patches[0], patches[1], patches[2], patches[3], patches[4], patches[5]:
                COLLATOR._build_package_from_verified_inputs(
                    root, contract, pairs, fingerprint, archive, conjfdr, {},
                    runner=lambda argv, cwd: subprocess.CompletedProcess(argv, 0, "", ""),
                )
                COLLATOR.publish_package(root, fingerprint)
                verify_patches = (
                    mock.patch.object(COLLATOR, "verify_pre_result_contract", return_value=contract),
                    mock.patch.object(COLLATOR, "_load_production_modules", return_value={"archive_guard": object()}),
                    mock.patch.object(COLLATOR, "verify_archive_state", return_value=archive),
                    mock.patch.object(COLLATOR, "validate_conjfdr_blocked", return_value=conjfdr),
                )
                with verify_patches[0], verify_patches[1], verify_patches[2], verify_patches[3]:
                    first = COLLATOR.verify_results(root)
                    self.assertTrue(first["zero_family"]["is_zero"])
                    for relative in COLLATOR.COMPETING_PATHS:
                        path = root / relative
                        path.parent.mkdir(parents=True, exist_ok=True)
                        path.write_text("authorized downstream\n", encoding="utf-8")
                    second = COLLATOR.verify_results(root)
                    self.assertEqual(second["fingerprint"], fingerprint)

    def test_special_files_and_package_symlinks_are_rejected(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            directory = root / "package"
            directory.mkdir()
            (directory / "real").write_text("x")
            (directory / "link").symlink_to(directory / "real")
            with self.assertRaisesRegex(COLLATOR.CollationError, "special file"):
                COLLATOR._inventory_files(directory)
            if hasattr(os, "mkfifo"):
                (directory / "link").unlink()
                os.mkfifo(directory / "fifo")
                with self.assertRaisesRegex(COLLATOR.CollationError, "special file"):
                    COLLATOR._inventory_files(directory)

    def test_same_size_named_inode_replacement_during_hash_is_detected(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            target = root / "large.bin"
            replacement = root / "replacement.bin"
            target.write_bytes(b"a" * (5 * 1024 * 1024))
            replacement.write_bytes(b"b" * (5 * 1024 * 1024))
            original_read = COLLATOR.os.read
            replaced = False

            def adversarial_read(descriptor, size):
                nonlocal replaced
                payload = original_read(descriptor, size)
                if payload and not replaced:
                    replaced = True
                    os.replace(replacement, target)
                return payload

            with (
                mock.patch.object(COLLATOR.os, "read", side_effect=adversarial_read),
                self.assertRaisesRegex(COLLATOR.CollationError, "changed while hashing"),
            ):
                COLLATOR.stable_identity(root, "large.bin")

    def test_runtime_reference_and_plink_are_frozen_as_exact_inodes(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            reference = root / "ref"
            reference.mkdir()
            contract = {"ld_reference": {}}
            for suffix, payload in (("bed", b"bed"), ("bim", b"1 rs1 0 1 A C\n"), ("fam", b"fam")):
                path = reference / f"g1000_eur.{suffix}"
                path.write_bytes(payload)
                contract["ld_reference"][suffix] = COLLATOR.portable_identity(
                    root, path.relative_to(root),
                )
            plink = reference / "plink"
            plink.write_text("#!/bin/sh\necho 'PLINK v1.9.0-b.7.11 synthetic'\n", encoding="utf-8")
            plink.chmod(0o755)
            contract["ld_reference"]["plink"] = {
                **COLLATOR.portable_identity(root, plink.relative_to(root)),
                "version": "PLINK v1.9.0-b.7.11",
            }
            attempt = root / "attempt"
            attempt.mkdir()
            frozen = COLLATOR.freeze_runtime_family(root, contract, attempt)
            self.assertEqual(
                os.stat(reference / "g1000_eur.bim").st_ino,
                os.stat(attempt / "runtime_reference/g1000_eur.bim").st_ino,
            )
            audit = COLLATOR.verify_and_release_runtime_family(root, frozen)
            self.assertEqual(audit["reference"]["bim"], contract["ld_reference"]["bim"])
            self.assertFalse((attempt / "runtime_reference").exists())
            self.assertFalse((attempt / "runtime_software").exists())

    def test_archive_guard_accepts_exactly_two_states_and_consumer_audit_passes(self) -> None:
        base = {
            "extracted_payload_count": 44, "archive_family_sha256": "a" * 64,
            "extracted_family_sha256": "b" * 64,
            "reference_provenance_identity": {"path": "x", "bytes": 1, "sha256": "c" * 64},
        }

        class Guard:
            def __init__(self, state):
                self.state = state

            def verify_reference_state(self, _root):
                return {**base, "state": self.state}

        for state in (COLLATOR.ARCHIVES_PRESENT, COLLATOR.ARCHIVES_EVICTED):
            self.assertEqual(COLLATOR.verify_archive_state(Guard(state), ROOT)["state"], state)
        with self.assertRaisesRegex(COLLATOR.CollationError, "two-state"):
            COLLATOR.verify_archive_state(Guard("PARTIAL_ARCHIVE_EVICTION_RESTART_REQUIRED"), ROOT)
        guard_path = ROOT / "scripts/146_manage_lava_reference_archives.py"
        specification = importlib.util.spec_from_file_location("archive_guard_for_147_test", guard_path)
        guard = importlib.util.module_from_spec(specification)
        specification.loader.exec_module(guard)
        self.assertTrue(
            guard._has_executable_guard_integration(
                Path("scripts/147_build_track_b_pleiotropy_results_v2.py"),
                SCRIPT.read_text(encoding="utf-8"),
            )
        )

    def test_schemas_paths_and_ram_namespace_are_additive(self) -> None:
        policy = json.loads((ROOT / "config/track_b_pleiotropy_policy.json").read_text())
        self.assertEqual(
            COLLATOR.SHARED_LOCI_08_FIELDS[:len(COLLATOR.SHARED_LOCI_08_PREFIX)],
            policy["required_future_result_schemas"]["shared_loci_08_required_prefix"],
        )
        self.assertEqual(
            COLLATOR.RAM_FIELDS,
            [
                "analysis", "pair", "locus", "chromosome", "n_snps", "peak_ram_gb",
                "runtime_sec", "exit_status", "output_hash",
            ],
        )
        self.assertNotEqual(COLLATOR.RAM_COMPONENT_PATH, Path("results/track_b/RAM_BENCHMARK.tsv"))
        source = SCRIPT.read_text(encoding="utf-8")
        self.assertNotIn("FINEMAPPING_GATE_SCRIPT", source)
        self.assertNotIn("scripts/145_build_track_b_finemapping_manifest.py", source)

    def test_disk_preflight_uses_exact_07_projection_and_explicit_allowance(self) -> None:
        expected = 1234 + 256 * 1024**2 + 7 * 4096
        observed = COLLATOR.disk_preflight(ROOT, 1234, 7, free_bytes=expected)
        self.assertEqual(observed["required_bytes"], expected)
        with self.assertRaisesRegex(COLLATOR.CollationError, "insufficient storage"):
            COLLATOR.disk_preflight(ROOT, 1234, 7, free_bytes=expected - 1)

    def test_rss_units_and_pgid_sampling_include_reparented_descendants(self) -> None:
        self.assertEqual(COLLATOR.normalize_ru_maxrss(123, platform="darwin"), 123)
        self.assertEqual(COLLATOR.normalize_ru_maxrss(123, platform="linux"), 123 * 1024)
        sample = COLLATOR.parse_ps_rss_snapshot(
            "10 5 100\n20 77 200\n21 77 300\n", process_group_id=77, parent_pid=10,
        )
        self.assertEqual(sample["process_group_bytes"], 500 * 1024)
        self.assertEqual(sample["parent_bytes"], 100 * 1024)
        self.assertEqual(sample["aggregate_bytes"], 600 * 1024)

    def test_default_runner_measures_live_parent_plus_owned_process_group(self) -> None:
        sample = {
            "process_group_bytes": 4_000_000, "parent_bytes": 5_000_000,
            "aggregate_bytes": 9_000_000,
        }
        with (
            tempfile.TemporaryDirectory() as temporary,
            mock.patch.object(COLLATOR, "process_group_and_parent_rss_bytes", return_value=sample),
            mock.patch.object(COLLATOR, "process_group_rss_bytes", return_value=0),
        ):
            result = COLLATOR._default_runner(
                [sys.executable, "-c", "import time; x=bytearray(4_000_000); time.sleep(.15)"],
                Path(temporary),
            )
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertGreater(result.peak_aggregate_rss_bytes, 0)
        self.assertIn(
            result.rss_measurement_method,
            {
                "LIVE_PARENT_PLUS_OWNED_FRESH_SESSION_PGID_PS_SAMPLED",
                "NORMALIZED_SELF_PLUS_CHILD_RUSAGE_CONSERVATIVE_FAST_EXIT_FALLBACK",
            },
        )

    def test_default_runner_rejects_owned_pgid_that_outlives_leader(self) -> None:
        sample = {
            "process_group_bytes": 1_000_000, "parent_bytes": 2_000_000,
            "aggregate_bytes": 3_000_000,
        }
        with (
            tempfile.TemporaryDirectory() as temporary,
            mock.patch.object(COLLATOR, "process_group_and_parent_rss_bytes", return_value=sample),
            mock.patch.object(COLLATOR, "process_group_rss_bytes", return_value=1_000_000),
            mock.patch.object(COLLATOR, "terminate_owned_process_group", return_value=True),
        ):
            result = COLLATOR._default_runner([sys.executable, "-c", "pass"], Path(temporary))
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("process group outlived its leader", result.stderr)


if __name__ == "__main__":
    unittest.main()
