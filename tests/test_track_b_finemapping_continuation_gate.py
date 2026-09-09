import copy
import importlib.util
import json
import tempfile
import unittest
from contextlib import redirect_stderr
from io import StringIO
from pathlib import Path
from unittest import mock


ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "scripts/145_build_track_b_finemapping_continuation_gate.py"


def load_module():
    specification = importlib.util.spec_from_file_location(
        "track_b_finemapping_continuation_gate_test", SCRIPT,
    )
    assert specification and specification.loader
    module = importlib.util.module_from_spec(specification)
    specification.loader.exec_module(module)
    return module


GATE = load_module()


def identity(path="x", size=1, digest="a" * 64):
    return {"path": path, "bytes": size, "sha256": digest}


def block_fixture():
    return GATE.parse_blocks_text(
        "\n".join([
            "LOC CHR START STOP",
            "1 1 1 100",
            "2 1 101 200",
            "3 1 301 400",
            "4 1 501 600",
            "5 1 701 800",
            "6 2 1 100",
            "",
        ]),
        expected_count=6,
    )


def make_lead(
    pair="A", snp="rs10", chromosome=1, position=10, p="1e-9", q="0.2",
    members=None, pair_block_leads="PENDING",
):
    _, index = block_fixture()
    block = GATE.map_to_block(chromosome, position, index)
    lead_evidence = f"PLACO:{pair}:{snp}:{chromosome}:{position}"
    members = list(members or [lead_evidence])
    trait1, trait2, _ = GATE.PAIR_IDENTITIES[pair]
    upstream_role = GATE.UPSTREAM_PAIR_ROLES[pair]
    p_float, q_float = float(p), float(q)
    if block is None:
        block_id, block_start, block_end = "UNMAPPED", "NA", "NA"
        lava_accounting = "TERMINAL_FAILED_QC_ACCOUNTING_ONLY;OFFICIAL_BLOCK=UNMAPPED"
        conditional = "BLOCKED_BY_OFFICIAL_LAVA_PARTITION_COVERAGE"
    else:
        block_id = f"LOC{block.locus}"
        block_start, block_end = str(block.start), str(block.stop)
        lava_accounting = (
            "TERMINAL_FAILED_QC_ACCOUNTING_ONLY;LOCUS_STATUS=FAILED_QC;"
            "BIVARIATE_STATUS=FAILED_QC;SCIENTIFIC_SUPPORT=FORBIDDEN"
        )
        conditional = "TERMINAL_FAILED_QC_ACCOUNTING_ONLY;FAILURE_STATUSES=MODEL_FAILED"
    return {
        "pair_id": pair,
        "lead_variant": snp,
        "chr": str(chromosome),
        "position": str(position),
        "PLACO_P": p,
        "FDR": q,
        "trait1_P": "0.01",
        "trait2_P": "0.02",
        "locus_start": "1",
        "locus_end": str(position + 1_000_000),
        "independent_signal": "TRUE",
        "annotations": (
            "STATISTICAL_CROSS_TRAIT_PLEIOTROPY_ONLY;"
            "NO_SHARED_CAUSAL_VARIANT_NO_MEDIATION_NO_CAUSAL_DIRECTION"
        ),
        "analysis_id": GATE.UPSTREAM_ANALYSIS_ID,
        "family_role": upstream_role,
        "trait1": trait1,
        "trait2": trait2,
        "lead_A1": "A",
        "lead_A2": "C",
        "lead_Z1": "3",
        "lead_Z2": "4",
        "evidence_labels": ";".join(GATE._expected_labels(pair, p_float, q_float)),
        "primary_headline": (
            "NOT_APPLICABLE_CONTROL" if pair == "CONTROL"
            else str(p_float <= 2.5e-8).upper()
        ),
        "eligibility_basis": ";".join(GATE._expected_eligibility(p_float, q_float)),
        "clump_id": f"PLACO:{pair}:{snp}",
        "clump_evidence_ids": ";".join(members),
        "lava_block_id": block_id,
        "lava_block_start": block_start,
        "lava_block_end": block_end,
        "pair_block_evidence_ids": pair_block_leads,
        "ld_reference_status": "EXACT_RSID_CHR_BP_UNORDERED_NONPALINDROMIC_ALLELES_MATCHED",
        "plink_clump_status": "INDEPENDENT_LEAD_REAL_PLINK_1_9",
        "pairwise_r2_status": "PASS_NO_WITHIN_WINDOW_LEAD_PAIR",
        "max_pairwise_r2_within_1mb": "NA",
        "lava_terminal_state": "TERMINAL_FAILED_QC_DIAGNOSTIC_ONLY",
        "lava_failed_qc_accounting": lava_accounting,
        "conditional_status": conditional,
        "method_availability": GATE.METHOD_AVAILABILITY,
        "claim_limit": GATE.CONTROL_CLAIM if pair == "CONTROL" else GATE.PRIMARY_CLAIM,
    }


def complete_family(*, has_hits=True):
    pairs = {}
    total = 0
    for offset, pair in enumerate(GATE.PAIR_ORDER, start=1):
        rows = 100 + offset
        total += rows
        terminal = "COMPLETE_WITH_HITS" if has_hits else "TESTED_NO_HIT"
        pairs[pair] = {
            "terminal_state": terminal,
            "rows": rows,
            "ledger": identity(f"upstream/{pair}.tsv.gz", rows, str(offset) * 64),
            "provenance": identity(f"upstream/{pair}.json", 10, str(offset + 3) * 64),
            "counts": {
                "rows": rows,
                "failures": 0,
                "primary": 1 if has_hits and pair != "CONTROL" else 0,
                "pairwise": 1 if has_hits else 0,
                "bh": 0,
                "minimum_p": 1e-9 if has_hits else 1.0,
                "terminal_status": terminal,
            },
        }
    return {
        "order": list(GATE.PAIR_ORDER),
        "rows": total,
        "canonical_07": identity(
            "results/track_b/07_placo_plus_variants.tsv", 1234, "f" * 64,
        ),
        "pairs": pairs,
    }


def upstream_payload(leads=(), state=GATE.UPSTREAM_COMPLETE, coverage_failures=0):
    leads = [copy.deepcopy(row) for row in leads]
    grouped = {}
    for row in leads:
        grouped.setdefault((row["pair_id"], row["lava_block_id"]), []).append(row)
    for family in grouped.values():
        lead_ids = ";".join(
            f"PLACO:{row['pair_id']}:{row['lead_variant']}:{row['chr']}:{row['position']}"
            for row in family
        )
        for row in family:
            row["pair_block_evidence_ids"] = lead_ids
    leads.sort(
        key=lambda row: (
            GATE.PAIR_RANK[row["pair_id"]], int(row["chr"]), int(row["position"]),
            row["lead_variant"],
        )
    )
    evidence_count = sum(len(row["clump_evidence_ids"].split(";")) for row in leads)
    if state == GATE.UPSTREAM_ZERO:
        eligible, counts, failures = 0, {}, 0
        family = complete_family(has_hits=False)
    elif state == GATE.UPSTREAM_BLOCKED:
        eligible, counts, failures = 2, {"NO_REFERENCE_MATCH": 2}, 2
        family = complete_family(has_hits=True)
    else:
        eligible = evidence_count + coverage_failures
        counts = {"MATCHED": evidence_count}
        if coverage_failures:
            counts["NO_REFERENCE_MATCH"] = coverage_failures
        failures = coverage_failures
        family = complete_family(has_hits=True)
        lead_pairs = {row["pair_id"] for row in leads}
        for pair in set(GATE.PAIR_ORDER) - lead_pairs:
            family["pairs"][pair]["terminal_state"] = "TESTED_NO_HIT"
            family["pairs"][pair]["counts"].update({
                "primary": 0, "pairwise": 0, "bh": 0, "minimum_p": 1.0,
                "terminal_status": "TESTED_NO_HIT",
            })
    archive = {
        "state": GATE.ARCHIVES_PRESENT,
        "extracted_payload_count": 44,
        "archive_family_sha256": "8" * 64,
        "extracted_family_sha256": "9" * 64,
    }
    return {
        "schema_version": GATE.UPSTREAM_VERIFIER_SCHEMA,
        "state": state,
        "fingerprint": "e" * 64,
        "zero_family": {"is_zero": state == GATE.UPSTREAM_ZERO, "state": state},
        "primary_leads": [row for row in leads if row["pair_id"] != "CONTROL"],
        "control_leads": [row for row in leads if row["pair_id"] == "CONTROL"],
        "complete_family": family,
        "archive_state": GATE.ARCHIVES_PRESENT,
        "archive_reference": archive,
        "ld_reference_coverage": {
            "eligible_rows": eligible,
            "status_counts": counts,
            "coverage_failure_rows": failures,
        },
        "method_availability": {
            "PLACO": "COMPLETE", "CONJFDR": GATE.CONJFDR_BLOCKED,
        },
    }


def validate(payload):
    _, index = block_fixture()
    return GATE.validate_upstream_payload(payload, index)


class TrackBFineMappingContinuationGateTests(unittest.TestCase):
    def test_exact_public_schemas_and_pinned_147_api(self):
        self.assertEqual(GATE.FAMILY_FIELDS, [
            "analysis_id", "locus_entry_id", "pair_id", "family_role", "trait1", "trait2",
            "CHR", "START", "STOP", "ld_block_id", "inclusion_sources", "priority_tier",
            "local_evidence_ids", "local_lead_variants", "pleiotropy_lead_evidence_ids",
            "pleiotropy_lead_variants", "pleiotropy_clump_evidence_ids",
            "method_availability", "partition_role", "selection_rule", "claim_limit",
        ])
        self.assertEqual(GATE.UNAVAILABLE_FIELDS, [
            "analysis_id", "unavailable_entry_id", "pair_id", "family_role", "trait1",
            "trait2", "CHR", "BP", "lead_variant", "ld_block_id", "evidence_ids",
            "inclusion_sources", "terminal_status", "blocker_reason", "claim_limit",
        ])
        self.assertEqual(
            GATE.UPSTREAM_VERIFIER_SHA256,
            "cc0cbb0b9b3f8c6045b02367a4d4d4ed2a8b21804f376318efce429ea26518a5",
        )
        source = SCRIPT.read_text(encoding="utf-8")
        self.assertIn("verifier.verify_results(root=root)", source)
        self.assertNotIn("143_run_track_b_placo", source)
        self.assertNotIn("placo_pair_manifest", source)

    def test_live_policy_and_official_partition_retain_whole_locus_contract(self):
        policy, _ = GATE.validate_policy(ROOT)
        self.assertIsNone(policy["dense_input_contract"]["maximum_locus_variants"])
        self.assertTrue(
            policy["ram_aware_execution_contract"]["full_locus_variant_universe_required"]
        )
        blocks, index, observed = GATE.load_official_blocks(ROOT)
        self.assertEqual(len(blocks), 2495)
        self.assertEqual(observed["sha256"], GATE.OFFICIAL_BLOCK_SHA256)
        self.assertEqual(GATE.map_to_block(1, blocks[0].start, index), blocks[0])

    def test_partition_rejects_overlap_and_returns_none_for_unmapped_gap(self):
        with self.assertRaisesRegex(GATE.GateError, "overlap"):
            GATE.parse_blocks_text(
                "LOC CHR START STOP\n1 1 1 100\n2 1 100 200\n", expected_count=2,
            )
        _, index = block_fixture()
        self.assertIsNone(GATE.map_to_block(1, 250, index))

    def test_pair_block_group_preserves_parallel_leads_and_full_clump_union(self):
        first_id = "PLACO:A:rs10:1:10"
        second_id = "PLACO:A:rs20:1:20"
        first = make_lead(
            snp="rs10", position=10,
            members=[first_id, "PLACO:A:rs11:1:11"],
        )
        second = make_lead(
            snp="rs20", position=20,
            members=[second_id, "PLACO:A:rs21:1:21", "PLACO:A:rs22:1:22"],
        )
        built = GATE.build_family_rows(validate(upstream_payload([first, second])))
        self.assertEqual(len(built["mapped_rows"]), 1)
        row = built["mapped_rows"][0]
        self.assertEqual(row["locus_entry_id"], "TBFM_A_LOC1")
        self.assertEqual(row["local_evidence_ids"], "NA")
        self.assertEqual(row["local_lead_variants"], "NA")
        self.assertEqual(row["pleiotropy_lead_evidence_ids"], f"{first_id};{second_id}")
        self.assertEqual(row["pleiotropy_lead_variants"], "rs10;rs20")
        self.assertEqual(
            row["pleiotropy_clump_evidence_ids"],
            ";".join([
                first_id, "PLACO:A:rs11:1:11", second_id,
                "PLACO:A:rs21:1:21", "PLACO:A:rs22:1:22",
            ]),
        )
        self.assertEqual(len(built["upstream_evidence_ids"]), 5)
        self.assertEqual(row["partition_role"], GATE.PARTITION_ROLE)
        self.assertNotIn("LAVA", row["inclusion_sources"])

    def test_same_official_block_remains_distinct_across_pairs_and_control(self):
        leads = [
            make_lead(pair="A", snp="rs10", position=10),
            make_lead(pair="B", snp="rs20", position=20),
            make_lead(pair="CONTROL", snp="rs30", position=30),
        ]
        built = GATE.build_family_rows(validate(upstream_payload(leads)))
        self.assertEqual(
            [row["locus_entry_id"] for row in built["mapped_rows"]],
            ["TBFM_A_LOC1", "TBFM_B_LOC1", "TBFM_CONTROL_LOC1"],
        )
        self.assertEqual(built["mapped_rows"][-1]["priority_tier"], "CONTROL_TIER")
        self.assertEqual(
            built["mapped_rows"][-1]["family_role"], "POSITIVE_CONTROL_NON_NOVELTY",
        )

    def test_priority_is_top_three_per_primary_pair_but_does_not_filter_or_reorder(self):
        leads = [
            make_lead(snp="rs10", position=10, p="4e-8", q="0.04"),
            make_lead(snp="rs110", position=110, p="2e-8", q="0.9"),
            make_lead(snp="rs310", position=310, p="1e-8", q="0.9"),
            make_lead(snp="rs510", position=510, p="3e-8", q="0.01"),
            make_lead(snp="rs710", position=710, p="2e-8", q="0.01"),
        ]
        built = GATE.build_family_rows(validate(upstream_payload(leads)))
        rows = built["mapped_rows"]
        self.assertEqual([int(row["START"]) for row in rows], [1, 101, 301, 501, 701])
        self.assertEqual(len(rows), 5)
        primary = {row["ld_block_id"] for row in rows if row["priority_tier"] == "PRIMARY_TIER"}
        self.assertEqual(primary, {"LOC2", "LOC3", "LOC5"})
        self.assertEqual(
            {row["ld_block_id"] for row in rows if row["priority_tier"] == "SECONDARY_TIER"},
            {"LOC1", "LOC4"},
        )
        self.assertTrue(all("NON_FILTERING" in row["selection_rule"] for row in rows))

    def test_unmapped_lead_is_blocked_by_data_and_retains_full_evidence(self):
        evidence = ["PLACO:B:rs250:1:250", "PLACO:B:rs251:1:251"]
        lead = make_lead(
            pair="B", snp="rs250", position=250, p="4e-8", q="0.2", members=evidence,
        )
        built = GATE.build_family_rows(validate(upstream_payload([lead])))
        self.assertEqual(built["mapped_rows"], [])
        self.assertEqual(len(built["unavailable_rows"]), 1)
        row = built["unavailable_rows"][0]
        self.assertEqual(row["ld_block_id"], "UNMAPPED")
        self.assertEqual(row["terminal_status"], "BLOCKED_BY_DATA")
        self.assertEqual(row["blocker_reason"], "NO_PREDECLARED_SIGNED_LD_BLOCK")
        self.assertEqual(row["evidence_ids"], ";".join(evidence))
        self.assertNotIn("NULL", str(row))
        self.assertEqual(built["upstream_evidence_ids"], evidence)

    def test_zero_and_ld_coverage_blocked_states_are_distinct(self):
        zero = validate(upstream_payload([], state=GATE.UPSTREAM_ZERO))
        built = GATE.build_family_rows(zero)
        self.assertEqual(built["state"], GATE.ZERO_STATE)
        self.assertEqual(built["mapped_rows"], [])
        blocked = validate(upstream_payload([], state=GATE.UPSTREAM_BLOCKED))
        with self.assertRaisesRegex(GATE.UpstreamBlocked, "never zero"):
            GATE.build_family_rows(blocked)

    def test_complete_state_with_any_reference_failure_fails_closed(self):
        lead = make_lead()
        with self.assertRaisesRegex(GATE.UpstreamBlocked, "coverage failures"):
            validate(upstream_payload([lead], coverage_failures=1))

    def test_exact_147_rows_method_and_lava_accounting_fail_closed(self):
        payload = upstream_payload([make_lead()])
        changed = copy.deepcopy(payload)
        changed["primary_leads"][0]["method_availability"] = (
            "PLACO:COMPLETE_WITH_HITS;CONJFDR:COMPLETE_WITH_HITS"
        )
        with self.assertRaisesRegex(GATE.GateError, "method, or claim"):
            validate(changed)
        changed = copy.deepcopy(payload)
        changed["primary_leads"][0]["lava_failed_qc_accounting"] = (
            "TERMINAL_FAILED_QC_ACCOUNTING_ONLY;SCIENTIFIC_SUPPORT=USED"
        )
        with self.assertRaisesRegex(GATE.GateError, "LAVA accounting"):
            validate(changed)
        changed = copy.deepcopy(payload)
        changed["primary_leads"][0]["unexpected"] = "field"
        with self.assertRaisesRegex(GATE.GateError, "exact verified schema"):
            validate(changed)

    def test_pair_block_and_disjoint_clump_evidence_are_verified(self):
        first = make_lead(
            snp="rs10", position=10,
            members=["PLACO:A:rs10:1:10", "PLACO:A:rs99:1:99"],
        )
        second = make_lead(
            snp="rs20", position=20,
            members=["PLACO:A:rs20:1:20", "PLACO:A:rs98:1:98"],
        )
        payload = upstream_payload([first, second])
        changed = copy.deepcopy(payload)
        changed["primary_leads"][0]["pair_block_evidence_ids"] = "PLACO:A:rs10:1:10"
        with self.assertRaisesRegex(GATE.GateError, "pair/block lead-evidence"):
            validate(changed)
        changed = copy.deepcopy(payload)
        changed["primary_leads"][1]["clump_evidence_ids"] += ";PLACO:A:rs99:1:99"
        changed["ld_reference_coverage"]["eligible_rows"] += 1
        changed["ld_reference_coverage"]["status_counts"]["MATCHED"] += 1
        with self.assertRaisesRegex(GATE.GateError, "overlap"):
            validate(changed)

    def test_no_replace_publish_and_stable_verify_api_allow_downstream_results(self):
        validated = validate(upstream_payload([make_lead()]))
        built = GATE.build_family_rows(validated)
        manifest = GATE.tsv_bytes(GATE.FAMILY_FIELDS, built["mapped_rows"])
        manifest_identity = GATE.identity_for_content(GATE.FAMILY_REL, manifest)
        lock = {
            "schema_version": GATE.SCHEMA,
            "analysis_id": GATE.ANALYSIS_ID,
            "artifact_role": GATE.ARTIFACT_ROLE,
            "state": GATE.READY_STATE,
            "results_accessed_before_lock": False,
            "pair_order": list(GATE.PAIR_ORDER),
            "manifest": manifest_identity,
            "locus_count": 1,
            "locus_entry_ids_in_order": ["TBFM_A_LOC1"],
            "evidence_count": 1,
            "evidence_ids_sha256": GATE.digest_json(built["mapped_evidence_ids"]),
            "policy": identity("policy", 1, "1" * 64),
            "official_blocks": identity("blocks", 1, "2" * 64),
            "upstream_family_provenance": {"post_placo_results": identity()},
            "generator_scripts": {"family_freezer": identity()},
            "selection_contract": dict(GATE.SELECTION_CONTRACT),
            "zero_family_provenance": None,
            "unavailable_manifest": None,
            "unavailable_count": 0,
            "unavailable_entry_ids_in_order": [],
            "unavailable_evidence_count": 0,
            "unavailable_evidence_ids_sha256": GATE.digest_json([]),
            "upstream_evidence_count": 1,
            "upstream_evidence_ids_sha256": GATE.digest_json(built["upstream_evidence_ids"]),
        }
        lock_content = json.dumps(lock, indent=2, sort_keys=True).encode() + b"\n"
        expected = {
            "publication_allowed": True,
            "upstream_state": GATE.UPSTREAM_COMPLETE,
            "upstream_fingerprint": "e" * 64,
            "state": GATE.READY_STATE,
            "mapped_rows": built["mapped_rows"],
            "unavailable_rows": [],
            "upstream_evidence_ids": built["upstream_evidence_ids"],
            "manifest_content": manifest,
            "unavailable_content": None,
            "zero_content": None,
            "lock": lock,
            "lock_content": lock_content,
        }
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            GATE.publish_no_replace(root, GATE.FAMILY_REL, manifest)
            GATE.publish_no_replace(root, GATE.FAMILY_LOCK_REL, lock_content)
            downstream = root / "results/track_b/10_finemap_trait1.tsv"
            downstream.parent.mkdir(parents=True, exist_ok=True)
            downstream.write_text("authorized downstream\n", encoding="utf-8")
            with mock.patch.object(GATE, "assemble_expected", return_value=expected):
                observed = GATE.verify_results(root)
            self.assertEqual(set(observed), {
                "state", "rows", "unavailable_rows", "family_path", "family_lock_path",
                "family_identity", "family_lock_identity", "unavailable_path",
                "unavailable_lock_path", "unavailable_identity", "unavailable_lock_identity",
                "zero_provenance_path", "zero_provenance_identity", "mapped_count",
                "unavailable_count", "upstream_evidence_denominator", "family_sha256",
                "unavailable_sha256",
            })
            self.assertEqual(observed["state"], GATE.READY_STATE)
            self.assertEqual(observed["rows"], built["mapped_rows"])
            self.assertIsNone(observed["unavailable_lock_path"])
            with self.assertRaisesRegex(GATE.GateError, "differs"):
                GATE.publish_no_replace(root, GATE.FAMILY_REL, b"different\n")

    def test_zero_receipt_is_separate_from_later_148_zero_receipt(self):
        self.assertEqual(
            GATE.ZERO_PROVENANCE_REL,
            Path("results/track_b/DEEP_LOCUS_ZERO_FAMILY.provenance.json"),
        )
        self.assertNotEqual(
            GATE.ZERO_PROVENANCE_REL,
            Path("results/track_b/finemapping/zero_family.provenance.json"),
        )

    def test_archive_guard_integration_is_executable_and_two_state(self):
        path = ROOT / "scripts/146_manage_lava_reference_archives.py"
        specification = importlib.util.spec_from_file_location("archive_guard_for_145_test", path)
        assert specification and specification.loader
        module = importlib.util.module_from_spec(specification)
        specification.loader.exec_module(module)
        self.assertTrue(
            module._has_executable_guard_integration(
                GATE.SCRIPT_REL, SCRIPT.read_text(encoding="utf-8"),
            )
        )

    def test_mutating_cli_requires_execute_before_any_work(self):
        error = StringIO()
        with mock.patch.object(GATE, "freeze_family") as freeze, redirect_stderr(error):
            status = GATE.main(["--freeze"])
        self.assertEqual(status, 1)
        freeze.assert_not_called()
        self.assertIn("explicit --execute", error.getvalue())


if __name__ == "__main__":
    unittest.main()
