import csv
import hashlib
import json
import subprocess
import tempfile
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "scripts/127_prepare_track_b_pleiotropy_ld.py"
FIXTURE = ROOT / "tests/fixtures/track_b_pleiotropy_ld/make_fixture.py"
POLICY = ROOT / "config/track_b_pleiotropy_policy.json"


class TrackBPleiotropyLDTests(unittest.TestCase):
    def make_fixture(self, directory: Path, mode: str = "valid") -> None:
        result = subprocess.run(
            ["python3", str(FIXTURE), str(directory), "--mode", mode],
            check=False, text=True, capture_output=True,
        )
        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)

    def run_ld(self, directory: Path, action: str) -> subprocess.CompletedProcess[str]:
        return subprocess.run(
            [
                "python3", str(SCRIPT), "--root", str(directory),
                "--policy", "fixture_policy.json", "--test-fixture", action,
            ],
            check=False, text=True, capture_output=True,
        )

    def test_tiny_zip_materializes_verifies_and_is_idempotent_no_replace(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            self.make_fixture(root)
            preflight = self.run_ld(root, "--preflight")
            self.assertEqual(preflight.returncode, 0, preflight.stdout + preflight.stderr)
            self.assertIn("TRACK_B_PLEIOTROPY_LD_PREFLIGHT_READY", preflight.stdout)
            materialized = self.run_ld(root, "--materialize")
            self.assertEqual(materialized.returncode, 0, materialized.stdout + materialized.stderr)
            self.assertIn("TRACK_B_PLEIOTROPY_LD_MATERIALIZED", materialized.stdout)
            destination = root / "ref/track_b/pleiotropy/g1000_eur"
            self.assertEqual({path.name for path in destination.iterdir()}, {
                "g1000_eur.bed", "g1000_eur.bim", "g1000_eur.fam",
                "g1000_eur.manifest.tsv", "g1000_eur.provenance.json",
            })
            provenance_path = destination / "g1000_eur.provenance.json"
            before = hashlib.sha256(provenance_path.read_bytes()).hexdigest()
            verified = self.run_ld(root, "--verify")
            self.assertEqual(verified.returncode, 0, verified.stdout + verified.stderr)
            self.assertIn("variants=4 autosomal=3", verified.stdout)
            repeated = self.run_ld(root, "--materialize")
            self.assertEqual(repeated.returncode, 0, repeated.stdout + repeated.stderr)
            self.assertEqual(hashlib.sha256(provenance_path.read_bytes()).hexdigest(), before)

            provenance = json.loads(provenance_path.read_text(encoding="utf-8"))
            self.assertEqual(provenance["reference"]["variant_count"], 4)
            self.assertEqual(provenance["reference"]["eligible_autosomal_variant_count"], 3)
            self.assertEqual(provenance["reference"]["excluded_source_fidelity_chromosomes"], [23])
            self.assertEqual(
                provenance["publication"]["mode"],
                "ATOMIC_DIRECTORY_RENAME_WITH_KERNEL_NO_REPLACE",
            )
            self.assertEqual(list((root / "ref/track_b/pleiotropy/.staging").iterdir()), [])
            self.assertEqual(
                provenance["validation"]["bim"]["chromosome_variant_counts"],
                {"1": 1, "2": 1, "22": 1, "23": 1},
            )
            self.assertIn("no LD-block splitting", provenance["materialization_transform_rule"])
            metrics = provenance["materialization_resource_metrics"]
            self.assertGreater(metrics["peak_rss_bytes"], 0)
            self.assertGreaterEqual(metrics["staging_extraction_validation_wall_seconds"], 0)
            self.assertEqual(metrics["streaming_chunk_bytes"], 8 * 1024 * 1024)
            fixture_policy = json.loads((root / "fixture_policy.json").read_text(encoding="utf-8"))
            fixture_members = fixture_policy["ld_reference_policy"]["source_members"]
            self.assertEqual(
                metrics["extracted_uncompressed_bytes"],
                sum(member["expected_bytes"] for member in fixture_members.values()),
            )
            with (destination / "g1000_eur.manifest.tsv").open(newline="", encoding="utf-8") as handle:
                rows = list(csv.DictReader(handle, delimiter="\t"))
            self.assertEqual([row["component_id"] for row in rows], ["BED", "BIM", "FAM"])
            self.assertTrue(all(row["path"].startswith("ref/track_b/pleiotropy/g1000_eur/") for row in rows))

    def test_active_lava_marker_blocks_before_staging(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            self.make_fixture(root)
            (root / "active_lava_download.tmp").write_text("active\n", encoding="utf-8")
            result = self.run_ld(root, "--materialize")
            self.assertNotEqual(result.returncode, 0)
            self.assertIn("BLOCKED_ACTIVE_LAVA_DOWNLOAD", result.stderr)
            self.assertFalse((root / "ref/track_b/pleiotropy/g1000_eur").exists())
            self.assertFalse((root / "ref/track_b/pleiotropy/.staging").exists())

    def test_invalid_bim_fails_closed_and_preserves_fingerprinted_stage(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            self.make_fixture(root, mode="invalid-bim")
            result = self.run_ld(root, "--materialize")
            self.assertNotEqual(result.returncode, 0)
            self.assertIn("BIM coordinate/ID/allele drift", result.stderr)
            self.assertFalse((root / "ref/track_b/pleiotropy/g1000_eur").exists())
            stages = list((root / "ref/track_b/pleiotropy/.staging").iterdir())
            self.assertEqual(len(stages), 1)
            self.assertRegex(stages[0].name, r"^g1000_eur_[0-9a-f]{64}$")
            repeated = self.run_ld(root, "--materialize")
            self.assertNotEqual(repeated.returncode, 0)
            self.assertIn("staging directory already exists", repeated.stderr)

    def test_missing_zip_member_and_tampered_archive_are_rejected(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            self.make_fixture(root, mode="missing-member")
            missing = self.run_ld(root, "--preflight")
            self.assertNotEqual(missing.returncode, 0)
            self.assertIn("ZIP member family drifted", missing.stderr)
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            self.make_fixture(root)
            archive = root / "fixture_source/g1000_eur.zip"
            archive.write_bytes(archive.read_bytes() + b"tamper")
            tampered = self.run_ld(root, "--preflight")
            self.assertNotEqual(tampered.returncode, 0)
            self.assertIn("archive size or SHA-256 drifted", tampered.stderr)

    def test_partial_or_tampered_destination_is_never_replaced(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            self.make_fixture(root)
            destination = root / "ref/track_b/pleiotropy/g1000_eur"
            destination.mkdir(parents=True)
            marker = destination / "g1000_eur.bed"
            marker.write_bytes(b"partial")
            result = self.run_ld(root, "--materialize")
            self.assertNotEqual(result.returncode, 0)
            self.assertIn("partial/preexisting", result.stderr)
            self.assertEqual(marker.read_bytes(), b"partial")
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            self.make_fixture(root)
            complete = self.run_ld(root, "--materialize")
            self.assertEqual(complete.returncode, 0, complete.stdout + complete.stderr)
            bed = root / "ref/track_b/pleiotropy/g1000_eur/g1000_eur.bed"
            payload = bytearray(bed.read_bytes())
            payload[-1] ^= 0x01
            bed.write_bytes(payload)
            verify = self.run_ld(root, "--verify")
            self.assertNotEqual(verify.returncode, 0)
            self.assertIn("manifest differs from live full hashes", verify.stderr)
            retry = self.run_ld(root, "--materialize")
            self.assertNotEqual(retry.returncode, 0)
            self.assertIn("manifest differs from live full hashes", retry.stderr)

    def test_disk_and_plink_identity_gates_fail_closed(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            self.make_fixture(root)
            policy_path = root / "fixture_policy.json"
            policy = json.loads(policy_path.read_text(encoding="utf-8"))
            policy["ld_reference_policy"]["minimum_free_storage_bytes"] = 10**18
            policy_path.write_text(json.dumps(policy, indent=2, sort_keys=True) + "\n", encoding="utf-8")
            disk = self.run_ld(root, "--preflight")
            self.assertNotEqual(disk.returncode, 0)
            self.assertIn("insufficient disk", disk.stderr)
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            self.make_fixture(root)
            plink = root / "fixture_source/plink"
            plink.write_text("#!/bin/sh\necho wrong\n", encoding="utf-8")
            plink.chmod(0o755)
            identity = self.run_ld(root, "--preflight")
            self.assertNotEqual(identity.returncode, 0)
            self.assertIn("PLINK SHA-256 drifted", identity.stderr)

    def test_production_policy_preserves_raw_x_but_forbids_it_from_analysis(self) -> None:
        policy = json.loads(POLICY.read_text(encoding="utf-8"))
        ld = policy["ld_reference_policy"]
        self.assertEqual(ld["materialized_directory"], "ref/track_b/pleiotropy/g1000_eur")
        self.assertEqual(ld["source_archive_expected_bytes"], 511_626_945)
        self.assertEqual(
            ld["source_archive_sha256"],
            "83a48fd9dcaa0b9a874b18c63143a4ede93f05505b215b0bd8790130a0d7a954",
        )
        self.assertEqual(ld["materializer_script"], "scripts/127_prepare_track_b_pleiotropy_ld.py")
        self.assertEqual(
            ld["plink_version"], "PLINK v1.9.0-b.7.11 64-bit (19 Aug 2025)",
        )
        self.assertEqual(ld["variant_count"], 22_665_064)
        self.assertEqual(ld["expected_autosomal_variant_count"], 22_132_657)
        self.assertEqual(ld["reference_chromosomes"], list(range(1, 24)))
        self.assertEqual(ld["analysis_autosomes"], list(range(1, 23)))
        self.assertEqual(ld["nonanalysis_chromosomes_retained_for_source_fidelity"], [23])
        self.assertEqual(ld["expected_reference_chromosome_variant_counts"]["23"], 532_407)
        self.assertIn("chr23/X is always forbidden", ld["analysis_chromosome_rule"])
        self.assertIn("no LD-block splitting", ld["materialization_transform_rule"])
        self.assertIn("SNP reduction", ld["materialization_transform_rule"])


if __name__ == "__main__":
    unittest.main()
