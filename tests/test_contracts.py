import csv
import gzip
import hashlib
import importlib.util
import io
import math
import json
import os
from pathlib import Path
import shutil
import sqlite3
import subprocess
import sys
import tarfile
import tempfile
import unittest
import zipfile


ROOT = Path(__file__).resolve().parents[1]
MANIFEST = ROOT / "config" / "analysis_panel.tsv"
LOCK = ROOT / "config" / "analysis_panel.lock.json"


def build_fixture_variant_map(directory):
    directory = Path(directory)
    reference = directory / "eur_w_ld_chr"
    reference.mkdir()
    (reference / "w_hm3.snplist").write_text(
        "SNP\tA1\tA2\n"
        "rs123\tA\tG\n"
        "rs456\tC\tT\n"
        "rs789\tG\tC\n",
        encoding="utf-8",
    )
    with gzip.open(reference / "1.l2.ldscore.gz", "wt", encoding="utf-8", newline="") as handle:
        handle.write(
            "CHR\tSNP\tBP\tCM\tMAF\tL2\n"
            "1\trs123\t1000000\t0\t0.2\t1.0\n"
            "1\trs456\t2000000\t0\t0.3\t1.0\n"
            "1\trs999\t3000000\t0\t0.2\t1.0\n"
        )
    output = directory / "hm3_grch37_variant_map.tsv.gz"
    result = subprocess.run(
        [
            sys.executable,
            str(ROOT / "scripts" / "19_build_hm3_variant_map.py"),
            "--reference-dir", str(reference),
            "--out", str(output),
            "--chromosomes", "1",
            "--expected-hm3-count", "3",
            "--min-coverage", "0.6",
        ],
        cwd=ROOT,
        capture_output=True,
        text=True,
    )
    if result.returncode:
        raise AssertionError(result.stderr + result.stdout)
    return output


def build_fixture_liftover_chain(directory):
    path = Path(directory) / "hg38ToHg19.over.chain.gz"
    payload = (
        "chain 1000 chr1 1000 + 100 200 chr1 1000 + 50 150 1\n"
        "100\n\n"
        "chain 900 chr2 1000 + 100 200 chr2 1000 - 700 800 2\n"
        "100\n"
    )
    path.write_bytes(gzip.compress(payload.encode("ascii"), mtime=0))
    return path


def load_collator():
    spec = importlib.util.spec_from_file_location(
        "atlas_collator", ROOT / "scripts" / "05_collate.py"
    )
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def load_sbp_materializer():
    scripts_dir = str(ROOT / "scripts")
    spec = importlib.util.spec_from_file_location(
        "sbp_materializer", ROOT / "scripts" / "18_materialize_opengwas_vcf.py"
    )
    module = importlib.util.module_from_spec(spec)
    sys.path.insert(0, scripts_dir)
    try:
        spec.loader.exec_module(module)
    finally:
        sys.path.remove(scripts_dir)
    return module


def load_liftover_module():
    spec = importlib.util.spec_from_file_location(
        "atlas_liftover", ROOT / "scripts" / "liftover_chain.py"
    )
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def load_prostate_materializer():
    scripts_dir = str(ROOT / "scripts")
    spec = importlib.util.spec_from_file_location(
        "prostate_materializer", ROOT / "scripts" / "20_materialize_practical_prostate.py"
    )
    module = importlib.util.module_from_spec(spec)
    sys.path.insert(0, scripts_dir)
    try:
        spec.loader.exec_module(module)
    finally:
        sys.path.remove(scripts_dir)
    return module


def load_phase0_audit():
    spec = importlib.util.spec_from_file_location(
        "phase0_audit", ROOT / "scripts" / "10_phase0_audit.py"
    )
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def load_ranged_downloader():
    spec = importlib.util.spec_from_file_location(
        "ranged_downloader", ROOT / "scripts" / "14_ranged_download.py"
    )
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def load_substitution_auditor():
    spec = importlib.util.spec_from_file_location(
        "substitution_auditor", ROOT / "scripts" / "22_audit_substitution_source.py"
    )
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def load_pleiotropy_materializer():
    spec = importlib.util.spec_from_file_location(
        "pleiotropy_materializer", ROOT / "scripts" / "44_materialize_pleiotropy_pair.py"
    )
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def load_finemapping_materializer():
    spec = importlib.util.spec_from_file_location(
        "finemapping_materializer", ROOT / "scripts" / "58_materialize_finemapping_locus.py"
    )
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def load_numbered_script(filename, module_name):
    scripts_dir = str(ROOT / "scripts")
    spec = importlib.util.spec_from_file_location(module_name, ROOT / "scripts" / filename)
    module = importlib.util.module_from_spec(spec)
    sys.path.insert(0, scripts_dir)
    try:
        spec.loader.exec_module(module)
    finally:
        sys.path.remove(scripts_dir)
    return module


class PanelContractTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        with MANIFEST.open(newline="", encoding="utf-8") as handle:
            cls.panel = list(csv.DictReader(handle, delimiter="\t"))

    def test_exact_panel_shape(self):
        self.assertEqual(len(self.panel), 45)
        self.assertEqual(sum(row["domain"] == "sleep" for row in self.panel), 12)
        self.assertEqual(sum(row["domain"] != "sleep" for row in self.panel), 33)
        self.assertEqual(len({row["trait_id"] for row in self.panel}), 45)

    def test_identity_lock_rejects_silent_swap(self):
        with tempfile.TemporaryDirectory() as directory:
            changed = Path(directory) / "analysis_panel.tsv"
            rows = [dict(row) for row in self.panel]
            rows[0]["trait_id"] = "silently_swapped_trait"
            with changed.open("w", newline="", encoding="utf-8") as handle:
                writer = csv.DictWriter(
                    handle, fieldnames=rows[0].keys(), delimiter="\t", lineterminator="\n"
                )
                writer.writeheader()
                writer.writerows(rows)
            result = subprocess.run(
                [
                    sys.executable,
                    str(ROOT / "scripts" / "00_validate_panel.py"),
                    "--manifest",
                    str(changed),
                    "--lock",
                    str(LOCK),
                ],
                cwd=ROOT,
                capture_output=True,
                text=True,
            )
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("identity lock changed", result.stderr + result.stdout)

    def test_source_registry_maps_only_locked_traits(self):
        panel_ids = {row["trait_id"] for row in self.panel}
        with (ROOT / "config" / "public_gwas_sources.tsv").open(
            newline="", encoding="utf-8"
        ) as handle:
            sources = csv.DictReader(handle, delimiter="\t")
            mapped = {
                trait.strip()
                for source in sources
                for trait in source["trait_ids"].split(",")
                if trait.strip()
            }
        self.assertEqual(mapped, panel_ids)

    def test_substitution_candidates_are_exactly_the_six_blockers(self):
        selected = {row["trait_id"]: row["source_id"] for row in self.panel}
        with (ROOT / "config" / "public_gwas_substitution_candidates.tsv").open(
            newline="", encoding="utf-8"
        ) as handle:
            candidates = list(csv.DictReader(handle, delimiter="\t"))
        self.assertEqual(
            {row["trait_id"] for row in candidates},
            {"ms", "asthma", "t2d", "cad", "melanoma", "telomere_length"},
        )
        self.assertEqual(len(candidates), 6)
        self.assertEqual(
            sum(int(row["archive_bytes"]) for row in candidates),
            4_437_046_355,
        )
        for row in candidates:
            self.assertEqual(row["candidate_source_id"], selected[row["trait_id"]])
            self.assertNotEqual(row["replaces_source_id"], selected[row["trait_id"]])
            self.assertEqual(
                row["review_status"],
                "PROMOTED_SOURCE_VERIFIED",
            )
            self.assertRegex(row["archive_sha256"], r"^[0-9a-f]{64}$")
            self.assertTrue(row["download_url"].startswith("https://"))
            self.assertTrue(row["source_page_url"].startswith("https://"))
            self.assertIn(row["build"], {"GRCh38/hg38"})

    def test_ranged_downloader_supports_md5_bootstrap(self):
        downloader = load_ranged_downloader()
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "payload"
            path.write_bytes(b"atlas-candidate")
            self.assertEqual(
                downloader.file_digest(path, "md5"),
                hashlib.md5(b"atlas-candidate").hexdigest(),
            )
            self.assertEqual(
                downloader.sha256(path),
                hashlib.sha256(b"atlas-candidate").hexdigest(),
            )

    def test_substitution_source_auditor_scans_literal_finngen_schema(self):
        auditor = load_substitution_auditor()
        header = "\t".join(auditor.FINNGEN_HEADER)
        row = "\t".join([
            "1", "1000000", "A", "G", "rs123", "GENE", "0.01", "2",
            "0.2", "0.1", "0.3", "0.4", "0.3",
        ])
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "finngen.gz"
            with gzip.open(path, "wt", encoding="utf-8", newline="") as handle:
                handle.write(f"{header}\n{row}\n")
            candidate = {
                "candidate_source_id": "finngen_r9_ms",
                "trait_id": "ms",
                "archive_name": path.name,
                "archive_bytes": str(path.stat().st_size),
                "archive_md5": hashlib.md5(path.read_bytes()).hexdigest(),
                "archive_sha256": hashlib.sha256(path.read_bytes()).hexdigest(),
            }
            report = auditor.audit_source(path, candidate)
        self.assertEqual(report["source_rows"], 1)
        self.assertEqual(report["autosomal_rows"], 1)
        self.assertEqual(report["explicit_rsid_rows"], 1)
        self.assertEqual(report["valid_effect_se_p_rows"], 1)

    def test_nonjournal_citation_is_exactly_the_neale_grip_release(self):
        audit = load_phase0_audit()
        selected = {row["trait_id"]: row["source_id"] for row in self.panel}
        citations = audit.load_nonjournal_citations(
            ROOT / "config" / "nonjournal_source_citations.tsv", selected
        )
        self.assertEqual(set(citations), {"grip_strength"})
        citation = citations["grip_strength"]
        self.assertEqual(citation["source_id"], "neale_2018_left_grip_strength")
        self.assertEqual(citation["citation_type"], "PUBLIC_DATA_RELEASE")

    def test_nonjournal_citation_rejects_wrong_selected_source(self):
        audit = load_phase0_audit()
        with tempfile.TemporaryDirectory() as directory:
            citation_path = Path(directory) / "citations.tsv"
            citation_path.write_text(
                "source_id\ttrait_id\tcitation_type\ttitle\tpublisher\t"
                "release_date\tsource_url\tnotes\n"
                "wrong_source\tgrip_strength\tPUBLIC_DATA_RELEASE\tRelease\t"
                "Neale Lab\t2018-08-01\thttps://example.org/release\tExact release\n",
                encoding="utf-8",
            )
            with self.assertRaises(SystemExit) as error:
                audit.load_nonjournal_citations(
                    citation_path,
                    {"grip_strength": "neale_2018_left_grip_strength"},
                )
        self.assertIn("selected_source_id", str(error.exception))

    def test_prostate_materializer_filters_only_wrong_width_rows(self):
        materializer = load_prostate_materializer()
        header = materializer.HEADER.decode("ascii")
        fields = [
            "marker", "tag", "rs123", "1", "100", "a", "g", "0.2",
            "0.01", "0.1", "0.3", "0.05", "0.01", "1e-6", "+", "0.99",
        ]
        valid = "\t".join(fields)
        malformed = "\t".join(fields[:-1])
        projected = "\t".join(fields[index] for index in materializer.OUTPUT_COLUMN_INDEXES)
        projected_header = materializer.OUTPUT_HEADER.decode("ascii")
        with tempfile.TemporaryDirectory() as directory:
            directory = Path(directory)
            source = directory / "source.tsv"
            output = directory / "prostate.tsv.gz"
            source.write_text(
                f"{header}\n{valid}\n{malformed}\n{valid}\n", encoding="ascii"
            )
            counts = materializer.scan_or_materialize(
                source,
                output,
                expected_data_rows=3,
                expected_malformed_rows=1,
                expected_field_count=16,
                expected_non_rsid_rows=0,
                expected_non_snp_rows=0,
                expected_nonautosomal_rows=0,
                expected_retained_rows=2,
            )
            with gzip.open(output, "rt", encoding="ascii") as handle:
                materialized = handle.read().splitlines()
            provenance = Path(f"{output}.provenance.json").read_text(encoding="utf-8")
        self.assertEqual(counts, {
            "source_rows": 3,
            "retained_rows": 2,
            "malformed_rows": 1,
            "non_rsid_rows": 0,
            "non_snp_rows": 0,
            "nonautosomal_rows": 0,
        })
        self.assertEqual(materialized, [projected_header, projected, projected])
        self.assertIn('"malformed_rows": 1', provenance)

    def test_source_schemas_match_selected_source_trait_pairs(self):
        with (ROOT / "config" / "analysis_panel.tsv").open(newline="") as handle:
            selected = {
                row["trait_id"]: row["source_id"]
                for row in csv.DictReader(handle, delimiter="\t")
            }
        with (ROOT / "config" / "gwas_schemas.tsv").open(newline="") as handle:
            schemas = list(csv.DictReader(handle, delimiter="\t"))
        pairs = [(row["source_id"], row["trait_id"]) for row in schemas]
        self.assertEqual(len(pairs), len(set(pairs)))
        self.assertEqual({row["trait_id"] for row in schemas}, set(selected))
        for source_id, trait_id in pairs:
            self.assertEqual(selected.get(trait_id), source_id)
        allowed_statuses = {
            "SCHEMA_PENDING", "SCHEMA_VERIFIED",
            "SCHEMA_VERIFIED_REQUIRES_VARIANT_MAPPING",
        }
        self.assertFalse({row["schema_status"] for row in schemas} - allowed_statuses)
        jones = [row for row in schemas if row["source_id"] == "jones_2019_accelerometer_sleep"]
        self.assertEqual({row["trait_id"] for row in jones}, {
            "sleep_efficiency", "accel_sleep_duration", "sleep_timing",
        })
        self.assertEqual(len({row["effect"] for row in jones}), 3)

    def test_variant_mapping_plans_exactly_cover_mapping_required_schemas(self):
        with (ROOT / "config" / "gwas_schemas.tsv").open(newline="") as handle:
            schemas = list(csv.DictReader(handle, delimiter="\t"))
        with (ROOT / "config" / "variant_mapping_plans.tsv").open(newline="") as handle:
            mappings = list(csv.DictReader(handle, delimiter="\t"))
        required = {
            row["trait_id"]
            for row in schemas
            if row["schema_status"] == "SCHEMA_VERIFIED_REQUIRES_VARIANT_MAPPING"
        }
        self.assertEqual({row["trait_id"] for row in mappings}, required)
        self.assertEqual(len(mappings), len({row["trait_id"] for row in mappings}))
        self.assertEqual(
            {row["strategy"] for row in mappings},
            {"BY_COORD_ALLELES", "BY_RSID_ALLELES"},
        )
        self.assertEqual({row["map_path"] for row in mappings}, {
            "ref/hm3_grch37_variant_map.tsv.gz"
        })
        self.assertEqual({row["map_bytes"] for row in mappings}, {"9948937"})
        self.assertEqual({row["map_sha256"] for row in mappings}, {
            "6775a7a0d3ca90dc74e472180b1d77103bc238129c4f969358f46307e5c306b4"
        })

    def test_liftover_plans_exactly_cover_verified_hg38_sources(self):
        verified_hg38 = {
            row["trait_id"]
            for row in self.panel
            if row["source_status"] == "SOURCE_VERIFIED"
            and row["build"] in {"hg38", "GRCh38"}
        }
        with (ROOT / "config" / "liftover_plans.tsv").open(newline="") as handle:
            plans = list(csv.DictReader(handle, delimiter="\t"))
        self.assertEqual({row["trait_id"] for row in plans}, verified_hg38)
        self.assertEqual({row["strategy"] for row in plans}, {"UCSC_CHAIN_POINT"})
        self.assertEqual({row["chain_bytes"] for row in plans}, {"1246411"})
        self.assertEqual({row["chain_md5"] for row in plans}, {
            "ff3031d93792f4cbb86af44055efd903"
        })
        self.assertEqual({row["chain_sha256"] for row in plans}, {
            "14a712e8e147d9fc8e9d87d51977b46f6f8ddb93efbe5d0843d86b6205f587b1"
        })

    def test_public_substitutions_use_exact_endpoint_sample_counts(self):
        panel = {row["trait_id"]: row for row in self.panel}
        self.assertEqual(
            (panel["ms"]["ncase"], panel["ms"]["ncontrol"], panel["ms"]["n_total"]),
            ("2182", "373987", "376169"),
        )
        self.assertEqual(panel["ms"]["source_status"], "SOURCE_VERIFIED")
        self.assertEqual(
            (
                panel["melanoma"]["ncase"],
                panel["melanoma"]["ncontrol"],
                panel["melanoma"]["n_total"],
                panel["melanoma"]["build"],
            ),
            ("2993", "287137", "290130", "hg38"),
        )
        self.assertEqual(panel["melanoma"]["source_status"], "SOURCE_VERIFIED")
        self.assertEqual(panel["t2d"]["source_status"], "SOURCE_VERIFIED")

    def test_hm3_map_assigns_rsid_from_coordinate_and_alleles(self):
        payload = (
            "MarkerName\tAllele1\tAllele2\tEffect\tStdErr\tP.value\n"
            "1:1000000_G_A\tG\tA\t0.2\t0.1\t0.01\n"
            "1:2000000_A_G\tA\tG\t-0.3\t0.1\t0.02\n"
            "1:3000001_A_C\tA\tC\t0.1\t0.1\t0.03\n"
        )
        with tempfile.TemporaryDirectory() as directory:
            directory = Path(directory)
            variant_map = build_fixture_variant_map(directory)
            variant_map_sha256 = hashlib.sha256(variant_map.read_bytes()).hexdigest()
            source = directory / "ibd.tsv.gz"
            outdir = directory / "out"
            with gzip.open(source, "wt", encoding="utf-8", newline="") as handle:
                handle.write(payload)
            result = subprocess.run(
                [
                    sys.executable,
                    str(ROOT / "scripts" / "01_harmonize.py"),
                    "--trait", "ibd",
                    "--config", str(MANIFEST),
                    "--infile", str(source),
                    "--outdir", str(outdir),
                    "--source-build", "hg19",
                    "--variant-map", str(variant_map),
                    "--variant-map-strategy", "BY_COORD_ALLELES",
                    "--expected-variant-map-bytes", str(variant_map.stat().st_size),
                    "--expected-variant-map-sha256", variant_map_sha256,
                ],
                cwd=ROOT,
                capture_output=True,
                text=True,
            )
            self.assertEqual(result.returncode, 0, result.stderr + result.stdout)
            with gzip.open(outdir / "ibd.harmonized.tsv.gz", "rt", newline="") as handle:
                rows = list(csv.DictReader(handle, delimiter="\t"))
            qc = (outdir / "ibd.qc.txt").read_text(encoding="utf-8")
        self.assertEqual([row["SNP"] for row in rows], ["rs123", "rs456"])
        self.assertEqual(
            [(row["CHR"], row["BP"], row["A1"], row["A2"]) for row in rows],
            [("1", "1000000", "G", "A"), ("1", "2000000", "A", "G")],
        )
        self.assertIn("not present in pinned GRCh37 EUR HapMap3 map\t1\t2", qc)
        self.assertIn("variant_map_strategy\tBY_COORD_ALLELES", qc)

    def test_per_snp_effective_n_gate_uses_authoritative_source_maximum(self):
        payload = (
            "SNP\tChr\tPosition\tEA\tNEA\tEAF\tBeta\tSE\tP-value\tEffective_N\n"
            "1:1000000\t1\t1000000\tG\tA\t0.2\t0.2\t0.1\t0.01\t5808\n"
            "1:2000000\t1\t2000000\tC\tT\t0.3\t-0.3\t0.1\t0.02\t11615\n"
        )
        with tempfile.TemporaryDirectory() as directory:
            directory = Path(directory)
            variant_map = build_fixture_variant_map(directory)
            source = directory / "longevity.tsv.gz"
            outdir = directory / "out"
            with gzip.open(source, "wt", encoding="utf-8", newline="") as handle:
                handle.write(payload)
            result = subprocess.run(
                [
                    sys.executable,
                    str(ROOT / "scripts" / "01_harmonize.py"),
                    "--trait", "longevity",
                    "--config", str(MANIFEST),
                    "--infile", str(source),
                    "--outdir", str(outdir),
                    "--source-build", "hg19",
                    "--variant-map", str(variant_map),
                    "--variant-map-strategy", "BY_COORD_ALLELES",
                    "--expected-variant-map-bytes", str(variant_map.stat().st_size),
                    "--expected-variant-map-sha256",
                    hashlib.sha256(variant_map.read_bytes()).hexdigest(),
                ],
                cwd=ROOT,
                capture_output=True,
                text=True,
            )
            self.assertEqual(result.returncode, 0, result.stderr + result.stdout)
            with gzip.open(
                outdir / "longevity.harmonized.tsv.gz", "rt", newline=""
            ) as handle:
                rows = list(csv.DictReader(handle, delimiter="\t"))
            qc = (outdir / "longevity.qc.txt").read_text(encoding="utf-8")
        self.assertEqual([float(row["N"]) for row in rows], [11616.0, 23230.0])
        self.assertIn(
            "sample size below 50% of source maximum (23,230.0) [CDG3]\t0\t2",
            qc,
        )
        self.assertIn(
            "sample-size mode: per-SNP N_eff derived as 2 * source Effective_N",
            qc,
        )

    def test_hm3_map_bytes_are_deterministic_across_output_paths(self):
        with tempfile.TemporaryDirectory() as directory:
            directory = Path(directory)
            first = build_fixture_variant_map(directory)
            second = directory / "map-with-a-different-name.tsv.gz"
            result = subprocess.run(
                [
                    sys.executable,
                    str(ROOT / "scripts" / "19_build_hm3_variant_map.py"),
                    "--reference-dir", str(directory / "eur_w_ld_chr"),
                    "--out", str(second),
                    "--chromosomes", "1",
                    "--expected-hm3-count", "3",
                    "--min-coverage", "0.6",
                ],
                cwd=ROOT,
                capture_output=True,
                text=True,
            )
            self.assertEqual(result.returncode, 0, result.stderr + result.stdout)
            self.assertEqual(first.read_bytes(), second.read_bytes())

    def test_hm3_map_checksum_mismatch_fails_closed(self):
        payload = (
            "SNP\tA1\tA2\tBETA\tSE\tP\tN\n"
            "rs123\tG\tA\t0.2\t0.1\t0.01\t380000\n"
        )
        with tempfile.TemporaryDirectory() as directory:
            directory = Path(directory)
            variant_map = build_fixture_variant_map(directory)
            with variant_map.open("ab") as handle:
                handle.write(b"corrupt")
            source = directory / "rsid-only.tsv"
            source.write_text(payload, encoding="utf-8")
            result = subprocess.run(
                [
                    sys.executable,
                    str(ROOT / "scripts" / "01_harmonize.py"),
                    "--trait", "insomnia",
                    "--config", str(MANIFEST),
                    "--infile", str(source),
                    "--outdir", str(directory / "out"),
                    "--variant-map", str(variant_map),
                    "--variant-map-strategy", "BY_RSID_ALLELES",
                ],
                cwd=ROOT,
                capture_output=True,
                text=True,
            )
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("variant-map SHA-256 mismatch", result.stderr + result.stdout)

    def test_hm3_map_registered_hash_mismatch_fails_closed(self):
        payload = (
            "SNP\tA1\tA2\tBETA\tSE\tP\tN\n"
            "rs123\tG\tA\t0.2\t0.1\t0.01\t380000\n"
        )
        with tempfile.TemporaryDirectory() as directory:
            directory = Path(directory)
            variant_map = build_fixture_variant_map(directory)
            source = directory / "rsid-only.tsv"
            source.write_text(payload, encoding="utf-8")
            result = subprocess.run(
                [
                    sys.executable,
                    str(ROOT / "scripts" / "01_harmonize.py"),
                    "--trait", "insomnia",
                    "--config", str(MANIFEST),
                    "--infile", str(source),
                    "--outdir", str(directory / "out"),
                    "--variant-map", str(variant_map),
                    "--variant-map-strategy", "BY_RSID_ALLELES",
                    "--expected-variant-map-sha256", "0" * 64,
                ],
                cwd=ROOT,
                capture_output=True,
                text=True,
            )
        self.assertNotEqual(result.returncode, 0)
        self.assertIn(
            "variant-map SHA-256 differs from the registered map",
            result.stderr + result.stdout,
        )

    def test_hm3_map_assigns_grch37_coordinate_from_rsid_and_alleles(self):
        payload = (
            "SNP\tA1\tA2\tBETA\tSE\tP\tN\n"
            "rs123\tG\tA\t0.2\t0.1\t0.01\t380000\n"
            "rs456\tA\tG\t-0.3\t0.1\t0.02\t380000\n"
            "rs789\tA\tC\t0.1\t0.1\t0.03\t380000\n"
        )
        with tempfile.TemporaryDirectory() as directory:
            directory = Path(directory)
            variant_map = build_fixture_variant_map(directory)
            source = directory / "rsid-only.tsv"
            outdir = directory / "out"
            source.write_text(payload, encoding="utf-8")
            result = subprocess.run(
                [
                    sys.executable,
                    str(ROOT / "scripts" / "01_harmonize.py"),
                    "--trait", "insomnia",
                    "--config", str(MANIFEST),
                    "--infile", str(source),
                    "--outdir", str(outdir),
                    "--variant-map", str(variant_map),
                    "--variant-map-strategy", "BY_RSID_ALLELES",
                ],
                cwd=ROOT,
                capture_output=True,
                text=True,
            )
            self.assertEqual(result.returncode, 0, result.stderr + result.stdout)
            with gzip.open(outdir / "insomnia.harmonized.tsv.gz", "rt", newline="") as handle:
                rows = list(csv.DictReader(handle, delimiter="\t"))
            qc = (outdir / "insomnia.qc.txt").read_text(encoding="utf-8")
        self.assertEqual(
            [(row["SNP"], row["CHR"], row["BP"]) for row in rows],
            [("rs123", "1", "1000000"), ("rs456", "1", "2000000")],
        )
        self.assertIn("not present in pinned GRCh37 EUR HapMap3 map\t1\t2", qc)
        self.assertIn("output_build\thg19", qc)

    def test_mdd_literal_stderrlogor_header_harmonizes(self):
        payload = (
            "MarkerName A1 A2 Freq LogOR StdErrLogOR P\n"
            "rs123 g a 0.2 0.2 0.1 0.01\n"
            "rs456 a g 0.3 -0.3 0.1 0.02\n"
        )
        with tempfile.TemporaryDirectory() as directory:
            directory = Path(directory)
            variant_map = build_fixture_variant_map(directory)
            source = directory / "mdd.txt"
            source.write_text(payload, encoding="utf-8")
            outdir = directory / "out"
            result = subprocess.run(
                [
                    sys.executable,
                    str(ROOT / "scripts" / "01_harmonize.py"),
                    "--trait", "mdd",
                    "--config", str(MANIFEST),
                    "--infile", str(source),
                    "--outdir", str(outdir),
                    "--variant-map", str(variant_map),
                    "--variant-map-strategy", "BY_RSID_ALLELES",
                    "--expected-variant-map-bytes", str(variant_map.stat().st_size),
                    "--expected-variant-map-sha256",
                    hashlib.sha256(variant_map.read_bytes()).hexdigest(),
                ],
                cwd=ROOT,
                capture_output=True,
                text=True,
            )
            self.assertEqual(result.returncode, 0, result.stderr + result.stdout)
            with gzip.open(outdir / "mdd.harmonized.tsv.gz", "rt", newline="") as handle:
                rows = list(csv.DictReader(handle, delimiter="\t"))
        self.assertEqual([row["SNP"] for row in rows], ["rs123", "rs456"])
        self.assertEqual([float(row["SE"]) for row in rows], [0.1, 0.1])

    def test_hg38_liftover_maps_points_and_orients_reverse_strand_alleles(self):
        payload = (
            "rsids\t#chrom\tpos\talt\tref\tbeta\tsebeta\tpval\taf_alt\n"
            "rs123\t1\t101\tA\tC\t0.2\t0.1\t0.01\t0.2\n"
            "rs124\t2\t101\tA\tC\t-0.3\t0.1\t0.02\t0.3\n"
            "rs125\t1\t500\tA\tC\t0.1\t0.1\t0.03\t0.4\n"
        )
        with tempfile.TemporaryDirectory() as directory:
            directory = Path(directory)
            chain = build_fixture_liftover_chain(directory)
            chain_sha256 = hashlib.sha256(chain.read_bytes()).hexdigest()
            source = directory / "sleep-apnea.tsv"
            source.write_text(payload, encoding="utf-8")
            outdir = directory / "out"
            result = subprocess.run(
                [
                    sys.executable,
                    str(ROOT / "scripts" / "01_harmonize.py"),
                    "--trait", "sleep_apnea",
                    "--config", str(MANIFEST),
                    "--infile", str(source),
                    "--outdir", str(outdir),
                    "--source-build", "hg38",
                    "--liftover-chain", str(chain),
                    "--expected-liftover-chain-bytes", str(chain.stat().st_size),
                    "--expected-liftover-chain-sha256", chain_sha256,
                ],
                cwd=ROOT,
                capture_output=True,
                text=True,
            )
            self.assertEqual(result.returncode, 0, result.stderr + result.stdout)
            with gzip.open(
                outdir / "sleep_apnea.harmonized.tsv.gz", "rt", newline=""
            ) as handle:
                rows = list(csv.DictReader(handle, delimiter="\t"))
            qc = (outdir / "sleep_apnea.qc.txt").read_text(encoding="utf-8")
        self.assertEqual(
            [(row["SNP"], row["CHR"], row["BP"], row["A1"], row["A2"]) for row in rows],
            [
                ("rs123", "1", "51", "A", "C"),
                ("rs124", "2", "300", "T", "G"),
            ],
        )
        self.assertIn("not mapped by pinned UCSC hg38-to-hg19 chain\t1\t2", qc)
        self.assertIn("reverse-strand liftover mappings with allele complements: 1", qc)
        self.assertIn(f"liftover_chain_sha256\t{chain_sha256}", qc)

    def test_hg38_liftover_rejects_multiple_target_points(self):
        liftover = load_liftover_module()
        payload = (
            "chain 1000 chr1 1000 + 0 100 chr1 1000 + 0 100 1\n"
            "100\n\n"
            "chain 900 chr1 1000 + 0 100 chr1 1000 + 100 200 2\n"
            "100\n"
        )
        with tempfile.TemporaryDirectory() as directory:
            chain = Path(directory) / "ambiguous.chain.gz"
            chain.write_bytes(gzip.compress(payload.encode("ascii"), mtime=0))
            index, _ = liftover.load_chain(chain)
            status, mapped = index.map_point(1, 1)
        self.assertEqual(status, "ambiguous")
        self.assertIsNone(mapped)

    def test_new_sleep_source_headers_harmonize_with_documented_effects(self):
        cases = {
            "insomnia": (
                "SNP\tUNIQUE_ID\tCHR\tBP\tA1\tA2\tMAF\tOR\tSE\tP\tN\tINFO\n"
                "rs123\t1:1000000:A_C\t1\t1000000\tA\tC\t0.2\t1.2\t0.1\t0.01\t380000\t0.99\n",
                math.log(1.2),
            ),
            "chronotype": (
                "SNP\tCHR\tBP\tALLELE1\tALLELE0\tA1FREQ\tINFO\tLOGOR\tLOGOR_SE\tP_BOLT_LMM\tHWE_P\n"
                "rs123\t1\t1000000\tA\tC\t0.2\t0.99\t0.2\t0.1\t0.01\t0.9\n",
                0.2,
            ),
        }
        with tempfile.TemporaryDirectory() as directory:
            directory = Path(directory)
            for trait, (payload, expected_beta) in cases.items():
                with self.subTest(trait=trait):
                    source = directory / f"{trait}.tsv"
                    outdir = directory / f"{trait}-out"
                    source.write_text(payload, encoding="utf-8")
                    result = subprocess.run(
                        [
                            sys.executable,
                            str(ROOT / "scripts" / "01_harmonize.py"),
                            "--trait", trait,
                            "--config", str(MANIFEST),
                            "--infile", str(source),
                            "--outdir", str(outdir),
                            "--source-build", "hg19",
                        ],
                        cwd=ROOT,
                        capture_output=True,
                        text=True,
                    )
                    self.assertEqual(result.returncode, 0, result.stderr + result.stdout)
                    with gzip.open(
                        outdir / f"{trait}.harmonized.tsv.gz", "rt", newline=""
                    ) as handle:
                        rows = list(csv.DictReader(handle, delimiter="\t"))
                    self.assertEqual(len(rows), 1)
                    self.assertAlmostEqual(float(rows[0]["BETA"]), expected_beta)

    def test_pgc_metadata_and_neffdiv2_harmonize_by_named_fields(self):
        payload = (
            '##fileFormat=PGCsumstatsVCFv1.0\n'
            '##genomeReference="GRCh37"\n'
            '#CHROM\tPOS\tID\tA1\tA2\tBETA\tSE\tPVAL\tFCON\tIMPINFO\tNEFFDIV2\tNCAS\tNCON\n'
            '1\t1000000\trs123\tA\tC\t0.2\t0.1\t0.01\t0.2\t0.99\t50981\t41917\t371549\n'
        )
        with tempfile.TemporaryDirectory() as directory:
            directory = Path(directory)
            source = directory / "bipolar.tsv.gz"
            outdir = directory / "out"
            with gzip.open(source, "wt", encoding="utf-8", newline="") as handle:
                handle.write(payload)
            result = subprocess.run(
                [
                    sys.executable,
                    str(ROOT / "scripts" / "01_harmonize.py"),
                    "--trait", "bipolar",
                    "--config", str(MANIFEST),
                    "--infile", str(source),
                    "--outdir", str(outdir),
                    "--source-build", "hg19",
                ],
                cwd=ROOT,
                capture_output=True,
                text=True,
            )
            self.assertEqual(result.returncode, 0, result.stderr + result.stdout)
            with gzip.open(outdir / "bipolar.harmonized.tsv.gz", "rt", newline="") as handle:
                rows = list(csv.DictReader(handle, delimiter="\t"))
        self.assertEqual(len(rows), 1)
        self.assertEqual(rows[0]["SNP"], "rs123")
        self.assertEqual(rows[0]["FRQ"], "0.2")
        self.assertEqual(float(rows[0]["N"]), 101962.0)

    def test_scz_release_neff_is_doubled_by_registered_schema(self):
        payload = (
            '##fileFormat=PGCsumstatsVCFv1.0\n'
            '##genomeReference="GRCh37"\n'
            'CHROM\tID\tPOS\tA1\tA2\tFCON\tIMPINFO\tBETA\tSE\tPVAL\tNCAS\tNCON\tNEFF\n'
            '1\trs123\t1000000\tA\tC\t0.2\t0.99\t0.2\t0.1\t0.01\t53386\t77258\t58749.13\n'
        )
        with tempfile.TemporaryDirectory() as directory:
            directory = Path(directory)
            source = directory / "scz.tsv.gz"
            outdir = directory / "out"
            with gzip.open(source, "wt", encoding="utf-8", newline="") as handle:
                handle.write(payload)
            result = subprocess.run(
                [
                    sys.executable,
                    str(ROOT / "scripts" / "01_harmonize.py"),
                    "--trait", "scz",
                    "--config", str(MANIFEST),
                    "--infile", str(source),
                    "--outdir", str(outdir),
                    "--source-build", "hg19",
                ],
                cwd=ROOT,
                capture_output=True,
                text=True,
            )
            self.assertEqual(result.returncode, 0, result.stderr + result.stdout)
            with gzip.open(outdir / "scz.harmonized.tsv.gz", "rt", newline="") as handle:
                rows = list(csv.DictReader(handle, delimiter="\t"))
            qc = (outdir / "scz.qc.txt").read_text(encoding="utf-8")
        self.assertEqual(len(rows), 1)
        self.assertAlmostEqual(float(rows[0]["N"]), 117498.26)
        self.assertIn("sample_size_schema_mapping\tDERIVED_N_EFF=2*NEFF", qc)
        self.assertIn("sample-size mode: per-SNP N_eff derived as 2 * source NEFF", qc)

    def test_opengwas_sbp_vcf_materializes_by_named_format_fields(self):
        materializer = load_sbp_materializer()
        vcf_url = "https://example.org/public/ieu-b-38/ieu-b-38.vcf.gz?token=short-lived"
        self.assertEqual(
            materializer.select_vcf_url(
                {"files": [vcf_url, vcf_url.replace(".vcf.gz?", ".vcf.gz.tbi?")]},
                "ieu-b-38",
            ),
            vcf_url,
        )
        metadata = [
            "##fileformat=VCFv4.2",
            "##contig=<ID=1,length=249250621,assembly=HG19/GRCh37>",
            "##SAMPLE=<ID=ieu-b-38,TotalVariants=3,VariantsNotRead=0,HarmonisedVariants=3,VariantsNotHarmonised=0,SwitchedAlleles=1,NormalisedVariants=0,StudyType=Continuous>",
        ]
        for format_id in ["ES", "SE", "LP", "AF", "SS", "ID"]:
            metadata.append(
                f'##FORMAT=<ID={format_id},Number=1,Type=String,Description="fixture, {format_id}">'
            )
        payload = "\n".join(metadata) + "\n" + (
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tieu-b-38\n"
            "1\t752566\trs3094315\tG\tA\t.\tPASS\tAF=0.8363\t"
            "ES:SE:LP:AF:SS:ID\t0.0687:0.0447:0.906578:0.8363:639410:rs3094315\n"
            "1\t1000000\trs123\tC\tT\t.\tPASS\tAF=0.25\t"
            "ID:SS:AF:LP:SE:ES\trs123:757601:0.25:350.5:0.02:-0.1\n"
            "1\t1000100\t.\tG\tA\t.\tPASS\tAF=0.2\t"
            "ES:SE:LP:AF:SS:ID\t0.1:0.02:2:0.2:757601:.\n"
        )
        with tempfile.TemporaryDirectory() as directory:
            directory = Path(directory)
            source = directory / "ieu-b-38.vcf.gz"
            first = directory / "sbp-first.txt.gz"
            second = directory / "sbp-second.txt.gz"
            with gzip.open(source, "wt", encoding="utf-8", newline="") as handle:
                handle.write(payload)
            counts = materializer.materialize_vcf(
                source, first, expected_variants=3
            )
            materializer.materialize_vcf(source, second, expected_variants=3)
            with gzip.open(first, "rt", newline="") as handle:
                rows = list(csv.DictReader(handle, delimiter="\t"))
            self.assertEqual(first.read_bytes(), second.read_bytes())
        self.assertEqual(counts, {
            "source_rows": 3,
            "materialized_rows": 2,
            "non_rsid_rows": 1,
            "p_values_floored_at_1e-300": 1,
        })
        self.assertEqual(rows[0]["SNP"], "rs3094315")
        self.assertEqual((rows[0]["A1"], rows[0]["A2"]), ("A", "G"))
        self.assertAlmostEqual(float(rows[0]["P"]), 10 ** -0.906578)
        self.assertEqual(rows[1]["P"], "1.000000000000000E-300")

    def test_shell_trait_lookup_is_header_aware(self):
        result = subprocess.run(
            [
                "bash",
                "-c",
                "source scripts/_common.sh && trait_field snoring source_id",
            ],
            cwd=ROOT,
            capture_output=True,
            text=True,
        )
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertEqual(result.stdout.strip(), "campos_2020_snoring")

    def test_munge_accepts_empty_optional_argument_arrays(self):
        """Empty arrays work on Bash 3.2 and an absent source FRQ is ignored."""
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            for name in ["scripts", "config", "data/raw", "ref/eur_w_ld_chr", "ldsc"]:
                (root / name).mkdir(parents=True, exist_ok=True)
            for name in ["02_munge.sh", "_common.sh"]:
                shutil.copy2(ROOT / "scripts" / name, root / "scripts" / name)
            (root / "config" / "analysis_panel.tsv").write_text(
                "atlas_version\ttrait_id\tsource_id\traw_file\tbuild\t"
                "source_status\tancestry\n"
                "atlas-v1.0\tfixture\tfixture_source\tfixture.txt.gz\thg19\t"
                "SOURCE_VERIFIED\tEUR\n"
                "atlas-v1.0\tmapped_fixture\tmapped_source\tmapped.txt.gz\t"
                "UNRESOLVED\tSOURCE_VERIFIED\tEUR\n",
                encoding="utf-8",
            )
            (root / "config" / "variant_mapping_plans.tsv").write_text(
                "trait_id\tstrategy\tmap_path\tmap_bytes\tmap_sha256\n"
                "mapped_fixture\tBY_RSID_ALLELES\tref/map.tsv.gz\t1\tabc\n",
                encoding="utf-8",
            )
            for path in [
                root / "config" / "analysis_panel.lock.json",
                root / "data" / "raw" / "fixture.txt.gz",
                root / "data" / "raw" / "mapped.txt.gz",
                root / "ref" / "w_hm3.snplist",
                root / "ref" / "eur_w_ld_chr" / "1.l2.ldscore.gz",
                root / "ref" / "map.tsv.gz",
                root / "ref" / "map.tsv.gz.provenance.json",
                root / "ldsc" / "ldsc.py",
                root / "ldsc" / "munge_sumstats.py",
                root / "scripts" / "00_validate_panel.py",
                root / "scripts" / "01_harmonize.py",
            ]:
                path.touch()
            runner = root / "fake-python"
            runner.write_text(
                "#!/usr/bin/env bash\n"
                "set -euo pipefail\n"
                "target=$1; shift\n"
                "case \"$target\" in\n"
                "  scripts/00_validate_panel.py) exit 0 ;;\n"
                "  scripts/01_harmonize.py)\n"
                "    trait=; outdir=\n"
                "    while [ $# -gt 0 ]; do\n"
                "      case \"$1\" in\n"
                "        --trait) trait=$2; shift 2 ;;\n"
                "        --outdir) outdir=$2; shift 2 ;;\n"
                "        *) shift ;;\n"
                "      esac\n"
                "    done\n"
                "    mkdir -p \"$outdir\"\n"
                "    : > \"$outdir/$trait.harmonized.tsv.gz\"\n"
                "    printf '%s\\n' 'FRQ column absent - source-level MAF QC must be documented' > \"$outdir/$trait.qc.txt\" ;;\n"
                "  ldsc/munge_sumstats.py)\n"
                "    out=; ignored=\n"
                "    while [ $# -gt 0 ]; do\n"
                "      case \"$1\" in\n"
                "        --out) out=$2; shift 2 ;;\n"
                "        --ignore) ignored=$2; shift 2 ;;\n"
                "        *) shift ;;\n"
                "      esac\n"
                "    done\n"
                "    [ \"$ignored\" = FRQ ] || exit 3\n"
                "    mkdir -p \"$(dirname \"$out\")\"\n"
                "    : > \"$out.sumstats.gz\"\n"
                "    : > \"$out.log\" ;;\n"
                "  *) exit 2 ;;\n"
                "esac\n",
                encoding="utf-8",
            )
            runner.chmod(0o755)
            result = subprocess.run(
                ["/bin/bash", "scripts/02_munge.sh", "fixture", "mapped_fixture"],
                cwd=root,
                env={
                    **os.environ,
                    "PYTHON_BIN": str(runner),
                    "LDSC_PYTHON": str(runner),
                    "LDSC_DIR": "ldsc",
                },
                capture_output=True,
                text=True,
            )
            self.assertEqual(result.returncode, 0, result.stderr + result.stdout)
            self.assertTrue((root / "data/munged/fixture.sumstats.gz").is_file())
            self.assertTrue((root / "data/munged/mapped_fixture.sumstats.gz").is_file())

    def test_harmonizer_does_not_load_unmapped_annotation_columns(self):
        script = (ROOT / "scripts" / "01_harmonize.py").read_text(encoding="utf-8")
        self.assertIn('read_kwargs["usecols"] = selected_source_columns', script)
        self.assertIn('raw.rename(columns=source_to_standard, inplace=True)', script)
        self.assertIn('"pos_b37"', script)

    def test_streaming_hm3_prefilter_preserves_selected_source_rows(self):
        with tempfile.TemporaryDirectory() as directory:
            directory = Path(directory)
            source = directory / "source.tsv.gz"
            allowlist = directory / "hm3.tsv"
            output = directory / "filtered.tsv.gz"
            provenance = directory / "filtered.provenance.json"
            payload = (
                b"rsID\tvalue\tnote\n"
                b"rs1\t1.00\tkeep exact text\n"
                b"rs2\t2e-3\tdrop\n"
                b"RS3\tNA\tkeep case-insensitively\n"
            )
            with gzip.open(source, "wb") as handle:
                handle.write(payload)
            allowlist.write_text("SNP\tA1\tA2\nrs1\tA\tG\nrs3\tC\tT\n", encoding="utf-8")
            source_bytes = source.stat().st_size
            source_sha = hashlib.sha256(source.read_bytes()).hexdigest()
            result = subprocess.run(
                [
                    sys.executable, str(ROOT / "scripts" / "21_prefilter_hm3.py"),
                    "--input", str(source), "--output", str(output),
                    "--provenance", str(provenance), "--snp-column", "rsID",
                    "--allowlist", str(allowlist), "--expected-input-bytes",
                    str(source_bytes), "--expected-input-sha256", source_sha,
                ],
                cwd=ROOT, capture_output=True, text=True,
            )
            self.assertEqual(result.returncode, 0, result.stderr + result.stdout)
            reused = subprocess.run(
                [
                    sys.executable, str(ROOT / "scripts" / "21_prefilter_hm3.py"),
                    "--input", str(source), "--output", str(output),
                    "--provenance", str(provenance), "--snp-column", "rsID",
                    "--allowlist", str(allowlist), "--expected-input-bytes",
                    str(source_bytes), "--expected-input-sha256", source_sha,
                ],
                cwd=ROOT, capture_output=True, text=True,
            )
            self.assertEqual(reused.returncode, 0, reused.stderr + reused.stdout)
            self.assertIn("Reusing verified prefilter", reused.stdout)
            with gzip.open(output, "rb") as handle:
                filtered = handle.read()
            ledger = json.loads(provenance.read_text(encoding="utf-8"))
        self.assertEqual(
            filtered,
            payload.splitlines(keepends=True)[0]
            + payload.splitlines(keepends=True)[1]
            + payload.splitlines(keepends=True)[3],
        )
        self.assertEqual(ledger["source_rows"], 3)
        self.assertEqual(ledger["retained_rows"], 2)

    def test_hm3_prefilter_plan_is_exactly_the_registered_large_traits(self):
        with (ROOT / "config" / "hm3_prefilter_plans.tsv").open(
            newline="", encoding="utf-8"
        ) as handle:
            rows = list(csv.DictReader(handle, delimiter="\t"))
        self.assertEqual(
            {row["trait_id"] for row in rows},
            {
                "ldl", "hdl", "triglycerides", "ms", "asthma", "t2d",
                "cad", "telomere_length", "melanoma",
            },
        )
        self.assertTrue(all(row["strategy"] == "HAPMAP3_RSID_ALLOWLIST" for row in rows))

    def test_retained_prefilter_binds_removed_raw_to_registry(self):
        audit = load_phase0_audit()
        with tempfile.TemporaryDirectory() as directory:
            directory = Path(directory)
            harmonized = directory / "harmonized"
            (harmonized / ".prefilter").mkdir(parents=True)
            allowlist = directory / "w_hm3.snplist"
            allowlist.write_text("SNP\nrs123\n", encoding="ascii")
            output = harmonized / ".prefilter" / "hdl.hm3.tsv.gz"
            output.write_bytes(gzip.compress(b"rsID\tBETA\nrs123\t0.1\n", mtime=0))
            source_sha = "a" * 64
            provenance = {
                "strategy": "HAPMAP3_RSID_ALLOWLIST",
                "input_bytes": 123456,
                "input_sha256": source_sha,
                "snp_column": "rsID",
                "allowlist_path": str(allowlist),
                "allowlist_sha256": hashlib.sha256(allowlist.read_bytes()).hexdigest(),
                "source_rows": 10,
                "retained_rows": 1,
                "output_bytes": output.stat().st_size,
                "output_sha256": hashlib.sha256(output.read_bytes()).hexdigest(),
            }
            (harmonized / "hdl.prefilter.provenance.json").write_text(
                json.dumps(provenance), encoding="utf-8"
            )
            plan = {
                "strategy": "HAPMAP3_RSID_ALLOWLIST",
                "snp_column": "rsID",
                "allowlist_path": str(allowlist),
            }
            source = {
                "archive_member": "DIRECT_GZIP",
                "raw_files": "hdl.txt.gz",
                "archive_bytes": "123456",
                "archive_sha256": source_sha,
            }
            trait = {"trait_id": "hdl", "raw_file": "hdl.txt.gz"}
            self.assertEqual(
                audit.retained_prefilter_issues(
                    plan, source, trait, str(harmonized)
                ),
                [],
            )
            source["archive_sha256"] = "b" * 64
            self.assertIn(
                "retained_prefilter_source_checksum_mismatch",
                audit.retained_prefilter_issues(
                    plan, source, trait, str(harmonized)
                ),
            )

    def test_production_code_has_no_historic_manifest_reference(self):
        paths = list((ROOT / "scripts").glob("*.py"))
        paths += list((ROOT / "scripts").glob("*.sh"))
        paths += [ROOT / "Snakefile"]
        offenders = [
            str(path.relative_to(ROOT))
            for path in paths
            if "config/traits.tsv" in path.read_text(encoding="utf-8")
        ]
        self.assertEqual(offenders, [])

    def test_mixer_proxy_is_not_called_mixer_pass(self):
        paths = list((ROOT / "scripts").glob("*.py"))
        offenders = [
            str(path.relative_to(ROOT))
            for path in paths
            if "mixer_pass" in path.read_text(encoding="utf-8")
        ]
        self.assertEqual(offenders, [])

    def test_h2_collator_rejects_out_of_panel_log(self):
        collator = load_collator()
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "h2_intruder.log"
            path.write_text(
                "Total Observed scale h2: 0.10 (0.01)\nIntercept: 1.01 (0.01)\n",
                encoding="utf-8",
            )
            with self.assertRaises(SystemExit):
                collator.parse_h2(directory, str(MANIFEST))

    def test_h2_wrapper_rebuilds_the_cumulative_table(self):
        script = (ROOT / "scripts" / "03_h2_qc.sh").read_text(encoding="utf-8")
        collate_call = script.split('scripts/05_collate.py --mode h2', 1)[1]
        self.assertNotIn('--traits "$@"', collate_call)

    def test_collator_fails_closed_without_manifest(self):
        collator = load_collator()
        with tempfile.TemporaryDirectory() as directory:
            with self.assertRaises(SystemExit):
                collator.parse_h2(directory, str(Path(directory) / "missing.tsv"))

    def test_rg_collator_rejects_out_of_panel_pair(self):
        collator = load_collator()
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "rg_intruder.log"
            path.write_text(
                "Summary of Genetic Correlation Results\n"
                "p1 p2 rg se z p\n"
                "data/munged/intruder.sumstats.gz data/munged/mdd.sumstats.gz "
                "0.1 0.02 5 1e-6\n\n",
                encoding="utf-8",
            )
            with self.assertRaises(SystemExit):
                collator.parse_rg(directory, str(MANIFEST))

    def test_rg_collator_preserves_scalar_p_values_for_multi_pair_log(self):
        collator = load_collator()
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "rg_snoring.log"
            path.write_text(
                "P: 1.0292e-90\n"
                "P: 0.4283\n"
                "Summary of Genetic Correlation Results\n"
                "p1 p2 rg se z p\n"
                "data/munged/snoring.sumstats.gz data/munged/bmi.sumstats.gz "
                "0.3779 0.0187 20.1975 0.0000\n"
                "data/munged/snoring.sumstats.gz data/munged/longevity.sumstats.gz "
                "-0.0444 0.0561 -0.7922 0.4283\n\n",
                encoding="utf-8",
            )
            rows = collator.parse_rg(directory, str(MANIFEST))
        by_trait = rows.set_index("disease_trait")
        self.assertEqual(by_trait.loc["bmi", "p"], 1.0292e-90)
        self.assertEqual(by_trait.loc["longevity", "p"], 0.4283)
        self.assertGreater(by_trait.loc["bmi", "fdr"], 0)

    def test_rg_matrix_collation_ignores_controlled_pair_logs(self):
        collator = load_collator()
        with tempfile.TemporaryDirectory() as directory:
            directory = Path(directory)
            inclusion = directory / "inclusion.tsv"
            inclusion.write_text(
                "trait_id\tdomain\tinclude_phase1\n"
                "snoring\tsleep\tTrue\n"
                "bmi\tmetabolic\tTrue\n",
                encoding="utf-8",
            )
            body = (
                "Summary of Genetic Correlation Results\n"
                "p1 p2 rg se z p\n"
                "data/munged/snoring.sumstats.gz data/munged/bmi.sumstats.gz "
                "0.3 0.02 15 1e-8\n\n"
            )
            (directory / "rg_snoring.log").write_text(body, encoding="utf-8")
            (directory / "rg_snoring__bmi.log").write_text(body, encoding="utf-8")
            rows = collator.parse_rg(directory, str(MANIFEST), str(inclusion))
        self.assertEqual(len(rows), 1)
        self.assertEqual(rows.iloc[0]["sleep_trait"], "snoring")
        self.assertEqual(rows.iloc[0]["disease_trait"], "bmi")

    def test_complete_rg_merge_labels_qc_failed_sensitivity(self):
        with MANIFEST.open(newline="", encoding="utf-8") as handle:
            panel = list(csv.DictReader(handle, delimiter="\t"))
        sleeps = [row["trait_id"] for row in panel if row["domain"] == "sleep"]
        diseases = [row["trait_id"] for row in panel if row["domain"] != "sleep"]
        failed = diseases[-1]
        with tempfile.TemporaryDirectory() as directory:
            directory = Path(directory)
            h2 = directory / "h2.tsv"
            primary = directory / "primary.tsv"
            sensitivity = directory / "sensitivity.tsv"
            output = directory / "out.tsv"
            with h2.open("w", newline="", encoding="utf-8") as handle:
                writer = csv.DictWriter(
                    handle, fieldnames=["trait", "verdict", "qc_reason"],
                    delimiter="\t",
                )
                writer.writeheader()
                for row in panel:
                    is_failed = row["trait_id"] == failed
                    writer.writerow({
                        "trait": row["trait_id"],
                        "verdict": "DROP" if is_failed else "PASS",
                        "qc_reason": "fixture_failure" if is_failed else "pass",
                    })
            fields = ["sleep_trait", "disease_trait", "rg", "se", "p", "fdr"]
            for path, selected in [
                (primary, [disease for disease in diseases if disease != failed]),
                (sensitivity, [failed]),
            ]:
                with path.open("w", newline="", encoding="utf-8") as handle:
                    writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
                    writer.writeheader()
                    for sleep in sleeps:
                        for disease in selected:
                            writer.writerow({
                                "sleep_trait": sleep, "disease_trait": disease,
                                "rg": 0.1, "se": 0.05, "p": 0.01, "fdr": 0.02,
                            })
            result = subprocess.run(
                [
                    sys.executable, str(ROOT / "scripts" / "23_merge_rg_families.py"),
                    "--config", str(MANIFEST), "--h2", str(h2),
                    "--primary", str(primary), "--sensitivity", str(sensitivity),
                    "--out", str(output),
                ],
                cwd=ROOT, capture_output=True, text=True,
            )
            self.assertEqual(result.returncode, 0, result.stderr + result.stdout)
            with output.open(newline="", encoding="utf-8") as handle:
                rows = list(csv.DictReader(handle, delimiter="\t"))
        self.assertEqual(len(rows), 396)
        failed_rows = [row for row in rows if row["disease_trait"] == failed]
        self.assertEqual(len(failed_rows), 12)
        self.assertTrue(all(
            row["analysis_tier"] == "QC_FAILED_SENSITIVITY"
            and row["interpretation_status"] == "EXCLUDED_FROM_PRIMARY_INFERENCE"
            and row["disease_h2_qc_reason"] == "fixture_failure"
            for row in failed_rows
        ))

    def test_acceptance_audit_distinguishes_contract_from_science(self):
        result = subprocess.run(
            [
                sys.executable,
                str(ROOT / "scripts" / "99_atlas_acceptance.py"),
                "--report-only",
            ],
            cwd=ROOT,
            capture_output=True,
            text=True,
        )
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertIn("PASS    locked_panel", result.stdout)
        self.assertIn("PASS    source_curation", result.stdout)
        self.assertIn("PASS    source_schemas", result.stdout)
        gate_lines = [
            line for line in result.stdout.splitlines()
            if line.startswith(("PASS", "BLOCKED"))
        ]
        self.assertEqual(len(gate_lines), 23)
        self.assertIn("Acceptance:", result.stdout)
        self.assertIn("/23 gates passed", result.stdout)

    def test_full_covariance_runtime_is_locked_to_real_45_trait_inputs(self):
        script = (ROOT / "scripts" / "25_genomicsem_covariance.R").read_text(encoding="utf-8")
        setup = (ROOT / "scripts" / "25_setup_genomicsem.sh").read_text(encoding="utf-8")
        env = (ROOT / "environment" / "genomicsem.yml").read_text(encoding="utf-8")
        self.assertIn('nrow(panel) != 45L', script)
        self.assertIn('sum(panel$domain == "sleep") != 12L', script)
        self.assertIn('unique(panel$panel_version) != "atlas-v1.0"', script)
        self.assertIn('getRversion() != "4.3.3"', script)
        self.assertIn('GenomicSEM::ldsc(', script)
        self.assertIn('stand = TRUE', script)
        self.assertIn('ldsc_sampling_covariance_1035x1035.tsv.gz', script)
        self.assertIn('genetic_correlation_se = se_correlation', script)
        self.assertIn('dim(covstruc$V_Stand)', script)
        self.assertIn('(length(trait_ids) + 2L) / 2L) + 1L', script)
        validator = (ROOT / "scripts" / "26_validate_covariance.py").read_text(encoding="utf-8")
        self.assertIn('"jackknife_blocks": "1082"', validator)
        self.assertIn('genetic_correlation_fdr_off_diagonal', validator)
        sem_qc = (ROOT / "scripts" / "27_genomicsem_trait_qc.py").read_text(encoding="utf-8")
        self.assertIn('INTERCEPT_MAX = 1.20', sem_qc)
        self.assertIn('phase1_vs_genomicsem_implementation_sensitivity', sem_qc)
        split = (ROOT / "scripts" / "28_chromosome_split_covariance.R").read_text(encoding="utf-8")
        self.assertIn('sum(selected) != 42L', split)
        self.assertIn('GenomicSEM::ldsc(', split)
        model = (ROOT / "scripts" / "29_genomicsem_model.R").read_text(encoding="utf-8")
        self.assertIn('cfi >= 0.90 && srmr <= 0.10', model)
        self.assertIn('discovery_heywood', model)
        acceptance = (ROOT / "scripts" / "99_atlas_acceptance.py").read_text(encoding="utf-8")
        self.assertIn('no candidate passed held-out validation', acceptance)
        self.assertIn('6b65ca5db39fdade08b0d811477be1cdd57b5039', setup)
        self.assertIn('simsalapar_1.0-13.tar.gz', setup)
        self.assertIn('r-base=4.3.3', env)

    def test_lava_policy_and_reference_plan_are_locked(self):
        policy = json.loads((ROOT / "config" / "lava_analysis_policy.json").read_text(encoding="utf-8"))
        self.assertEqual(policy["lava_version"], "0.1.5")
        self.assertEqual(policy["lava_commit"], "e729a245f7b6923967a96804fbf5246eadf2d6c6")
        self.assertEqual(policy["expected_loci"], 2495)
        self.assertEqual(policy["expected_traits"], 45)
        self.assertEqual(policy["expected_sleep_non_sleep_pairs"], 396)
        self.assertEqual(policy["planned_univariate_tests"], 112275)
        self.assertAlmostEqual(
            policy["univariate_p_threshold"],
            policy["univariate_alpha"] / policy["planned_univariate_tests"],
        )
        with (ROOT / "config" / "lava_reference_sources.tsv").open(newline="", encoding="utf-8") as handle:
            sources = list(csv.DictReader(handle, delimiter="\t"))
        self.assertEqual(len(sources), 7)
        self.assertEqual(sum(int(row["archive_bytes"]) for row in sources), 14110596095)

    def test_lava_runtime_fails_closed_and_does_not_filter_on_global_rg(self):
        setup = (ROOT / "scripts" / "30_setup_lava.sh").read_text(encoding="utf-8")
        downloader = (ROOT / "scripts" / "32_download_lava_reference.sh").read_text(encoding="utf-8")
        prepare = (ROOT / "scripts" / "31_prepare_lava.py").read_text(encoding="utf-8")
        runtime = (ROOT / "scripts" / "33_run_lava.R").read_text(encoding="utf-8")
        validator = (ROOT / "scripts" / "34_validate_lava.py").read_text(encoding="utf-8")
        acceptance = (ROOT / "scripts" / "99_atlas_acceptance.py").read_text(encoding="utf-8")
        self.assertIn("15477c9547c3533d681cbfb6491508a29342fe99b725845054887b9d29dc2bec", setup)
        self.assertIn("The 15 GiB UK Biobank LD reference is deliberately not downloaded", setup)
        self.assertIn("DOWNLOAD=false", downloader)
        self.assertIn('if [ "$DOWNLOAD" != true ]', downloader)
        self.assertIn('len(rg) != 396 or observed != expected', prepare)
        self.assertNotIn('global_rg_p <=', runtime)
        self.assertIn('univ$p <= numeric_policy("univariate_p_threshold")', runtime)
        self.assertIn('run.bivar(locus', runtime)
        self.assertIn('expected_bivar', validator)
        self.assertIn('scripts/34_validate_lava.py', acceptance)

    def test_mixer_policy_is_version_and_scope_locked(self):
        policy = json.loads(
            (ROOT / "config" / "mixer_analysis_policy.json").read_text(encoding="utf-8")
        )
        self.assertEqual(policy["mixer_release"], "2.2.1")
        self.assertEqual(
            policy["mixer_release_tag_commit"],
            "cf65c57d5d1ad76597db1d4fa3907f1d711d39e7",
        )
        self.assertEqual(
            policy["container_amd64_manifest_digest"],
            "sha256:5bf54ddd6f7f81b93eeb5450b1a6dc809d1fcf926b3a815d6d514f36ce04f51c",
        )
        self.assertEqual(policy["expected_traits"], 45)
        self.assertEqual(policy["expected_sleep_non_sleep_pairs"], 396)
        self.assertEqual(policy["fit_replicates"], 20)
        self.assertEqual(policy["univariate_aic_threshold"], 0)
        self.assertIn("HapMap3-prefiltered inputs are forbidden", policy["input_requirement"])

    def test_mixer_runtime_is_fail_closed_and_has_no_implicit_pull(self):
        preflight = (ROOT / "scripts" / "35_mixer_preflight.py").read_text(encoding="utf-8")
        prepare = (ROOT / "scripts" / "36_prepare_mixer_inputs.py").read_text(encoding="utf-8")
        pull = (ROOT / "scripts" / "37_pull_mixer_image.sh").read_text(encoding="utf-8")
        runtime = (ROOT / "scripts" / "38_run_mixer_task.sh").read_text(encoding="utf-8")
        collate = (ROOT / "scripts" / "39_collate_mixer.py").read_text(encoding="utf-8")
        validator = (ROOT / "scripts" / "40_validate_mixer.py").read_text(encoding="utf-8")
        acceptance = (ROOT / "scripts" / "99_atlas_acceptance.py").read_text(encoding="utf-8")
        self.assertIn('machine in policy["supported_architectures"]', preflight)
        self.assertIn('memory_bytes >= policy["minimum_memory_bytes"]', preflight)
        self.assertIn('physical_cores >= policy["recommended_physical_cores"]', preflight)
        self.assertIn('prefilter == "not supplied"', preflight)
        self.assertIn('prefilter != "not supplied"', prepare)
        self.assertIn('if [ "$PULL" != true ]', pull)
        self.assertNotIn("docker pull", runtime)
        self.assertIn("for rep in $(seq 1 20)", runtime)
        self.assertNotIn("global_rg", runtime)
        self.assertIn('policy["univariate_aic_threshold"]', collate)
        self.assertIn("expected_pairs", validator)
        self.assertIn("scripts/40_validate_mixer.py", acceptance)

    def test_mixer_input_conversion_is_atomic_and_deterministic(self):
        spec = importlib.util.spec_from_file_location(
            "mixer_prepare", ROOT / "scripts" / "36_prepare_mixer_inputs.py"
        )
        module = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(module)
        with tempfile.TemporaryDirectory() as directory:
            directory = Path(directory)
            source = directory / "source.tsv.gz"
            header = "SNP\tCHR\tBP\tA1\tA2\tN\tBETA\tSE\n"
            body = "rs1\t1\t101\tA\tG\t1000\t0.2\t0.1\nrs2\t2\t202\tC\tT\t800\t-0.3\t0.2\n"
            with gzip.open(source, "wt", encoding="utf-8", newline="") as handle:
                handle.write(header + body)
            first = directory / "first.sumstats.gz"
            second = directory / "second.sumstats.gz"
            metrics = module.convert(source, first)
            module.convert(source, second)
            self.assertEqual(metrics["rows"], 2)
            self.assertEqual(first.read_bytes(), second.read_bytes())
            with gzip.open(first, "rt", encoding="utf-8", newline="") as handle:
                rows = list(csv.DictReader(handle, delimiter="\t"))
            self.assertEqual(list(rows[0]), module.OUTPUT_COLUMNS)
            self.assertEqual([float(row["Z"]) for row in rows], [2.0, -1.5])

            invalid = directory / "invalid.tsv.gz"
            with gzip.open(invalid, "wt", encoding="utf-8", newline="") as handle:
                handle.write(header + "rs3\t3\t303\tA\tC\t900\t0.1\t0\n")
            failed = directory / "failed.sumstats.gz"
            with self.assertRaises(SystemExit):
                module.convert(invalid, failed)
            self.assertFalse(failed.exists())
            self.assertFalse(Path(str(failed) + ".tmp").exists())

    def test_pleiotropy_policy_locks_both_methods_and_all_pairs(self):
        policy = json.loads(
            (ROOT / "config" / "pleiotropy_analysis_policy.json").read_text(
                encoding="utf-8"
            )
        )
        self.assertEqual(policy["expected_sleep_non_sleep_pairs"], 396)
        self.assertIn("All locked 396", policy["pair_scope"])
        self.assertEqual(policy["placo_version"], "0.2.0")
        self.assertEqual(
            policy["placo_source_sha256"],
            "fb684a8ed88f27dd138f5c2e8613b092904e84f364030d036058f17db95b7124",
        )
        self.assertEqual(policy["pleiofdr_mode"], "conjfdr")
        self.assertEqual(policy["pleiofdr_random_prune_iterations"], 500)
        self.assertEqual(policy["pleiofdr_reference_bytes"], 2383912974)
        self.assertEqual(policy["pleiofdr_variant_template_variants"], 9545380)
        self.assertEqual(policy["pleiofdr_variant_template_bytes"], 274423819)
        self.assertEqual(
            policy["pleiofdr_variant_template_sha256"],
            "06268420a0ec04e4529e832e1d4f4a53231b078cc5a3741a3eb215a5a4e1a9d5",
        )
        self.assertTrue(policy["pleiofdr_correct_sample_overlap"])
        self.assertEqual(
            policy["pleiofdr_overlap_patch_sha256"],
            "a5ea51cafd08b783903733f6485272fd8d3563664b060ad31ae8a9266dced200",
        )
        self.assertAlmostEqual(
            policy["placo_locked_pair_family_threshold"], 5e-8 / 396
        )
        self.assertIn("PLACO+ association and conjFDR", policy["canonical_rule"])

    def test_pleiotropy_setup_and_manifest_are_fail_closed(self):
        setup = (ROOT / "scripts" / "41_setup_pleiotropy.sh").read_text(
            encoding="utf-8"
        )
        preflight = (ROOT / "scripts" / "42_pleiotropy_preflight.py").read_text(
            encoding="utf-8"
        )
        prepare = (ROOT / "scripts" / "43_prepare_pleiotropy_pairs.py").read_text(
            encoding="utf-8"
        )
        self.assertIn("DOWNLOAD_REFERENCE=false", setup)
        self.assertIn('elif [ "$DOWNLOAD_REFERENCE" = true ]', setup)
        self.assertIn("REFERENCE_BYTES=2383912974", setup)
        self.assertIn("TEMPLATE_BYTES=274423819", setup)
        self.assertIn("TEMPLATE_SHA=06268420a0ec04e4529e832e1d4f4a53231b078cc5a3741a3eb215a5a4e1a9d5", setup)
        self.assertIn('prefilter == "not supplied"', preflight)
        self.assertIn('matlab = shutil.which("matlab")', preflight)
        self.assertIn("ready_pair_scans", preflight)
        self.assertIn("set(rg_by_pair) != set(expected)", prepare)
        self.assertIn('"global_rg_filters_pair_eligibility": False', prepare)
        self.assertNotIn("global_rg_p <=", prepare)

    def test_pleiotropy_pair_materializer_aligns_and_rejects_bad_rows(self):
        module = load_pleiotropy_materializer()
        self.assertEqual(module.align("A", "C", "A", "C", 2.0), (2.0, "DIRECT"))
        self.assertEqual(module.align("A", "C", "C", "A", 2.0), (-2.0, "SWAPPED"))
        self.assertEqual(module.align("A", "C", "T", "G", 2.0), (2.0, "COMPLEMENT"))
        self.assertEqual(
            module.align("A", "C", "G", "T", 2.0),
            (-2.0, "COMPLEMENT_SWAPPED"),
        )
        self.assertIsNone(module.align("A", "C", "A", "G", 2.0))
        with tempfile.TemporaryDirectory() as directory:
            source = Path(directory) / "source.tsv.gz"
            with gzip.open(source, "wt", encoding="utf-8", newline="") as handle:
                handle.write(
                    "SNP\tCHR\tBP\tA1\tA2\tBETA\tSE\tP\tN\n"
                    "rs1\t1\t101\tA\tC\t0.2\t0.1\t0.05\t1000\n"
                    "rs2\t2\t202\tG\tT\t-0.3\t0.2\t0.01\t900\n"
                )
            rows = list(module.variant_rows(source))
            self.assertEqual([row[0] for row in rows], ["rs1", "rs2"])
            self.assertAlmostEqual(rows[0][5], 2.0)
            self.assertAlmostEqual(rows[1][5], -1.5)
            duplicate = Path(directory) / "duplicate.tsv.gz"
            with gzip.open(duplicate, "wt", encoding="utf-8", newline="") as handle:
                handle.write(
                    "SNP\tCHR\tBP\tA1\tA2\tBETA\tSE\tP\tN\n"
                    "rs1\t1\t101\tA\tC\t0.2\t0.1\t0.05\t1000\n"
                    "rs1\t1\t101\tA\tC\t0.2\t0.1\t0.05\t1000\n"
                )
            connection = sqlite3.connect(":memory:")
            with self.assertRaises(SystemExit):
                module.load_table(connection, "variants", duplicate)
            connection.close()

    def test_pleiotropy_materializer_has_no_implicit_large_run(self):
        materializer = (ROOT / "scripts" / "44_materialize_pleiotropy_pair.py").read_text(
            encoding="utf-8"
        )
        self.assertIn('parser.add_argument("--materialize", action="store_true")', materializer)
        self.assertIn('if not args.materialize:', materializer)
        self.assertIn('lock.get("policy_sha256") != sha256(policy_path)', materializer)
        self.assertIn('policy["placo_z_squared_maximum"]', materializer)

    def test_placo_task_and_runtime_are_checksum_locked(self):
        task = (ROOT / "scripts" / "45_prepare_placo_task.py").read_text(
            encoding="utf-8"
        )
        runtime = (ROOT / "scripts" / "46_run_placo_pair.R").read_text(
            encoding="utf-8"
        )
        self.assertIn('manifest_lock.get("policy_sha256") != sha256(policy_path)', task)
        self.assertIn('provenance.get("output_sha256") == pair_hash', task)
        self.assertIn('sha256(placo_source) != policy["placo_source_sha256"]', task)
        self.assertIn('args[[3L]] != "--execute"', runtime)
        self.assertIn('lock$task_sha256 != sha256(task_path)', runtime)
        self.assertIn('placo.plus(', runtime)
        self.assertIn('failure_fraction > maximum_failure_fraction', runtime)
        self.assertIn('family_threshold >= conventional_threshold', runtime)
        self.assertNotIn("global_rg", runtime)

    def test_pleiofdr_trait_materializer_is_reference_ordered_and_opt_in(self):
        script = (ROOT / "scripts" / "47_prepare_pleiofdr_trait.py").read_text(
            encoding="utf-8"
        )
        self.assertIn('parser.add_argument("--materialize", action="store_true")', script)
        self.assertIn('if not args.materialize:', script)
        self.assertIn('policy["pleiofdr_variant_template_sha256"]', script)
        self.assertIn('reader.fieldnames != TEMPLATE_FIELDS', script)
        self.assertIn('np.full(variant_count, np.nan', script)
        self.assertIn('do_compression=True', script)
        self.assertIn('normalise_mat_header(temporary)', script)

    def test_conjfdr_task_and_runtime_require_pinned_overlap_correction(self):
        task = (ROOT / "scripts" / "48_prepare_conjfdr_task.py").read_text(
            encoding="utf-8"
        )
        runtime = (ROOT / "scripts" / "49_run_conjfdr_pair.py").read_text(
            encoding="utf-8"
        )
        self.assertIn('dirty:', task)
        self.assertIn('policy["pleiofdr_overlap_patch_sha256"]', task)
        self.assertIn('"randprune_n={policy[\'pleiofdr_random_prune_iterations\']}"', task)
        self.assertIn('"exclude_from_discovery=false"', task)
        self.assertIn('"manh_plot=false"', task)
        self.assertIn('if not args.execute:', runtime)
        self.assertIn('task.get("correct_sample_overlap") != "TRUE"', runtime)
        self.assertIn('shutil.which("matlab")', runtime)
        self.assertIn('shutil.copytree(code, runtime', runtime)
        self.assertIn('["patch", "-p1", "-i", str(patch)]', runtime)

    def test_pleiotropy_collator_requires_all_pairs_and_same_block_consensus(self):
        collator = (ROOT / "scripts" / "50_collate_pleiotropy.py").read_text(
            encoding="utf-8"
        )
        acceptance = (ROOT / "scripts" / "99_atlas_acceptance.py").read_text(
            encoding="utf-8"
        )
        self.assertIn('len(manifest) != policy["expected_sleep_non_sleep_pairs"]', collator)
        self.assertIn('set(placo_index) & set(conj_index)', collator)
        self.assertIn('PLACO_PLUS_AND_CONJFDR_SAME_LOCKED_LD_BLOCK', collator)
        self.assertIn('"effect_direction"', collator)
        self.assertIn('hit["Z1"]', collator)
        self.assertIn('hit["Z2"]', collator)
        self.assertIn('"--validate-only", "--quiet"', acceptance)
        self.assertIn('396 locked PLACO+/conjunction-FDR pair scans', acceptance)
        workflow = (ROOT / "Snakefile").read_text(encoding="utf-8")
        self.assertIn("rule placo_pair:", workflow)
        self.assertIn("rule conjfdr_pair:", workflow)
        self.assertIn("rule pleiotropy:", workflow)

    def test_downstream_policy_and_atlas_core_are_locked_without_placeholders(self):
        policy = json.loads(
            (ROOT / "config/downstream_analysis_policy.json").read_text(encoding="utf-8")
        )
        self.assertEqual(policy["expected_traits"], 45)
        self.assertEqual(policy["expected_sleep_non_sleep_pairs"], 396)
        self.assertEqual(policy["fine_mapping"]["susieR_version"], "0.14.2")
        self.assertFalse(policy["fine_mapping"]["estimate_residual_variance"])
        self.assertEqual(policy["colocalization"]["coloc_version"], "5.2.3")
        self.assertEqual(len(policy["robustness"]["required_families"]), 9)
        self.assertEqual(len(policy["integrated_atlas"]["tables"]), 10)
        builder = (ROOT / "scripts/51_build_atlas_core.py").read_text(encoding="utf-8")
        self.assertIn('len(rg_rows) != 396 or observed_pairs != expected_pairs', builder)
        self.assertIn('"CORE_ATLAS_COMPLETE_DOWNSTREAM_LAYERS_PENDING"', builder)
        self.assertNotIn("PLACEHOLDER", builder)

    def test_integrated_atlas_schema_is_exact_and_currently_fails_closed(self):
        schema = json.loads(
            (ROOT / "config/atlas_table_schema.json").read_text(encoding="utf-8")
        )
        policy = json.loads(
            (ROOT / "config/downstream_analysis_policy.json").read_text(encoding="utf-8")
        )
        self.assertEqual(list(schema["tables"]), policy["integrated_atlas"]["tables"])
        self.assertEqual(
            schema["tables"]["edges.tsv"]["fields"],
            policy["integrated_atlas"]["edge_fields"],
        )
        self.assertEqual(schema["tables"]["traits.tsv"]["expected_rows"], 45)
        self.assertEqual(schema["tables"]["trait_pairs.tsv"]["expected_rows"], 396)
        for name in policy["integrated_atlas"]["tables"][2:]:
            expected_minimum = 1 if name in {"loci.tsv", "variants.tsv"} else 0
            self.assertEqual(schema["tables"][name]["minimum_rows"], expected_minimum)
            self.assertIn("provenance_id", schema["tables"][name]["fields"])
        result = subprocess.run(
            [sys.executable, str(ROOT / "scripts/52_validate_integrated_atlas.py")],
            cwd=ROOT, capture_output=True, text=True,
        )
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("loci.tsv", result.stdout + result.stderr)

    def test_robustness_and_release_are_content_validated_not_presence_gated(self):
        policy = json.loads(
            (ROOT / "config/downstream_analysis_policy.json").read_text(encoding="utf-8")
        )
        self.assertEqual(len(policy["robustness"]["required_families"]), 9)
        self.assertEqual(
            policy["robustness"]["table_fields"][0:3],
            ["conclusion_id", "conclusion", "robustness_family"],
        )
        robustness = subprocess.run(
            [sys.executable, str(ROOT / "scripts/53_validate_robustness.py")],
            cwd=ROOT, capture_output=True, text=True,
        )
        self.assertNotEqual(robustness.returncode, 0)
        self.assertIn("robustness_summary.tsv", robustness.stdout + robustness.stderr)
        acceptance = (ROOT / "scripts/99_atlas_acceptance.py").read_text(encoding="utf-8")
        self.assertIn("scripts/52_validate_integrated_atlas.py", acceptance)
        self.assertIn("scripts/53_validate_robustness.py", acceptance)
        self.assertIn("scripts/55_validate_release.py", acceptance)
        release_builder = (ROOT / "scripts/54_build_release.py").read_text(encoding="utf-8")
        self.assertIn("non-release gates remain blocked", release_builder)
        self.assertIn("releases are never overwritten", release_builder)
        self.assertIn('parser.add_argument("--execute", action="store_true"', release_builder)
        for directory in ("results/molecular", "results/interpretation", "results/robustness"):
            self.assertIn(f'"{directory}"', release_builder)
        release_validator = (ROOT / "scripts/55_validate_release.py").read_text(
            encoding="utf-8"
        )
        self.assertIn("release contains unmanifested files", release_validator)
        self.assertIn("frozen pre-release acceptance evidence is incomplete", release_validator)

    def test_finemapping_policy_is_all_primary_loci_and_signed_lava_ld(self):
        policy = json.loads(
            (ROOT / "config/fine_mapping_analysis_policy.json").read_text(encoding="utf-8")
        )
        downstream = json.loads(
            (ROOT / "config/downstream_analysis_policy.json").read_text(encoding="utf-8")
        )
        self.assertIn("Every cross-method shared locus", policy["entry_rule"])
        self.assertIn("no result-ranked maximum", policy["entry_rule"])
        self.assertEqual(policy["reference"]["mode"], "LAVA_BCOR_SIGNED_CORRELATION")
        self.assertEqual(policy["reference"]["build"], "GRCh37")
        self.assertEqual(policy["reference"]["ancestry"], "EUR")
        self.assertEqual(policy["fine_mapping"]["version"], downstream["fine_mapping"]["susieR_version"])
        self.assertEqual(policy["colocalization"]["version"], downstream["colocalization"]["coloc_version"])
        self.assertFalse(policy["fine_mapping"]["estimate_residual_variance"])
        self.assertEqual(policy["colocalization"]["p12_sensitivity"], [1e-6, 5e-6, 1e-5, 5e-5])
        engine = ROOT / policy["shared_engine"]["path"]
        self.assertEqual(hashlib.sha256(engine.read_bytes()).hexdigest(), policy["shared_engine"]["sha256"])

    def test_finemapping_allele_alignment_is_signed_and_conservative(self):
        module = load_finemapping_materializer()
        self.assertEqual(module.aligned_sign("A", "C", "A", "C"), 1)
        self.assertEqual(module.aligned_sign("C", "A", "A", "C"), -1)
        self.assertEqual(module.aligned_sign("T", "G", "A", "C"), 1)
        self.assertEqual(module.aligned_sign("G", "T", "A", "C"), -1)
        self.assertIsNone(module.aligned_sign("A", "T", "A", "T"))
        self.assertIsNone(module.aligned_sign("A", "G", "A", "C"))

    def test_finemapping_workflow_fails_closed_and_does_not_claim_molecular_coloc(self):
        preflight = subprocess.run(
            [sys.executable, str(ROOT / "scripts/56_finemapping_preflight.py"), "--report-only"],
            cwd=ROOT, capture_output=True, text=True,
        )
        self.assertEqual(preflight.returncode, 0, preflight.stdout + preflight.stderr)
        report = json.loads((ROOT / "results/tables/fine_mapping_preflight.json").read_text(encoding="utf-8"))
        self.assertTrue(all(check["archive_present_and_pinned"] for check in report["code_checks"]))
        self.assertTrue(report["runtime"]["pass"])
        self.assertEqual(report["full_input_traits"], 36)
        self.assertEqual(report["lava_reference_files"], 0)
        self.assertFalse(report["ready"])
        prepare = (ROOT / "scripts/57_prepare_finemapping_loci.py").read_text(encoding="utf-8")
        self.assertIn('row["analysis_tier"] == "PRIMARY_PHASE1"', prepare)
        self.assertNotIn("global_rg", prepare)
        materialize = (ROOT / "scripts/58_materialize_finemapping_locus.py").read_text(encoding="utf-8")
        self.assertIn("full non-HapMap3 harmonized input is absent", prepare)
        self.assertIn('"results_accessed_before_lock": False', materialize)
        self.assertIn('"single_signal_fallback_justification": "NOT_JUSTIFIED"', materialize)
        collator = (ROOT / "scripts/60_collate_finemapping.py").read_text(encoding="utf-8")
        self.assertIn("results/tables/trait_trait_colocalization.tsv", collator)
        acceptance = (ROOT / "scripts/99_atlas_acceptance.py").read_text(encoding="utf-8")
        self.assertIn("trait-trait and molecular-QTL signal-level colocalization", acceptance)
        downloader = (ROOT / "scripts/32_download_lava_reference.sh").read_text(encoding="utf-8")
        self.assertIn("extracted_manifest.tsv", downloader)
        workflow = (ROOT / "Snakefile").read_text(encoding="utf-8")
        self.assertIn("checkpoint fine_mapping_loci:", workflow)
        for rule in ("fine_mapping_preflight", "fine_mapping_input", "fine_mapping_locus", "fine_mapping"):
            self.assertIn(f"rule {rule}:", workflow)

    def test_molecular_source_alignment_is_oriented_to_locked_ld_alleles(self):
        module = load_numbered_script("63_query_eqtl_catalogue.py", "molecular_query")
        self.assertEqual(module.orient_to_locked_alleles("A", "G", 0.2, "A", "G"), ("A", "G", 0.2))
        self.assertEqual(module.orient_to_locked_alleles("G", "A", 0.2, "A", "G"), ("A", "G", -0.2))
        with self.assertRaises(ValueError):
            module.orient_to_locked_alleles("A", "C", 0.2, "A", "G")

    def test_twas_eligibility_rejects_invalid_h2_without_capping(self):
        module = load_numbered_script("69_prepare_twas_manifest.py", "twas_manifest")
        trait = {"trait_id": "trait", "domain": "disease", "type": "binary"}
        eligible = module.trait_eligibility(trait, {"n_total": "1000", "h2": "0.2", "h2_scale": "liability"})
        self.assertEqual(eligible["analysis_status"], "ELIGIBLE")
        invalid = module.trait_eligibility(trait, {"n_total": "1000", "h2": "1.2", "h2_scale": "liability"})
        self.assertEqual(invalid["analysis_status"], "NOT_APPLICABLE")
        self.assertEqual(invalid["reason"], "H2_OUTSIDE_OPEN_CLOSED_0_1")
        self.assertEqual(invalid["gwas_h2"], "1.2")

    def test_phi_model_inventory_snapshot_has_exact_49_context_family(self):
        module = load_numbered_script("67_lock_twas_model_inventory.py", "phi_inventory")
        files = []
        for page in sorted((ROOT / "results/sources/twas_phi_model_inventory").glob("page_*.html")):
            parser = module.BoxListingParser(); parser.feed(page.read_text(encoding="utf-8"))
            files.extend(parser.files)
        self.assertEqual(len(files), 98)
        contexts = {}
        for row in files:
            name = row["filename"]
            context = name[3:-3] if name.endswith(".db") else name[3:-7]
            role = "db" if name.endswith(".db") else "cov"
            contexts.setdefault(context, set()).add(role)
        self.assertEqual(len(contexts), 49)
        self.assertTrue(all(roles == {"db", "cov"} for roles in contexts.values()))
        self.assertEqual(sum(int(row["bytes"]) for row in files), 3135665776)

    def test_phi_model_index_locks_missing_phi_as_pre_result_exclusion(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            (root / "config").mkdir(); (root / "results/tables").mkdir(parents=True)
            model_dir = root / "ref/models"; model_dir.mkdir(parents=True)
            model = model_dir / "en_Whole_Blood.db"
            connection = sqlite3.connect(model)
            connection.execute("CREATE TABLE weights (gene TEXT, rsid TEXT, varID TEXT, ref_allele TEXT, eff_allele TEXT, weight REAL)")
            connection.execute("CREATE TABLE extra (gene TEXT, phi REAL)")
            connection.executemany(
                "INSERT INTO weights VALUES (?,?,?,?,?,?)",
                [
                    ("ENSG1.1", "rs1", "chr1_100_A_G_b38", "A", "G", 0.1),
                    ("ENSG2.1", "rs2", "chr1_200_C_T_b38", "C", "T", 0.2),
                ],
            )
            connection.execute("INSERT INTO extra VALUES (?,?)", ("ENSG1.1", 0.0001))
            connection.commit(); connection.close()
            covariance = model_dir / "en_Whole_Blood.txt.gz"
            with gzip.open(covariance, "wt", encoding="utf-8") as handle:
                handle.write("GENE RSID1 RSID2 VALUE\n")
            policy = {
                "twas": {
                    "model_families": ["GTEx_v8_ELASTIC_NET_PHI_eQTL"],
                    "phi_model_source": {
                        "inventory_path": "results/tables/inventory.tsv",
                        "inventory_lock_path": "results/tables/inventory.lock.json",
                        "download_lock_path": "results/tables/download.lock.json",
                        "install_dir": "ref/models", "expected_context_count": 1,
                        "release_name": "fixture",
                    },
                },
                "claim_limit": "fixture claim limit",
            }
            policy_path = root / "config/policy.json"
            policy_path.write_text(json.dumps(policy) + "\n", encoding="utf-8")
            fields = ["file_id", "filename", "context", "file_role", "bytes", "download_url", "results_accessed_before_lock"]
            inventory = root / "results/tables/inventory.tsv"
            rows = [
                ["f_1", covariance.name, "Whole_Blood", "COVARIANCE", covariance.stat().st_size, "https://example.test/cov", "NO"],
                ["f_2", model.name, "Whole_Blood", "MODEL_DB", model.stat().st_size, "https://example.test/db", "NO"],
            ]
            with inventory.open("w", encoding="utf-8", newline="") as handle:
                writer = csv.writer(handle, delimiter="\t", lineterminator="\n"); writer.writerow(fields); writer.writerows(rows)
            inventory_lock = root / "results/tables/inventory.lock.json"
            inventory_lock.write_text(json.dumps({
                "inventory_sha256": hashlib.sha256(inventory.read_bytes()).hexdigest(),
                "policy_sha256": hashlib.sha256(policy_path.read_bytes()).hexdigest(),
                "file_ids_in_locked_order": ["f_1", "f_2"],
            }) + "\n", encoding="utf-8")
            download_lock = root / "results/tables/download.lock.json"
            download_lock.write_text(json.dumps({
                "inventory_sha256": hashlib.sha256(inventory.read_bytes()).hexdigest(),
                "inventory_lock_sha256": hashlib.sha256(inventory_lock.read_bytes()).hexdigest(),
                "policy_sha256": hashlib.sha256(policy_path.read_bytes()).hexdigest(),
                "files": {
                    covariance.name: {"file_id": "f_1", "bytes": covariance.stat().st_size, "sha256": hashlib.sha256(covariance.read_bytes()).hexdigest()},
                    model.name: {"file_id": "f_2", "bytes": model.stat().st_size, "sha256": hashlib.sha256(model.read_bytes()).hexdigest()},
                },
            }) + "\n", encoding="utf-8")
            result = subprocess.run([
                sys.executable, str(ROOT / "scripts/68_index_twas_models.py"), "--root", str(root),
                "--policy", "config/policy.json", "--materialize",
            ], cwd=ROOT, capture_output=True, text=True)
            self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
            with (root / "results/tables/twas_model_registry.tsv").open(encoding="utf-8", newline="") as handle:
                registry = list(csv.DictReader(handle, delimiter="\t"))
            self.assertEqual(registry[0]["phi_gene_count"], "1")
            self.assertEqual(registry[0]["phi_missing_gene_count"], "1")
            with (root / "results/tables/twas_model_phi_exclusions.tsv").open(encoding="utf-8", newline="") as handle:
                exclusions = list(csv.DictReader(handle, delimiter="\t"))
            self.assertEqual(exclusions[0]["gene_id"], "ENSG2.1")
            self.assertEqual(exclusions[0]["exclusion_reason"], "MODEL_FEATURE_MISSING_PHI")

    def test_molecular_coloc_prior_robustness_never_pools_signal_pairs(self):
        module = load_numbered_script("73_collate_molecular.py", "molecular_collator")
        grid = [1e-6, 5e-6, 1e-5, 5e-5]
        rows = [{
            "comparison_id": "MQC", "p12": str(value), "signal1": "s1", "signal2": "q1",
            "PP_H4": "0.9", "PP_H4_over_PP_H3": "6", "coloc_method": "COLOC_SUSIE",
            "fine_mapping_qc": "PASS",
        } for value in grid]
        self.assertEqual(module.robust_signal_pairs(rows, grid, 0.8, 5), {("s1", "q1")})
        rows[-1]["signal2"] = "q2"
        self.assertEqual(module.robust_signal_pairs(rows, grid, 0.8, 5), set())

    def test_twas_gene_to_locus_mapping_uses_exact_grch38_overlap(self):
        module = load_numbered_script("73_collate_molecular.py", "molecular_gene_mapping")

        class IdentityChain:
            def map_point(self, chromosome, position):
                return "mapped", (int(chromosome), int(position) + 100, "+")

        genes = {
            "ENSG1": {"gene_id": "ENSG1", "gene_symbol": "G1", "chromosome": 1, "start": 1090, "end": 1200},
            "ENSG2": {"gene_id": "ENSG2", "gene_symbol": "G2", "chromosome": 1, "start": 1301, "end": 1400},
            "ENSG3": {"gene_id": "ENSG3", "gene_symbol": "G3", "chromosome": 2, "start": 1090, "end": 1200},
        }
        loci = [{"locus_id": "L1", "chromosome": "1", "start_bp": "1000", "end_bp": "1200"}]
        overlap, mapped = module.overlap_gene_loci(genes, loci, IdentityChain())
        self.assertEqual(overlap, {"L1": ["ENSG1"]})
        self.assertEqual(mapped, {"L1": (1, 1100, 1300)})

    def test_molecular_workflow_pins_annotation_and_guards_large_download(self):
        policy = json.loads((ROOT / "config/molecular_analysis_policy.json").read_text(encoding="utf-8"))
        annotation = next(asset for asset in policy["metadata_assets"] if asset["id"] == "GTEX_V8_GENCODE_V26_COLLAPSED_GENES")
        self.assertEqual(annotation["bytes"], 134502408)
        self.assertEqual(annotation["sha256"], "2a7ae0883b35e7ff66a8f71a1779198acc6609d648ff55d7b670d010a9f243f6")
        self.assertEqual(policy["twas"]["gene_annotation_asset_id"], annotation["id"])
        fetcher = (ROOT / "scripts/67_fetch_twas_resources.sh").read_text(encoding="utf-8")
        self.assertIn("--acknowledge-large-download", fetcher)
        self.assertIn("warn the user first", fetcher)
        collator = (ROOT / "scripts/73_collate_molecular.py").read_text(encoding="utf-8")
        self.assertIn("rather than a fabricated gene", policy["twas"]["locus_integration_rule"])
        self.assertIn('"molecular_locus_coverage.tsv"', (ROOT / "scripts/52_validate_integrated_atlas.py").read_text(encoding="utf-8"))
        self.assertIn("supported_gene_keys", collator)
        tool_versions = (ROOT / "environment/tool_versions.tsv").read_text(encoding="utf-8")
        self.assertIn("MetaXcan\t0.8.1\tPINNED", tool_versions)
        release_validator = (ROOT / "scripts/55_validate_release.py").read_text(encoding="utf-8")
        self.assertIn('"susieR", "coloc", "MetaXcan"', release_validator)
        workflow = (ROOT / "Snakefile").read_text(encoding="utf-8")
        for rule in (
            "molecular_metadata", "molecular_preflight", "molecular_source_search",
            "molecular_coloc_run", "twas_model_index", "twas_run", "molecular_integration",
        ):
            self.assertIn(f"rule {rule}:", workflow)

    def test_interpretation_task_family_is_exact_cross_product(self):
        module = load_numbered_script("75_prepare_interpretation_tasks.py", "interpretation_tasks")
        policy = json.loads((ROOT / "config/interpretation_analysis_policy.json").read_text(encoding="utf-8"))
        with (ROOT / "config/interpretation_source_registry.tsv").open(encoding="utf-8", newline="") as handle:
            sources = list(csv.DictReader(handle, delimiter="\t"))
        panel = [{"trait_id": "sleep", "domain": "sleep"}, {"trait_id": "disease", "domain": "cardio"}]
        pairs = [{"pair_id": "sleep__disease", "sleep_trait": "sleep", "non_sleep_trait": "disease"}]
        loci = [{"locus_id": "L1", "sleep_trait": "sleep", "non_sleep_trait": "disease"}]
        variants = [{"locus_id": "L1", "variant_id": "rs1"}]
        genes = [{"locus_id": "L1", "gene_id": "ENSG1", "evidence_level": "HIGH_CONFIDENCE_CONVERGENT_GENE"}]
        tasks = module.build_tasks(policy, panel, pairs, loci, variants, genes, sources)
        counts = {family: sum(row["analysis_family"] == family for row in tasks) for family in (
            "regulatory", "cell_type", "pathway", "causal",
        )}
        self.assertEqual(counts, {"regulatory": 24, "cell_type": 40, "pathway": 4, "causal": 10})
        self.assertEqual(len({row["task_id"] for row in tasks}), len(tasks))

    def test_interpretation_bh_and_null_evidence_contract(self):
        module = load_numbered_script("77_collate_interpretation.py", "interpretation_collator")
        self.assertEqual(module.bh([0.01, 0.04, 0.03]), [0.03, 0.04, 0.04])
        schema = json.loads((ROOT / "config/atlas_table_schema.json").read_text(encoding="utf-8"))
        for name in ("genes.tsv", "regulatory_elements.tsv", "cell_types.tsv", "pathways.tsv", "causal_tests.tsv"):
            self.assertEqual(schema["tables"][name]["minimum_rows"], 0)
        validator = (ROOT / "scripts/52_validate_integrated_atlas.py").read_text(encoding="utf-8")
        self.assertIn("interpretation coverage is not the exact locked task family", validator)
        self.assertIn("access-blocked interpretation task remains unresolved", validator)
        release_validator = (ROOT / "scripts/55_validate_release.py").read_text(encoding="utf-8")
        for component in ("MAGMA", "TwoSampleMR", "MR-PRESSO", "CAUSE_or_LHC_MR"):
            self.assertIn(f'"{component}"', release_validator)

    def test_regulatory_adapters_preserve_unlinked_screen_and_supported_promoters(self):
        module = load_numbered_script("82_run_regulatory_task.py", "regulatory_adapter")
        policy = json.loads((ROOT / "config/interpretation_analysis_policy.json").read_text(encoding="utf-8"))
        variant = {"variant_id": "rs1"}
        screen_task = {
            "task_id": "REG__L1__ENCODE_CCRE_V4_ENHANCER__brain",
            "source_id": "ENCODE_CCRE_V4_ENHANCER", "source_release": "SCREEN_TEST",
            "locus_id": "L1", "domain": "brain",
        }
        promoter_task = {
            "task_id": "REG__L1__GENCODE_V26_PROMOTERS__brain",
            "source_id": "GENCODE_V26_PROMOTERS", "source_release": "GENCODE_TEST",
            "locus_id": "L1", "domain": "brain",
        }
        with tempfile.TemporaryDirectory() as temporary:
            directory = Path(temporary)
            coordinates = directory / "ccres.bed"
            coordinates.write_text("chr1\t100\t110\tEH38D1\tEH38E1\tpELS\n", encoding="ascii")
            matrix = directory / "matrix.tsv.gz"
            with gzip.open(matrix, "wt", encoding="ascii", newline="") as handle:
                handle.write("cCRE\tENCSR1\nEH38E1\tpELS\n")
            metadata = directory / "metadata.js"
            metadata.write_text(
                '{"name":"astrocyte","collection":"Core","ontology":"brain",'
                '"lifeStage":"adult","sampleType":"primary cell","displayName":"astrocyte",'
                '"assays":[{"id":"dnase-ENCFF1","assay":"dnase","url":"x",'
                '"experimentAccession":"ENCSR1","fileAccession":"ENCFF1"}]}',
                encoding="utf-8",
            )
            mapped = [{"variant": variant, "chromosome_grch38": 1, "position_grch38": 101, "target_strand": "+"}]
            screen = module.screen_rows(screen_task, policy, mapped, coordinates, matrix, metadata)
            self.assertEqual(len(screen), 1)
            self.assertEqual(screen[0]["target_gene_id"], "NA")
            self.assertEqual(screen[0]["cell_type"], "astrocyte")
            self.assertEqual(list(screen[0]), policy["regulatory_mapping"]["canonical_fields"])

            gtf = directory / "genes.gtf"
            gtf.write_text(
                'chr1\tGENCODE\tgene\t1000\t2000\t.\t+\t.\tgene_id "ENSG1.5"; gene_name "G1";\n',
                encoding="utf-8",
            )
            promoter_mapped = [{"variant": variant, "chromosome_grch38": 1, "position_grch38": 900, "target_strand": "+"}]
            promoters = module.promoter_rows(
                promoter_task, policy, promoter_mapped, [{"locus_id": "L1", "gene_id": "ENSG1"}], gtf,
            )
            self.assertEqual(len(promoters), 1)
            self.assertEqual(promoters[0]["target_gene_id"], "ENSG1")
            self.assertEqual(promoters[0]["element_type"], "promoter")

        schema = json.loads((ROOT / "config/atlas_table_schema.json").read_text(encoding="utf-8"))
        self.assertEqual(schema["tables"]["regulatory_elements.tsv"]["primary_key"], ["regulatory_evidence_id"])
        edge_builder = (ROOT / "scripts/78_build_atlas_edges.py").read_text(encoding="utf-8")
        self.assertIn('if regulatory["target_gene_id"] not in {"", "NA"}', edge_builder)

    def test_hocomoco_allele_scanner_uses_pinned_score_thresholds(self):
        module = load_numbered_script("83_run_motif_task.py", "motif_adapter")
        motif = {
            "pwm": [[1.0, 0.0, 0.0, 0.0]],
            "threshold_scores": [0.0, 1.0],
            "threshold_p": [0.5, 0.01],
        }
        hit = module.best_allele_hit("A", "C", 0, motif, 0.05)
        self.assertIsNotNone(hit)
        self.assertEqual(hit["ref_p"], 0.01)
        self.assertEqual(hit["alt_p"], 0.5)
        self.assertEqual(hit["delta"], -1.0)
        self.assertIsNone(module.best_allele_hit("A", "C", 0, motif, 0.001))

        policy = json.loads((ROOT / "config/interpretation_analysis_policy.json").read_text(encoding="utf-8"))
        bundle_path = ROOT / policy["hocomoco_v14"]["component_manifest"]
        self.assertEqual(hashlib.sha256(bundle_path.read_bytes()).hexdigest(), policy["hocomoco_v14"]["component_manifest_sha256"])
        bundle = json.loads(bundle_path.read_text(encoding="utf-8"))
        self.assertEqual({row["component_id"] for row in bundle["components"]}, {
            "H14CORE_MEME", "H14CORE_ANNOTATION", "H14CORE_PWM", "H14CORE_THRESHOLDS",
        })

    def test_abc_adapter_requires_nonself_source_link_and_supported_gene(self):
        module = load_numbered_script("84_prepare_abc_overlap_cache.py", "abc_adapter")
        policy = json.loads((ROOT / "config/interpretation_analysis_policy.json").read_text(encoding="utf-8"))
        manifest_path = ROOT / policy["abc_2021"]["component_manifest"]
        self.assertEqual(
            hashlib.sha256(manifest_path.read_bytes()).hexdigest(),
            policy["abc_2021"]["component_manifest_sha256"],
        )
        manifest = {
            "release": "ABC_TEST", "row_count": 2, "biosample_count": 1,
            "domain_biosamples": {"brain": ["astrocyte"], "immune": [], "metabolic": [], "vascular": []},
        }
        variants = [{
            "locus_id": "L1", "variant_id": "rs1", "chromosome": "1", "position_bp": "101",
            "qc_status": "PASS", "credible_set_sleep": "CS1", "credible_set_non_sleep": "NA",
            "shared_signal_posterior": "NA",
        }]
        genes = [{"locus_id": "L1", "gene_id": "ENSG1", "gene_symbol": "G1"}]
        with tempfile.TemporaryDirectory() as temporary:
            source = Path(temporary) / "abc.tsv.gz"
            base = {field: "NA" for field in module.ABC_FIELDS}
            base.update({
                "chr": "chr1", "start": "100", "end": "110", "name": "element1",
                "class": "genic", "activity_base": "1", "TargetGene": "G1",
                "TargetGeneTSS": "1000", "TargetGeneExpression": "1",
                "TargetGenePromoterActivityQuantile": "1", "TargetGeneIsExpressed": "True",
                "distance": "899", "isSelfPromoter": "False", "ABC.Score": "0.025",
                "CellType": "astrocyte",
            })
            self_promoter = dict(base, name="promoter", isSelfPromoter="True")
            with gzip.open(source, "wt", encoding="utf-8", newline="") as handle:
                writer = csv.DictWriter(handle, fieldnames=module.ABC_FIELDS, delimiter="\t", lineterminator="\n")
                writer.writeheader()
                writer.writerows([base, self_promoter])
            rows, counters = module.prepare_rows(policy, manifest, source, variants, genes)
        self.assertEqual(len(rows), 1)
        self.assertEqual(rows[0]["target_gene_id"], "ENSG1")
        self.assertEqual(rows[0]["context_domain"], "brain")
        self.assertEqual(rows[0]["effect"], "0.025")
        self.assertEqual(rows[0]["evidence_level"], "ABC_PRIMARY_GE_0.02")
        self.assertEqual(counters["self_promoter_rows_excluded"], 1)
        self.assertEqual(list(rows[0]), policy["regulatory_mapping"]["canonical_fields"])

    def test_pchic_adapter_maps_both_contact_directions_at_chicago_five(self):
        module = load_numbered_script("86_prepare_pchic_overlap_cache.py", "pchic_adapter")
        policy = json.loads((ROOT / "config/interpretation_analysis_policy.json").read_text(encoding="utf-8"))
        manifest_path = ROOT / policy["pchic_2016"]["component_manifest"]
        self.assertEqual(
            hashlib.sha256(manifest_path.read_bytes()).hexdigest(),
            policy["pchic_2016"]["component_manifest_sha256"],
        )
        cells = {"Mon": "monocytes", "nB": "naive_B_cells"}
        manifest = {"release": "PCHIC_TEST", "row_count": 1, "cell_types": cells}
        variants = [
            {
                "locus_id": "L1", "variant_id": "rs_bait", "chromosome": "1", "position_bp": "101",
                "qc_status": "PASS", "credible_set_sleep": "CS1", "credible_set_non_sleep": "NA",
                "shared_signal_posterior": "NA",
            },
            {
                "locus_id": "L1", "variant_id": "rs_oe", "chromosome": "1", "position_bp": "201",
                "qc_status": "PASS", "credible_set_sleep": "CS1", "credible_set_non_sleep": "NA",
                "shared_signal_posterior": "NA",
            },
        ]
        genes = [
            {"locus_id": "L1", "gene_id": "ENSG1", "gene_symbol": "G1"},
            {"locus_id": "L1", "gene_id": "ENSG2", "gene_symbol": "G2"},
        ]
        fields = module.BASE_FIELDS + list(cells) + ["clusterID", "clusterPostProb"]
        source_row = {field: "NA" for field in fields}
        source_row.update({
            "baitChr": "1", "baitStart": "100", "baitEnd": "110", "baitID": "1",
            "baitName": "G1", "oeChr": "1", "oeStart": "200", "oeEnd": "210",
            "oeID": "2", "oeName": "G2", "dist": "100", "Mon": "6", "nB": "4.9",
            "clusterID": "1", "clusterPostProb": "0.9",
        })
        with tempfile.TemporaryDirectory() as temporary:
            matrix = Path(temporary) / "pchic.tsv.gz"
            with gzip.open(matrix, "wt", encoding="utf-8", newline="") as handle:
                writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
                writer.writeheader()
                writer.writerow(source_row)
            rows, counters = module.prepare_rows(policy, manifest, matrix, variants, genes)
        self.assertEqual(len(rows), 2)
        observed = {(row["variant_id"], row["target_gene_id"]) for row in rows}
        self.assertEqual(observed, {("rs_bait", "ENSG2"), ("rs_oe", "ENSG1")})
        self.assertTrue(all(row["biosample"] == "Mon" for row in rows))
        self.assertTrue(all(row["effect"] == "6" for row in rows))
        self.assertEqual(counters["source_rows"], 1)
        self.assertTrue(all(list(row) == policy["regulatory_mapping"]["canonical_fields"] for row in rows))

    def test_fuma_scrna_source_and_magma_method_are_pre_result_locked(self):
        policy = json.loads((ROOT / "config/interpretation_analysis_policy.json").read_text(encoding="utf-8"))
        spec = policy["fuma_scrna"]
        manifest_path = ROOT / spec["component_manifest"]
        self.assertEqual(hashlib.sha256(manifest_path.read_bytes()).hexdigest(), spec["component_manifest_sha256"])
        manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
        self.assertEqual(manifest["commit"], "dd526163ea80af1a80a6cdc80db167144500694b")
        self.assertEqual(manifest["processed_matrix_count"], 179)
        self.assertEqual(manifest["gene_count"], 20260)
        self.assertEqual(manifest["magma_version"], "1.10")
        self.assertEqual(manifest["gene_window_kb"], [1, 1])
        self.assertEqual(manifest["reference_sample_count"], 503)
        self.assertEqual(manifest["reference_archive_uncompressed_bytes"], 3600697827)
        domains = {domain: 0 for domain in policy["cell_types"]["required_domains"]}
        for dataset in manifest["datasets"]:
            domains[dataset["domain"]] += 1
            self.assertEqual(len(dataset["compressed_sha256"]), 64)
            self.assertEqual(len(dataset["text_sha256"]), 64)
            self.assertGreaterEqual(dataset["cell_type_count"], 4)
        self.assertEqual(domains, {domain: 2 for domain in domains})
        self.assertEqual(spec["gene_property_model"], "condition-hide=Average direction=greater")
        self.assertEqual(spec["gene_model"], "snp-wise=mean")
        self.assertEqual(spec["p_value_floor"], 1e-300)
        with (ROOT / policy["source_registry"]).open(encoding="utf-8", newline="") as handle:
            source = next(row for row in csv.DictReader(handle, delimiter="\t") if row["source_id"] == "FUMA_SCRNA")
        self.assertEqual(source["source_status"], "SOURCE_VERIFIED")
        self.assertEqual(source["exact_release"], manifest["release"])
        self.assertEqual(int(source["expected_bytes"]), 438582802)
        workflow = (ROOT / "Snakefile").read_text(encoding="utf-8")
        for rule in ("fuma_scrna_matrices", "magma_eur_reference", "magma_gene_annotation", "magma_gene_results"):
            self.assertIn(f"rule {rule}:", workflow)
        self.assertIn("91_run_fuma_scrna_task.py", (ROOT / "scripts/76_run_interpretation_task.py").read_text())

    def test_fuma_matrix_materializer_preserves_exact_text(self):
        module = load_numbered_script("88_materialize_fuma_resources.py", "fuma_resources")
        text = b"GENE\tA\tB\tAverage\nENSG00000000001\t1\t3\t2\n"
        compressed = gzip.compress(text, mtime=0)
        dataset = {
            "dataset_id": "TEST", "member": "processed_data/TEST.txt.gz",
            "domain": "brain", "species": "human", "tissue": "brain",
            "compressed_bytes": len(compressed),
            "compressed_sha256": hashlib.sha256(compressed).hexdigest(),
            "text_bytes": len(text), "text_sha256": hashlib.sha256(text).hexdigest(),
            "gene_rows": 1, "cell_type_count": 2,
        }
        with tempfile.TemporaryDirectory() as temporary:
            temporary_path = Path(temporary)
            archive = temporary_path / "source.tar.gz"
            member_name = "commit/processed_data/TEST.txt.gz"
            with tarfile.open(archive, "w:gz") as handle:
                info = tarfile.TarInfo(member_name)
                info.size = len(compressed)
                handle.addfile(info, io.BytesIO(compressed))
            records = module.materialize_matrices(archive, "commit/", [dataset], temporary_path / "out")
            output = temporary_path / "out/TEST.txt"
            self.assertEqual(output.read_bytes(), text)
            self.assertEqual(records[0]["sha256"], hashlib.sha256(text).hexdigest())

    def test_catlas_adult_source_and_enrichment_design_are_pre_result_locked(self):
        policy = json.loads((ROOT / "config/interpretation_analysis_policy.json").read_text(encoding="utf-8"))
        spec = policy["catlas_adult_v4"]
        manifest_path = ROOT / spec["component_manifest"]
        self.assertEqual(hashlib.sha256(manifest_path.read_bytes()).hexdigest(), spec["component_manifest_sha256"])
        manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
        self.assertEqual(manifest["source_release"], "Mendeley_Data_yv4fzv6cnm_v4_2022-01-24")
        self.assertEqual(manifest["data_doi"], "10.17632/yv4fzv6cnm.4")
        self.assertEqual(
            manifest["source_counts"],
            {
                "adult_nuclei": 615998,
                "adult_cell_types": 111,
                "adult_tissues": 30,
                "all_life_stage_ccres": 1154611,
                "adult_present_ccres": 890130,
                "cell_type_restricted_ccres": 435142,
                "cell_type_restricted_assignments": 1091805,
                "peak_archive_members": 111,
            },
        )
        selected = manifest["selected_cells"]
        self.assertEqual(len(selected), 43)
        observed_domains = {
            domain: sum(row["domain"] == domain for row in selected)
            for domain in policy["cell_types"]["required_domains"]
        }
        self.assertEqual(observed_domains, {"brain": 10, "immune": 9, "metabolic": 9, "vascular": 15})
        self.assertEqual(observed_domains, spec["selected_cells_by_domain"])
        for field in ("cell_type", "metadata_cell_type", "archive_member"):
            values = [row[field] for row in selected]
            self.assertEqual(len(values), len(set(values)))
        design = manifest["enrichment_design"]
        self.assertEqual(spec["instrument_p_threshold"], 5e-8)
        self.assertEqual(spec["minimum_maf"], 0.01)
        self.assertEqual(spec["minimum_trait_variant_coverage"], 0.5)
        self.assertEqual(spec["ld_pruning_r2"], 0.1)
        self.assertEqual(spec["ld_pruning_window_kb"], 1000)
        self.assertEqual(design["test"], "one-sided hypergeometric overrepresentation")
        self.assertEqual(design["effect"], "observed_over_expected_overlap_ratio")
        self.assertIn("BH FDR across all 43 selected cells", design["multiple_testing"])
        self.assertIn("not evidence", spec["claim_limit"])
        with (ROOT / policy["source_registry"]).open(encoding="utf-8", newline="") as handle:
            source = next(
                row for row in csv.DictReader(handle, delimiter="\t")
                if row["source_id"] == spec["source_id"]
            )
        peak_component = next(
            row for row in manifest["components"] if row["component_id"] == "CELL_TYPE_RESTRICTED_PEAKS"
        )
        self.assertEqual(source["source_status"], "SOURCE_VERIFIED")
        self.assertEqual(source["exact_release"], manifest["source_release"])
        self.assertEqual(int(source["expected_bytes"]), peak_component["bytes"])
        self.assertEqual(source["expected_sha256"], peak_component["sha256"])
        preflight = (ROOT / "scripts/74_interpretation_preflight.py").read_text(encoding="utf-8")
        self.assertIn("def validate_catlas_bundle", preflight)
        self.assertIn("catlas_ready, catlas_blocker, catlas_validation", preflight)
        dispatcher = (ROOT / "scripts/76_run_interpretation_task.py").read_text(encoding="utf-8")
        self.assertIn("97_run_catlas_task.py", dispatcher)
        workflow = (ROOT / "Snakefile").read_text(encoding="utf-8")
        self.assertIn("rule catlas_fixed_variant_universe:", workflow)
        self.assertIn("rule catlas_trait_cache:", workflow)

    def test_catlas_reference_annotation_uses_half_open_adult_intervals(self):
        module = load_numbered_script("95_prepare_catlas_reference.py", "catlas_reference")
        with tempfile.TemporaryDirectory() as temporary:
            temporary_path = Path(temporary)
            ccre = temporary_path / "ccre.tsv.gz"
            with gzip.open(ccre, "wt", encoding="utf-8", newline="") as handle:
                handle.write(
                    "#Chromosome\tStart\tEnd\tClass\tPresent in fetal tissues\t"
                    "Present in adult tissues\tCRE module\n"
                    "chr1\t99\t102\tDistal\tno\tyes\t1\n"
                    "chr1\t102\t104\tDistal\tyes\tno\t1\n"
                )
            peaks = temporary_path / "peaks.zip"
            member = "2E_Cell_type_restricted_peaks/A.bed.gz"
            with zipfile.ZipFile(peaks, "w") as archive:
                archive.writestr(member, gzip.compress(b"chr1\t99\t102\n", mtime=0))
            selected = [{"metadata_cell_type": "Cell A", "archive_member": member}]
            adult, membership, counts = module.annotate_catlas(
                ccre, peaks, selected, {"rs1": (1, 100), "rs2": (1, 103)}, "CATLAS_TEST",
            )
        self.assertEqual(adult, {"rs1"})
        self.assertEqual(membership["rs1"], {"CATLAS_TEST::Cell A"})
        self.assertEqual(membership["rs2"], set())
        self.assertEqual(counts, {"CATLAS_TEST::Cell A": 1})

    def test_catlas_trait_cache_and_enrichment_retain_complete_null_cells(self):
        trait_module = load_numbered_script("96_prepare_catlas_trait_cache.py", "catlas_trait")
        task_module = load_numbered_script("97_run_catlas_task.py", "catlas_task")
        with tempfile.TemporaryDirectory() as temporary:
            gwas = Path(temporary) / "trait.tsv.gz"
            with gzip.open(gwas, "wt", encoding="utf-8", newline="") as handle:
                handle.write(
                    "SNP\tCHR\tBP\tP\n"
                    "rs3\t1\t300\t0.5\n"
                    "rs1\t1\t100\t1e-9\n"
                    "outside\t1\t999\t0.2\n"
                    "rs2\t1\t200\t0.8\n"
                    "rs4\t1\t401\t0.2\n"
                )
            trait_rows, counts = trait_module.materialize_trait(
                gwas, ["rs1", "rs2", "rs3", "rs4"],
                {"rs1": (1, 100), "rs2": (1, 200), "rs3": (1, 300), "rs4": (1, 400)},
            )
        self.assertEqual([row["SNP"] for row in trait_rows], ["rs1", "rs2", "rs3"])
        self.assertEqual(counts["matched_variants"], 3)
        self.assertEqual(counts["coordinate_mismatch_variants"], 1)

        policy = json.loads((ROOT / "config/interpretation_analysis_policy.json").read_text(encoding="utf-8"))
        policy = json.loads(json.dumps(policy))
        policy["catlas_adult_v4"]["selected_cells_by_domain"]["brain"] = 2
        source_manifest = {"selected_cells": [
            {
                "domain": "brain", "cell_type": "Cell A", "tissue": "brain",
                "metadata_cell_type": "Cell A",
            },
            {
                "domain": "brain", "cell_type": "Cell B", "tissue": "brain",
                "metadata_cell_type": "Cell B",
            },
        ]}
        task = {
            "domain": "brain", "trait_id": "trait", "method": "scATAC_enrichment",
            "source_id": "SCATAC_MULTI_DOMAIN", "source_release": "TEST",
        }
        universe_rows = [
            {"SNP": "rs1", "IN_ADULT_CCRE": "1", "CELL_TYPE_IDS": "CATLAS_ADULT_V4::Cell A"},
            {"SNP": "rs2", "IN_ADULT_CCRE": "1", "CELL_TYPE_IDS": "CATLAS_ADULT_V4::Cell B"},
            {"SNP": "rs3", "IN_ADULT_CCRE": "1", "CELL_TYPE_IDS": "CATLAS_ADULT_V4::Cell A;CATLAS_ADULT_V4::Cell B"},
            {"SNP": "rs4", "IN_ADULT_CCRE": "0", "CELL_TYPE_IDS": "CATLAS_ADULT_V4::Cell A"},
        ]
        rows, result_counts = task_module.compute_enrichment_rows(
            policy, source_manifest, task, universe_rows,
            [*trait_rows, {"SNP": "rs4", "P": "1e-12"}],
        )
        self.assertEqual(len(rows), 2)
        self.assertAlmostEqual(float(rows[0]["effect"]), 1.5)
        self.assertAlmostEqual(float(rows[0]["p_value"]), 2 / 3)
        self.assertEqual(rows[1]["effect"], "0")
        self.assertEqual(rows[1]["p_value"], "1")
        self.assertTrue(all(row["fdr"] == "NA" and row["se"] == "NA" for row in rows))
        self.assertTrue(all(list(row) == policy["cell_types"]["canonical_fields"] for row in rows))
        self.assertEqual(result_counts["adult_ccre_background_variants"], 3)
        self.assertEqual(result_counts["genome_wide_significant_background_variants"], 1)

    def test_magma_pvalue_materializer_and_gsa_parser_are_fail_closed(self):
        gene_module = load_numbered_script("90_prepare_magma_gene_results.py", "magma_genes")
        task_module = load_numbered_script("91_run_fuma_scrna_task.py", "magma_cell_task")
        with tempfile.TemporaryDirectory() as temporary:
            temporary_path = Path(temporary)
            source = temporary_path / "trait.tsv.gz"
            with gzip.open(source, "wt", encoding="utf-8", newline="") as handle:
                handle.write("SNP\tP\tN\nrs1\t0\t1000\nrs2\t0.5\t900\n")
            pvalues = temporary_path / "trait.pval.tsv"
            details = gene_module.materialize_pvalues(
                source, pvalues,
                {"gwas_snp_field": "SNP", "gwas_p_field": "P", "gwas_n_field": "N", "p_value_floor": 1e-300},
            )
            self.assertEqual(details["rows"], 2)
            self.assertEqual(details["p_values_floored"], 1)
            self.assertIn("rs1\t1e-300\t1000", pvalues.read_text(encoding="utf-8"))
            gsa = temporary_path / "test.gsa.out"
            gsa.write_text(
                "# TOTAL_GENES = 12000\nVARIABLE TYPE NGENES BETA BETA_STD SE P\n"
                "A COVAR 12000 0.1 0.2 0.03 0.001\nB COVAR 12000 -0.2 -0.3 0.04 0.9\n",
                encoding="utf-8",
            )
            rows = task_module.parse_gsa(gsa, ["A", "B"])
            self.assertEqual([row["VARIABLE"] for row in rows], ["A", "B"])

    def test_four_pathway_sources_are_exactly_locked_before_results(self):
        policy = json.loads((ROOT / "config/interpretation_analysis_policy.json").read_text(encoding="utf-8"))
        spec = policy["public_pathway_sources"]
        manifest_path = ROOT / spec["component_manifest"]
        self.assertEqual(hashlib.sha256(manifest_path.read_bytes()).hexdigest(), spec["component_manifest_sha256"])
        bundle = json.loads(manifest_path.read_text(encoding="utf-8"))
        self.assertEqual(bundle["set_size_range"], [10, 1000])
        self.assertEqual(bundle["resources"]["REACTOME"]["release"], "Reactome_v97_2026-06-30")
        self.assertEqual(bundle["resources"]["REACTOME"]["eligible_sets"], 1744)
        self.assertEqual(bundle["resources"]["GO"]["release"], "GO_pipeline_2026-08-05_ontology_2026-07-26")
        self.assertEqual(bundle["resources"]["GO"]["eligible_sets"], 7719)
        self.assertEqual(bundle["resources"]["GO"]["propagation_relations"], ["is_a", "part_of"])
        self.assertEqual(bundle["resources"]["MSIGDB"]["release"], "MSigDB_v2026.1.Hs_2026-01-29")
        self.assertEqual(bundle["resources"]["MSIGDB"]["source_sets"], 35361)
        self.assertEqual(bundle["resources"]["MSIGDB"]["eligible_sets"], 29004)
        self.assertEqual(bundle["resources"]["MAGMA_GENE_SETS"]["source_sets"], 17023)
        self.assertEqual(bundle["resources"]["MAGMA_GENE_SETS"]["eligible_sets"], 12960)
        with (ROOT / policy["source_registry"]).open(encoding="utf-8", newline="") as handle:
            sources = {row["source_id"]: row for row in csv.DictReader(handle, delimiter="\t")}
        for source_id in spec["source_ids"]:
            self.assertEqual(sources[source_id]["source_status"], "SOURCE_VERIFIED")
            self.assertEqual(sources[source_id]["exact_release"], bundle["resources"][source_id]["release"])
        dispatcher = (ROOT / "scripts/76_run_interpretation_task.py").read_text(encoding="utf-8")
        self.assertIn("92_run_pathway_task.py", dispatcher)
        workflow = (ROOT / "Snakefile").read_text(encoding="utf-8")
        self.assertIn("PATHWAY_SOURCE_PATHS", workflow)

    def test_pathway_parsers_and_hypergeometric_tail_are_fail_closed(self):
        module = load_numbered_script("92_run_pathway_task.py", "pathway_task")
        with tempfile.TemporaryDirectory() as temporary:
            directory = Path(temporary)
            gtf = directory / "genes.gtf"
            gtf.write_text(
                "chr1\tX\tgene\t1\t2\t.\t+\t.\tgene_id \"ENSG00000000001.1\"; gene_name \"A\";\n"
                "chr1\tX\tgene\t3\t4\t.\t+\t.\tgene_id \"ENSG00000000002.1\"; gene_name \"B\";\n"
                "chr1\tX\tgene\t5\t6\t.\t+\t.\tgene_id \"ENSG00000000003.1\"; gene_name \"DUP\";\n"
                "chr1\tX\tgene\t7\t8\t.\t+\t.\tgene_id \"ENSG00000000004.1\"; gene_name \"DUP\";\n",
                encoding="utf-8",
            )
            symbols, stats = module.load_unique_gencode_symbols(gtf)
            self.assertEqual(symbols, {"A": "ENSG00000000001", "B": "ENSG00000000002"})
            self.assertEqual(stats["ambiguous_symbols"], 1)
            archive = directory / "reactome.zip"
            with zipfile.ZipFile(archive, "w") as handle:
                handle.writestr("ReactomePathways.gmt", "Path\tR-HSA-1\tA\tB\tDUP\n")
            sets, observed = module.load_reactome_sets(archive, symbols, 1, 10)
            self.assertEqual(sets["R-HSA-1"][1], {"ENSG00000000001", "ENSG00000000002"})
            self.assertEqual(observed["eligible_sets"], 1)
            msigdb = directory / "msigdb.gmt"
            msigdb.write_text(
                "SET_A\thttps://www.gsea-msigdb.org/gsea/msigdb/human/geneset/SET_A\tA\tB\tDUP\n",
                encoding="utf-8",
            )
            sets, observed = module.load_msigdb_symbol_sets(msigdb, symbols, 1, 10)
            self.assertEqual(sets["SET_A"][1], {"ENSG00000000001", "ENSG00000000002"})
            self.assertEqual(observed["mapped_genes"], 2)
            magma = directory / "magma.gmt"
            magma.write_text(
                "SET_B\thttp://www.gsea-msigdb.org/gsea/msigdb/human/geneset/SET_B\t"
                "ENSG00000000001\tENSG00000000002\tENSG00000000009\n",
                encoding="utf-8",
            )
            sets, observed = module.load_magma_gene_sets(
                magma, {"ENSG00000000001", "ENSG00000000002"}, 1, 10,
            )
            self.assertEqual(sets["SET_B"][1], {"ENSG00000000001", "ENSG00000000002"})
            self.assertEqual(observed["source_gene_ids"], 3)
        expected = sum(
            math.comb(3, value) * math.comb(7, 4 - value) / math.comb(10, 4)
            for value in range(2, 4)
        )
        self.assertAlmostEqual(module.hypergeometric_right_tail(2, 10, 3, 4), expected)

    def test_causal_runtime_and_estimator_family_are_exactly_pre_result_locked(self):
        policy = json.loads((ROOT / "config/interpretation_analysis_policy.json").read_text(encoding="utf-8"))
        spec = policy["causal_inference"]
        manifest_path = ROOT / spec["component_manifest"]
        self.assertEqual(hashlib.sha256(manifest_path.read_bytes()).hexdigest(), spec["component_manifest_sha256"])
        manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
        self.assertEqual(manifest["runtime"]["r_version"], "4.3.3")
        self.assertEqual(manifest["runtime"]["architecture"], "aarch64")
        self.assertEqual(manifest["runtime"]["dependency_source_count"], 63)
        self.assertEqual(manifest["runtime"]["installed_dependency_closure_count"], 153)
        self.assertEqual({row["component_id"] for row in manifest["components"]}, {
            "TWOSAMPLEMR_SOURCE", "MRPRESSO_SOURCE", "CAUSE_SOURCE", "LHCMR_SOURCE",
            "PLINK_1_9_STABLE_MAC_ARCHIVE", "PLINK_1_9_STABLE_UNIVERSAL_BINARY",
        })
        self.assertEqual(manifest["ld_reference"]["sample_count"], 503)
        self.assertEqual(manifest["ld_reference"]["bed_record_bytes"], 126)
        self.assertTrue(manifest["determinism"]["local_ld_only"])
        self.assertEqual(manifest["determinism"]["harmonise_action"], 3)
        self.assertEqual(manifest["estimators"]["MR_PRESSO"]["minimum_instruments"], 4)
        self.assertEqual(manifest["estimators"]["MR_PRESSO"]["nb_distribution"], 100000)
        robust = manifest["estimators"]["CAUSE_or_LHC_MR_where_appropriate"]
        self.assertEqual(robust["primary"], "CAUSE")
        self.assertEqual(robust["forbidden_lhc_mr_fallbacks"], ["run_ldsc=FALSE", "run_MR=FALSE"])

        dependency_path = ROOT / manifest["runtime"]["dependency_source_manifest"]
        installed_path = ROOT / manifest["runtime"]["installed_package_manifest"]
        self.assertEqual(hashlib.sha256(dependency_path.read_bytes()).hexdigest(), manifest["runtime"]["dependency_source_manifest_sha256"])
        self.assertEqual(hashlib.sha256(installed_path.read_bytes()).hexdigest(), manifest["runtime"]["installed_package_manifest_sha256"])
        with dependency_path.open(encoding="utf-8", newline="") as handle:
            dependencies = list(csv.DictReader(handle, delimiter="\t"))
        with installed_path.open(encoding="utf-8", newline="") as handle:
            installed = list(csv.DictReader(handle, delimiter="\t"))
        self.assertEqual(len(dependencies), 63)
        self.assertEqual(len(installed), 153)
        self.assertEqual({row["package"] for row in installed if row["library"] == "causal"}, {
            row["package"] for row in dependencies
        } | {"TwoSampleMR", "MRPRESSO", "cause", "lhcMR"})

        with (ROOT / policy["source_registry"]).open(encoding="utf-8", newline="") as handle:
            sources = [row for row in csv.DictReader(handle, delimiter="\t") if row["analysis_family"] == "causal"]
        self.assertEqual({row["method"] for row in sources}, set(spec["methods"]))
        self.assertTrue(all(row["source_status"] == "SOURCE_VERIFIED" for row in sources))
        setup = (ROOT / "scripts/93_setup_causal_runtime.sh").read_text(encoding="utf-8")
        self.assertIn("CAUSAL_RUNTIME_READY", setup)
        self.assertIn("plink_mac_20250819.zip", setup)
        self.assertNotIn("install_github", setup)
        preflight = (ROOT / "scripts/74_interpretation_preflight.py").read_text(encoding="utf-8")
        self.assertIn("validate_causal_runtime_bundle", preflight)
        self.assertIn("causal installed dependency closure differs", preflight)

    def test_causal_reference_builder_streams_exact_snp_major_records(self):
        module = load_numbered_script("94_prepare_causal_reference.py", "causal_reference")
        bim = (
            b"1 rs1 0 1 A G\n"
            b"1 rs2 0 2 C T\n"
            b"1 rs3 0 3 C A\n"
        )
        bed = b"\x6c\x1b\x01" + b"ab" + b"cd" + b"ef"
        with tempfile.TemporaryDirectory() as temporary:
            directory = Path(temporary)
            source = directory / "reference.zip"
            with zipfile.ZipFile(source, "w") as archive:
                archive.writestr("g1000_eur.bim", bim)
                archive.writestr("g1000_eur.bed", bed)
            with zipfile.ZipFile(source) as archive:
                indices, counters, bim_hash = module.stream_bim(
                    archive, archive.getinfo("g1000_eur.bim"),
                    {"bytes": len(bim), "sha256": hashlib.sha256(bim).hexdigest()},
                    {"rs1": frozenset(("A", "G")), "rs3": frozenset(("A", "C"))},
                    directory / "subset.bim", 3, 2,
                )
                bed_hash = module.stream_bed(
                    archive, archive.getinfo("g1000_eur.bed"),
                    {"bytes": len(bed), "sha256": hashlib.sha256(bed).hexdigest()},
                    indices, directory / "subset.bed", 2, 3,
                )
            self.assertEqual(indices, [0, 2])
            self.assertEqual(counters["allele_matched_variants"], 2)
            self.assertEqual((directory / "subset.bim").read_bytes(), b"1 rs1 0 1 A G\n1 rs3 0 3 C A\n")
            self.assertEqual((directory / "subset.bed").read_bytes(), b"\x6c\x1b\x01abef")
            self.assertEqual(bim_hash, hashlib.sha256(bim).hexdigest())
            self.assertEqual(bed_hash, hashlib.sha256(bed).hexdigest())

    def test_robustness_applicability_is_locked_before_results(self):
        policy = json.loads((ROOT / "config/interpretation_analysis_policy.json").read_text(encoding="utf-8"))
        families = policy["robustness"]["required_families"]
        for conclusion_type, matrix in policy["robustness"]["applicability_by_conclusion_type"].items():
            self.assertIn(conclusion_type, policy["robustness"]["major_conclusion_edge_levels"])
            self.assertEqual(list(matrix), families)
            self.assertTrue(all(value == "APPLICABLE" or len(value) > 15 for value in matrix.values()))
        workflow = (ROOT / "Snakefile").read_text(encoding="utf-8")
        for rule in ("interpretation_preflight", "interpretation_task", "interpretation", "atlas_edges", "robustness_task", "robustness"):
            self.assertIn(f"rule {rule}:", workflow)


if __name__ == "__main__":
    unittest.main()
