import csv
import gzip
import hashlib
import importlib.util
import math
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest


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

    def test_pending_sources_do_not_reuse_stage_or_headline_sample_counts(self):
        panel = {row["trait_id"]: row for row in self.panel}
        self.assertEqual(
            (panel["ms"]["ncase"], panel["ms"]["ncontrol"], panel["ms"]["n_total"]),
            ("14802", "26703", "41505"),
        )
        self.assertEqual(panel["ms"]["source_status"], "SOURCE_PENDING")
        self.assertEqual(
            (
                panel["melanoma"]["ncase"],
                panel["melanoma"]["ncontrol"],
                panel["melanoma"]["n_total"],
                panel["melanoma"]["build"],
            ),
            ("30134", "81415", "111549", "hg38"),
        )
        self.assertEqual(panel["melanoma"]["source_status"], "SOURCE_PENDING")
        self.assertEqual(panel["t2d"]["source_status"], "SOURCE_PENDING")

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
            "##contig=<ID=1,length=249250621,assembly=GRCh37>",
            "##SAMPLE=<ID=ieu-b-38,TotalVariants=2,VariantsNotRead=0,HarmonisedVariants=2,VariantsNotHarmonised=0,SwitchedAlleles=1,NormalisedVariants=0,StudyType=Continuous>",
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
        )
        with tempfile.TemporaryDirectory() as directory:
            directory = Path(directory)
            source = directory / "ieu-b-38.vcf.gz"
            first = directory / "sbp-first.txt.gz"
            second = directory / "sbp-second.txt.gz"
            with gzip.open(source, "wt", encoding="utf-8", newline="") as handle:
                handle.write(payload)
            counts = materializer.materialize_vcf(
                source, first, expected_variants=2
            )
            materializer.materialize_vcf(source, second, expected_variants=2)
            with gzip.open(first, "rt", newline="") as handle:
                rows = list(csv.DictReader(handle, delimiter="\t"))
            self.assertEqual(first.read_bytes(), second.read_bytes())
        self.assertEqual(counts, {
            "source_rows": 2,
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
        self.assertIn("BLOCKED source_curation", result.stdout)
        self.assertIn("BLOCKED source_schemas", result.stdout)
        self.assertIn("Acceptance: 1/23 gates passed", result.stdout)


if __name__ == "__main__":
    unittest.main()
