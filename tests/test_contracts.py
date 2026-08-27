import csv
import gzip
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
