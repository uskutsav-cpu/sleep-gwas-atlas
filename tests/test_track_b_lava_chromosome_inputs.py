import csv
import gzip
import hashlib
import io
import json
import os
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "scripts/129_prepare_track_b_lava_chromosome_inputs.py"
TRAITS = [
    "snoring",
    "parental_lifespan",
    "insomnia",
    "adhd",
    "frailty",
    "bmi",
    "sleep_apnea",
    "mdd",
]
INFO_FIELDS = ["phenotype", "cases", "controls", "prevalence", "filename"]
PROVENANCE_FIELDS = [
    "phenotype",
    "filename",
    "bytes",
    "sha256",
    "observed_columns",
    "validation_status",
]


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def write_tsv(path: Path, fields: list[str], rows: list[dict[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def write_gzip(path: Path, lines: list[bytes]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("wb") as raw:
        with gzip.GzipFile(filename="", fileobj=raw, mode="wb", compresslevel=6, mtime=0) as stream:
            for line in lines:
                stream.write(line)


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def build_fixture(root: Path) -> dict[str, list[bytes]]:
    policy_path = root / "config/track_b_local_analysis_policy.json"
    policy_path.parent.mkdir(parents=True, exist_ok=True)
    policy = {
        "analysis_id": "track-b-v1.0-local",
        "expected_analysis_traits": 8,
        "trait_order": TRAITS,
        "reference_prefix": "ref/lava/test-ref",
    }
    policy_path.write_text(json.dumps(policy, sort_keys=True) + "\n", encoding="utf-8")

    reference_records = []
    for chromosome in range(1, 23):
        path = root / f"ref/lava/test-ref_chr{chromosome}.info"
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(
            "SNP\tCHR\tPOS\tA1\tA2\tNOBS\tNMISS\tFREQ\tNCORRS\n"
            f"RS{chromosome}A\t{chromosome}\t100\tA\tG\t1000\t0\t0.2\t999\n"
            f"rs{chromosome}b\t{chromosome}\t200\tC\tT\t1000\t0\t0.3\t999\n",
            encoding="utf-8",
        )
        reference_records.append({
            "path": path.relative_to(root).as_posix(),
            "bytes": path.stat().st_size,
            "sha256": sha256(path),
        })
    reference_provenance = root / "ref/lava/reference.provenance.json"
    reference_provenance.write_text(
        json.dumps(
            {
                "schema_version": "sleep-atlas-lava-reference.1",
                "extracted_file_count": len(reference_records),
                "extracted_files": reference_records,
                "verification": "SHA-256 verified after official HTTPS acquisition",
            },
            indent=2,
            sort_keys=True,
        )
        + "\n",
        encoding="utf-8",
    )

    header = b"SNP\tA1\tA2\tZ\tN\n"
    source_lines: dict[str, list[bytes]] = {}
    input_rows = []
    provenance_rows = []
    for trait_index, phenotype in enumerate(TRAITS, start=1):
        lines = [header, f"rs_not_reference\tA\tG\t{trait_index}.01\t1000\n".encode()]
        # Reverse chromosome order and put B before A so order preservation is observable.
        for chromosome in range(22, 0, -1):
            # A reference-matching row with absent analysis fields must still be
            # retained: this tool performs no scientific completeness filtering.
            lines.append(
                b"RS1B\t\t\t\t\n"
                if chromosome == 1
                else f"RS{chromosome}B\tT\tC\t-{trait_index}.{chromosome:02d}\t1000\n".encode()
            )
            lines.append(f"rs{chromosome}a\tA\tG\t{trait_index}.{chromosome:02d}\t1000\n".encode())
        source = root / f"data/munged/{phenotype}.sumstats.gz"
        write_gzip(source, lines)
        source_lines[phenotype] = lines
        relative = source.relative_to(root).as_posix()
        input_rows.append({
            "phenotype": phenotype,
            "cases": "NA",
            "controls": "NA",
            "prevalence": "NA",
            "filename": relative,
        })
        provenance_rows.append({
            "phenotype": phenotype,
            "filename": relative,
            "bytes": source.stat().st_size,
            "sha256": sha256(source),
            "observed_columns": "SNP,A1,A2,Z,N",
            "validation_status": "VALIDATED_FOR_LAVA_LOCAL_NOT_FINE_MAPPING",
        })

    input_info = root / "results/track_b/lava_input_info.tsv"
    input_provenance = root / "results/track_b/lava_input_provenance.tsv"
    write_tsv(input_info, INFO_FIELDS, input_rows)
    write_tsv(input_provenance, PROVENANCE_FIELDS, provenance_rows)
    input_lock = {
        "analysis_id": policy["analysis_id"],
        "selection_timing": "BEFORE_LOCAL_RESULT_ACCESS",
        "local_results_accessed_before_input_freeze": False,
        "policy_sha256": sha256(policy_path),
        "input_artifact_sha256": {
            input_info.relative_to(root).as_posix(): sha256(input_info),
            input_provenance.relative_to(root).as_posix(): sha256(input_provenance),
        },
        "required_reference_provenance": reference_provenance.relative_to(root).as_posix(),
    }
    lock_path = root / "results/track_b/local_analysis_input.lock.json"
    lock_path.write_text(json.dumps(input_lock, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return source_lines


def run_script(root: Path, *arguments: str) -> subprocess.CompletedProcess[str]:
    environment = dict(os.environ)
    environment["PYTHONDONTWRITEBYTECODE"] = "1"
    return subprocess.run(
        [sys.executable, str(SCRIPT), "--root", str(root), *arguments],
        cwd=ROOT,
        env=environment,
        text=True,
        capture_output=True,
        check=False,
    )


class TrackBLAVAChromosomeInputTests(unittest.TestCase):
    def test_builds_exact_full_family_and_verify_is_read_only(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            source_lines = build_fixture(root)
            built = run_script(root)
            self.assertEqual(built.returncode, 0, built.stdout + built.stderr)
            self.assertIn("chromosomes=22 traits=8 shards=176", built.stdout)

            lock_path = root / "results/track_b/lava_chromosome_inputs.provenance.json"
            payload = json.loads(lock_path.read_text(encoding="utf-8"))
            self.assertEqual(payload["schema_version"], "track-b-lava-chromosome-inputs.1")
            self.assertEqual(payload["chromosomes"], list(range(1, 23)))
            self.assertEqual(payload["trait_order"], TRAITS)
            self.assertEqual(payload["totals"]["shard_count"], 22 * 8)
            self.assertEqual(len(payload["reference_info"]), 22)
            self.assertEqual(len(payload["chromosome_input_info"]), 22)

            shard = root / "results/track_b/lava_chromosome_inputs/chr01/snoring.sumstats.gz"
            with gzip.open(shard, "rb") as handle:
                observed = handle.readlines()
            expected = [source_lines["snoring"][0]] + [
                line
                for line in source_lines["snoring"][1:]
                if line.split(b"\t", 1)[0].decode().strip().lower() in {"rs1a", "rs1b"}
            ]
            self.assertEqual(observed, expected)
            self.assertEqual([line.split(b"\t", 1)[0] for line in observed[1:]], [b"RS1B", b"rs1a"])
            self.assertEqual(observed[1], b"RS1B\t\t\t\t\n")

            info = read_tsv(root / "results/track_b/lava_chromosome_inputs/chr01/lava_input_info.tsv")
            self.assertEqual([row["phenotype"] for row in info], TRAITS)
            self.assertEqual(
                info[0]["filename"],
                "results/track_b/lava_chromosome_inputs/chr01/snoring.sumstats.gz",
            )

            protected = [lock_path, shard, root / "results/track_b/lava_chromosome_inputs/chr01/lava_input_info.tsv"]
            before = {path: (path.stat().st_ino, path.stat().st_mtime_ns, sha256(path)) for path in protected}
            verified = run_script(root, "--verify")
            self.assertEqual(verified.returncode, 0, verified.stdout + verified.stderr)
            after = {path: (path.stat().st_ino, path.stat().st_mtime_ns, sha256(path)) for path in protected}
            self.assertEqual(after, before)

    def test_partial_family_resumes_without_replacing_completed_shards(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            build_fixture(root)
            first = run_script(root)
            self.assertEqual(first.returncode, 0, first.stdout + first.stderr)
            lock_path = root / "results/track_b/lava_chromosome_inputs.provenance.json"
            preserved = root / "results/track_b/lava_chromosome_inputs/chr01/snoring.sumstats.gz"
            missing = root / "results/track_b/lava_chromosome_inputs/chr22/mdd.sumstats.gz"
            preserved_identity = (preserved.stat().st_ino, preserved.stat().st_mtime_ns, sha256(preserved))
            missing_bytes = missing.read_bytes()
            lock_path.unlink()
            missing.unlink()

            resumed = run_script(root)
            self.assertEqual(resumed.returncode, 0, resumed.stdout + resumed.stderr)
            self.assertEqual(
                (preserved.stat().st_ino, preserved.stat().st_mtime_ns, sha256(preserved)),
                preserved_identity,
            )
            self.assertEqual(missing.read_bytes(), missing_bytes)
            self.assertTrue(lock_path.is_file())
            verified = run_script(root, "--verify")
            self.assertEqual(verified.returncode, 0, verified.stdout + verified.stderr)

    def test_tampered_shard_is_rejected_and_never_replaced(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            build_fixture(root)
            first = run_script(root)
            self.assertEqual(first.returncode, 0, first.stdout + first.stderr)
            shard = root / "results/track_b/lava_chromosome_inputs/chr07/adhd.sumstats.gz"
            write_gzip(
                shard,
                [
                    b"SNP\tA1\tA2\tZ\tN\n",
                    b"rs_not_reference\tA\tG\t99\t1000\n",
                ],
            )
            tampered = shard.read_bytes()

            verified = run_script(root, "--verify")
            self.assertNotEqual(verified.returncode, 0)
            self.assertIn("outside its chromosome reference", verified.stdout + verified.stderr)
            rebuilt = run_script(root)
            self.assertNotEqual(rebuilt.returncode, 0)
            self.assertEqual(shard.read_bytes(), tampered)

    def test_source_or_reference_binding_drift_fails_closed(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            build_fixture(root)
            first = run_script(root)
            self.assertEqual(first.returncode, 0, first.stdout + first.stderr)
            source = root / "data/munged/snoring.sumstats.gz"
            with source.open("ab") as handle:
                handle.write(b"tamper")
            verified = run_script(root, "--verify")
            self.assertNotEqual(verified.returncode, 0)
            self.assertIn("source differs from frozen provenance", verified.stdout + verified.stderr)

        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            build_fixture(root)
            first = run_script(root)
            self.assertEqual(first.returncode, 0, first.stdout + first.stderr)
            reference = root / "ref/lava/test-ref_chr12.info"
            with reference.open("a", encoding="utf-8") as handle:
                handle.write("rs12tamper\t12\t300\tA\tC\t1000\t0\t0.4\t999\n")
            verified = run_script(root, "--verify")
            self.assertNotEqual(verified.returncode, 0)
            self.assertIn("reference differs from sealed provenance", verified.stdout + verified.stderr)

    def test_script_parses_without_import_side_effects(self) -> None:
        result = subprocess.run(
            [
                sys.executable,
                "-c",
                "import ast,pathlib; ast.parse(pathlib.Path('scripts/129_prepare_track_b_lava_chromosome_inputs.py').read_text())",
            ],
            cwd=ROOT,
            env={**os.environ, "PYTHONDONTWRITEBYTECODE": "1"},
            text=True,
            capture_output=True,
            check=False,
        )
        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)


if __name__ == "__main__":
    unittest.main()
