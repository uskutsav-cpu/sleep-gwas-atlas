import csv
import gzip
import importlib.util
import io
import json
import math
import shutil
import subprocess
import tempfile
import unittest
from pathlib import Path
from types import SimpleNamespace
from unittest import mock


ROOT = Path(__file__).resolve().parents[1]


def load_materializer():
    spec = importlib.util.spec_from_file_location(
        "track_b_placo_materializer_test", ROOT / "scripts/125_materialize_track_b_placo_pair.py"
    )
    if spec is None or spec.loader is None:
        raise RuntimeError("could not load Track B PLACO+ materializer")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def read_tsv(path: Path, *, compressed: bool = False) -> list[dict[str, str]]:
    opener = gzip.open if compressed else path.open
    if compressed:
        handle = opener(path, "rt", encoding="utf-8", newline="")
    else:
        handle = opener(newline="", encoding="utf-8")
    with handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_gzip_tsv(path: Path, fields: list[str], rows: list[dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("wb") as raw:
        with gzip.GzipFile(filename="", mode="wb", fileobj=raw, mtime=0) as compressed:
            with io.TextIOWrapper(compressed, encoding="utf-8", newline="") as text:
                writer = csv.DictWriter(text, fieldnames=fields, delimiter="\t", lineterminator="\n")
                writer.writeheader()
                writer.writerows(rows)


def write_tsv(path: Path, fields: list[str], rows: list[dict[str, str]]) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def bh(values: list[float]) -> list[float]:
    count = len(values)
    order = sorted(range(count), key=values.__getitem__)
    adjusted = [0.0] * count
    running = 1.0
    for rank0 in range(count - 1, -1, -1):
        index = order[rank0]
        running = min(running, min(1.0, values[index] * count / (rank0 + 1)))
        adjusted[index] = running
    return adjusted


def dense_rows(*, second: bool) -> list[dict[str, str]]:
    rows: list[dict[str, str]] = []
    next_id = 1
    for chromosome in range(1, 23):
        z = 7.0 if chromosome == 1 and not second else ((chromosome * (3 if second else 2)) % 11 - 5) / 2
        z2 = ((chromosome * 5) % 13 - 6) / 2
        observed_z = z2 if second else z
        a1, a2 = ("C", "A") if second and chromosome == 1 else ("A", "C")
        rows.append({
            "SNP": f"rs{next_id}", "CHR": str(chromosome), "BP": str(chromosome * 1000),
            "A1": a1, "A2": a2, "FRQ": "0.4", "BETA": f"{observed_z:.12g}",
            "SE": "1", "P": "0.5", "N": "10000",
        })
        next_id += 1
    for extra in range(8):
        observed_z = ((extra * (4 if second else 3)) % 9 - 4) / 2
        rows.append({
            "SNP": f"rs{next_id}", "CHR": "1", "BP": str(2000 + extra),
            "A1": "A", "A2": "C", "FRQ": "0.35", "BETA": f"{observed_z:.12g}",
            "SE": "1", "P": "0.4", "N": "10000",
        })
        next_id += 1
    rows.append({
        "SNP": "rs90001", "CHR": "1", "BP": "5000", "A1": "A",
        "A2": "G" if second else "C", "FRQ": "0.3", "BETA": "0.5", "SE": "1", "P": "0.5", "N": "10000",
    })
    rows.append({
        "SNP": "rs90002", "CHR": "2", "BP": "6001" if second else "6000", "A1": "A",
        "A2": "C", "FRQ": "0.3", "BETA": "0.5", "SE": "1", "P": "0.5", "N": "10000",
    })
    duplicate = {
        "SNP": "rs90003", "CHR": "3", "BP": "7000", "A1": "A", "A2": "C", "FRQ": "0.3",
        "BETA": "0.5", "SE": "1", "P": "0.5", "N": "10000",
    }
    rows.extend([dict(duplicate), dict(duplicate)])
    rows.append({
        "SNP": "rs90004", "CHR": "4", "BP": "8000", "A1": "A", "A2": "C", "FRQ": "0.3",
        "BETA": "0.5", "SE": "0", "P": "0.5", "N": "10000",
    })
    rows.append({
        "SNP": "rs90005", "CHR": "5", "BP": "9000", "A1": "A", "A2": "C", "FRQ": "0.3",
        "BETA": "9", "SE": "1", "P": "1e-10", "N": "10000",
    })
    rows.append({
        "SNP": "rs91000", "CHR": "23", "BP": "10000", "A1": "A", "A2": "C", "FRQ": "0.3",
        "BETA": "0.5", "SE": "1", "P": "0.5", "N": "10000",
    })
    return rows


def fixture_source_text(force_failure: bool, force_hits: bool = False) -> str:
    if not force_failure and not force_hits:
        return (ROOT / ".r-env/share/placo/PLACO_v0.2.0.R").read_text(encoding="utf-8")
    forced_p = "1e-12" if force_hits else "min(1, 0.2 + exp(-abs(statistic)))"
    return """
var.placo <- function(Z.matrix, P.matrix, p.threshold=1e-4) diag(var(Z.matrix))
cor.pearson <- function(Z.matrix, P.matrix, p.threshold=1e-4, returnMatrix=TRUE) {
  value <- cor(Z.matrix)[1,2]
  if (returnMatrix) cor(Z.matrix) else value
}
placo.plus <- function(Z, VarZ, CorZ, AbsTol=1e-8) {
  if (FORCE_FAILURE && abs(Z[1] - 7) < 1e-12) stop("synthetic integration failure")
  statistic <- prod(Z)
  list(T.placo.plus=statistic, p.placo.plus=FORCED_P)
}
""".replace("FORCE_FAILURE", "TRUE" if force_failure else "FALSE").replace("FORCED_P", forced_p).lstrip()


def make_context(
    module, root: Path, *, pair_id: str = "A", force_failure: bool = False, force_hits: bool = False,
    maximum_failure_fraction: float | None = None,
) -> dict:
    (root / "scripts").mkdir(parents=True, exist_ok=True)
    shutil.copy2(ROOT / "scripts/125_materialize_track_b_placo_pair.py", root / "scripts/125_materialize_track_b_placo_pair.py")
    shutil.copy2(ROOT / "scripts/126_run_track_b_placo_pair.R", root / "scripts/126_run_track_b_placo_pair.R")
    placo_source = root / ".r-env/share/placo/PLACO_v0.2.0.R"
    placo_source.parent.mkdir(parents=True, exist_ok=True)
    placo_source.write_text(fixture_source_text(force_failure, force_hits), encoding="utf-8")
    rscript = root / ".r-env/bin/Rscript"
    rscript.parent.mkdir(parents=True, exist_ok=True)
    rscript.symlink_to(ROOT / ".r-env/bin/Rscript")
    first = root / "data/first.tsv.gz"
    second = root / "data/second.tsv.gz"
    first_rows, second_rows = dense_rows(second=False), dense_rows(second=True)
    write_gzip_tsv(first, module.SOURCE_FIELDS, first_rows)
    write_gzip_tsv(second, module.SOURCE_FIELDS, second_rows)
    policy = {
        "analysis_id": "track-b-placo-test",
        "dense_input_contract": {
            "alignment_rule": "exact test alignment",
            "pair_materialization_plausibility": {
                "minimum_aligned_eligible_variants": 30,
                "minimum_fraction_of_smaller_dense_input": 0.8,
                "required_autosomes": list(range(1, 23)),
                "nonzero_eligible_variants_on_every_required_autosome": True,
                "per_autosome_provenance_required": True,
            },
        },
        "placo_plus": {
            "source_sha256": module.sha256(placo_source),
            "marginal_p_threshold_for_nuisance_estimation": 0.0001,
            "z_squared_maximum": 80,
            "absolute_tolerance": 1e-8,
            "maximum_numerical_failure_fraction": (
                maximum_failure_fraction
                if maximum_failure_fraction is not None else (0.1 if force_failure else 0.0001)
            ),
            "primary_pair_family_headline_threshold": 2.5e-8,
            "per_pair_genome_wide_threshold": 5e-8,
            "control_genome_wide_threshold": 5e-8,
            "within_pair_bh_alpha": 0.05,
        },
    }
    (root / "config").mkdir(parents=True, exist_ok=True)
    (root / module.POLICY).write_text(json.dumps(policy) + "\n", encoding="utf-8")
    (root / module.CONTRACT_LOCK).parent.mkdir(parents=True, exist_ok=True)
    (root / module.CONTRACT_LOCK).write_text("synthetic contract lock\n", encoding="utf-8")
    (root / module.INPUT_LOCK).write_text("synthetic input lock\n", encoding="utf-8")
    trait1, trait2, family_role = module.PAIR_IDENTITIES[pair_id]
    sources = []
    for trait, path, rows in ((trait1, first, first_rows), (trait2, second, second_rows)):
        size, digest = module.file_identity(path)
        sources.append({
            "trait_id": trait, "path": str(path.relative_to(root)), "absolute_path": path,
            "bytes": size, "rows": len(rows), "sha256": digest, "schema": ",".join(module.SOURCE_FIELDS),
        })
    return {
        "root": root, "pair_id": pair_id, "trait1": trait1, "trait2": trait2,
        "family_role": family_role, "policy": policy, "contract_lock": {}, "input_lock": {},
        "sources": sources, "placo_source": placo_source, "fingerprint": ("e" if force_failure else "d") * 64,
        "fingerprint_payload": {},
    }


class TrackBPLACORunnerTests(unittest.TestCase):
    def test_preflight_has_storage_gate_but_no_assumed_whole_pipeline_ram_gate(self) -> None:
        module = load_materializer()
        with tempfile.TemporaryDirectory() as directory:
            context = make_context(module, Path(directory))
            with (
                mock.patch.object(module, "physical_memory_bytes", return_value=512 * 1024**2),
                mock.patch.object(module.shutil, "disk_usage", return_value=SimpleNamespace(free=32 * 1024**3)),
            ):
                resources = module.resource_preflight(context)
            self.assertEqual(resources["physical_memory_bytes"], 512 * 1024**2)
            self.assertLess(resources["estimated_work_bytes"], resources["free_bytes"])

    def test_all_three_frozen_pair_identities_materialize_exact_tasks(self) -> None:
        module = load_materializer()
        for pair_id, identity in module.PAIR_IDENTITIES.items():
            with self.subTest(pair_id=pair_id), tempfile.TemporaryDirectory() as directory:
                root = Path(directory)
                context = make_context(module, root, pair_id=pair_id)
                module.materialize_pair(context)
                task = read_tsv(module.materialized_paths(context)["task"])[0]
                self.assertEqual(
                    (task["trait1"], task["trait2"], task["family_role"]),
                    identity,
                )
                self.assertEqual(float(task["primary_headline_threshold"]), 2.5e-8)
                self.assertEqual(float(task["pairwise_gws_threshold"]), 5e-8)
                self.assertEqual(task["maximum_workers"], "4")
                self.assertEqual(task["pair_concurrency_limit"], "1")

    def test_upstream_block_stops_before_materialization(self) -> None:
        module = load_materializer()
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            row = {
                "dense_input_gate": "PASS_ALL_THREE_PAIRS",
                "replication_gate": "PASS_TERMINAL_PRIMARY_REPLICATION_FAMILY",
                "local_analysis_gate": "BLOCKED_LOCAL_ANALYSIS_NOT_TERMINAL",
                "software_gate": "READY_PINNED_PLACO_SOURCE",
            }
            with self.assertRaisesRegex(SystemExit, "BLOCKED_UPSTREAM.*LOCAL_ANALYSIS"):
                module.assert_upstream_ready(row)
            self.assertEqual(list(root.iterdir()), [])

    def test_materializer_exact_alignment_rejections_and_floors(self) -> None:
        module = load_materializer()
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            context = make_context(module, root)
            provenance = module.materialize_pair(context)
            counts = provenance["alignment_counts"]
            self.assertEqual(counts["trait1_rows"], 37)
            self.assertEqual(counts["trait2_rows"], 37)
            self.assertEqual(counts["duplicate_ids_dropped"], 4)
            self.assertEqual(counts["invalid_z_or_p_dropped"], 2)
            self.assertEqual(counts["allele_mismatches_dropped"], 1)
            self.assertEqual(counts["allele_flips"], 1)
            self.assertEqual(counts["z_squared_exclusions"], 1)
            self.assertEqual(counts["eligible_written"], 30)
            self.assertEqual(provenance["additional_rejection_counts"]["invalid_identity_or_frequency_dropped"], 2)
            self.assertGreaterEqual(provenance["coverage_fraction_of_smaller_frozen_dense_input"], 0.8)
            self.assertEqual([row["CHR"] for row in provenance["per_autosome"]], list(range(1, 23)))
            self.assertTrue(all(row["eligible_written"] > 0 for row in provenance["per_autosome"]))
            paths = module.materialized_paths(context)
            self.assertTrue(all(paths[name].is_file() for name in ("aligned", "provenance", "task")))
            aligned = read_tsv(paths["aligned"], compressed=True)
            self.assertEqual(len(aligned), 30)
            self.assertEqual(list(aligned[0]), module.ALIGNED_FIELDS)
            self.assertEqual(aligned[0]["Z2"], "0.5")
            self.assertEqual(aligned[0]["FRQ2"], "0.6")
            self.assertFalse(any(row["CHR"] == "23" for row in aligned))
            self.assertEqual(module.verify_materialized(context)["qc_status"], "PASS")

    def test_materialization_floor_failure_publishes_no_pair_or_task(self) -> None:
        module = load_materializer()
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            context = make_context(module, root)
            context["policy"]["dense_input_contract"]["pair_materialization_plausibility"]["minimum_aligned_eligible_variants"] = 31
            with self.assertRaisesRegex(SystemExit, "FAILED_QC.*plausibility floor"):
                module.materialize_pair(context)
            paths = module.materialized_paths(context)
            self.assertFalse(any(paths[name].exists() for name in ("aligned", "provenance", "task")))

    def test_fraction_and_autosome_floors_fail_closed(self) -> None:
        module = load_materializer()
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            context = make_context(module, root)
            floors = context["policy"]["dense_input_contract"]["pair_materialization_plausibility"]
            floors["minimum_aligned_eligible_variants"] = 1
            floors["minimum_fraction_of_smaller_dense_input"] = 0.82
            with self.assertRaisesRegex(SystemExit, "FAILED_QC.*coverage"):
                module.materialize_pair(context)
            paths = module.materialized_paths(context)
            self.assertFalse(any(paths[name].exists() for name in ("aligned", "provenance", "task")))

        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            context = make_context(module, root)
            first_rows = dense_rows(second=False)
            for row in first_rows:
                if row["CHR"] == "22":
                    row["BETA"] = "9"
                    break
            source = context["sources"][0]
            write_gzip_tsv(source["absolute_path"], module.SOURCE_FIELDS, first_rows)
            source["bytes"], source["sha256"] = module.file_identity(source["absolute_path"])
            source["rows"] = len(first_rows)
            floors = context["policy"]["dense_input_contract"]["pair_materialization_plausibility"]
            floors["minimum_aligned_eligible_variants"] = 1
            floors["minimum_fraction_of_smaller_dense_input"] = 0
            with self.assertRaisesRegex(SystemExit, r"missing_autosomes=\[22\]"):
                module.materialize_pair(context)
            paths = module.materialized_paths(context)
            self.assertFalse(any(paths[name].exists() for name in ("aligned", "provenance", "task")))

    def test_exact_frq_bearing_source_schema_is_required(self) -> None:
        module = load_materializer()
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            context = make_context(module, root)
            source = context["sources"][0]
            fields = [field for field in module.SOURCE_FIELDS if field != "FRQ"]
            rows = [{key: value for key, value in row.items() if key != "FRQ"} for row in dense_rows(second=False)]
            write_gzip_tsv(source["absolute_path"], fields, rows)
            source["bytes"], source["sha256"] = module.file_identity(source["absolute_path"])
            source["rows"] = len(rows)
            with self.assertRaisesRegex(SystemExit, "exact FRQ-bearing contract"):
                module.materialize_pair(context)
            paths = module.materialized_paths(context)
            self.assertFalse(any(paths[name].exists() for name in ("aligned", "provenance", "task")))

    def test_materialization_tamper_is_not_replaced(self) -> None:
        module = load_materializer()
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            context = make_context(module, root)
            module.materialize_pair(context)
            paths = module.materialized_paths(context)
            original = paths["aligned"].read_bytes()
            paths["aligned"].write_bytes(original + b"tamper")
            with self.assertRaisesRegex(SystemExit, "differs from its immutable provenance"):
                module.materialize_pair(context)
            self.assertEqual(paths["aligned"].read_bytes(), original + b"tamper")

    def test_official_runner_benchmark_preserves_global_nuisance_then_shards(self) -> None:
        module = load_materializer()
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            context = make_context(module, root)
            module.materialize_pair(context)
            paths = module.materialized_paths(context)
            measured = module.benchmark_pair(context, 5, 2)
            raw = read_tsv(paths["benchmark_raw"])[0]
            final = read_tsv(paths["benchmark"])[0]
            self.assertEqual(raw["input_rows"], "30")
            self.assertEqual(raw["benchmark_variants"], "5")
            self.assertEqual(raw["workers"], "2")
            self.assertEqual(raw["global_nuisance_input_rows"], "30")
            self.assertEqual(raw["nuisance_null_rows_variance"], "30")
            self.assertEqual(raw["nuisance_null_rows_correlation"], "30")
            self.assertIn("GLOBAL_OFFICIAL_NUISANCE_ESTIMATED_ON_ALL_VALID_VARIANTS", raw["scientific_equivalence"])
            self.assertGreater(float(raw["variants_per_second"]), 0)
            self.assertGreater(float(raw["projected_full_family_seconds"]), 0)
            self.assertGreater(int(final["peak_process_tree_rss_bytes"]), 0)
            self.assertGreater(float(final["wrapper_wall_seconds"]), 0)
            self.assertEqual(final["input_sha256"], module.sha256(paths["aligned"]))
            self.assertEqual(measured["peak_process_tree_rss_bytes"], final["peak_process_tree_rss_bytes"])
            self.assertFalse(paths["staged_ledger"].exists())
            with self.assertRaisesRegex(SystemExit, "benchmark already exists"):
                module.benchmark_pair(context, 5, 2)

    def test_full_failure_ledger_exact_bh_tamper_and_no_overwrite(self) -> None:
        module = load_materializer()
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            context = make_context(module, root, force_failure=True)
            module.materialize_pair(context)
            paths = module.materialized_paths(context)
            result = subprocess.run(
                [str(ROOT / ".r-env/bin/Rscript"), str(root / module.RUNNER), str(paths["task"].relative_to(root)),
                 "--root", str(root), "--execute", "--workers", "1", "--stage-only"],
                cwd=root, check=False, text=True, capture_output=True,
            )
            self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
            checkpoint = next(paths["checkpoint_dir"].glob("shard_*.rds"))
            checkpoint_original = checkpoint.read_bytes()
            checkpoint.write_bytes(checkpoint_original + b"tamper")
            tampered_resume = subprocess.run(
                [str(ROOT / ".r-env/bin/Rscript"), str(root / module.RUNNER), str(paths["task"].relative_to(root)),
                 "--root", str(root), "--execute", "--workers", "1", "--stage-only"],
                cwd=root, check=False, text=True, capture_output=True,
            )
            self.assertNotEqual(tampered_resume.returncode, 0)
            self.assertIn("tampered PLACO+ shard checkpoint", tampered_resume.stdout + tampered_resume.stderr)
            checkpoint.write_bytes(checkpoint_original)
            staged_hashes = (module.sha256(paths["staged_ledger"]), module.sha256(paths["run_summary"]))
            resumed = subprocess.run(
                [str(ROOT / ".r-env/bin/Rscript"), str(root / module.RUNNER), str(paths["task"].relative_to(root)),
                 "--root", str(root), "--execute", "--workers", "1", "--stage-only"],
                cwd=root, check=False, text=True, capture_output=True,
            )
            self.assertEqual(resumed.returncode, 0, resumed.stdout + resumed.stderr)
            self.assertEqual(
                staged_hashes,
                (module.sha256(paths["staged_ledger"]), module.sha256(paths["run_summary"])),
            )
            ledger = read_tsv(paths["staged_ledger"], compressed=True)
            self.assertEqual(len(ledger), 30)
            self.assertEqual({row["within_pair_family_n"] for row in ledger}, {"30"})
            failures = [row for row in ledger if row["analysis_status"] == "NUMERICAL_FAILURE_P_SET_TO_ONE"]
            self.assertEqual(len(failures), 1)
            self.assertEqual(failures[0]["P_PLACO_PLUS"], "1")
            self.assertIn("synthetic integration failure", failures[0]["numerical_error"])
            p_values = [float(row["P_PLACO_PLUS"]) for row in ledger]
            expected_q = bh(p_values)
            for row, wanted in zip(ledger, expected_q, strict=True):
                self.assertTrue(math.isclose(float(row["PLACO_BH_Q"]), wanted, rel_tol=1e-9, abs_tol=1e-12))

            original = paths["staged_ledger"].read_bytes()
            paths["staged_ledger"].write_bytes(original + b"tamper")
            with self.assertRaisesRegex(SystemExit, "ledger differs from its run summary"):
                module.publish_run(root, paths["task"], paths["run_summary"])
            self.assertFalse(paths["canonical_ledger"].exists())
            paths["staged_ledger"].write_bytes(original)

            original_summary = paths["run_summary"].read_bytes()
            changed = read_tsv(paths["staged_ledger"], compressed=True)
            changed[0]["SNP"] = "rs999999999"
            write_gzip_tsv(paths["staged_ledger"], module.LEDGER_FIELDS, changed)
            summary = read_tsv(paths["run_summary"])[0]
            summary["ledger_bytes"] = str(paths["staged_ledger"].stat().st_size)
            summary["ledger_sha256"] = module.sha256(paths["staged_ledger"])
            write_tsv(paths["run_summary"], module.SUMMARY_FIELDS, [summary])
            with self.assertRaisesRegex(SystemExit, "does not exactly preserve its aligned hypothesis family"):
                module.publish_run(root, paths["task"], paths["run_summary"])
            self.assertFalse(paths["canonical_ledger"].exists())
            paths["staged_ledger"].write_bytes(original)
            paths["run_summary"].write_bytes(original_summary)

            paths["staged_provenance"].write_text("occupied\n", encoding="utf-8")
            with self.assertRaisesRegex(SystemExit, "already exists.*overwrite is forbidden"):
                module.publish_run(root, paths["task"], paths["run_summary"])
            self.assertFalse(paths["canonical_ledger"].exists())
            self.assertFalse(paths["canonical_provenance"].exists())
            paths["staged_provenance"].unlink()

            provenance = module.publish_run(root, paths["task"], paths["run_summary"])
            self.assertEqual(provenance["output_rows"], 30)
            self.assertEqual(provenance["complete_family_counts"]["failures"], 1)
            self.assertEqual(provenance["terminal_result_state"], "TESTED_NO_HIT")
            self.assertTrue(paths["canonical_ledger"].is_file())
            self.assertTrue(paths["canonical_provenance"].is_file())
            with self.assertRaisesRegex(SystemExit, "already exists.*overwrite is forbidden"):
                module.publish_run(root, paths["task"], paths["run_summary"])

    def test_failure_fraction_gate_stops_before_scientific_output(self) -> None:
        module = load_materializer()
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            context = make_context(module, root, force_failure=True, maximum_failure_fraction=0.0001)
            module.materialize_pair(context)
            paths = module.materialized_paths(context)
            result = subprocess.run(
                [str(ROOT / ".r-env/bin/Rscript"), str(root / module.RUNNER), str(paths["task"].relative_to(root)),
                 "--root", str(root), "--execute", "--workers", "1", "--stage-only"],
                cwd=root, check=False, text=True, capture_output=True,
            )
            self.assertNotEqual(result.returncode, 0)
            self.assertIn("numerical failure fraction", result.stdout + result.stderr)
            self.assertFalse(any(paths[name].exists() for name in (
                "staged_ledger", "run_summary", "staged_provenance", "canonical_ledger", "canonical_provenance",
            )))

    def test_control_hits_never_enter_primary_ab_headline_count(self) -> None:
        module = load_materializer()
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            context = make_context(module, root, pair_id="CONTROL", force_hits=True)
            module.materialize_pair(context)
            paths = module.materialized_paths(context)
            result = subprocess.run(
                [str(ROOT / ".r-env/bin/Rscript"), str(root / module.RUNNER), str(paths["task"].relative_to(root)),
                 "--root", str(root), "--execute", "--workers", "1", "--stage-only"],
                cwd=root, check=False, text=True, capture_output=True,
            )
            self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
            summary = read_tsv(paths["run_summary"])[0]
            self.assertEqual(summary["primary_headline_count"], "0")
            self.assertEqual(summary["pairwise_gws_count"], "30")
            provenance = module.publish_run(root, paths["task"], paths["run_summary"])
            self.assertEqual(provenance["complete_family_counts"]["primary"], 0)
            self.assertEqual(provenance["complete_family_counts"]["pairwise"], 30)
            self.assertEqual(provenance["terminal_result_state"], "COMPLETE_WITH_HITS")


if __name__ == "__main__":
    unittest.main()
