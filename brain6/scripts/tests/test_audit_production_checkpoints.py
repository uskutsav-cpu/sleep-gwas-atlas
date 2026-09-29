import csv
import gzip
import hashlib
import json
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT / "brain6/scripts"))
from audit_external_source_archives import audit_archives
from audit_production_checkpoints import audit_lava, audit_placo_code_lock, audit_placo_collation, audit_placo_pair


def seal_chunk(path: Path, *, scientific_status="PASS"):
    path.mkdir(parents=True)
    payload = path / "result.tsv"
    payload.write_text("one\n", encoding="utf-8")
    data = payload.read_bytes()
    (path / "receipt.json").write_text(json.dumps({
        "stage": "native_job",
        "status": "COMPLETE",
        "scientific_status": scientific_status,
        "outputs": [{"path": "result.tsv", "bytes": len(data),
                     "sha256": hashlib.sha256(data).hexdigest()}],
    }), encoding="utf-8")


def test_placo_audit_counts_verified_contiguous_chunk_prefix(tmp_path):
    seal_chunk(tmp_path / "chunk_000000")
    seal_chunk(tmp_path / "chunk_000001")
    assert audit_placo_pair(tmp_path, "pair_x", 3) == {
        "verified_chunks": 2, "expected_chunks": 3, "complete": False,
        "first_missing_chunk": 2,
    }


def test_placo_audit_rejects_receipt_gaps_and_failed_scientific_status(tmp_path):
    seal_chunk(tmp_path / "chunk_000000")
    seal_chunk(tmp_path / "chunk_000002")
    with pytest.raises(ValueError, match="contiguous prefix"):
        audit_placo_pair(tmp_path, "pair_x", 3)

    failed = tmp_path / "failed"
    seal_chunk(failed / "chunk_000000", scientific_status="FAILED_QC")
    with pytest.raises(ValueError, match="failed its receipt status"):
        audit_placo_pair(failed, "pair_x", 1)


def test_placo_collation_audit_distinguishes_in_progress_and_not_started(tmp_path):
    assert audit_placo_collation(tmp_path)["status"] == "NOT_STARTED"
    (tmp_path / "collated.lock").write_text('{"pid": 123}', encoding="utf-8")
    (tmp_path / ".collated.test.partial").mkdir()
    result = audit_placo_collation(tmp_path)
    assert result == {
        "status": "COLLATION_IN_PROGRESS",
        "receipt_verified": False,
        "lock_present": True,
        "partial_directories": [".collated.test.partial"],
    }


def test_placo_collation_audit_verifies_final_pair_qc(tmp_path):
    final = tmp_path / "collated"
    final.mkdir()
    status = final / "status.json"
    status.write_text('{"status":"PASS"}', encoding="utf-8")
    (final / "receipt.json").write_text(json.dumps({
        "stage": "collate_placo",
        "status": "COMPLETE",
        "outputs": [{"path": "status.json", "bytes": status.stat().st_size,
                     "sha256": hashlib.sha256(status.read_bytes()).hexdigest()}],
    }), encoding="utf-8")
    result = audit_placo_collation(tmp_path)
    assert result == {
        "status": "COMPLETE_QC_PASS",
        "receipt_verified": True,
        "lock_present": False,
        "partial_directories": [],
    }


def test_placo_code_lock_verifies_inventory_and_adapters(tmp_path):
    files = {
        "extension.py": b"frozen extension\n",
        "adapter.R": b"frozen adapter\n",
        "bundled.R": b"frozen adapter\n",
    }
    for name, data in files.items():
        (tmp_path / name).write_bytes(data)
    lock = {
        "extension_code_inventory": [{"path": "extension.py",
                                       "sha256": hashlib.sha256(files["extension.py"]).hexdigest()}],
        "adapter_path": "adapter.R",
        "adapter_sha256": hashlib.sha256(files["adapter.R"]).hexdigest(),
        "bundled_adapter_path": "bundled.R",
        "bundled_adapter_sha256": hashlib.sha256(files["bundled.R"]).hexdigest(),
        "placo_source_path": str(tmp_path / "source.R"),
        "placo_source_sha256": hashlib.sha256(b"upstream source\n").hexdigest(),
    }
    (tmp_path / "source.R").write_bytes(b"upstream source\n")

    assert audit_placo_code_lock(tmp_path, lock) == {"status": "PASS", "files_checked": 4}


def test_placo_code_lock_rejects_modified_or_missing_files(tmp_path):
    source = tmp_path / "source.R"
    source.write_text("changed\n", encoding="utf-8")
    adapter = tmp_path / "adapter.R"
    adapter.write_text("adapter\n", encoding="utf-8")
    bundled = tmp_path / "bundled.R"
    bundled.write_text("bundled\n", encoding="utf-8")
    lock = {
        "extension_code_inventory": [],
        "adapter_path": "adapter.R",
        "adapter_sha256": hashlib.sha256(adapter.read_bytes()).hexdigest(),
        "bundled_adapter_path": "bundled.R",
        "bundled_adapter_sha256": hashlib.sha256(bundled.read_bytes()).hexdigest(),
        "placo_source_path": str(source),
        "placo_source_sha256": "f" * 64,
    }
    with pytest.raises(ValueError, match="hash mismatch"):
        audit_placo_code_lock(tmp_path, lock)
    lock["placo_source_path"] = str(tmp_path / "absent.R")
    with pytest.raises(ValueError, match="file is missing"):
        audit_placo_code_lock(tmp_path, lock)


def test_placo_code_lock_rejects_repository_path_escape(tmp_path):
    outside = tmp_path.parent / "outside.R"
    outside.write_text("outside\n", encoding="utf-8")
    lock = {
        "extension_code_inventory": [{"path": "../outside.R",
                                       "sha256": hashlib.sha256(outside.read_bytes()).hexdigest()}],
        "adapter_path": "unused.R", "adapter_sha256": "0" * 64,
        "bundled_adapter_path": "unused.R", "bundled_adapter_sha256": "0" * 64,
        "placo_source_path": str(outside),
        "placo_source_sha256": hashlib.sha256(outside.read_bytes()).hexdigest(),
    }
    with pytest.raises(ValueError, match="escapes repository"):
        audit_placo_code_lock(tmp_path, lock)


def test_external_archive_audit_binds_size_and_sha_and_keeps_missing_explicit(tmp_path):
    traits = ("insomnia", "longsleep", "adhd", "mdd", "scz", "bipolar", "parkinson")
    archive_root = tmp_path / "archives"
    archive_root.mkdir()
    manifest = tmp_path / "gwas.tsv"
    candidates = tmp_path / "candidates.tsv"
    manifest_rows = []
    candidate_rows = []
    for trait in traits:
        name = f"{trait}.gz"
        data = f"source-{trait}\n".encode()
        if trait == "parkinson":
            with gzip.open(archive_root / name, "wb") as stream:
                stream.write(data)
        else:
            (archive_root / name).write_bytes(data)
        manifest_rows.append({"trait": trait, "study": f"study-{trait}",
                              "source_URL": f"https://example.test/{trait}",
                              "expected_source_archive_bytes": str(len(data)),
                              "expected_source_archive_SHA256": hashlib.sha256(data).hexdigest()})
        candidate_rows.append({"trait": trait, "archive_path": name,
                               "content_transform": "gzip_decompress" if trait == "parkinson" else "identity"})
    for path, rows in ((manifest, manifest_rows), (candidates, candidate_rows)):
        with path.open("w", newline="", encoding="utf-8") as stream:
            writer = csv.DictWriter(stream, fieldnames=rows[0], delimiter="\t", lineterminator="\n")
            writer.writeheader()
            writer.writerows(rows)

    records = audit_archives(manifest, candidates, archive_root)
    statuses = {row["trait"]: row["status"] for row in records}
    assert statuses["insomnia"] == "VERIFIED_EXACT_ARCHIVE"
    assert statuses["parkinson"] == "VERIFIED_DECOMPRESSED_SOURCE"

    (archive_root / "mdd.gz").write_text("altered\n", encoding="utf-8")
    records = audit_archives(manifest, candidates, archive_root)
    assert next(row for row in records if row["trait"] == "mdd")["status"] == "MISMATCH"


def seal_lava_unit(path: Path, locus_id: str, run_id: str, result_rows: list[dict]) -> None:
    path.mkdir(parents=True)
    config = path / "worker_config.json"
    config.write_text('{"frozen":true}\n', encoding="utf-8")
    output_dir = path / "worker_output"
    output_dir.mkdir()
    fields = ["pair_id", "locus_id", "status", "p", "local_rg", "reason"]
    with (output_dir / "pair_results.tsv").open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader(); writer.writerows(result_rows)
    for name in ("univariate.tsv", "status.json"):
        (output_dir / name).write_text("stub\n", encoding="utf-8")
    outputs = []
    for item in sorted(output_dir.iterdir()):
        data = item.read_bytes()
        outputs.append({"path": f"worker_output/{item.name}", "bytes": len(data),
                        "sha256": hashlib.sha256(data).hexdigest()})
    (path / "receipt.json").write_text(json.dumps({
        "locus_id": locus_id, "run_id": run_id, "status": "COMPLETE",
        "worker_config_sha256": hashlib.sha256(config.read_bytes()).hexdigest(),
        "outputs": outputs,
    }), encoding="utf-8")


def test_lava_audit_explains_no_overlap_slots_from_verified_receipts(tmp_path):
    run_id = "run1"
    loci = tmp_path / "loci.txt"
    loci.write_text("LOC CHR START STOP\n10 1 100 200\n11 1 200 300\n", encoding="utf-8")
    rows10 = [
        {"pair_id": "p1", "locus_id": "10", "status": "NO_OVERLAP", "p": "NA", "local_rg": "NA",
         "reason": "FEWER_THAN_MIN_K_SHARED_VARIANTS_IN_SEALED_REFERENCE"},
        {"pair_id": "p2", "locus_id": "10", "status": "UNIVARIATE_UNDERPOWERED", "p": "NA",
         "local_rg": "NA", "reason": "LOCAL_H2_BELOW_THRESHOLD"},
    ]
    rows11 = [
        {"pair_id": "p1", "locus_id": "11", "status": "TESTED", "p": "0.5", "local_rg": "0.1", "reason": ""},
        {"pair_id": "p2", "locus_id": "11", "status": "NO_OVERLAP", "p": "NA", "local_rg": "NA",
         "reason": "FEWER_THAN_MIN_K_SHARED_VARIANTS_IN_SEALED_REFERENCE"},
    ]
    seal_lava_unit(tmp_path / "loci/10", "10", run_id, rows10)
    seal_lava_unit(tmp_path / "loci/11", "11", run_id, rows11)
    result = audit_lava(tmp_path, run_id, loci, {"p1", "p2"})
    assert result["verified_loci"] == 2
    assert result["pair_locus_status_counts"] == {"NO_OVERLAP": 2, "UNIVARIATE_UNDERPOWERED": 1, "TESTED": 1}
    assert result["pair_locus_reason_counts"]["NO_OVERLAP:FEWER_THAN_MIN_K_SHARED_VARIANTS_IN_SEALED_REFERENCE"] == 2
    assert result["no_overlap_unique_loci"] == 2
    assert result["no_overlap_slots_by_pair"] == {"p1": 1, "p2": 1}


def test_lava_parallel_receipt_audit_matches_serial_result(tmp_path):
    run_id = "run1"
    loci = tmp_path / "loci.txt"
    loci.write_text("LOC CHR START STOP\n10 1 100 200\n11 1 200 300\n", encoding="utf-8")
    rows = [
        {"pair_id": "p1", "locus_id": "10", "status": "NO_OVERLAP", "p": "NA",
         "local_rg": "NA", "reason": "FEWER_THAN_MIN_K_SHARED_VARIANTS_IN_SEALED_REFERENCE"},
        {"pair_id": "p2", "locus_id": "10", "status": "UNIVARIATE_UNDERPOWERED", "p": "NA",
         "local_rg": "NA", "reason": "LOCAL_H2_BELOW_THRESHOLD"},
        {"pair_id": "p1", "locus_id": "11", "status": "FAILED", "p": "NA",
         "local_rg": "NA", "reason": "PRESERVED_FAILURE"},
        {"pair_id": "p2", "locus_id": "11", "status": "TESTED", "p": "0.5",
         "local_rg": "0.1", "reason": ""},
    ]
    seal_lava_unit(tmp_path / "loci/10", "10", run_id, rows[:2])
    seal_lava_unit(tmp_path / "loci/11", "11", run_id, rows[2:])
    serial = audit_lava(tmp_path, run_id, loci, {"p1", "p2"}, max_workers=1)
    parallel = audit_lava(tmp_path, run_id, loci, {"p1", "p2"}, max_workers=4)
    assert parallel == serial


def test_lava_parallel_receipt_audit_requires_positive_worker_count(tmp_path):
    loci = tmp_path / "loci.txt"
    loci.write_text("LOC CHR START STOP\n10 1 100 200\n", encoding="utf-8")
    with pytest.raises(ValueError, match="max_workers must be positive"):
        audit_lava(tmp_path, "run1", loci, {"p1", "p2"}, max_workers=0)


def test_lava_audit_lists_failures_and_keeps_incomplete_ceiling_provisional(tmp_path):
    run_id = "run1"
    loci = tmp_path / "loci.txt"
    loci.write_text("LOC CHR START STOP\n10 1 100 200\n11 1 200 300\n12 1 300 400\n", encoding="utf-8")
    rows10 = [
        {"pair_id": "p1", "locus_id": "10", "status": "FAILED", "p": "NA", "local_rg": "NA",
         "reason": "Less than 3 SNPs shared across data sets"},
        {"pair_id": "p2", "locus_id": "10", "status": "UNIVARIATE_UNDERPOWERED", "p": "NA",
         "local_rg": "NA", "reason": "LOCAL_H2_BELOW_THRESHOLD"},
    ]
    rows11 = [
        {"pair_id": "p1", "locus_id": "11", "status": "TESTED", "p": "0.5", "local_rg": "0.1", "reason": ""},
        {"pair_id": "p2", "locus_id": "11", "status": "UNIVARIATE_UNDERPOWERED", "p": "NA",
         "local_rg": "NA", "reason": "LOCAL_H2_BELOW_THRESHOLD"},
    ]
    seal_lava_unit(tmp_path / "loci/10", "10", run_id, rows10)
    seal_lava_unit(tmp_path / "loci/11", "11", run_id, rows11)
    result = audit_lava(tmp_path, run_id, loci, {"p1", "p2"}, maximum_failure_rate=0.01)
    assert result["failed_slots"] == [{
        "locus_id": "10", "pair_id": "p1", "reason": "Less than 3 SNPs shared across data sets"
    }]
    assert result["observed_failed_slot_fraction_of_planned_family"] == pytest.approx(1 / 6)
    assert result["failure_rate_ceiling_status"] == "PROVISIONALLY_OVER_CEILING"


def test_lava_audit_passes_failure_ceiling_only_for_complete_family(tmp_path):
    run_id = "run1"
    loci = tmp_path / "loci.txt"
    loci.write_text("LOC CHR START STOP\n10 1 100 200\n", encoding="utf-8")
    rows10 = [
        {"pair_id": "p1", "locus_id": "10", "status": "FAILED", "p": "NA", "local_rg": "NA",
         "reason": "worker error"},
        {"pair_id": "p2", "locus_id": "10", "status": "TESTED", "p": "0.5", "local_rg": "0.1", "reason": ""},
    ]
    seal_lava_unit(tmp_path / "loci/10", "10", run_id, rows10)
    result = audit_lava(tmp_path, run_id, loci, {"p1", "p2"}, maximum_failure_rate=0.5)
    assert result["complete"] is True
    assert result["failure_rate_ceiling_status"] == "PASS"


def test_lava_audit_rejects_a_drifted_frozen_slot_denominator(tmp_path):
    loci = tmp_path / "loci.txt"
    loci.write_text("LOC CHR START STOP\n10 1 100 200\n", encoding="utf-8")
    with pytest.raises(ValueError, match="planned slot denominator drifted"):
        audit_lava(tmp_path, "run1", loci, {"p1", "p2"}, expected_slot_count=3)


def test_lava_audit_rejects_invalid_failure_ceiling(tmp_path):
    loci = tmp_path / "loci.txt"
    loci.write_text("LOC CHR START STOP\n10 1 100 200\n", encoding="utf-8")
    with pytest.raises(ValueError, match="maximum failure rate"):
        audit_lava(tmp_path, "run1", loci, {"p1", "p2"}, maximum_failure_rate=1.0)
