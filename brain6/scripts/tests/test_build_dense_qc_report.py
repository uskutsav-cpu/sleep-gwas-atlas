from copy import deepcopy
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT / "extensions/brain6"))
sys.path.insert(0, str(ROOT / "brain6/scripts"))
from build_dense_qc_report import compose_availability_rows


def test_availability_overlay_separates_raw_archive_from_harmonized_dense_input():
    manifest = [
        {"trait": "insomnia", "raw_file_path": "MISSING_FROM_CURRENT_WORKTREE",
         "raw_file_SHA256": "NOT_AVAILABLE", "current_readiness": "BLOCKED"},
        {"trait": "adhd_replication_finngen_r13", "raw_file_path": "NOT_RETAINED",
         "raw_file_SHA256": "KNOWN_HASH", "current_readiness": "ARCHIVED"},
    ]
    original = deepcopy(manifest)
    dense = {"insomnia": {
        "retained_rows": "123", "source_header": "SNP,CHR,BP,A1,A2,BETA,SE,P,N,FRQ",
        "source_path": "/archive/harmonized/insomnia.tsv.gz", "source_sha256": "dense-hash",
    }}
    archives = {"insomnia": {
        "archive_path": "/archive/raw/insomnia.txt.gz", "container_bytes": "456",
        "container_SHA256": "container-hash", "content_transform": "identity",
        "observed_source_SHA256": "source-hash", "status": "VERIFIED_EXACT_ARCHIVE",
    }}

    result = compose_availability_rows(manifest, dense, archives)
    primary, replication = result
    assert manifest == original
    assert primary["raw_file_path"] == "/archive/raw/insomnia.txt.gz"
    assert primary["raw_file_SHA256"] == "container-hash"
    assert primary["source_archive_content_SHA256"] == "source-hash"
    assert primary["harmonized_input_path"] == "/archive/harmonized/insomnia.tsv.gz"
    assert primary["harmonized_input_SHA256"] == "dense-hash"
    assert primary["current_readiness"] == (
        "SOURCE_CONTENT_HASH_VERIFIED; DENSE_INPUT_HASH_VERIFIED; ACQUISITION_DATE_UNKNOWN"
    )
    assert replication["raw_file_path"] == "NOT_RETAINED"
    assert replication["source_archive_status"] == "NA"


def test_availability_overlay_reports_missing_archive_without_fabricating_path():
    manifest = [{"trait": "parkinson", "raw_file_path": "MISSING", "raw_file_SHA256": "NA",
                 "current_readiness": "BLOCKED"}]
    archives = {"parkinson": {
        "archive_path": "/expected/not-found.gz", "container_bytes": "", "container_SHA256": "",
        "content_transform": "gzip_decompress", "observed_source_SHA256": "", "status": "MISSING",
    }}
    result = compose_availability_rows(manifest, {}, archives)[0]
    assert result["raw_file_path"] == "MISSING"
    assert result["current_readiness"] == (
        "SOURCE_ARCHIVE_MISSING; DENSE_INPUT_UNAVAILABLE; ACQUISITION_DATE_UNKNOWN"
    )
