import csv
import hashlib
import json
import shutil
from pathlib import Path

import pytest

from validate_finngen_r9_replication_candidates import validate


ROOT = Path(__file__).resolve().parents[2]
REPO_ROOT = ROOT.parent
FILES = (
    "qc/replication_sources/finngen_r9_manifest.tsv",
    "qc/finngen_r9_replication_candidates.tsv",
    "qc/finngen_r9_replication_candidates.provenance.json",
)


def _copy_fixture(destination: Path) -> Path:
    for rel in FILES:
        target = destination / "brain6" / rel
        target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(REPO_ROOT / "brain6" / rel, target)
    return destination


def test_r9_candidates_validate_as_metadata_only(tmp_path):
    result = validate(REPO_ROOT)
    assert result == {
        "status": "PASS_R9_METADATA_ONLY_CANDIDATES",
        "candidate_endpoints": 4,
        "summary_statistics_acquired": 0,
        "replication_estimates_admitted": 0,
        "official_access_form_required": True,
        "manifest_sha256": "13f3866a224e35471abc0290cba5ed151aad7a4e270ae70b3eacf3a92cf124a2",
    }


@pytest.mark.parametrize("corruption", ["candidate_bytes", "provenance_hash", "missing_candidate", "head_url"])
def test_r9_candidate_integrity_errors_fail_closed(tmp_path, corruption):
    root = _copy_fixture(tmp_path)
    candidate = root / "brain6/qc/finngen_r9_replication_candidates.tsv"
    provenance_path = root / "brain6/qc/finngen_r9_replication_candidates.provenance.json"
    provenance = json.loads(provenance_path.read_text())
    if corruption == "candidate_bytes":
        candidate.write_text(candidate.read_text() + "\n", encoding="utf-8")
    elif corruption == "provenance_hash":
        provenance["candidate_table"]["sha256"] = "0" * 64
        provenance_path.write_text(json.dumps(provenance), encoding="utf-8")
    elif corruption == "head_url":
        provenance["endpoint_head_metadata"]["F5_BIPO"]["url"] = "https://example.invalid/unrelated.gz"
        provenance_path.write_text(json.dumps(provenance), encoding="utf-8")
    else:
        with candidate.open(newline="", encoding="utf-8") as stream:
            rows = list(csv.DictReader(stream, delimiter="\t"))[:-1]
        with candidate.open("w", newline="", encoding="utf-8") as stream:
            writer = csv.DictWriter(stream, fieldnames=rows[0].keys(), delimiter="\t", lineterminator="\n")
            writer.writeheader()
            writer.writerows(rows)
        # Preserve the old receipt so the row-count and hash checks both reject the mutation.
    with pytest.raises((ValueError, KeyError)):
        validate(root)


def test_r9_manifest_endpoint_counts_are_bound_to_official_manifest(tmp_path):
    root = _copy_fixture(tmp_path)
    manifest = root / "brain6/qc/replication_sources/finngen_r9_manifest.tsv"
    manifest.write_text(manifest.read_text().replace("43280\t329192", "43281\t329192", 1), encoding="utf-8")
    with pytest.raises(ValueError, match="manifest does not match"):
        validate(root)


def test_r9_manifest_cannot_be_replaced_by_coordinated_sidecar_update(tmp_path):
    root = _copy_fixture(tmp_path)
    manifest = root / "brain6/qc/replication_sources/finngen_r9_manifest.tsv"
    provenance_path = root / "brain6/qc/finngen_r9_replication_candidates.provenance.json"
    manifest.write_text(manifest.read_text().replace("43280\t329192", "43281\t329192", 1), encoding="utf-8")
    provenance = json.loads(provenance_path.read_text())
    data = manifest.read_bytes()
    provenance["sources"]["manifest_bytes"] = len(data)
    provenance["sources"]["manifest_sha256"] = hashlib.sha256(data).hexdigest()
    provenance_path.write_text(json.dumps(provenance), encoding="utf-8")
    with pytest.raises(ValueError, match="frozen, independently verified"):
        validate(root)
