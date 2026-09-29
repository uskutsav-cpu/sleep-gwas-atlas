import csv
import gzip
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT / "extensions/brain6"))
sys.path.insert(0, str(ROOT / "brain6/scripts"))

from materialize_lava_family_details import canonical_univariate, write_gzip_tsv


def tsv(path, fields, rows):
    with path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def setup_univariate(tmp_path, *, aggregate_status="TESTED", raw_rows=None):
    run = tmp_path / "run"
    checkpoint = run / "loci" / "1" / "worker_output"
    checkpoint.mkdir(parents=True)
    tsv(run / "univariate_results.tsv", ["trait_id", "locus_id", "status", "p", "reason"], [
        {"trait_id": "trait_a", "locus_id": "1", "status": aggregate_status, "p": "0.01", "reason": ""},
        {"trait_id": "trait_b", "locus_id": "1", "status": "NOT_RUN", "p": "NA", "reason": "missing"},
    ])
    tsv(checkpoint / "univariate.tsv",
        ["phen", "h2.obs", "h2.latent", "ascertained", "p", "pair_id", "locus_id"],
        raw_rows if raw_rows is not None else [{
            "phen": "trait_a", "h2.obs": "0.02", "h2.latent": "0.03",
            "ascertained": "FALSE", "p": "0.01", "pair_id": "pair_a", "locus_id": "1",
        }])
    return run


def test_materializes_h2_values_and_explicit_unavailable_se(tmp_path):
    run = setup_univariate(tmp_path)
    out = canonical_univariate(run, [{"LOC": "1"}], ["trait_a", "trait_b"])
    assert out[0]["trait_id"] == "trait_a"
    assert out[0]["local_h2_observed"] == "0.02"
    assert out[0]["local_h2_latent"] == "0.03"
    assert out[0]["p_h2"] == "0.01"
    assert out[0]["se_local_h2"] == "NA"
    assert out[1]["status"] == "NOT_RUN"
    assert out[1]["local_h2_observed"] == "NA"


def test_rejects_disagreement_between_pairwise_univariate_outputs(tmp_path):
    rows = [
        {"phen": "trait_a", "h2.obs": "0.02", "h2.latent": "0.03", "ascertained": "FALSE",
         "p": "0.01", "pair_id": "pair_a", "locus_id": "1"},
        {"phen": "trait_a", "h2.obs": "0.04", "h2.latent": "0.03", "ascertained": "FALSE",
         "p": "0.01", "pair_id": "pair_b", "locus_id": "1"},
    ]
    run = setup_univariate(tmp_path, raw_rows=rows)
    with pytest.raises(ValueError, match="details disagree"):
        canonical_univariate(run, [{"LOC": "1"}], ["trait_a", "trait_b"])


def test_rejects_aggregate_marked_tested_without_receipt_bound_h2(tmp_path):
    run = setup_univariate(tmp_path, aggregate_status="TESTED", raw_rows=[])
    with pytest.raises(ValueError, match="status disagree"):
        canonical_univariate(run, [{"LOC": "1"}], ["trait_a", "trait_b"])


def test_gzip_tsv_is_deterministic_and_immutable(tmp_path):
    path = tmp_path / "details.tsv.gz"
    fields = ["trait", "value"]
    rows = [{"trait": "trait_a", "value": 1.25}]
    write_gzip_tsv(path, fields, rows)
    first = path.read_bytes()
    write_gzip_tsv(path, fields, rows)
    assert path.read_bytes() == first
    with gzip.open(path, "rt", newline="") as stream:
        assert list(csv.DictReader(stream, delimiter="\t")) == [{"trait": "trait_a", "value": "1.25"}]
    with pytest.raises(FileExistsError, match="Refusing to replace"):
        write_gzip_tsv(path, fields, [{"trait": "trait_b", "value": 2}])
