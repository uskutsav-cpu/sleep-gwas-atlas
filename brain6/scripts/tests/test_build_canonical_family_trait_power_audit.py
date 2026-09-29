import gzip
import hashlib
import importlib.util
import sys
from pathlib import Path

import pytest

SCRIPT = Path(__file__).resolve().parents[1] / "build_canonical_family_trait_power_audit.py"
sys.path.insert(0, str(SCRIPT.parent))
spec = importlib.util.spec_from_file_location("canonical_trait_power_audit", SCRIPT)
audit = importlib.util.module_from_spec(spec)
assert spec.loader is not None
spec.loader.exec_module(audit)


def write_info(path: Path, chromosome: int, snps: list[str]) -> None:
    with path.open("w", encoding="utf-8") as stream:
        stream.write("SNP\tCHR\tPOS\tA1\tA2\tNOBS\tMISS\tFREQ\tNCORRS\n")
        for index, snp in enumerate(snps, 1):
            stream.write(f"{snp}\t{chromosome}\t{index}\tA\tG\t100\t0\t0.2\t1\n")


def test_exact_id_overlap_handles_unordered_chromosomes(tmp_path, monkeypatch):
    ref_root = tmp_path / "ref"
    ref_root.mkdir()
    write_info(ref_root / "lava-ukb-v1.1_chr1.info", 1, ["rs1", "rs2"])
    write_info(ref_root / "lava-ukb-v1.1_chr2.info", 2, ["rs3", "rs4"])
    monkeypatch.setattr(audit, "ROOT", tmp_path)
    monkeypatch.setattr(audit, "REF_ROOT", ref_root)
    source = tmp_path / "source.tsv.gz"
    with gzip.open(source, "wt", encoding="utf-8", newline="") as stream:
        stream.write("SNP\tCHR\tBP\nrs3\t2\t1\nrs1\t1\t2\nrs9\t2\t3\n")
    digest = hashlib.sha256(source.read_bytes()).hexdigest()
    expected_refs = {
        "ref/lava-ukb-v1.1_chr1.info": audit.sha256(ref_root / "lava-ukb-v1.1_chr1.info"),
        "ref/lava-ukb-v1.1_chr2.info": audit.sha256(ref_root / "lava-ukb-v1.1_chr2.info"),
    }

    total, overlap, ref_hashes = audit.source_reference_overlap(
        source, source.stat().st_size, digest, expected_refs)

    assert (total, overlap) == (3, 2)
    assert ref_hashes == {"1": expected_refs["ref/lava-ukb-v1.1_chr1.info"],
                          "2": expected_refs["ref/lava-ukb-v1.1_chr2.info"]}


def test_exact_id_overlap_fails_on_duplicate_source_ids(tmp_path, monkeypatch):
    ref_root = tmp_path / "ref"
    ref_root.mkdir()
    write_info(ref_root / "lava-ukb-v1.1_chr1.info", 1, ["rs1"])
    monkeypatch.setattr(audit, "ROOT", tmp_path)
    monkeypatch.setattr(audit, "REF_ROOT", ref_root)
    source = tmp_path / "source.tsv.gz"
    with gzip.open(source, "wt", encoding="utf-8", newline="") as stream:
        stream.write("SNP\tCHR\tBP\nrs1\t1\t1\nrs1\t1\t2\n")
    refs = {"ref/lava-ukb-v1.1_chr1.info": audit.sha256(ref_root / "lava-ukb-v1.1_chr1.info")}

    with pytest.raises(ValueError, match="duplicate source SNP ID"):
        audit.source_reference_overlap(source, source.stat().st_size,
                                       hashlib.sha256(source.read_bytes()).hexdigest(), refs)


def test_quantiles_use_linear_interpolation():
    assert audit.q([1.0, 2.0, 3.0, 4.0], 0.25) == 1.75
    assert audit.q([1.0, 2.0, 3.0, 4.0], 0.5) == 2.5
