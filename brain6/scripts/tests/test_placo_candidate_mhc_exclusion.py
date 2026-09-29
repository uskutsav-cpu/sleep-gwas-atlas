import importlib.util
import json
from pathlib import Path
from tempfile import TemporaryDirectory

ROOT = Path(__file__).resolve().parents[3]
SCRIPT = ROOT / "brain6/scripts/analyze_placo_candidate_mhc_exclusion.py"
spec = importlib.util.spec_from_file_location("placo_candidate_mhc_exclusion", SCRIPT)
mhc = importlib.util.module_from_spec(spec)
assert spec.loader is not None
spec.loader.exec_module(mhc)


def test_mhc_interval_uses_chromosome_six_and_half_open_coordinates():
    assert mhc.is_mhc("6", 25_000_000, chr6="6", start=25_000_000, end=34_000_000)
    assert mhc.is_mhc("chr6", 33_999_999, chr6="6", start=25_000_000, end=34_000_000)
    assert not mhc.is_mhc("6", 24_999_999, chr6="6", start=25_000_000, end=34_000_000)
    assert not mhc.is_mhc("6", 34_000_000, chr6="6", start=25_000_000, end=34_000_000)
    assert not mhc.is_mhc("5", 30_000_000, chr6="6", start=25_000_000, end=34_000_000)


def test_current_partial_placo_candidate_sensitivity_is_checksum_bound(monkeypatch):
    with TemporaryDirectory(dir=ROOT) as temporary_directory:
        temporary_root = Path(temporary_directory)
        monkeypatch.setattr(mhc, "OUT", temporary_root / "mhc.tsv")
        monkeypatch.setattr(mhc, "PROVENANCE", temporary_root / "mhc.provenance.json")
        result = mhc.analyze()
        rows = mhc.read_tsv(mhc.OUT)
        provenance = json.loads(mhc.PROVENANCE.read_text(encoding="utf-8"))
    assert result["status"] == "PASS_PARTIAL_FAMILY_DESCRIPTIVE_COUNTS"
    assert result["family_complete"] is False
    assert rows[-1]["pair_id"] == "ALL_PARTIAL_PAIRS"
    assert rows[-1]["candidate_variants_all"] == "2297"
    assert rows[-1]["candidate_variants_inside_mhc"] == "0"
    assert rows[-1]["candidate_variants_outside_mhc"] == "2273"
    assert rows[-1]["candidate_variants_unpositioned_or_mismatched"] == "24"
    assert rows[-1]["candidate_intervals_overlapping_mhc"] == "0"
    assert provenance["output_path"].endswith("/mhc.tsv")
