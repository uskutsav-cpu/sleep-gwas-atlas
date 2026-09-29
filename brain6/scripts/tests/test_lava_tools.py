import json
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT / "extensions" / "brain6"))
sys.path.insert(0, str(ROOT / "brain6" / "scripts"))

from materialize_lava_inputs import load_loci, load_reference_ids
from run_lava_family import (EXPECTED_PAIRS, remove_appledouble_sidecars,
                             require_preflight_success, tsv_bytes, verify_unit,
                             write_bytes_immutable)
from brain6.io import read_json, sha256, write_json


def test_frozen_lava_policies_bind_exact_family():
    family_path = ROOT / "brain6/config/lava_family_v2.json"
    execution_path = ROOT / "brain6/config/lava_execution_v1.json"
    family = json.loads(family_path.read_text())
    execution = json.loads(execution_path.read_text())
    assert sha256(family_path) == execution["family_lock_sha256"]
    assert family["pairs"] == list(EXPECTED_PAIRS)
    assert family["locus_definition"]["n_loci"] == 2495
    assert family["bivariate_family"]["n_slots"] == len(EXPECTED_PAIRS) * 2495
    assert family["univariate_gate"]["n_tests"] == len(family["trait_ids"]) * 2495
    assert execution["worker_protocol"] == "one_frozen_locus_per_R_process; all_five_selected_pairs_per_locus"


def test_locus_reader_verifies_locked_2495_family():
    family = read_json(ROOT / "brain6/config/lava_family_v2.json")
    path = ROOT / family["locus_definition"]["path"]
    loci = load_loci(path, family)
    assert len(loci) == 2495
    assert len({row["LOC"] for row in loci}) == 2495


def test_reference_id_reader_fails_closed_on_duplicate(tmp_path):
    path = tmp_path / "chr1.info"
    path.write_text(
        "SNP\tCHR\tPOS\tA1\tA2\tNOBS\tMISS\tFREQ\tNCORRS\n"
        "rs1\t1\t100\tA\tC\t100\t0\t0.2\t100\n"
        "rs2\t1\t200\tG\tA\t100\t0\t0.3\t100\n"
    )
    assert load_reference_ids(path, 1) == {"rs1", "rs2"}
    path.write_text(path.read_text() + "rs1\t1\t300\tA\tG\t100\t0\t0.4\t100\n")
    with pytest.raises(ValueError, match="duplicate reference SNP"):
        load_reference_ids(path, 1)


def test_immutable_tsv_bytes_and_write_refusal(tmp_path):
    content = tsv_bytes(["pair_id", "locus_id"], [dict(pair_id="p", locus_id="L1")])
    assert content == b"pair_id\tlocus_id\np\tL1\n"
    path = tmp_path / "plan.tsv"
    write_bytes_immutable(path, content)
    write_bytes_immutable(path, content)
    with pytest.raises(FileExistsError, match="Refusing to replace"):
        write_bytes_immutable(path, content + b"extra\n")


def test_external_filesystem_appledouble_sidecars_are_removed(tmp_path):
    (tmp_path / "pair_results.tsv").write_text("pair_id\n")
    sidecar = tmp_path / "._pair_results.tsv"
    sidecar.write_bytes(b"filesystem resource-fork metadata")
    remove_appledouble_sidecars(tmp_path)
    assert (tmp_path / "pair_results.tsv").is_file()
    assert not sidecar.exists()


def test_failed_first_locus_preflight_is_fatal():
    with pytest.raises(RuntimeError, match="preflight failed for locus L1"):
        require_preflight_success("L1", 1)
    require_preflight_success("L1", 0)


def test_lava_unit_checkpoint_checks_identity_and_output_hashes(tmp_path):
    unit = tmp_path / "L1"
    output = unit / "worker_output"
    output.mkdir(parents=True)
    pairs = list(EXPECTED_PAIRS)
    result_path = output / "pair_results.tsv"
    result_path.write_text(
        "pair_id\tlocus_id\tstatus\tp\tlocal_rg\treason\n"
        + "".join(f"{pair}\tL1\tUNIVARIATE_UNDERPOWERED\tNA\tNA\tlow\n" for pair in pairs)
    )
    (output / "univariate.tsv").write_text("phen\tp\n")
    (output / "status.json").write_text('{"status":"PASS"}\n')
    config_path = unit / "worker_config.json"
    config_path.write_text('{"locus_id":"L1"}\n')
    outputs = []
    for item in sorted(output.iterdir()):
        outputs.append({"path": str(item.relative_to(unit)), "bytes": item.stat().st_size,
                        "sha256": sha256(item)})
    write_json(unit / "receipt.json", {"locus_id": "L1", "run_id": "run-a",
        "worker_config_sha256": sha256(config_path), "outputs": outputs, "status": "COMPLETE"})
    assert verify_unit(unit, "L1", set(pairs), "run-a") is not None
    with pytest.raises(ValueError, match="receipt identity mismatch"):
        verify_unit(unit, "L1", set(pairs), "run-b")
    result_path.write_text(result_path.read_text().replace("low", "changed"))
    with pytest.raises(ValueError, match="Corrupt LAVA locus checkpoint"):
        verify_unit(unit, "L1", set(pairs), "run-a")
