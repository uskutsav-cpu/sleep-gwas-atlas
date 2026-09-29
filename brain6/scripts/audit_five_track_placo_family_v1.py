#!/usr/bin/env python3
"""Source-bound five-row PLACO family audit without changing frozen outputs."""

from __future__ import annotations

import csv
import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
SOURCE_MANIFEST = ROOT / "brain6/manifests/gwas_external_availability.tsv"
FAMILY_LOCK = ROOT / "brain6/config/placo_family_v3/family_lock.json"
RUN_MANIFEST = ROOT / "brain6/manifests/placo_family_v3_run.json"
QC = ROOT / "brain6/results/placo/placo_v3_pair_qc_validation.json"
PAIR_LOCK = ROOT / "extensions/brain6/work/overnight-v03/brain6-pairs.lock.json"
CARD_ROOT = ROOT / "extensions/brain6/work/overnight-v03/source-cards"
PREP = Path("/Volumes/Extreme SSD/brain6-work/preparation-v1")
EXEC = Path("/Volumes/Extreme SSD/brain6-work/placo-family-v3-results")
ARCHIVE = Path("/Volumes/Extreme SSD/Utsav-Research-Archive/Sleep-GWAS/FAILED-PARTIAL-2026-09-01-sleep_gwas_atlas")
OUT = ROOT / "brain6/results/placo/five_track_completeness_v1"
PAIRS = (
    ("insomnia__adhd", "insomnia", "adhd"),
    ("insomnia__mdd", "insomnia", "mdd"),
    ("longsleep__scz", "longsleep", "scz"),
    ("longsleep__bipolar", "longsleep", "bipolar"),
    ("longsleep__parkinson", "longsleep", "parkinson"),
)


def digest(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(8 << 20), b""):
            h.update(block)
    return h.hexdigest()


def tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as stream:
        return list(csv.DictReader(stream, delimiter="\t"))


def main() -> None:
    manifest = {x["trait"]: x for x in tsv(SOURCE_MANIFEST)}
    lock = json.loads(FAMILY_LOCK.read_text())
    run = json.loads(RUN_MANIFEST.read_text())
    qc = json.loads(QC.read_text())
    selection = json.loads(PAIR_LOCK.read_text())
    assert lock["n_selected_tracks"] == len(PAIRS) == selection["n_pairs"] == 5
    assert lock["protected_legacy_pair"] == "insomnia__adhd"
    assert run["status"] == "PAIR_OUTPUTS_COMPLETE_PROTECTED_TRACK_B_PENDING"
    assert qc["full_family_complete"] is False
    assert {x[0] for x in PAIRS} == {x["pair_id"] for x in selection["pairs"]}
    assert digest(FAMILY_LOCK) == qc["family_lock_sha256"] == run["family_lock_sha256"]
    assert digest(Path(lock["placo_source_path"])) == lock["placo_source_sha256"]

    audits = []
    for pair, left, right in PAIRS:
        cards = [json.loads((CARD_ROOT / f"{name}.source.json").read_text()) for name in (left, right)]
        entries = [manifest[name] for name in (left, right)]
        for card, entry in zip(cards, entries):
            assert card["access_permitted"] is True and entry["license"] == "PUBLIC"
            assert card["genome_build"] == "GRCh37" and card["ancestry"] == "EUR"
            assert digest(Path(entry["source_archive_path"])) == entry["source_archive_container_SHA256"]
            assert digest(Path(card["path"])) == card["sha256"] == entry["harmonized_input_SHA256"]
        paired = PREP / f"pair_{pair}" / "pair.tsv.gz"
        prep_receipt = PREP / f"pair_{pair}" / "receipt.json"
        declared = next(x["sha256"] for x in json.loads(prep_receipt.read_text())["outputs"] if x["path"] == "pair.tsv.gz")
        assert digest(paired) == declared
        prep_qc = json.loads((paired.parent / "qc.json").read_text())
        item = {
            "track": "Track B / Pair B" if pair == "insomnia__adhd" else "Brain6 v3",
            "pair_id": pair,
            "source_gwas_identity": ";".join(card["study_id"] for card in cards),
            "source_public_urls": ";".join(card["source_uri"] for card in cards),
            "source_archive_sha256": ";".join(x["source_archive_container_SHA256"] for x in entries),
            "source_access": "PUBLIC_AGGREGATE_GWAS_BOTH",
            "sample_cohort_metadata": ";".join(f"{n}:N={x['sample_size']},cases={x['cases']},controls={x['controls']},cohorts={x['cohorts']}" for n, x in zip((left, right), entries)),
            "genome_build": "GRCh37;GRCh37",
            "harmonized_allele_effect_semantics": ";".join(f"{n}:A1_effect,BETA_{c['effect_scale']},N_{c['n_semantics']}" for n, c in zip((left, right), cards)),
            "harmonized_variant_count": ";".join(f"{n}:{x['variant_count']}" for n, x in zip((left, right), entries)),
            "harmonized_sha256": ";".join(c["sha256"] for c in cards),
            "pair_input_path": str(paired),
            "pair_input_rows": prep_qc["retained_rows"],
            "pair_input_sha256": declared,
            "genomewide_parameter_receipt_path": "",
            "genomewide_parameter_receipt_sha256": "",
            "placo_method": "PLACO_PLUS",
            "placo_implementation": "official PLACO 0.2.0",
            "placo_source_sha256": lock["placo_source_sha256"],
            "candidate_output_path": "",
            "candidate_output_sha256": "",
            "output_rows": "",
            "completed_chunks": "",
            "completeness_status": "",
            "admissibility_status": "",
        }
        if pair == "insomnia__adhd":
            item["candidate_output_path"] = "results/track_b/pleiotropy/results/placo/B.full.tsv.gz"
            item["completeness_status"] = "MISSING_TERMINAL_VALID_PROTECTED_RESULT"
            item["admissibility_status"] = "BLOCKED_ARCHIVED_RESULT_REJECTED_AND_V3_RERUN_FORBIDDEN"
        else:
            parameter = EXEC / pair / "parameters" / "receipt.json"
            output = ROOT / "brain6/results/placo" / pair / "variants.tsv.gz"
            expected = qc["published_output_sha256"][pair]
            assert digest(output) == expected
            assert run["pairs"][pair]["status"] == "COMPLETE"
            assert json.loads((output.parent / "status.json").read_text())["status"] == "PASS"
            item["genomewide_parameter_receipt_path"] = str(parameter)
            item["genomewide_parameter_receipt_sha256"] = digest(parameter)
            item["candidate_output_path"] = str(output)
            item["candidate_output_sha256"] = expected
            item["output_rows"] = run["pairs"][pair]["expected_rows"]
            item["completed_chunks"] = run["pairs"][pair]["completed_chunks"]
            item["completeness_status"] = "PAIR_QC_PASS_COMPLETE"
            item["admissibility_status"] = "ADMITTED_PAIR_RESULT_FULL_FAMILY_INCOMPLETE"
        audits.append(item)

    historical = ARCHIVE / "results/track_b/pleiotropy/results/placo/B.provenance.json"
    old = json.loads(historical.read_text())
    current_contract = ROOT / "results/track_b/pleiotropy/contract.lock.json"
    current_gate = ROOT / "results/track_b/pleiotropy/input_gate.lock.json"
    archived_gate = ARCHIVE / "results/track_b/pleiotropy/input_gate.lock.json"
    assert old["output_sha256"] == digest(ARCHIVE / old["output"])
    assert old["contract_lock_sha256"] == digest(current_contract)
    assert old["input_gate_lock_sha256"] != digest(archived_gate)
    assert old["input_gate_lock_sha256"] != digest(current_gate)
    assert not (ARCHIVE / old["run_summary"]["ledger"]).exists() or old["input_sha256"] != digest(ARCHIVE / old["run_summary"]["ledger"])
    amendment = ARCHIVE / "results/track_b/pleiotropy/continuations/post_lava_terminal_v2/implementation_compatibility_v4/AMENDMENT.md"
    assert "reject shards without successful full-scan monitor evidence" in amendment.read_text()

    OUT.mkdir(parents=True, exist_ok=False)
    table = OUT / "five_track_completeness.tsv"
    with table.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(audits[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(audits)
    report = OUT / "track_b_unblock_requirements.md"
    report.write_text("\n".join([
        "# Brain6 Track B: exact unblock requirements", "",
        "Track B is the legacy **Pair B: insomnia–ADHD** analysis. Its locked sources are the publicly released Jansen 2019 UK Biobank-only insomnia GWAS and Demontis 2023 iPSYCH/deCODE/PGC European ADHD GWAS. Both original aggregates and harmonized GRCh37 files are present and SHA256-verified; their prepared 5,514,402-row pair input is present and verified. The access restriction is on admission of a protected historical result under the frozen Brain6 v3 design, not on obtaining these public GWAS summaries.", "",
        "The frozen Brain6 v3 lock names this pair as `protected_legacy_pair` and the checked-in `plan-placo` and `execute-placo` functions explicitly forbid rerunning it into v3. The four other selected pairs have complete, receipt-bound QC-passing outputs. Their five-track Bonferroni multiplier is already applied, but the family remains incomplete because Track B has no admissible terminal output.", "",
        "## Required legacy output package", "",
        "- `results/track_b/pleiotropy/results/placo/B.full.tsv.gz` and `B.provenance.json`, with full genome-wide rows, per-row P and within-pair BH q, exact allele/coordinate schema, full family counts and QC status. The canonical expected directory is fixed by `config/track_b_pleiotropy_policy.json`.",
        "- Matching `contract.lock.json` and `input_gate.lock.json`; SHA256 references in the sidecar must match the supplied bytes. The current contract lock is `" + digest(current_contract) + "`; the current gate lock is `" + digest(current_gate) + "`.",
        "- The exact aligned input and its checksum, pinned PLACO 0.2.0 source, nuisance parameter checkpoint, runner and materializer source files, complete shard receipts, and successful full-scan process/RAM monitor and terminal exit receipts. The sidecar must bind them and the result output hash. See `result_provenance_required` in the frozen Track B policy.",
        "", "## Why the mounted failed-partial archive is inadmissible", "",
        "The mounted archive's `B.full.tsv.gz` is intact (SHA256 `" + old["output_sha256"] + "`, " + str(old["output_rows"]) + " rows), but its sidecar states `input_gate_lock_sha256=" + old["input_gate_lock_sha256"] + "`; the archived gate file hashes to `" + digest(archived_gate) + "`, and the current gate to `" + digest(current_gate) + "`. The aligned input named in that sidecar (SHA256 `" + old["input_sha256"] + "`) is not supplied with an independently verifiable terminal package. Its runner/materializer hashes differ from the available sources. The archive's V4 amendment explicitly rejects old shards without successful full-scan monitor evidence and requires recomputation under a new implementation fingerprint. The old result cannot be imported into the protected slot.",
        "", "## Acquisition route and next computation", "",
        "The source GWAS summaries themselves are public: [Jansen insomnia release](https://cncr.nl/research/summary_statistics/) and [Demontis ADHD release](https://doi.org/10.6084/m9.figshare.22564390.v1). No login, license, or manual GWAS download is needed for the verified local copies. To complete the frozen protected slot, provide a different terminal-valid legacy Track B result package from the authorized run with the exact files and matching receipts above. If none exists, the frozen v3 family cannot be completed as designed. A newly executed public-data result can be analyzed as a separately named sensitivity, but cannot retroactively become the protected legacy output.",
        "", "Once an admissible protected package exists, verify all hashes and terminal receipts, import without changing frozen outputs, rebuild complete five-track family accounting under the pinned BH and Bonferroni policy, then recompute candidate significance and independent loci before downstream work.", "",
    ]))
    sidecar = {"analysis_id": "brain6_five_track_completeness_v1",
               "created_utc": datetime.now(timezone.utc).isoformat(),
               "track_count": len(audits), "complete_admitted": 4, "blocked": 1,
               "inputs_sha256": {str(x): digest(x) for x in (SOURCE_MANIFEST, FAMILY_LOCK, RUN_MANIFEST, QC, PAIR_LOCK, current_contract, current_gate, archived_gate, historical, amendment, Path(__file__))},
               "outputs_sha256": {x.name: digest(x) for x in (table, report)}}
    (OUT / "provenance.json").write_text(json.dumps(sidecar, indent=2) + "\n")
    print(json.dumps({"complete_admitted": 4, "blocked": 1, "path": str(table)}))


if __name__ == "__main__":
    main()
