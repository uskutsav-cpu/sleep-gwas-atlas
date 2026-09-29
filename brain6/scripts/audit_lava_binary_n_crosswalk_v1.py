#!/usr/bin/env python3
"""Receipt-bound, read-only audit of the binary N supplied to canonical LAVA.

The atlas harmonizer's N is an LDSC effective N for binary traits. This audit
checks the actual one-locus LAVA inputs and source headers without altering any
historical artifact or inferring that a different N would pass a LAVA gate.
"""
from __future__ import annotations

import csv
import gzip
import hashlib
import json
import math
import zipfile
from datetime import datetime, timezone
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "brain6/results/lava_binary_n_crosswalk_v1"
TRAITS = ("adhd", "bipolar", "insomnia", "longsleep", "mdd", "parkinson", "scz")
SOURCE_N_FIELDS = {
    "adhd": ("Nca", "Nco"),
    "bipolar": ("NCAS", "NCON", "NEFFDIV2"),
    "insomnia": ("N",),
    "longsleep": (),
    "mdd": (),
    "parkinson": ("N_cases", "N_controls"),
    "scz": ("NCAS", "NCON", "NEFF"),
}
INPUT_ROOT = Path("/Volumes/Extreme SSD/brain6-work/preparation-v1/lava-inputs-v1")
CANONICAL_DECISION = ROOT / ("work/lava-canonical-v3-production/"
    "d730debf45266d298401564f3260bdecb14739d1c1f1835a5aebd615c83fa60b/"
    "canonical_family_decision.json")
EXPECTED_CANONICAL_SHA = "3569ce81d33479321f9b3adec7bd09c9b54a10a01f1fa4abd15888d9d2e66b22"


def sha(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for block in iter(lambda: f.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as f:
        return list(csv.DictReader(f, delimiter="\t"))


def first_n(path: Path) -> float:
    with gzip.open(path, "rt", encoding="utf-8") as f:
        reader = csv.DictReader(f, delimiter="\t")
        assert reader.fieldnames and "N" in reader.fieldnames
        row = next(reader)
        return float(row["N"])


def source_header(path: Path) -> set[str]:
    if path.suffix == ".zip":
        with zipfile.ZipFile(path) as z:
            members = [n for n in z.namelist() if not n.startswith("__MACOSX/")]
            assert len(members) == 1, path
            with z.open(members[0]) as f:
                for raw in f:
                    line = raw.decode("utf-8", "replace").strip()
                    if line and not line.startswith("##"):
                        return set(line.lstrip("#").split())
    else:
        opener = gzip.open if path.suffix == ".gz" else open
        with opener(path, "rt", encoding="utf-8", errors="replace") as f:
            for line in f:
                line = line.strip()
                if line and not line.startswith("##"):
                    return set(line.lstrip("#").split())
    raise ValueError(f"No header in {path}")


def main() -> None:
    dense_audit = ROOT / "brain6/manifests/locked_dense_input_audit.tsv"
    external = ROOT / "brain6/manifests/gwas_external_availability.tsv"
    failure = ROOT / "brain6/results/lava_confirmatory_rescue_plan_v1/trait_failure_matrix.tsv"
    harmonizer = ROOT / "scripts/01_harmonize.py"
    runner = ROOT / "brain6/scripts/run_lava_canonical_batch_v3.R"
    lava_binary = ROOT / "brain6/results/lava_confirmatory_protocol_v1/source_docs/LAVA_binary_processing_e729a245.R"
    lava_input = ROOT / "brain6/results/lava_confirmatory_protocol_v1/source_docs/LAVA_input_processing_e729a245.R"
    info = INPUT_ROOT / "locus_1/input_info.tsv"
    if sha(CANONICAL_DECISION) != EXPECTED_CANONICAL_SHA:
        raise ValueError("Canonical decision hash drifted")
    dense = {r["trait_id"]: r for r in tsv(dense_audit)}
    availability = {r["trait"]: r for r in tsv(external)}
    reasons = {r["trait_id"]: r for r in tsv(failure)}
    counts = {r["phenotype"]: r for r in tsv(info)}
    assert set(TRAITS).issubset(dense) and set(TRAITS).issubset(availability)
    assert set(TRAITS).issubset(reasons) and set(TRAITS) == set(counts)
    code = harmonizer.read_text()
    if not all(s in code for s in ('n_mode = "constant N_eff derived from config ncase/ncontrol"',
                                  'n_mode = "per-SNP N_eff from source"',
                                  'n_mode = "per-SNP N_eff derived from NCASE/NCONTROL"')):
        raise ValueError("Harmonizer's effective-N branches changed")
    if 'info$prop_cases <- info$cases / info$N' not in runner.read_text():
        raise ValueError("LAVA runner's case-fraction assignment changed")
    if 'N.case = N * case.prop' not in lava_binary.read_text():
        raise ValueError("Pinned LAVA binary model changed")
    rows: list[dict[str, str]] = []
    bound: dict[str, str] = {}
    for trait in TRAITS:
        card_path = Path(dense[trait]["source_card"])
        card = json.loads(card_path.read_text())
        archive = Path(availability[trait]["source_archive_path"])
        dense_path = Path(dense[trait]["dense_path"])
        shard = INPUT_ROOT / f"locus_1/{trait}.sumstats.tsv.gz"
        assert card["trait_id"] == trait and card["path"] == str(dense_path)
        assert dense[trait]["status"] == "VERIFIED"
        assert card["sha256"] == dense[trait]["expected_sha256"] == dense[trait]["observed_sha256"]
        assert sha(card_path) == dense[trait]["source_card_sha256"]
        assert dense_path.stat().st_size == int(dense[trait]["dense_bytes"])
        assert archive.is_file() and archive.stat().st_size > 0
        assert availability[trait]["source_archive_content_SHA256"]
        header = source_header(archive)
        fields = SOURCE_N_FIELDS[trait]
        assert set(fields).issubset(header), (trait, fields, header)
        n = first_n(shard)
        dense_n = first_n(dense_path)
        cases = int(counts[trait]["cases"])
        controls = int(counts[trait]["controls"])
        total = cases + controls
        pooled_neff = 4 * cases * controls / total
        fraction = cases / total
        assert math.isfinite(n) and n > 0 and math.isfinite(dense_n) and dense_n > 0
        rows.append({
            "trait_id": trait,
            "canonical_NOT_RUN": reasons[trait]["NOT_RUN"],
            "source_card_n_semantics": card["n_semantics"],
            "source_sample_fields": ";".join(fields) if fields else "NONE",
            "source_study_id": card["study_id"],
            "input_info_cases": str(cases),
            "input_info_controls": str(controls),
            "input_info_total": str(total),
            "input_info_case_fraction": f"{fraction:.12g}",
            "pooled_count_neff": f"{pooled_neff:.12g}",
            "dense_first_n": f"{dense_n:.12g}",
            "locus_1_first_n": f"{n:.12g}",
            "locus_1_n_over_total": f"{n / total:.12g}",
            "modeled_cases_at_locus_1_first_n": f"{n * fraction:.12g}",
            "source_archive_path": str(archive),
            "source_archive_content_sha256_from_manifest": availability[trait]["source_archive_content_SHA256"],
            "dense_sha256_from_locked_audit": card["sha256"],
            "locus_1_sumstats_sha256": sha(shard),
            "interpretation": "EFFECTIVE_N_WITH_OBSERVED_AGGREGATE_CASE_FRACTION_REQUIRES_REVIEW",
        })
        bound[f"source_card_{trait}"] = sha(card_path)
        bound[f"locus_1_sumstats_{trait}"] = rows[-1]["locus_1_sumstats_sha256"]
    OUT.mkdir(parents=True, exist_ok=True)
    table = OUT / "binary_n_crosswalk.tsv"
    with table.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)
    report = OUT / "binary_n_interpretation.md"
    report.write_text(
        "# Canonical LAVA binary-N crosswalk — read-only audit\n\n"
        "The seven canonical trait inputs came from `scripts/01_harmonize.py`, whose binary-trait branches write conventional effective `N` for LDSC. The pinned LAVA 0.1.5 binary path instead consumes that per-SNP `N` with a trait-wide observed case fraction from `input_info.tsv`, setting reconstructed cases to `N * case.prop`. The first locus's exact files and all seven source-card hashes are bound in `provenance.json`; the locked dense-input audit certifies the full harmonized file hashes.\n\n"
        "| Trait | NOT_RUN | source fields relevant to N | study cases / controls | first locus N | N / total | modeled cases |\n"
        "|---|---:|---|---:|---:|---:|---:|\n" +
        "\n".join(f"| {r['trait_id']} | {r['canonical_NOT_RUN']} | {r['source_sample_fields']} | {r['input_info_cases']} / {r['input_info_controls']} | {float(r['locus_1_first_n']):,.1f} | {float(r['locus_1_n_over_total']):.3f} | {float(r['modeled_cases_at_locus_1_first_n']):,.1f} |" for r in rows) +
        "\n\nFive source cards label the harmonized `N` as `total`; inspection of their inputs and the harmonizer shows effective `N` instead. The two cards labeled `effective` are internally consistent but still coupled to an observed aggregate case fraction in LAVA. Source archives expose variant-level case/control counts for ADHD, bipolar disorder, Parkinson disease, and schizophrenia, and variant-level total `N` for insomnia. The Dashti long-sleep and Howard MDD releases inspected here lack variant-level sample counts.\n\n"
        "This is a **definite input-model semantic mismatch**, not a proven explanation for any particular low-local-h² `NOT_RUN` cell. Association Z, LD, and true local signal may dominate testability. No canonical output, input, or frozen threshold has been altered. A separate source-specific method decision and outcome-blinded pilot lock must precede any repaired run; do not substitute total `N` for a heterogeneous meta-analysis solely to increase testability.\n",
        encoding="utf-8",
    )
    provenance = {
        "schema_version": 1,
        "audit_id": "brain6_lava_binary_n_crosswalk_v1",
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "read_only_sources": {str(p): sha(p) for p in (dense_audit, external, failure,
            harmonizer, runner, lava_binary, lava_input, info, CANONICAL_DECISION)},
        "bound_sample_sha256": bound,
        "output_sha256": {str(table.relative_to(ROOT)): sha(table),
                          str(report.relative_to(ROOT)): sha(report)},
        "canonical_modified": False,
        "source_replacement_admitted": False,
        "new_lava_run_started": False,
    }
    (OUT / "provenance.json").write_text(json.dumps(provenance, indent=2, sort_keys=True) + "\n")
    print(json.dumps({"traits": len(rows), "not_run": sum(int(r["canonical_NOT_RUN"]) for r in rows),
                      "table_sha256": sha(table), "report_sha256": sha(report)}, sort_keys=True))


if __name__ == "__main__":
    main()
