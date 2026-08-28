#!/usr/bin/env python3
"""Build the complete real-data trait and trait-pair core of the evidence atlas."""
from __future__ import annotations

import argparse
import csv
import hashlib
import io
import json
import math
from pathlib import Path


TRAIT_FIELDS = [
    "trait_id", "label", "domain", "type", "phenotype_definition", "source_id",
    "dataset_version", "pmid", "doi", "ancestry", "build", "ncase", "ncontrol",
    "n_total", "population_prevalence", "source_status", "harmonized_path",
    "harmonized_sha256", "harmonized_rows", "munged_path", "munged_sha256",
    "h2_scale", "h2", "h2_se", "h2_z", "ldsc_intercept", "ldsc_ratio",
    "phase1_verdict", "phase1_qc_reason",
]
PAIR_FIELDS = [
    "pair_id", "sleep_trait", "non_sleep_trait", "analysis_tier",
    "interpretation_status", "global_rg", "global_rg_se", "global_rg_z",
    "global_rg_p", "global_rg_fdr_all_396", "global_rg_fdr_primary_372",
    "global_rg_primary_significant", "effect_direction", "ldsc_input_log",
]


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def table_text(fields: list[str], rows: list[dict[str, object]]) -> str:
    output = io.StringIO(newline="")
    writer = csv.DictWriter(output, fieldnames=fields, delimiter="\t", lineterminator="\n")
    writer.writeheader()
    writer.writerows(rows)
    return output.getvalue()


def atomic_text(path: Path, payload: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(payload, encoding="utf-8")
    temporary.replace(path)


def qc_value(path: Path, key: str) -> str:
    matches = []
    for line in path.read_text(encoding="utf-8").splitlines():
        fields = line.split("\t")
        if len(fields) == 2 and fields[0] == key:
            matches.append(fields[1])
    if len(matches) != 1:
        raise SystemExit(f"ERROR: expected one {key} value in {path}")
    return matches[0]


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", default=".")
    parser.add_argument("--trait-out", default="results/atlas/traits.tsv")
    parser.add_argument("--pair-out", default="results/atlas/trait_pairs.tsv")
    parser.add_argument("--provenance-out", default="results/atlas/core.provenance.json")
    parser.add_argument("--validate-only", action="store_true")
    args = parser.parse_args()
    root = Path(args.root).resolve()
    panel_path = root / "config/analysis_panel.tsv"
    panel_lock_path = root / "config/analysis_panel.lock.json"
    h2_path = root / "results/tables/h2_summary.tsv"
    rg_path = root / "results/tables/rg_matrix.tsv"
    policy_path = root / "config/downstream_analysis_policy.json"
    panel = read_tsv(panel_path)
    panel_lock = json.loads(panel_lock_path.read_text(encoding="utf-8"))
    policy = json.loads(policy_path.read_text(encoding="utf-8"))
    trait_ids = [row["trait_id"] for row in panel]
    ordered_hash = hashlib.sha256("".join(f"{trait}\n" for trait in trait_ids).encode()).hexdigest()
    if (
        len(panel) != policy["expected_traits"]
        or len(set(trait_ids)) != len(trait_ids)
        or ordered_hash != panel_lock["ordered_trait_ids_sha256"]
    ):
        raise SystemExit("ERROR: atlas core does not match the locked panel")
    h2_rows = read_tsv(h2_path)
    h2 = {row["trait"]: row for row in h2_rows}
    if len(h2_rows) != 45 or set(h2) != set(trait_ids):
        raise SystemExit("ERROR: h2 summary is not the exact 45-trait family")
    rg_rows = read_tsv(rg_path)
    sleeps = [row["trait_id"] for row in panel if row["domain"] == "sleep"]
    non_sleep = [row["trait_id"] for row in panel if row["domain"] != "sleep"]
    expected_pairs = [(sleep, other) for sleep in sleeps for other in non_sleep]
    observed_pairs = [(row["sleep_trait"], row["disease_trait"]) for row in rg_rows]
    if len(rg_rows) != 396 or observed_pairs != expected_pairs:
        raise SystemExit("ERROR: genetic-correlation table is not in exact locked pair order")

    traits = []
    for row in panel:
        trait = row["trait_id"]
        harmonized = root / "data/harmonized" / f"{trait}.harmonized.tsv.gz"
        qc = root / "data/harmonized" / f"{trait}.qc.txt"
        munged = root / "data/munged" / f"{trait}.sumstats.gz"
        if not harmonized.is_file() or not qc.is_file() or not munged.is_file():
            raise SystemExit(f"ERROR: canonical data are absent for {trait}")
        estimate = h2[trait]
        traits.append({
            "trait_id": trait, "label": row["label"], "domain": row["domain"],
            "type": row["type"], "phenotype_definition": row["phenotype_definition"],
            "source_id": row["source_id"], "dataset_version": row["dataset_version"],
            "pmid": row["pmid"], "doi": row["doi"], "ancestry": row["ancestry"],
            "build": row["build"], "ncase": row["ncase"], "ncontrol": row["ncontrol"],
            "n_total": row["n_total"], "population_prevalence": row["pop_prev"],
            "source_status": row["source_status"],
            "harmonized_path": str(harmonized.relative_to(root)),
            "harmonized_sha256": sha256(harmonized), "harmonized_rows": qc_value(qc, "rows_out"),
            "munged_path": str(munged.relative_to(root)), "munged_sha256": sha256(munged),
            "h2_scale": estimate["scale"], "h2": estimate["h2"], "h2_se": estimate["se"],
            "h2_z": estimate["z"], "ldsc_intercept": estimate["intercept"],
            "ldsc_ratio": estimate["ratio"], "phase1_verdict": estimate["verdict"],
            "phase1_qc_reason": estimate["qc_reason"],
        })
    pairs = []
    for row in rg_rows:
        value = float(row["rg"])
        if not math.isfinite(value):
            raise SystemExit(f"ERROR: invalid rg for {row['sleep_trait']} x {row['disease_trait']}")
        primary_fdr = float(row["fdr_primary_phase1"]) if row["fdr_primary_phase1"] not in {"", "NA"} else math.nan
        significant = row["analysis_tier"] == "PRIMARY_PHASE1" and math.isfinite(primary_fdr) and primary_fdr <= 0.05
        pairs.append({
            "pair_id": f"{row['sleep_trait']}__{row['disease_trait']}",
            "sleep_trait": row["sleep_trait"], "non_sleep_trait": row["disease_trait"],
            "analysis_tier": row["analysis_tier"], "interpretation_status": row["interpretation_status"],
            "global_rg": row["rg"], "global_rg_se": row["se"], "global_rg_z": row["z"],
            "global_rg_p": row["p"], "global_rg_fdr_all_396": row["fdr"],
            "global_rg_fdr_primary_372": row["fdr_primary_phase1"] or "NA",
            "global_rg_primary_significant": str(significant).upper(),
            "effect_direction": "POSITIVE" if value > 0 else ("NEGATIVE" if value < 0 else "ZERO"),
            "ldsc_input_log": row["input_log"],
        })
    payloads = {
        root / args.trait_out: table_text(TRAIT_FIELDS, traits),
        root / args.pair_out: table_text(PAIR_FIELDS, pairs),
    }
    provenance = {
        "analysis_id": policy["analysis_id"], "panel_sha256": sha256(panel_path),
        "panel_lock_sha256": sha256(panel_lock_path), "h2_summary_sha256": sha256(h2_path),
        "rg_matrix_sha256": sha256(rg_path), "downstream_policy_sha256": sha256(policy_path),
        "trait_count": len(traits), "pair_count": len(pairs),
        "primary_pair_count": sum(row["analysis_tier"] == "PRIMARY_PHASE1" for row in pairs),
        "sensitivity_pair_count": sum(row["analysis_tier"] != "PRIMARY_PHASE1" for row in pairs),
        "primary_fdr_significant_pair_count": sum(row["global_rg_primary_significant"] == "TRUE" for row in pairs),
        "outputs": {str(path.relative_to(root)): hashlib.sha256(payload.encode()).hexdigest() for path, payload in payloads.items()},
        "status": "CORE_ATLAS_COMPLETE_DOWNSTREAM_LAYERS_PENDING",
    }
    provenance_text = json.dumps(provenance, indent=2, sort_keys=True) + "\n"
    provenance_path = root / args.provenance_out
    if args.validate_only:
        for path, payload in payloads.items():
            if not path.is_file() or path.read_text(encoding="utf-8") != payload:
                raise SystemExit(f"ERROR: atlas core output drifted: {path}")
        if not provenance_path.is_file() or provenance_path.read_text(encoding="utf-8") != provenance_text:
            raise SystemExit("ERROR: atlas core provenance drifted")
    else:
        for path, payload in payloads.items():
            atomic_text(path, payload)
        atomic_text(provenance_path, provenance_text)
    print(f"Atlas core complete: {len(traits)} traits; {len(pairs)} trait pairs")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
