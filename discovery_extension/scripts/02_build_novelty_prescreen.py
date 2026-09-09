#!/usr/bin/env python3
"""Build a source-grounded prior sleep-screen inventory and annotate Pan-UKB traits."""

from __future__ import annotations

import argparse
import csv
import glob
import hashlib
import json
import re
import unicodedata
from pathlib import Path
from typing import Iterable


def fail(message: str) -> None:
    raise SystemExit(f"ERROR: {message}")


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def normalize(value: str) -> str:
    value = unicodedata.normalize("NFKD", value or "").encode("ascii", "ignore").decode()
    return " ".join(re.sub(r"[^a-z0-9]+", " ", value.lower()).split())


def clean_key(value: str) -> str:
    value = (value or "").strip()
    return "" if value.upper() == "NA" else value


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def read_embedded_header(path: Path, required: str) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.reader(handle, delimiter="\t"))
    header_index = next((i for i, row in enumerate(rows) if required in row), None)
    if header_index is None:
        fail(f"could not find header field {required!r} in {path}")
    header = [cell.strip() for cell in rows[header_index]]
    output = []
    for row in rows[header_index + 1 :]:
        padded = row + [""] * (len(header) - len(row))
        record = dict(zip(header, padded))
        if any(value.strip() for value in record.values()):
            output.append(record)
    return output


def write_tsv(path: Path, fields: list[str], rows: Iterable[dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=fields, lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--universe",
        type=Path,
        default=Path("discovery_extension/results/panukbb_eur_eligible_universe.tsv"),
    )
    parser.add_argument(
        "--prior-dir",
        type=Path,
        default=Path("discovery_extension/provenance/prior_screens"),
    )
    parser.add_argument(
        "--annotated-out",
        type=Path,
        default=Path("discovery_extension/results/panukbb_novelty_prescreen.tsv"),
    )
    parser.add_argument(
        "--inventory-out",
        type=Path,
        default=Path("discovery_extension/results/prior_sleep_screen_inventory.tsv"),
    )
    parser.add_argument(
        "--registry-out",
        type=Path,
        default=Path("discovery_extension/provenance/prior_screens/prior_screen_registry.tsv"),
    )
    parser.add_argument(
        "--run-metadata",
        type=Path,
        default=Path("discovery_extension/provenance/prior_screens/novelty_prescreen_build.json"),
    )
    args = parser.parse_args()

    extracted = args.prior_dir / "extracted"
    goodman_path = extracted / "goodman_2025_SD24_ldsc.tsv"
    dashti_path = extracted / "dashti_2019_supplementary_data_18.tsv"
    lizhao_path = extracted / "li_zhao_2020_s3_table.tsv"
    morrison_paths = sorted(
        Path(path)
        for path in glob.glob(str(extracted / "morrison_2024_S1[0-5]_*_LDSC.tsv"))
    )
    required_paths = [args.universe, goodman_path, dashti_path, lizhao_path, *morrison_paths]
    for path in required_paths:
        if not path.is_file():
            fail(f"required source/extract is missing: {path}")
    if len(morrison_paths) != 6:
        fail(f"expected six Morrison LDSC sheets, found {len(morrison_paths)}")

    universe = read_tsv(args.universe)
    goodman = read_embedded_header(goodman_path, "phenocode")
    dashti = read_embedded_header(dashti_path, "Trait")
    lizhao_rows = read_embedded_header(lizhao_path, "Category Details")
    morrison_by_sleep: dict[str, list[dict[str, str]]] = {}
    morrison_path_by_sleep: dict[str, Path] = {}
    for path in morrison_paths:
        rows = read_tsv(path)
        if not rows or "phenotype" not in rows[0] or "sleep" not in rows[0]:
            fail(f"unexpected Morrison schema: {path}")
        sleep_names = {row["sleep"] for row in rows if row["sleep"]}
        if len(sleep_names) != 1:
            fail(f"expected one sleep construct in {path}, found {sorted(sleep_names)}")
        sleep_name = next(iter(sleep_names))
        morrison_by_sleep[sleep_name] = rows
        morrison_path_by_sleep[sleep_name] = path

    # Li & Zhao's first three columns are Category, Source, and external Trait.
    if not lizhao_rows or "Category Details" not in lizhao_rows[0]:
        fail(f"unexpected Li & Zhao S3 extract schema: {lizhao_path}")

    inventory: list[dict[str, str]] = []
    for row in goodman:
        trait = row.get("description", "").strip()
        if trait:
            inventory.append(
                {
                    "screen_id": "goodman_2025",
                    "sleep_construct": "six_composite_sleep_health_scores",
                    "external_trait": trait,
                    "external_trait_normalized": normalize(trait),
                    "phenocode": clean_key(row.get("phenocode", "")),
                    "coding": clean_key(row.get("coding", "")),
                    "source_extract": str(goodman_path),
                }
            )
    for sleep, rows in sorted(morrison_by_sleep.items()):
        for row in rows:
            trait = row.get("phenotype", "").strip()
            if trait:
                inventory.append(
                    {
                        "screen_id": "morrison_2024",
                        "sleep_construct": sleep,
                        "external_trait": trait,
                        "external_trait_normalized": normalize(trait),
                        "phenocode": "",
                        "coding": "",
                        "source_extract": str(morrison_path_by_sleep[sleep]),
                    }
                )
    for row in dashti:
        trait = row.get("Trait", "").strip()
        if trait:
            inventory.append(
                {
                    "screen_id": "dashti_2019",
                    "sleep_construct": "sleep_duration_short_sleep_long_sleep",
                    "external_trait": trait,
                    "external_trait_normalized": normalize(trait),
                    "phenocode": "",
                    "coding": "",
                    "source_extract": str(dashti_path),
                }
            )
    for row in lizhao_rows:
        trait = row.get("Category Details", "").strip()
        if trait:
            inventory.append(
                {
                    "screen_id": "li_zhao_2020",
                    "sleep_construct": "12_wearable_sleep_and_circadian_traits",
                    "external_trait": trait,
                    "external_trait_normalized": normalize(trait),
                    "phenocode": "",
                    "coding": "",
                    "source_extract": str(lizhao_path),
                }
            )

    write_tsv(
        args.inventory_out,
        [
            "screen_id",
            "sleep_construct",
            "external_trait",
            "external_trait_normalized",
            "phenocode",
            "coding",
            "source_extract",
        ],
        inventory,
    )

    registry = [
        {
            "screen_id": "morrison_2024",
            "citation": "Morrison et al. Multivariate genome-wide association study of sleep health demonstrates unity and diversity. Sleep. 2024.",
            "PMID": "38109788",
            "DOI": "10.1093/sleep/zsad320",
            "source_url": "https://academic.oup.com/sleep/article/47/2/zsad320/7477860",
            "sleep_traits_or_constructs": "six latent sleep-health factors derived from the atlas source sleep traits",
            "reported_external_trait_scope": "1403 phenotypes per factor; 8418 LDSC rows in six supplement sheets",
            "local_inventory_rows": str(sum(len(rows) for rows in morrison_by_sleep.values())),
            "scope_note": "Complete S10-S15 external-trait inventories extracted locally.",
        },
        {
            "screen_id": "goodman_2025",
            "citation": "Goodman et al. Genome-wide association analysis of composite sleep health scores in 413,904 individuals. Communications Biology. 2025.",
            "PMID": "39856408",
            "DOI": "10.1038/s42003-025-07514-0",
            "source_url": "https://www.nature.com/articles/s42003-025-07514-0",
            "sleep_traits_or_constructs": "six composite sleep-health scores",
            "reported_external_trait_scope": "375 representative heritable Pan-UKB traits; 2250 tests",
            "local_inventory_rows": str(sum(1 for row in goodman if row.get("description", "").strip())),
            "scope_note": "Publication reports 375 representative traits and 2250 tests; Supplementary Data 24 has 380 nonblank phenotype rows, including six sleep-related reference rows. All 380 are retained for conservative coverage matching.",
        },
        {
            "screen_id": "dashti_2019",
            "citation": "Dashti et al. Genome-wide association study identifies genetic loci for self-reported habitual sleep duration. Nature Communications. 2019.",
            "PMID": "30846698",
            "DOI": "10.1038/s41467-019-08917-4",
            "source_url": "https://www.nature.com/articles/s41467-019-08917-4",
            "sleep_traits_or_constructs": "sleep duration, short sleep, and long sleep",
            "reported_external_trait_scope": "multiplicity threshold accounts for 224 tested traits; Data 18 displays 45 trait rows",
            "local_inventory_rows": str(sum(1 for row in dashti if row.get("Trait", "").strip())),
            "scope_note": "Do not interpret absence from the 45 displayed rows as absence from all 224 tests.",
        },
        {
            "screen_id": "li_zhao_2020",
            "citation": "Li and Zhao. Automated feature extraction from population wearable device data identified novel loci associated with sleep and circadian rhythms. PLOS Genetics. 2020.",
            "PMID": "33075057",
            "DOI": "10.1371/journal.pgen.1009089",
            "source_url": "https://journals.plos.org/plosgenetics/article?id=10.1371/journal.pgen.1009089",
            "sleep_traits_or_constructs": "12 wearable sleep and circadian traits",
            "reported_external_trait_scope": "53 external traits; 636 tests",
            "local_inventory_rows": str(
                sum(1 for row in lizhao_rows if row.get("Category Details", "").strip())
            ),
            "scope_note": "Complete S3 Table external-trait inventory extracted locally.",
        },
    ]
    write_tsv(args.registry_out, list(registry[0]), registry)

    goodman_keys = {
        (clean_key(row.get("phenocode", "")), clean_key(row.get("coding", "")))
        for row in goodman
        if clean_key(row.get("phenocode", ""))
    }
    study_names: dict[str, set[str]] = {}
    for row in inventory:
        study_names.setdefault(row["screen_id"], set()).add(row["external_trait_normalized"])

    annotated: list[dict[str, str]] = []
    category_counts: dict[str, int] = {}
    for row in universe:
        hits: list[str] = []
        pan_key = (clean_key(row.get("phenocode", "")), clean_key(row.get("coding", "")))
        if pan_key[0] and pan_key in goodman_keys:
            hits.append("goodman_2025:exact_panukbb_key")
        name = normalize(row.get("phenotype_name", ""))
        for study in ("morrison_2024", "dashti_2019", "li_zhao_2020"):
            if name and name in study_names.get(study, set()):
                hits.append(f"{study}:exact_normalized_label")
        distinct_studies = sorted({hit.split(":", 1)[0] for hit in hits})
        if len(distinct_studies) >= 2:
            novelty = "HEAVILY_STUDIED"
        elif len(distinct_studies) == 1:
            novelty = "PREVIOUSLY_SCREENED"
        else:
            novelty = "UNDEREXPLORED"
        category_counts[novelty] = category_counts.get(novelty, 0) + 1
        enriched = dict(row)
        enriched.update(
            {
                "prior_sleep_screen_coverage": ";".join(hits) or "NO_EXACT_LOCAL_INVENTORY_MATCH",
                "prior_direct_sleep_genetics_evidence": (
                    f"EXACT_MATCH_IN_{len(distinct_studies)}_PRIOR_SCREEN_STUDIES"
                    if distinct_studies
                    else "NONE_IN_EXTRACTED_INVENTORIES__NOT_PROOF_OF_NOVELTY"
                ),
                "novelty_priority": novelty,
                "prior_screen_study_match_count": str(len(distinct_studies)),
            }
        )
        annotated.append(enriched)
    fields = list(annotated[0]) if annotated else []
    write_tsv(args.annotated_out, fields, annotated)

    metadata = {
        "schema_version": "1.0.0",
        "matching_policy": {
            "Goodman_2025": "exact Pan-UKB phenocode plus coding",
            "Morrison_2024": "exact normalized external-trait label",
            "Dashti_2019": "exact normalized displayed external-trait label only",
            "Li_Zhao_2020": "exact normalized S3 external-trait label",
            "warning": "No-match means only no exact match in the locally extracted inventories; it is not proof that a pair is novel.",
        },
        "source_sha256": {str(path): sha256(path) for path in required_paths},
        "universe_rows": len(universe),
        "inventory_rows": len(inventory),
        "screen_inventory_counts": {
            key: sum(row["screen_id"] == key for row in inventory)
            for key in sorted({row["screen_id"] for row in inventory})
        },
        "novelty_category_counts": category_counts,
        "annotated_output_sha256": sha256(args.annotated_out),
        "inventory_output_sha256": sha256(args.inventory_out),
        "registry_output_sha256": sha256(args.registry_out),
    }
    args.run_metadata.parent.mkdir(parents=True, exist_ok=True)
    args.run_metadata.write_text(json.dumps(metadata, indent=2, sort_keys=True) + "\n")
    print(
        f"NOVELTY_PRESCREEN_OK universe={len(universe)} inventory={len(inventory)} "
        f"categories={json.dumps(category_counts, sort_keys=True)}"
    )


if __name__ == "__main__":
    main()
