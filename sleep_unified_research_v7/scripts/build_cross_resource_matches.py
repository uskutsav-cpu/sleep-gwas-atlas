#!/usr/bin/env python3
"""Build source-aware V7 match annotations from frozen V4 tables.

Comparator source tables are read in place. The output stores row references,
metadata, and adjudication classes; it does not redistribute comparator rows.
"""
from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
import re
import unicodedata
from collections import Counter, defaultdict
from pathlib import Path
from typing import Iterable


ROOT = Path(__file__).resolve().parents[2]
V4 = ROOT / "sleep_unified_research_v4"
SOURCES = ROOT / "discovery_extension" / "provenance" / "prior_screens" / "extracted"
OUT = ROOT / "sleep_unified_research_v7" / "tables"

RESULTS = (
    V4 / "tables" / "core_native_full_precision_rg_v4.tsv",
    V4 / "tables" / "extension_native_full_precision_rg_v4.tsv",
    V4 / "tables" / "validation_native_full_precision_rg_v4.tsv",
)
PHENOTYPES = V4 / "cell_genomics_readiness_v1" / "release_candidate" / "outcome_phenotypes.tsv"
SLEEP_CONSTRUCTS = V4 / "cell_genomics_readiness_v1" / "release_candidate" / "sleep_constructs.tsv"

OUTCOME_SYNONYMS = {
    "alz": ("alzheimer's disease", "alzheimer disease"),
    "parkinson": ("parkinson's disease", "parkinson disease"),
    "mdd": ("major depressive disorder", "depression"),
    "scz": ("schizophrenia",),
    "bipolar": ("bipolar disorder",),
    "adhd": ("attention deficit hyperactivity disorder", "adhd"),
    "ibd": ("inflammatory bowel disease", "ibd"),
    "crohn": ("crohn's disease", "crohn disease"),
    "uc": ("ulcerative colitis",),
    "ra": ("rheumatoid arthritis",),
    "ms": ("multiple sclerosis",),
    "asthma": ("asthma",),
    "bmi": ("body mass index", "bmi"),
    "t2d": ("type 2 diabetes", "type 2 diabetes mellitus"),
    "ldl": ("ldl cholesterol", "low density lipoprotein cholesterol"),
    "hdl": ("hdl cholesterol", "high density lipoprotein cholesterol"),
    "triglycerides": ("triglyceride", "triglycerides"),
    "cad": ("coronary artery disease", "ischemic heart disease", "ischaemic heart disease"),
    "stroke": ("stroke", "any stroke"),
    "atrial_fibrillation": ("atrial fibrillation",),
    "sbp": ("systolic blood pressure",),
    "longevity": ("longevity",),
    "parental_lifespan": ("parental lifespan",),
    "healthspan": ("healthspan", "health span"),
    "frailty": ("frailty index", "frailty"),
    "telomere_length": ("leukocyte telomere length", "telomere length"),
    "grip_strength": ("hand grip strength", "grip strength"),
    "breast_cancer": ("breast cancer",),
    "prostate_cancer": ("prostate cancer",),
    "colorectal_cancer": ("colorectal cancer",),
    "lung_cancer": ("lung cancer",),
    "melanoma": ("cutaneous melanoma", "melanoma"),
}

SLEEP_SYNONYMS = {
    "insomnia": ("insomnia", "sleeplessness"),
    "sleepdur": ("sleep duration", "total sleep duration"),
    "shortsleep": ("short sleep", "short sleep duration"),
    "longsleep": ("long sleep", "long sleep duration"),
    "chronotype": ("chronotype", "morningness", "morning evening person"),
    "sleepiness": ("daytime sleepiness", "sleepiness", "daytime dozing"),
    "napping": ("daytime napping", "napping", "nap"),
    "snoring": ("snoring",),
    "sleep_apnea": ("sleep apnea", "sleep apnoea", "obstructive sleep apnea", "obstructive sleep apnoea"),
    "sleep_efficiency": ("sleep efficiency",),
    "accel_sleep_duration": ("sleep duration mean", "mean sleep duration", "sleep duration"),
    "sleep_timing": ("sleep midpoint", "sleep timing", "l5 timing", "m10 timing", "circadian preference"),
}

MORRISON_DOMAIN = {
    "insomnia": ("Insomnia", "related construct; Morrison sleep factor and atlas source are not treated as identical"),
    "sleepdur": ("Duration", "related construct; factor/source differs from the atlas continuous self-report measure"),
    "shortsleep": ("Duration", "related construct; categorical short-sleep estimand differs"),
    "longsleep": ("Duration", "related construct; categorical long-sleep estimand differs"),
    "chronotype": ("Circadian Preference", "related construct; source and coding may differ"),
    "sleep_timing": ("Circadian Preference", "related construct; sleep timing is not assumed equivalent to chronotype"),
    "sleepiness": ("Alertness", "related construct; source and coding may differ"),
    "napping": ("Alertness", "related construct; napping is not assumed equivalent to alertness"),
    "sleep_efficiency": ("Efficiency", "related construct; source and coding may differ"),
    "accel_sleep_duration": ("Duration", "related construct; device-derived source differs"),
}


def tsv_rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8-sig") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def write_tsv(path: Path, rows: list[dict[str, object]], fields: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=fields, delimiter="\t", extrasaction="ignore", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def norm(text: object) -> str:
    value = unicodedata.normalize("NFKC", str(text or "")).casefold()
    value = re.sub(r"[^\w]+", " ", value, flags=re.UNICODE)
    return " ".join(value.split())


def match_type(candidate: object, aliases: Iterable[str]) -> str:
    value = norm(candidate)
    if not value:
        return "none"
    normalized = [norm(x) for x in aliases if norm(x)]
    if value in normalized:
        return "exact"
    for alias in normalized:
        if re.search(r"(?:^| )" + re.escape(alias) + r"(?: |$)", value):
            return "synonym_or_qualified"
    return "none"


def aliases_for(row: dict[str, str], domain: str) -> tuple[str, ...]:
    trait = row.get("trait_id", "")
    base = [row.get("label", ""), row.get("trait_id", "")]
    table = SLEEP_SYNONYMS if domain == "sleep" else OUTCOME_SYNONYMS
    base.extend(table.get(trait, ()))
    return tuple(dict.fromkeys(x for x in base if x and x != "NA"))


def numeric(value: object) -> float | None:
    try:
        result = float(str(value))
    except (TypeError, ValueError):
        return None
    return result if result == result else None


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as fh:
        for chunk in iter(lambda: fh.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def read_human_atlas(directory: Path) -> tuple[list[dict[str, str]], set[tuple[str, str]], dict[str, dict[str, str]], dict[str, int]]:
    meta_path = directory / "gwasATLAS_v20191115.txt.gz"
    pair_path = directory / "gwasATLAS_v20191115_GC.txt.gz"
    if not meta_path.is_file() or not pair_path.is_file():
        raise FileNotFoundError(f"Expected both official Human GWAS ATLAS files in {directory}")
    with gzip.open(meta_path, "rt", encoding="utf-8") as fh:
        metadata = list(csv.DictReader(fh, delimiter="\t"))
    by_id = {row["id"]: row for row in metadata}
    pairs: set[tuple[str, str]] = set()
    pair_rows = 0
    with gzip.open(pair_path, "rt", encoding="utf-8") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            pair_rows += 1
            if row["id1"] not in by_id or row["id2"] not in by_id:
                raise ValueError(f"Pair table references missing metadata ID: {row['id1']}|{row['id2']}")
            a, b = sorted((row["id1"], row["id2"]))
            pairs.add((a, b))
    if pair_rows != len(pairs):
        raise ValueError(f"Pair table contains duplicate keys: rows={pair_rows}, unique={len(pairs)}")
    file_info = {
        "metadata_sha256": sha256(meta_path),
        "pair_table_sha256": sha256(pair_path),
        "metadata_bytes": meta_path.stat().st_size,
        "pair_table_bytes": pair_path.stat().st_size,
        "metadata_rows": len(metadata),
        "pair_table_rows": pair_rows,
        "unique_pair_keys": len(pairs),
    }
    return metadata, pairs, by_id, file_info


def read_publisher_sheet(path: Path, header_excel_row: int = 2) -> list[dict[str, str]]:
    """Read a publisher XLSX sheet in place; return source rows without copying result values."""
    try:
        from openpyxl import load_workbook
    except ImportError as exc:
        raise RuntimeError("Fan supplement extraction requires openpyxl in the selected Python environment") from exc
    workbook = load_workbook(path, read_only=True, data_only=True)
    try:
        sheet = workbook.active
        header_values = next(sheet.iter_rows(min_row=header_excel_row, max_row=header_excel_row, values_only=True))
        headers = [str(x).strip() if x is not None else "" for x in header_values]
        result = []
        for excel_row, values in enumerate(sheet.iter_rows(min_row=header_excel_row + 1, values_only=True), start=header_excel_row + 1):
            if not any(value is not None for value in values):
                continue
            row = {headers[i]: str(values[i]).strip() if i < len(values) and values[i] is not None else "" for i in range(len(headers)) if headers[i]}
            row["_excel_row"] = str(excel_row)
            result.append(row)
        return result
    finally:
        workbook.close()


def fan_match(record: dict[str, str], fan_sleep_sources: list[dict[str, str]],
              fan_outcome_sources: list[dict[str, str]], fan_results: list[dict[str, str]],
              phenotype_rows: dict[str, dict[str, str]], sleep_rows: dict[str, dict[str, str]]) -> dict[str, object]:
    """Conservative source-aware join to Fan SD11; never exports third-party statistics."""
    sleep_target = sleep_rows[record["sleep_trait"]]
    outcome_target = phenotype_rows[record["outcome_trait"]]
    sleep_aliases = aliases_for(sleep_target, "sleep")
    outcome_aliases = aliases_for(outcome_target, "outcome")
    candidates = []
    for row in fan_results:
        sleep_match = match_type(row.get("Sleep trait", ""), sleep_aliases)
        outcome_match = match_type(row.get("Disease/trait", ""), outcome_aliases)
        if sleep_match != "none" and outcome_match != "none":
            candidates.append((row, sleep_match, outcome_match))

    sleep_meta = [r for r in fan_sleep_sources if match_type(r.get("Phenotype", ""), sleep_aliases) != "none"]
    outcome_meta = [r for r in fan_outcome_sources if match_type(r.get("Disorder/trait", ""), outcome_aliases) != "none"]
    atlas_sleep_pmid = str(record.get("sleep_pmid", ""))
    atlas_outcome_pmid = str(record.get("outcome_pmid", ""))
    same_outcome_pmid = [r for r in outcome_meta if str(r.get("PMID", "")) == atlas_outcome_pmid and atlas_outcome_pmid not in {"", "NA", "UNRESOLVED"}]
    atlas_n = numeric(record.get("sleep_n"))
    same_sleep_author_n = []
    for r in sleep_meta:
        author = norm(r.get("Author", "")).split()
        source = norm(record.get("sleep_source_id", ""))
        sample = numeric(r.get("Sample size"))
        if author and author[0] in source and atlas_n is not None and sample is not None and abs(atlas_n - sample) <= max(1.0, 0.005 * atlas_n):
            same_sleep_author_n.append(r)

    if candidates:
        candidates.sort(key=lambda x: (x[1] != "exact", x[2] != "exact", int(x[0]["_excel_row"])))
        refs = []
        for row, _, _ in candidates:
            # Pointers only: Fan numerical estimates and source labels are not redistributed.
            refs.append(f"Fan2026:SD11:row{row['_excel_row']}")
        source_refs = [f"Fan2026:SD2:row{r['_excel_row']}" for r in sleep_meta]
        source_refs.extend(f"Fan2026:SD4:row{r['_excel_row']}" for r in outcome_meta)
        if same_outcome_pmid and same_sleep_author_n:
            classification = "E"
            note = "Both primary-source identities are supported by the source supplements, but the Fan outcome sample/ancestry and exact analysis release are not sufficiently specified for an exact source/estimand classification."
            sleep_state, outcome_state = "author_and_N_match; definition/release not established", "PMID_match; sample/ancestry not established"
        elif sleep_meta and outcome_meta:
            classification = "B"
            note = "A label-level pair exists, but one or both cited GWAS sources or phenotype/sample definitions differ from the atlas source; it is related coverage, not an exact-source match."
            sleep_state = "related_source_or_definition" if not same_sleep_author_n else "author_and_N_match; definition/release not established"
            outcome_state = "related_source" if not same_outcome_pmid else "PMID_match; sample/ancestry not established"
        else:
            classification = "E"
            note = "A label-level result exists, but the cited source identity or phenotype comparability cannot be established from SD2/SD4."
            sleep_state = "insufficient_source_evidence"
            outcome_state = "insufficient_source_evidence"
        label_state = "exact" if candidates[0][1] == "exact" and candidates[0][2] == "exact" else "synonym_or_qualified"
        return {
            "match_class": classification, "label_match_type": label_state,
            "candidate_result_count": len(candidates), "candidate_result_refs": ";".join(refs + source_refs),
            "sleep_source_match": sleep_state, "outcome_source_match": outcome_state,
            "candidate_sleep_trait": "", "candidate_outcome_trait": "",
            "candidate_sleep_pmid": "", "candidate_outcome_pmid": "",
            "adjudication_note": note,
        }
    # SD11 is the complete 34 x 119 released result matrix. Absence is specific to that table.
    return {
        "match_class": "D", "label_match_type": "none", "candidate_result_count": 0,
        "candidate_result_refs": "", "sleep_source_match": "no label-pair row in complete SD11 scope",
        "outcome_source_match": "no label-pair row in complete SD11 scope",
        "candidate_sleep_trait": "", "candidate_outcome_trait": "",
        "candidate_sleep_pmid": "", "candidate_outcome_pmid": "",
        "adjudication_note": "No label-pair counterpart in the complete 2026 SD11 table; this is release/table-specific absence only, not evidence of novelty.",
    }


def get_comparator_traits(metadata: list[dict[str, str]], target: dict[str, str], side: str,
                          cache: dict[tuple[str, ...], list[tuple[dict[str, str], str]]]) -> list[tuple[dict[str, str], str]]:
    aliases = aliases_for(target, side)
    key = tuple(norm(x) for x in aliases)
    if key in cache:
        return cache[key]
    hits = []
    for row in metadata:
        hit = max((match_type(row.get(k, ""), aliases) for k in ("Trait", "uniqTrait")), key=lambda x: {"none": 0, "synonym_or_qualified": 1, "exact": 2}[x])
        if hit != "none":
            hits.append((row, hit))
    cache[key] = hits
    return hits


def ancestry_compatible(atlas: dict[str, str], comparator: dict[str, str]) -> bool:
    target = norm(atlas.get("ancestry", ""))
    observed = norm(comparator.get("Population", ""))
    if not target or target in {"na", "unknown", "unresolved"}:
        return False
    if "eur" in target:
        return "eur" in observed
    return target in observed


def sample_compatible(atlas: dict[str, str], comparator: dict[str, str]) -> bool:
    left = numeric(atlas.get("sample_size") or atlas.get("N"))
    right = numeric(comparator.get("N"))
    return left is not None and right is not None and abs(left - right) <= max(1.0, 0.005 * left)


def human_match(atlas_row: dict[str, str], metadata: list[dict[str, str]], pairs: set[tuple[str, str]],
                label_cache: dict[tuple[str, ...], list[tuple[dict[str, str], str]]]) -> dict[str, object]:
    # Result rows carry separate sleep and outcome definitions; identify the sleep target by its role in the frozen sleep dictionary.
    # `phenotype` is the outcome row for lookup, while sleep IDs are passed in atlas_row.
    sleep_source_id = atlas_row["sleep_source_id"]
    sleep_source_pmid = atlas_row["sleep_pmid"]
    outcome_source_pmid = atlas_row["outcome_pmid"]
    sleep_aliases = tuple(x for x in atlas_row["sleep_aliases"].split("||") if x)
    outcome_aliases = tuple(x for x in atlas_row["outcome_aliases"].split("||") if x)

    # Search metadata over declared trait labels/synonyms, then resolve only pairs present in the released pair table.
    sleep_target = {"trait_id": atlas_row["sleep_trait"], "label": atlas_row["sleep_trait"]}
    outcome_target = {"trait_id": atlas_row["outcome_trait"], "label": atlas_row["outcome_trait"]}
    sleep_hits = get_comparator_traits(metadata, sleep_target, "sleep", label_cache)
    outcome_hits = get_comparator_traits(metadata, outcome_target, "outcome", label_cache)

    candidate_rows = []
    for srow, smatch in sleep_hits:
        for orow, omatch in outcome_hits:
            if srow["id"] == orow["id"]:
                continue
            key = tuple(sorted((srow["id"], orow["id"])))
            if key not in pairs:
                continue
            s_exact = str(srow.get("PMID", "")) == sleep_source_pmid and sleep_source_pmid not in {"", "NA", "UNRESOLVED"}
            o_exact = str(orow.get("PMID", "")) == outcome_source_pmid and outcome_source_pmid not in {"", "NA", "UNRESOLVED"}
            s_source_fit = s_exact and sample_compatible({"sample_size": atlas_row["sleep_n"]}, srow) and ancestry_compatible({"ancestry": atlas_row["ancestry"]}, srow)
            o_source_fit = o_exact and sample_compatible({"sample_size": atlas_row["outcome_n"]}, orow) and ancestry_compatible({"ancestry": atlas_row["ancestry"]}, orow)
            candidate_rows.append({
                "sleep_id": srow["id"], "outcome_id": orow["id"], "sleep_match": smatch,
                "outcome_match": omatch, "sleep_source_match": s_exact, "outcome_source_match": o_exact,
                "sleep_source_fit": s_source_fit, "outcome_source_fit": o_source_fit,
                "sleep_pmid": srow.get("PMID", ""), "outcome_pmid": orow.get("PMID", ""),
                "sleep_trait": srow.get("Trait", ""), "outcome_trait": orow.get("Trait", ""),
                "population": srow.get("Population", ""), "outcome_population": orow.get("Population", ""),
            })
    candidate_rows.sort(key=lambda r: (-(int(r["sleep_source_fit"]) + int(r["outcome_source_fit"])), r["sleep_id"], r["outcome_id"]))

    if candidate_rows:
        best = candidate_rows[0]
        if (best["sleep_source_fit"] and best["outcome_source_fit"]
                and best["sleep_match"] == "exact" and best["outcome_match"] == "exact"):
            # PMID/N/ancestry/label agreement does not establish release,
            # ascertainment/coding, effect direction, or estimator compatibility.
            classification = "E"
            note = "Bibliographic candidate only: both PMIDs, approximate sample sizes, ancestry labels, and phenotype labels align, but source release, ascertainment/coding, effect direction, and estimator/reference compatibility are not established; exact source/estimand equivalence is unresolved."
        elif best["sleep_source_match"] and best["outcome_source_match"]:
            classification = "B"
            note = "Same source publications are identified, but sample-size, ancestry, or release compatibility is not established to the frozen exact-match rule."
        elif best["sleep_match"] in {"exact", "synonym_or_qualified"} and best["outcome_match"] in {"exact", "synonym_or_qualified"}:
            classification = "B"
            note = "A related label-level result exists, but at least one source differs; this is not an exact-source comparison."
        else:
            classification = "E"
            note = "Candidate row found but phenotype/source comparability cannot be established."
        refs = [f"HumanGWASATLAS:v20191115:{r['sleep_id']}|{r['outcome_id']}" for r in candidate_rows]
        source_state = "both_pmids_match_metadata_only" if best["sleep_source_match"] and best["outcome_source_match"] else ("one_pmid_matches_metadata_only" if best["sleep_source_match"] or best["outcome_source_match"] else "different_or_unresolved_sources")
        label_state = "exact" if best["sleep_match"] == "exact" and best["outcome_match"] == "exact" else "synonym_or_qualified"
    else:
        unresolved = any(x in {"", "NA", "UNRESOLVED"} for x in (sleep_source_pmid, outcome_source_pmid))
        classification = "E" if unresolved else "D"
        note = "No counterpart in the complete v20191115 pair table for the searched source/label rules; release-specific absence only." if classification == "D" else "No counterpart identified and at least one source PMID is unresolved."
        refs = []
        source_state = "unresolved" if unresolved else "not_present_in_release"
        label_state = "none"
        best = {}
    return {
        "match_class": classification, "label_match_type": label_state,
        "candidate_result_count": len(candidate_rows), "candidate_result_refs": ";".join(refs),
        "sleep_source_match": source_state, "outcome_source_match": ("pmid_match_metadata_only" if best.get("outcome_source_match") else "different_or_absent") if best else "not_present_in_release",
        "candidate_sleep_trait": best.get("sleep_trait", ""), "candidate_outcome_trait": best.get("outcome_trait", ""),
        "candidate_sleep_pmid": best.get("sleep_pmid", ""), "candidate_outcome_pmid": best.get("outcome_pmid", ""),
        "adjudication_note": note,
    }


def parse_goodman_rows(path: Path) -> list[dict[str, str]]:
    raw = list(csv.reader(path.open(newline="", encoding="utf-8-sig"), delimiter="\t"))
    header_idx = next(i for i, row in enumerate(raw) if len(row) >= 20 and row[1:5] == ["phenocode", "coding", "description", "coding_description"])
    header = raw[header_idx]
    score_names = raw[header_idx - 1]
    score_cols = [(score_names[i], i, i + 1) for i in range(5, 17, 2)]
    out = []
    for line_num, row in enumerate(raw[header_idx + 1:], header_idx + 2):
        if len(row) < len(header) or not row[1] or row[1] in {"phenocode"}:
            continue
        for score, rg_col, p_col in score_cols:
            if row[rg_col] and row[p_col] and row[rg_col] != "NA":
                out.append({"line": str(line_num), "phenocode": row[1], "coding": row[2], "description": row[3], "score": score, "rg": row[rg_col], "p": row[p_col]})
    return out


def panukbb_code(trait_id: str) -> tuple[str, str] | None:
    parts = trait_id.split("__")
    if len(parts) < 5 or not trait_id.startswith("panukbb_"):
        return None
    kind = parts[0].replace("panukbb_", "")
    if kind not in {"continuous", "categorical"}:
        return None
    return parts[1], parts[3]


def build() -> dict[str, object]:
    parser = argparse.ArgumentParser()
    parser.add_argument("--gwas-atlas-dir", type=Path, required=True, help="Directory containing the two verified Human GWAS ATLAS release 3 gzip files")
    parser.add_argument("--fan-data2", type=Path, required=True, help="Verified Fan 2026 SD2 source-definition workbook on the SSD")
    parser.add_argument("--fan-data4", type=Path, required=True, help="Verified Fan 2026 SD4 source-definition workbook on the SSD")
    parser.add_argument("--fan-data11", type=Path, required=True, help="Verified Fan 2026 SD11 result workbook on the SSD")
    args = parser.parse_args()

    metadata_rows = tsv_rows(PHENOTYPES)
    sleep_rows = tsv_rows(SLEEP_CONSTRUCTS)
    pheno_by_id = {r["trait_id"]: r for r in metadata_rows}
    sleep_by_id = {r["trait_id"]: r for r in sleep_rows}

    atlas_records: list[dict[str, str]] = []
    for path in RESULTS:
        for result in tsv_rows(path):
            sleep = sleep_by_id[result["sleep_trait"]]
            outcome = pheno_by_id[result["outcome_trait"]]
            record_id = f"{result['stage']}|{result['sleep_trait']}|{result['outcome_trait']}"
            atlas_records.append({
                **result,
                "atlas_record_id": record_id,
            "sleep_source_id": sleep["source_id"],
            "sleep_source_release": sleep["source_release"],
                "sleep_pmid": pheno_by_id[result["sleep_trait"]].get("PMID", "UNRESOLVED"),
                "outcome_source_id": outcome["source_id"],
                "outcome_source_release": outcome["source_release_label"],
                "outcome_pmid": outcome.get("PMID", "UNRESOLVED"),
                "ancestry": outcome.get("ancestry", ""),
                "sleep_n": sleep.get("N", ""),
                "outcome_n": outcome.get("sample_size", ""),
                "sleep_aliases": "||".join(aliases_for(pheno_by_id[result["sleep_trait"]], "sleep")),
                "outcome_aliases": "||".join(aliases_for(outcome, "outcome")),
                "sleep_definition": sleep.get("wording_device_and_coding", sleep.get("construct", "")),
                "outcome_definition": outcome.get("phenotype_definition", ""),
                "h2_z_pass": result.get("pairwise_h2_z_ge_4_both", ""),
                "intercept_pass": result.get("pairwise_h2_intercept_le_1_2_both", ""),
            })
    if len(atlas_records) != 1637 or len({r["atlas_record_id"] for r in atlas_records}) != 1637:
        raise ValueError(f"Expected 1,637 unique frozen records, observed {len(atlas_records)}")

    metadata, pairs, by_id, atlas_info = read_human_atlas(args.gwas_atlas_dir)
    human_label_cache: dict[tuple[str, ...], list[tuple[dict[str, str], str]]] = {}
    all_rows: list[dict[str, object]] = []
    for record in atlas_records:
        match = human_match(record, metadata, pairs, human_label_cache)
        all_rows.append({
            **{k: record[k] for k in ("atlas_record_id", "stage", "sleep_trait", "outcome_trait", "sleep_source_id", "sleep_source_release", "outcome_source_id", "outcome_source_release", "rg", "se", "p", "reproduction_status", "h2_z_pass", "intercept_pass", "sleep_definition", "outcome_definition", "ancestry")},
            "comparator": "human_gwas_atlas", "comparator_release": "v20191115",
            "source_table": "gwasATLAS_v20191115_GC.txt.gz; metadata gwasATLAS_v20191115.txt.gz",
            **match,
        })

    # Morrison tables contain complete per-domain source/phenome result tables.
    morrison_files = sorted(SOURCES.glob("morrison_2024_S1[0-5]_*.tsv"))
    morrison_index: dict[str, list[dict[str, str]]] = defaultdict(list)
    for path in morrison_files:
        rows = tsv_rows(path)
        for line, row in enumerate(rows, start=2):
            morrison_index[row.get("sleep", "")].append({**row, "_line": str(line), "_file": path.name})
    for record in atlas_records:
        outcome = pheno_by_id[record["outcome_trait"]]
        target_aliases = aliases_for(outcome, "outcome")
        domain, relation = MORRISON_DOMAIN.get(record["sleep_trait"], ("", ""))
        candidates = []
        if domain:
            for row in morrison_index.get(domain, []):
                mt = match_type(row.get("phenotype", ""), target_aliases)
                if mt != "none":
                    candidates.append((row, mt))
        candidates.sort(key=lambda item: (item[1] != "exact", item[0]["_file"], int(item[0]["_line"])))
        if candidates:
            cls = "B"
            mt = candidates[0][1]
            refs = ";".join(f"{x['_file']}:line{x['_line']}" for x, _ in candidates)
            note = relation + "; exact source equivalence is not asserted."
        else:
            cls = "D"
            mt = "none"
            refs = ""
            note = "No outcome-label counterpart for this sleep-domain table in the complete extracted S10-S15 results; absence is limited to these tables."
        all_rows.append({
            **{k: record[k] for k in ("atlas_record_id", "stage", "sleep_trait", "outcome_trait", "sleep_source_id", "sleep_source_release", "outcome_source_id", "outcome_source_release", "rg", "se", "p", "reproduction_status", "h2_z_pass", "intercept_pass", "sleep_definition", "outcome_definition", "ancestry")},
            "comparator": "morrison_2024", "comparator_release": "SLEEP 2024; S10-S15 LDSC tables",
            "source_table": "discovery_extension/provenance/prior_screens/extracted/",
            "match_class": cls, "label_match_type": mt, "candidate_result_count": len(candidates),
            "candidate_result_refs": refs, "sleep_source_match": "related_factor_source",
            "outcome_source_match": "source_not_identical_or_not_verified",
            "candidate_sleep_trait": domain, "candidate_outcome_trait": candidates[0][0].get("phenotype", "") if candidates else "",
            "candidate_sleep_pmid": "", "candidate_outcome_pmid": "",
            "adjudication_note": note,
        })

    # Goodman SD24 uses six composite sleep-health scores and UK Biobank field-coded outcomes.
    goodman_path = SOURCES / "goodman_2025_SD24_ldsc.tsv"
    goodman_rows = parse_goodman_rows(goodman_path)
    goodman_index: dict[tuple[str, str], list[dict[str, str]]] = defaultdict(list)
    goodman_by_label: dict[str, list[dict[str, str]]] = defaultdict(list)
    for row in goodman_rows:
        goodman_index[(row["phenocode"], norm(row["coding"]))].append(row)
        goodman_by_label[norm(row["description"])].append(row)
    for record in atlas_records:
        outcome = pheno_by_id[record["outcome_trait"]]
        code = panukbb_code(record["outcome_trait"])
        candidates: list[dict[str, str]] = []
        label_type = "none"
        if code:
            candidates = goodman_index.get((code[0], norm(code[1])), [])
            if candidates:
                label_type = "source_code_match"
        if not candidates:
            aliases = aliases_for(outcome, "outcome")
            for row in goodman_rows:
                mt = match_type(row["description"], aliases)
                if mt != "none":
                    candidates.append(row)
            if candidates:
                label_type = "synonym_or_qualified"
        candidates = sorted({(r["line"], r["score"]): r for r in candidates}.values(), key=lambda x: (x["line"], x["score"]))
        if candidates:
            cls = "B" if label_type == "source_code_match" else "E"
            note = "Outcome field/code match to a published UK Biobank result, but the comparator exposure is a composite sleep-health score, not the atlas's individual sleep phenotype; source versions are not treated as identical." if cls == "B" else "Label candidate exists, but outcome coding/source equivalence is insufficiently established."
            refs = ";".join(f"goodman_2025_SD24_ldsc.tsv:line{x['line']}:{x['score']}" for x in candidates)
        else:
            cls = "D"
            label_type = "none"
            refs = ""
            note = "No outcome phenotype counterpart in the complete extracted SD24 table for the searched code/label rules; this is limited to that table's representative-trait scope."
        all_rows.append({
            **{k: record[k] for k in ("atlas_record_id", "stage", "sleep_trait", "outcome_trait", "sleep_source_id", "sleep_source_release", "outcome_source_id", "outcome_source_release", "rg", "se", "p", "reproduction_status", "h2_z_pass", "intercept_pass", "sleep_definition", "outcome_definition", "ancestry")},
            "comparator": "goodman_2025", "comparator_release": "Communications Biology 2025; Supplementary Data 24",
            "source_table": "discovery_extension/provenance/prior_screens/extracted/goodman_2025_SD24_ldsc.tsv",
            "match_class": cls, "label_match_type": label_type, "candidate_result_count": len(candidates),
            "candidate_result_refs": refs, "sleep_source_match": "composite_sleep_health_score",
            "outcome_source_match": "UKB_field_code" if label_type == "source_code_match" else "not_verified",
            "candidate_sleep_trait": ";".join(sorted(set(x["score"] for x in candidates))),
            "candidate_outcome_trait": candidates[0]["description"] if candidates else "",
            "candidate_sleep_pmid": "", "candidate_outcome_pmid": "",
            "adjudication_note": note,
        })

    # Fan et al. 2026: source definitions are read from SD2/SD4; SD11 numeric cells are never exported.
    fan_sleep_sources = read_publisher_sheet(args.fan_data2)
    fan_outcome_sources = read_publisher_sheet(args.fan_data4)
    fan_results = read_publisher_sheet(args.fan_data11)
    if len(fan_sleep_sources) != 34 or len(fan_outcome_sources) != 113 or len(fan_results) != 4046:
        raise ValueError(f"Unexpected Fan supplement dimensions: SD2={len(fan_sleep_sources)}, SD4={len(fan_outcome_sources)}, SD11={len(fan_results)}")
    fan_pair_keys = [(x.get("Sleep trait", ""), x.get("Disease/trait", "")) for x in fan_results]
    if len(set(fan_pair_keys)) != len(fan_pair_keys):
        raise ValueError("Fan SD11 has duplicate sleep/outcome label keys")
    for record in atlas_records:
        match = fan_match(record, fan_sleep_sources, fan_outcome_sources, fan_results, pheno_by_id, sleep_by_id)
        all_rows.append({
            **{k: record[k] for k in ("atlas_record_id", "stage", "sleep_trait", "outcome_trait", "sleep_source_id", "sleep_source_release", "outcome_source_id", "outcome_source_release", "rg", "se", "p", "reproduction_status", "h2_z_pass", "intercept_pass", "sleep_definition", "outcome_definition", "ancestry")},
            "comparator": "fan_2026", "comparator_release": "Communications Medicine 2026; SD2/SD4/SD11",
            "source_table": "Fan 2026 publisher supplements; SSD source hashes in build manifest",
            **match,
        })

    fields = ["atlas_record_id", "stage", "sleep_trait", "outcome_trait", "sleep_source_id", "sleep_source_release", "outcome_source_id", "outcome_source_release", "sleep_definition", "outcome_definition", "ancestry", "rg", "se", "p", "reproduction_status", "h2_z_pass", "intercept_pass", "comparator", "comparator_release", "source_table", "match_class", "label_match_type", "candidate_result_count", "candidate_result_refs", "sleep_source_match", "outcome_source_match", "candidate_sleep_trait", "candidate_outcome_trait", "candidate_sleep_pmid", "candidate_outcome_pmid", "adjudication_note"]
    write_tsv(OUT / "cross_resource_match_annotations.tsv", all_rows, fields)

    summaries = []
    for comparator in ("human_gwas_atlas", "morrison_2024", "goodman_2025", "fan_2026"):
        rows = [r for r in all_rows if r["comparator"] == comparator]
        counts = Counter(str(r["match_class"]) for r in rows)
        candidates = [r for r in rows if int(r["candidate_result_count"]) > 0]
        failed = sum(str(r["match_class"]) in {"B", "C", "E"} for r in candidates)
        qualified = sum(
            r["stage"] in {"core", "extension"}
            and r["h2_z_pass"] == "True" and r["intercept_pass"] == "True"
            and str(r["sleep_source_id"]) not in {"", "UNKNOWN", "UNRESOLVED"}
            and str(r["outcome_source_id"]) not in {"", "UNKNOWN", "UNRESOLVED"}
            and bool(str(r.get("sleep_definition", "")).strip())
            and bool(str(r.get("outcome_definition", "")).strip())
            and str(r.get("sleep_source_release", "")) not in {"", "UNKNOWN", "UNRESOLVED"}
            and str(r.get("outcome_source_release", "")) not in {"", "UNKNOWN", "UNRESOLVED"}
            and str(r["match_class"]) in {"C", "D"}
            for r in rows
        )
        summaries.append({
            "comparator": comparator, "atlas_records": len(rows),
            "A_exact_source_estimand": counts["A"], "B_related_nonidentical": counts["B"],
            "C_incompatible": counts["C"], "D_absent_in_snapshot_scope": counts["D"],
            "E_insufficient_evidence": counts["E"], "exact_match_coverage_percent": f"{100*counts['A']/len(rows):.4f}",
            "apparent_label_or_code_candidates": len(candidates),
            "strict_match_failure_candidates_B_C_E": failed,
            "strict_match_failure_percent": f"{100*failed/len(candidates):.4f}" if candidates else "NA",
            "eligible_core_extension_without_A_B_E_in_snapshot": qualified,
            "qualified_match_absent_core_extension_across_four_snapshots": "NA",
        })
    # Union classification is intentionally conservative: an unresolved comparator prevents an atlas-only statement.
    by_record: dict[str, list[dict[str, object]]] = defaultdict(list)
    for row in all_rows:
        by_record[str(row["atlas_record_id"])].append(row)
    union_rows = []
    qualified_rows = []
    record_meta = {r["atlas_record_id"]: r for r in atlas_records}
    for record_id, group in by_record.items():
        classes = {str(x["match_class"]) for x in group}
        if "A" in classes:
            cls = "A"
        elif "B" in classes:
            cls = "B"
        elif "E" in classes:
            cls = "E"
        elif "C" in classes:
            cls = "C"
        else:
            cls = "D"
        rmeta = record_meta[record_id]
        qao = bool(
            rmeta["stage"] in {"core", "extension"}
            and rmeta["h2_z_pass"] == "True" and rmeta["intercept_pass"] == "True"
            and rmeta["sleep_source_id"] not in {"", "UNKNOWN", "UNRESOLVED"}
            and rmeta["outcome_source_id"] not in {"", "UNKNOWN", "UNRESOLVED"}
            and bool(str(rmeta.get("sleep_definition", "")).strip())
            and bool(str(rmeta.get("outcome_definition", "")).strip())
            and rmeta["sleep_source_release"] not in {"", "UNKNOWN", "UNRESOLVED"}
            and rmeta["outcome_source_release"] not in {"", "UNKNOWN", "UNRESOLVED"}
            and cls in {"C", "D"}
        )
        union_rows.append({"atlas_record_id": record_id, "match_class": cls, "candidate_result_count": sum(int(x["candidate_result_count"]) for x in group), "qualified_match_absent": str(qao)})
        if qao:
            by_comp = {str(x["comparator"]): str(x["match_class"]) for x in group}
            qualified_rows.append({
                "atlas_record_id": record_id, "stage": rmeta["stage"], "sleep_trait": rmeta["sleep_trait"],
                "outcome_trait": rmeta["outcome_trait"], "sleep_source_id": rmeta["sleep_source_id"],
                "outcome_source_id": rmeta["outcome_source_id"], "rg": rmeta["rg"], "se": rmeta["se"], "p": rmeta["p"],
                "human_gwas_atlas_class": by_comp.get("human_gwas_atlas", ""),
                "morrison_2024_class": by_comp.get("morrison_2024", ""),
                "goodman_2025_class": by_comp.get("goodman_2025", ""),
                "fan_2026_class": by_comp.get("fan_2026", ""),
                "interpretation": "Qualified match-absent coverage within the four searched row-level snapshots; not a global absence, novelty, or biological-discovery claim.",
            })
    union_counts = Counter(x["match_class"] for x in union_rows)
    union_candidates = [x for x in union_rows if x["candidate_result_count"] > 0]
    union_failed = sum(x["match_class"] in {"B", "C", "E"} for x in union_candidates)
    summaries.append({
        "comparator": "union_accessible_snapshots", "atlas_records": len(union_rows),
        "A_exact_source_estimand": union_counts["A"], "B_related_nonidentical": union_counts["B"],
        "C_incompatible": union_counts["C"], "D_absent_in_snapshot_scope": union_counts["D"],
        "E_insufficient_evidence": union_counts["E"], "exact_match_coverage_percent": f"{100*union_counts['A']/len(union_rows):.4f}",
        "apparent_label_or_code_candidates": len(union_candidates),
        "strict_match_failure_candidates_B_C_E": union_failed,
        "strict_match_failure_percent": f"{100*union_failed/len(union_candidates):.4f}" if union_candidates else "NA",
        "eligible_core_extension_without_A_B_E_in_snapshot": "NA",
        "qualified_match_absent_core_extension_across_four_snapshots": sum(x["qualified_match_absent"] == "True" for x in union_rows),
    })
    sum_fields = list(summaries[0])
    write_tsv(OUT / "cross_resource_match_summary.tsv", summaries, sum_fields)
    qualified_path = ROOT / "sleep_unified_research_v7" / "tables" / "qualified_match_absent_across_four_snapshots.tsv"
    write_tsv(qualified_path, qualified_rows, list(qualified_rows[0]) if qualified_rows else ["atlas_record_id", "interpretation"])
    old_qualified_path = ROOT / "sleep_unified_research_v7" / "tables" / "qualified_atlas_only.tsv"
    if old_qualified_path.exists():
        old_qualified_path.unlink()

    availability = [
        {"comparator": "human_gwas_atlas", "snapshot": "release 3 v20191115", "status": "ACCESSIBLE_COMPLETE_PAIR_TABLE", "table_data_rows": atlas_info["metadata_rows"], "result_cells": atlas_info["pair_table_rows"], "source_sha256": atlas_info["pair_table_sha256"], "source_reference": "https://atlas.ctglab.nl/; https://atlas.ctglab.nl/documentation", "included_in_union": "yes", "scope_note": "Complete released pair table; absence is release-specific, not novelty."},
        {"comparator": "morrison_2024", "snapshot": "SLEEP 2024 S10-S15", "status": "ACCESSIBLE_COMPLETE_SUPPLEMENT_TABLES", "table_data_rows": sum(len(v) for v in morrison_index.values()), "result_cells": sum(len(v) for v in morrison_index.values()), "source_sha256": "see build manifest by source table", "source_reference": "https://academic.oup.com/sleep/article/47/2/zsad320/7477860", "included_in_union": "yes", "scope_note": "Six published sleep-factor tables; outcome absence is confined to this table scope."},
        {"comparator": "goodman_2025", "snapshot": "Communications Biology 2025 Supplementary Data 24", "status": "ACCESSIBLE_COMPLETE_SUPPLEMENT_TABLE", "table_data_rows": len({(x["line"], x["phenocode"], x["coding"]) for x in goodman_rows}), "result_cells": len(goodman_rows), "source_sha256": "see build manifest by source table", "source_reference": "https://www.nature.com/articles/s42003-025-07514-0", "included_in_union": "yes", "scope_note": "Six composite sleep-health scores and representative UK Biobank traits; not individual-sleep-trait equivalents."},
        {"comparator": "CTG-VL", "snapshot": "live platform export", "status": "NO_DIRECT_EXPORT_CAPTURED", "table_data_rows": "NA", "result_cells": "NA", "source_sha256": "NA", "source_reference": "https://ctg.cncr.nl/", "included_in_union": "no", "scope_note": "Morrison S10-S15 are included as published result tables; no additional direct query/export is counted."},
        {"comparator": "SleepChart 2026", "snapshot": "Sleep chart of biological ageing clocks", "status": "PRIOR_V4_AUDIT_CONTEXT_ONLY", "table_data_rows": "1050 populated rg/SE/P estimates in previously inspected SF5", "result_cells": "NA", "source_sha256": "prior inspection hash not present in active checkout", "source_reference": "https://www.nature.com/articles/s41586-026-10524-5; V4 benchmark_matrix.tsv", "included_in_union": "no", "scope_note": "V4 audit documents partial/related short/long-duration disease coverage. Verified workbook is absent from the active checkout and was not redownloaded; not a negative search and not wholly incommensurate."},
        {"comparator": "fan_2026", "snapshot": "Communications Medicine 2026 SD2/SD4/SD11", "status": "ACCESSIBLE_COMPLETE_PUBLISHED_PAIR_TABLE", "table_data_rows": len(fan_results), "result_cells": len(fan_results), "source_sha256": "see build manifest for SD2/SD4/SD11", "source_reference": "https://www.nature.com/articles/s43856-026-01656-w", "included_in_union": "yes", "scope_note": "Complete 34 x 119 SD11 source-pair matrix; SD2/SD4 provide sleep/outcome metadata. V7 outputs source row references and classifications only."},
    ]
    write_tsv(OUT / "comparator_access_status.tsv", availability, list(availability[0]))

    manifest = {
        "protocol_sha256": sha256(ROOT / "sleep_unified_research_v7" / "protocol" / "empirical_benchmark_preregistration_v1.md"),
        "human_gwas_atlas": atlas_info,
        "frozen_atlas_records": len(atlas_records),
        "annotation_rows": len(all_rows),
        "comparator_rows": {k: sum(x["comparator"] == k for x in all_rows) for k in ("human_gwas_atlas", "morrison_2024", "goodman_2025", "fan_2026")},
        "input_files": {str(p.relative_to(ROOT)): sha256(p) for p in (*RESULTS, PHENOTYPES, SLEEP_CONSTRUCTS)},
        "comparator_source_sha256": {str(p.relative_to(ROOT)): sha256(p) for p in (*morrison_files, goodman_path)},
        "fan_publisher_inputs": {
            "SD2": {"sha256": sha256(args.fan_data2), "bytes": args.fan_data2.stat().st_size, "data_rows": len(fan_sleep_sources), "path": str(args.fan_data2)},
            "SD4": {"sha256": sha256(args.fan_data4), "bytes": args.fan_data4.stat().st_size, "data_rows": len(fan_outcome_sources), "path": str(args.fan_data4)},
            "SD11": {"sha256": sha256(args.fan_data11), "bytes": args.fan_data11.stat().st_size, "data_rows": len(fan_results), "unique_sleep_outcome_pairs": len(set(fan_pair_keys)), "path": str(args.fan_data11)},
        },
        "protocol_deviations_sha256": sha256(ROOT / "sleep_unified_research_v7" / "protocol" / "protocol_deviations.md"),
        "build_script_sha256": sha256(Path(__file__).resolve()),
        "generated_table_sha256": {
            str(path.relative_to(ROOT)): sha256(path)
            for path in (
                OUT / "cross_resource_match_annotations.tsv",
                OUT / "cross_resource_match_summary.tsv",
                OUT / "comparator_access_status.tsv",
                qualified_path,
            )
        },
        "generated_utc": "2026-10-10",
    }
    (ROOT / "sleep_unified_research_v7" / "manifests").mkdir(parents=True, exist_ok=True)
    (ROOT / "sleep_unified_research_v7" / "manifests" / "cross_resource_build_manifest.json").write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps({"summary": summaries, "manifest": manifest}, indent=2))
    return manifest


if __name__ == "__main__":
    build()
