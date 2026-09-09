#!/usr/bin/env python3
"""Build and validate the frozen six-candidate deep novelty reaudit."""
from __future__ import annotations

import argparse
import csv
import hashlib
import io
import json
import re
import zipfile
from pathlib import Path
from urllib.parse import urlencode
from xml.etree import ElementTree as ET


EXPECTED_CANDIDATES = [
    "longsleep__parental_lifespan",
    "shortsleep__parental_lifespan",
    "sleep_apnea__healthspan",
    "sleep_apnea__parental_lifespan",
    "sleep_efficiency__frailty",
    "snoring__healthspan",
]
ALLOWED_CLASSES = {
    "DIRECT_RG_ALREADY_REPORTED",
    "RELATED_GENETIC_EVIDENCE",
    "SHARED_LOCUS_ONLY",
    "MR_ONLY",
    "OBSERVATIONAL_ONLY",
    "NO_DIRECT_GENETIC_RESULT_FOUND",
    "APPARENTLY_UNREPORTED",
    "UNCERTAIN",
}
AS_OF_DATE = "2026-08-28"
MAIN_NS = "http://schemas.openxmlformats.org/spreadsheetml/2006/main"
REL_NS = "http://schemas.openxmlformats.org/officeDocument/2006/relationships"
PKG_REL_NS = "http://schemas.openxmlformats.org/package/2006/relationships"
CELL = re.compile(r"([A-Z]+)([0-9]+)$")


BROAD_QUERIES = {
    "longsleep__parental_lifespan": {
        "query": 'TITLE_ABS:("long sleep" AND lifespan AND (genetic OR GWAS OR LDSC OR MiXeR OR conjFDR OR "Mendelian randomization"))',
        "hit_count": 2,
        "inspected": "38388528;25651884",
        "accepted": "38388528",
    },
    "shortsleep__parental_lifespan": {
        "query": 'TITLE_ABS:("short sleep" AND lifespan AND (genetic OR GWAS OR LDSC OR MiXeR OR conjFDR OR "Mendelian randomization"))',
        "hit_count": 10,
        "inspected": "38388528;40989023;29240764;25651884;22197809;20429945",
        "accepted": "38388528",
    },
    "sleep_apnea__healthspan": {
        "query": 'TITLE_ABS:(("sleep apnea" OR "sleep apnoea" OR OSA) AND (healthspan OR "healthy lifespan") AND (genetic OR GWAS OR LDSC OR "Mendelian randomization"))',
        "hit_count": 0,
        "inspected": "NA",
        "accepted": "NA",
    },
    "sleep_apnea__parental_lifespan": {
        "query": 'TITLE_ABS:(("sleep apnea" OR "sleep apnoea" OR OSA) AND (lifespan OR longevity OR "parental survival" OR "age at death") AND (genetic OR GWAS OR LDSC OR "Mendelian randomization"))',
        "hit_count": 10,
        "inspected": "38711061;41695109;39288580;41740872;40003319;40299556;29067653;17975198;26411274;21143160",
        "accepted": "NA",
    },
    "sleep_efficiency__frailty": {
        "query": 'TITLE_ABS:("sleep efficiency" AND (frailty OR "frailty index"))',
        "hit_count": 20,
        "inspected": "41628923;41527929;40179140;38043367;39002328;37438957;35804185;32820329;31622445;31162039;19793160;22705247",
        "accepted": "19793160;22705247",
    },
    "snoring__healthspan": {
        "query": 'TITLE_ABS:(snoring AND (healthspan OR "healthy lifespan" OR lifespan OR longevity))',
        "hit_count": 12,
        "inspected": "42554603;41682175;39172898;41695109;31562442;34211497;27402049;12619308;17035465;15006946;7713177",
        "accepted": "34211497",
    },
}


SUMMARY = {
    "longsleep__parental_lifespan": {
        "updated_novelty_classification": "DIRECT_RG_ALREADY_REPORTED",
        "classification_confidence": "HIGH",
        "representative_pmids": "38388528",
        "representative_dois": "10.1038/s41398-024-02826-x",
        "phenotype_match": "EXACT_SLEEP_TAIL_AND_EXACT_TIMMERS_PARENTAL_SURVIVAL_GWAS",
        "prior_rg": "-0.24", "prior_se": "NA", "prior_p": "5.50e-9",
        "prior_method": "LDSC",
        "evidence_strength": "DIRECT_NUMERIC_SAME_INPUT_PAIR",
        "dataset_independence": "NOT_INDEPENDENT_SAME_GWAS_PAIR",
        "conclusion": "The Phase-1 estimate reproduces a directly reported LDSC result; this is not an apparently unreported connection.",
        "classification_change": "DOWNGRADED_FROM_APPARENTLY_NOVEL",
        "change_reason": "The original query required parental-lifespan wording and missed a lifespan-titled paper whose Methods identify the exact Timmers parental-survival GWAS.",
    },
    "shortsleep__parental_lifespan": {
        "updated_novelty_classification": "DIRECT_RG_ALREADY_REPORTED",
        "classification_confidence": "HIGH",
        "representative_pmids": "38388528",
        "representative_dois": "10.1038/s41398-024-02826-x",
        "phenotype_match": "EXACT_SLEEP_TAIL_AND_EXACT_TIMMERS_PARENTAL_SURVIVAL_GWAS",
        "prior_rg": "-0.30", "prior_se": "NA", "prior_p": "1.27e-15",
        "prior_method": "LDSC",
        "evidence_strength": "DIRECT_NUMERIC_SAME_INPUT_PAIR",
        "dataset_independence": "NOT_INDEPENDENT_SAME_GWAS_PAIR",
        "conclusion": "The Phase-1 estimate reproduces a directly reported LDSC result; this is not an apparently unreported connection.",
        "classification_change": "DOWNGRADED_FROM_APPARENTLY_NOVEL",
        "change_reason": "The original query required parental-lifespan wording and missed a lifespan-titled paper whose Methods identify the exact Timmers parental-survival GWAS.",
    },
    "sleep_apnea__healthspan": {
        "updated_novelty_classification": "APPARENTLY_UNREPORTED",
        "classification_confidence": "MODERATE",
        "representative_pmids": "36525587;36989840;41292649",
        "representative_dois": "10.1093/sleep/zsac310;10.1016/j.ebiom.2023.104536;10.1101/2025.11.08.25339824",
        "phenotype_match": "NO_EXACT_PAIR_PAPER_FOUND",
        "prior_rg": "NA", "prior_se": "NA", "prior_p": "NA",
        "prior_method": "NO_DIRECT_PAIR_METHOD_FOUND",
        "evidence_strength": "NEGATIVE_SEARCH_PLUS_CHECKSUM_PINNED_114_TRAIT_OSA_SCAN",
        "dataset_independence": "NOT_APPLICABLE_NO_PRIOR_PAIR_ESTIMATE",
        "conclusion": "No direct or pair-specific genetic result was found; apparently unreported remains a search-qualified statement, not a first-ever claim.",
        "classification_change": "UNCHANGED_BUT_WORDING_HARDENED",
        "change_reason": "A broader synonym search returned zero records and the newest 114-trait OSA LDSC supplement did not include healthspan or any lifespan synonym.",
    },
    "sleep_apnea__parental_lifespan": {
        "updated_novelty_classification": "APPARENTLY_UNREPORTED",
        "classification_confidence": "MODERATE",
        "representative_pmids": "36525587;36989840;41292649",
        "representative_dois": "10.1093/sleep/zsac310;10.1016/j.ebiom.2023.104536;10.1101/2025.11.08.25339824",
        "phenotype_match": "NO_EXACT_PAIR_PAPER_FOUND",
        "prior_rg": "NA", "prior_se": "NA", "prior_p": "NA",
        "prior_method": "NO_DIRECT_PAIR_METHOD_FOUND",
        "evidence_strength": "NEGATIVE_SEARCH_PLUS_CHECKSUM_PINNED_114_TRAIT_OSA_SCAN",
        "dataset_independence": "NOT_APPLICABLE_NO_PRIOR_PAIR_ESTIMATE",
        "conclusion": "No direct or pair-specific genetic result was found; apparently unreported remains a search-qualified statement, not a first-ever claim.",
        "classification_change": "UNCHANGED_BUT_WORDING_HARDENED",
        "change_reason": "All ten broad-query hits were false-positive co-occurrences and the newest 114-trait OSA LDSC supplement omitted parental survival, lifespan, longevity, and age at death.",
    },
    "sleep_efficiency__frailty": {
        "updated_novelty_classification": "OBSERVATIONAL_ONLY",
        "classification_confidence": "HIGH",
        "representative_pmids": "19793160;22705247",
        "representative_dois": "10.1111/j.1532-5415.2009.02490.x;10.1016/j.sleep.2012.04.010",
        "phenotype_match": "EXACT_OBJECTIVE_SLEEP_EFFICIENCY_RELATED_FRAILTY_PHENOTYPE",
        "prior_rg": "NA", "prior_se": "NA", "prior_p": "NA",
        "prior_method": "CROSS_SECTIONAL_AND_PROSPECTIVE_OBSERVATIONAL",
        "evidence_strength": "EXACT_PAIR_OBSERVATIONAL_NO_DIRECT_GENETIC_RESULT",
        "dataset_independence": "INDEPENDENT_OBSERVATIONAL_MROS_COHORT",
        "conclusion": "Exact observational evidence exists, but no direct rg, shared-locus, colocalization, or pair-specific MR result was found.",
        "classification_change": "DOWNGRADED_FROM_APPARENTLY_NOVEL_TO_OBSERVATIONAL_ONLY",
        "change_reason": "The broad search recovered an exact actigraphy sleep-efficiency/frailty study (MOR 1.37), while genetic searches and the source GWAS supplement yielded no frailty result.",
    },
    "snoring__healthspan": {
        "updated_novelty_classification": "RELATED_GENETIC_EVIDENCE",
        "classification_confidence": "HIGH",
        "representative_pmids": "32060260;34211497",
        "representative_dois": "10.1038/s41467-020-14625-1;10.3389/fgene.2021.663449",
        "phenotype_match": "RELATED_PARENTAL_AGE_AT_DEATH_GENETIC_RESULT_PLUS_EXACT_HEALTHSPAN_OBSERVATION",
        "prior_rg": "-0.1279", "prior_se": "0.0633", "prior_p": "0.0432",
        "prior_method": "LDSC_RELATED_AGING_TRAIT",
        "evidence_strength": "RELATED_DIRECT_RG_AND_EXACT_OBSERVATIONAL_NULL",
        "dataset_independence": "SAME_SNORING_GWAS_RELATED_EXTERNAL_TRAIT;OVERLAPPING_UKB_OBSERVATION",
        "conclusion": "No prior exact snoring-healthspan genetic result was found; related aging rg exists, while the exact observational healthspan association was null after adjustment.",
        "classification_change": "DOWNGRADED_FROM_APPARENTLY_NOVEL_TO_RELATED_GENETIC_EVIDENCE",
        "change_reason": "The checksum-pinned snoring supplement reports parents' age at death, and a targeted search found a UKB healthspan analysis that included snoring.",
    },
}


def fail(message: str) -> None:
    raise SystemExit(f"ERROR: {message}")


def read_tsv(path: Path) -> list[dict[str, str]]:
    if not path.is_file() or path.stat().st_size == 0:
        fail(f"missing non-empty input: {path}")
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def table_text(fields: list[str], rows: list[dict[str, object]]) -> str:
    buffer = io.StringIO(newline="")
    writer = csv.DictWriter(buffer, fieldnames=fields, delimiter="\t", lineterminator="\n")
    writer.writeheader()
    for row in rows:
        writer.writerow({field: row.get(field, "NA") for field in fields})
    return buffer.getvalue()


def atomic_text(path: Path, payload: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(payload, encoding="utf-8")
    temporary.replace(path)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def column_index(reference: str) -> int:
    match = CELL.match(reference)
    if not match:
        fail(f"invalid XLSX cell reference: {reference}")
    value = 0
    for character in match.group(1):
        value = value * 26 + ord(character) - 64
    return value - 1


def xlsx_sheet(path: Path, sheet_name: str) -> list[list[object]]:
    with zipfile.ZipFile(path) as archive:
        shared = []
        if "xl/sharedStrings.xml" in archive.namelist():
            root = ET.fromstring(archive.read("xl/sharedStrings.xml"))
            for item in root.findall(f"{{{MAIN_NS}}}si"):
                shared.append("".join(node.text or "" for node in item.iter(f"{{{MAIN_NS}}}t")))
        workbook = ET.fromstring(archive.read("xl/workbook.xml"))
        relationships = ET.fromstring(archive.read("xl/_rels/workbook.xml.rels"))
        targets = {node.attrib["Id"]: node.attrib["Target"] for node in relationships.findall(f"{{{PKG_REL_NS}}}Relationship")}
        target = None
        for node in workbook.find(f"{{{MAIN_NS}}}sheets"):
            if node.attrib["name"] == sheet_name:
                target = targets[node.attrib[f"{{{REL_NS}}}id"]]
                break
        if target is None:
            fail(f"worksheet not found: {path}/{sheet_name}")
        member = target.lstrip("/") if target.startswith("/") else "xl/" + target.lstrip("./")
        root = ET.fromstring(archive.read(member))
        rows = []
        for row_node in root.findall(f".//{{{MAIN_NS}}}sheetData/{{{MAIN_NS}}}row"):
            values: dict[int, object] = {}
            for cell in row_node.findall(f"{{{MAIN_NS}}}c"):
                index = column_index(cell.attrib["r"])
                kind = cell.attrib.get("t")
                value_node = cell.find(f"{{{MAIN_NS}}}v")
                if kind == "inlineStr":
                    value: object = "".join(node.text or "" for node in cell.iter(f"{{{MAIN_NS}}}t"))
                elif value_node is None:
                    value = ""
                elif kind == "s":
                    value = shared[int(value_node.text)]
                else:
                    raw = value_node.text or ""
                    try:
                        number = float(raw)
                        value = int(number) if number.is_integer() else number
                    except ValueError:
                        value = raw
                values[index] = value
            rows.append([values.get(index, "") for index in range(max(values, default=-1) + 1)])
        return rows


def normalized(value: object) -> str:
    return re.sub(r"[^a-z0-9]+", " ", str(value).casefold()).strip()


def verify_frozen_family(root: Path) -> list[dict[str, str]]:
    manifest_path = root / "results/validation/candidate_validation_manifest.tsv"
    lock_path = root / "results/validation/candidate_validation_manifest.lock.json"
    manifest = read_tsv(manifest_path)
    if [row["candidate_id"] for row in manifest] != EXPECTED_CANDIDATES:
        fail("candidate family/order differs from the immutable six")
    lock = json.loads(lock_path.read_text(encoding="utf-8"))
    if lock.get("manifest_sha256") != sha256(manifest_path):
        fail("candidate manifest differs from its lock")
    if lock.get("validation_results_accessed_before_lock") is not False:
        fail("candidate family was not frozen before validation")
    return manifest


def verify_local_supplements(root: Path) -> tuple[dict[str, object], dict[str, object]]:
    registry = {row["evidence_source_id"]: row for row in read_tsv(root / "config/published_rg_source_workbooks.tsv")}
    checks = {
        "campos_2020_snoring_supp_data_2": ("rG", "Parents age at death"),
        "jones_2019_actigraphy_supp_data_7": ("Supp Data 7", "frailty"),
    }
    inspected: dict[str, dict[str, object]] = {}
    for source_id, (sheet_name, term) in checks.items():
        source = registry[source_id]
        path = root / "results/_literature_cache/published_rg_workbooks" / source["filename"]
        if sha256(path) != source["sha256"]:
            fail(f"checksum mismatch for {source_id}")
        rows = xlsx_sheet(path, sheet_name)
        matches = [(index + 1, row) for index, row in enumerate(rows) if term.casefold() in " ".join(map(str, row)).casefold()]
        inspected[source_id] = {"path": path, "sha256": source["sha256"], "matches": matches, "row_count": len(rows)}
    parent_matches = inspected["campos_2020_snoring_supp_data_2"]["matches"]
    if len(parent_matches) != 1:
        fail("expected one parents-age-at-death row in the snoring supplement")
    _, snoring_row = parent_matches[0]
    expected = ("parents age at death", 27015805, -0.1279, 0.0633, 0.0432)
    observed = (normalized(snoring_row[1]), int(snoring_row[2]), float(snoring_row[5]), float(snoring_row[6]), float(snoring_row[7]))
    if observed != expected:
        fail(f"snoring related-aging row drifted: {observed}")
    if inspected["jones_2019_actigraphy_supp_data_7"]["matches"]:
        fail("the source sleep-efficiency supplement unexpectedly contains a frailty row")
    return inspected["campos_2020_snoring_supp_data_2"], inspected["jones_2019_actigraphy_supp_data_7"]


def build_queries(root: Path) -> list[dict[str, object]]:
    original = read_tsv(root / "results/analysis/literature_search_queries.tsv")
    selected = [row for row in original if f"{row['sleep_trait']}__{row['external_trait']}" in EXPECTED_CANDIDATES]
    if len(selected) != 24 or any(row["query_status"] != "COMPLETE" for row in selected):
        fail("expected 24 complete original fixed queries for the six candidates")
    rows: list[dict[str, object]] = []
    for row in selected:
        rows.append({
            "candidate_id": f"{row['sleep_trait']}__{row['external_trait']}",
            "query_id": row["query_id"], "database": "Europe_PMC", "search_date": row["as_of_date"],
            "query_scope": "PHASE1_FIXED_NARROW", "query": row["europe_pmc_query"],
            "query_url": row["europe_pmc_url"], "hit_count": row["hit_count"],
            "returned_result_count": row["returned_result_count"], "papers_inspected": "NA",
            "accepted_evidence_pmids": "NA", "status": "COMPLETE",
        })
    for candidate_id in EXPECTED_CANDIDATES:
        item = BROAD_QUERIES[candidate_id]
        rows.append({
            "candidate_id": candidate_id, "query_id": f"{candidate_id}__DEEP_SYNONYM_REAUDIT",
            "database": "Europe_PMC", "search_date": AS_OF_DATE, "query_scope": "DEEP_SYNONYM_REAUDIT",
            "query": item["query"],
            "query_url": "https://www.ebi.ac.uk/europepmc/webservices/rest/search?" + urlencode({
                "query": item["query"], "format": "json", "resultType": "core", "pageSize": 25,
            }),
            "hit_count": item["hit_count"], "returned_result_count": min(int(item["hit_count"]), 25),
            "papers_inspected": item["inspected"], "accepted_evidence_pmids": item["accepted"],
            "status": "COMPLETE_CHECKED_2026-08-28",
        })
    return rows


def evidence_rows(campos: dict[str, object], jones: dict[str, object]) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    for candidate_id, rg, p, mixer, mr in (
        ("longsleep__parental_lifespan", "-0.24", "5.50e-9", "2100 shared variants; 3 conjFDR loci", "IVW beta=0.24, 95% CI=-0.52 to 1.00, p=0.54; null"),
        ("shortsleep__parental_lifespan", "-0.30", "1.27e-15", "3000 shared variants; 3 conjFDR loci", "IVW beta=-0.60, 95% CI=-0.95 to -0.26, p=5.78e-4"),
    ):
        rows.append({
            "candidate_id": candidate_id, "evidence_id": f"wu_2024_{candidate_id}",
            "evidence_type": "DIRECT_RG_ALREADY_REPORTED", "pmid": "38388528", "doi": "10.1038/s41398-024-02826-x",
            "title": "Shared genetic architecture and causal relationship between sleep behaviors and lifespan",
            "source_location": "Main text Methods and Results; Figure 1; Table 1", "phenotype_match": "EXACT_SAME_GWAS_PAIR",
            "dataset_relation": "SAME_SLEEP_GWAS;SAME_TIMMERS_PARENTAL_SURVIVAL_GWAS", "method": "LDSC;MiXeR;conjFDR;two-sample_MR",
            "rg": rg, "se": "NA", "p": p, "other_numeric_result": f"{mixer}; {mr}",
            "evidence_conclusion": "Direct prior report falsifies the original apparently-novel classification.",
            "source_url": "https://pmc.ncbi.nlm.nih.gov/articles/PMC10883970/", "source_sha256": "NA",
            "verification_status": "PRIMARY_FULL_TEXT_NUMERIC_REVIEWED",
        })
    for candidate_id in ("sleep_apnea__healthspan", "sleep_apnea__parental_lifespan"):
        rows.append({
            "candidate_id": candidate_id, "evidence_id": f"garcia_marin_2025_table29_coverage_{candidate_id}",
            "evidence_type": "BROAD_GENETIC_SCAN_NO_EXACT_TRAIT", "pmid": "41292649", "doi": "10.1101/2025.11.08.25339824",
            "title": "Genome-wide analysis in over 1.6 million participants uncovers 147 loci associated with obstructive sleep apnoea",
            "source_location": "Supplementary Table 29; 228 nonblank rows, 114 phenotypes for OSA and OSA-no-BMI",
            "phenotype_match": "NO_HEALTHSPAN_LIFESPAN_PARENTAL_SURVIVAL_LONGEVITY_FRAILTY_OR_AGE_AT_DEATH_ROW",
            "dataset_relation": "EXPANDED_DIFFERENT_OSA_GWAS;NO_TARGET_EXTERNAL_GWAS", "method": "LDSC_TABLE_COVERAGE_AUDIT",
            "rg": "NA", "se": "NA", "p": "NA", "other_numeric_result": "114 distinct comparison phenotypes; zero target-synonym matches",
            "evidence_conclusion": "A broad recent OSA genetic-correlation scan did not test the candidate aging phenotype.",
            "source_url": "https://www.ebi.ac.uk/europepmc/webservices/rest/PMC12642705/supplementaryFiles",
            "source_sha256": "f77c6b82076f11204128bbca9560c330ed78eb33b6c3b35f78506048e4a2a07d",
            "verification_status": "CHECKSUM_PINNED_XLSX_TABLE_REVIEWED",
        })
    rows.extend([
        {
            "candidate_id": "sleep_efficiency__frailty", "evidence_id": "ensrud_2009_sleep_efficiency_frailty",
            "evidence_type": "OBSERVATIONAL_ONLY", "pmid": "19793160", "doi": "10.1111/j.1532-5415.2009.02490.x",
            "title": "Sleep disturbances and frailty status in older community-dwelling men",
            "source_location": "Abstract Results", "phenotype_match": "OBJECTIVE_ACTIGRAPHY_SLEEP_EFFICIENCY;RELATED_CLINICAL_FRAILTY_INDEX",
            "dataset_relation": "INDEPENDENT_MROS_OBSERVATIONAL_COHORT", "method": "CROSS_SECTIONAL_MULTIVARIABLE_ORDINAL_ASSOCIATION",
            "rg": "NA", "se": "NA", "p": "NA", "other_numeric_result": "sleep efficiency <70%: MOR=1.37, 95% CI=1.12-1.67",
            "evidence_conclusion": "Lower objective sleep efficiency was associated with greater frailty, but this is not genetic evidence.",
            "source_url": "https://pubmed.ncbi.nlm.nih.gov/19793160/", "source_sha256": "NA",
            "verification_status": "PRIMARY_ABSTRACT_NUMERIC_REVIEWED",
        },
        {
            "candidate_id": "sleep_efficiency__frailty", "evidence_id": "jones_2019_supp7_no_frailty",
            "evidence_type": "SOURCE_GWAS_SUPPLEMENT_NO_EXACT_TRAIT", "pmid": "30952852", "doi": "10.1038/s41467-019-09576-1",
            "title": "Genetic studies of accelerometer-based sleep measures yield new insights into human sleep behaviour",
            "source_location": f"Supp Data 7; {jones['row_count']} XML worksheet rows", "phenotype_match": "NO_FRAILTY_ROW",
            "dataset_relation": "SAME_SLEEP_EFFICIENCY_GWAS;NO_FRAILTY_COMPARISON", "method": "PUBLISHED_LDSC_SUPPLEMENT_COVERAGE_AUDIT",
            "rg": "NA", "se": "NA", "p": "NA", "other_numeric_result": "zero case-insensitive frailty matches",
            "evidence_conclusion": "The discovery sleep GWAS supplement did not report a frailty comparison.",
            "source_url": "https://pmc.ncbi.nlm.nih.gov/articles/PMC6451011/", "source_sha256": jones["sha256"],
            "verification_status": "LOCAL_CHECKSUM_PINNED_XLSX_REVALIDATED",
        },
        {
            "candidate_id": "snoring__healthspan", "evidence_id": "campos_2020_parents_age_at_death",
            "evidence_type": "RELATED_GENETIC_EVIDENCE", "pmid": "32060260", "doi": "10.1038/s41467-020-14625-1",
            "title": "Insights into the aetiology of snoring from observational and genetic investigations in the UK Biobank",
            "source_location": "Supplementary Data 2 rG row 412", "phenotype_match": "RELATED_PARENTAL_AGE_AT_DEATH_NOT_HEALTHSPAN",
            "dataset_relation": "SAME_SNORING_GWAS;RELATED_AGING_EXTERNAL_GWAS_PMID_27015805", "method": "LDSC",
            "rg": "-0.1279", "se": "0.0633", "p": "0.0432", "other_numeric_result": "FDR=0.081459854",
            "evidence_conclusion": "A nominal related-aging genetic correlation exists, but it is not an exact healthspan result and did not pass its source FDR.",
            "source_url": "https://pmc.ncbi.nlm.nih.gov/articles/PMC7021827/", "source_sha256": campos["sha256"],
            "verification_status": "LOCAL_CHECKSUM_PINNED_XLSX_NUMERIC_REVALIDATED",
        },
        {
            "candidate_id": "snoring__healthspan", "evidence_id": "liu_2021_snoring_healthspan",
            "evidence_type": "OBSERVATIONAL_ONLY", "pmid": "34211497", "doi": "10.3389/fgene.2021.663449",
            "title": "Associations Between Sleep Quality and Health Span: A Prospective Cohort Study Based on 328,850 UK Biobank Participants",
            "source_location": "Table 3", "phenotype_match": "EXACT_SNORING_COMPONENT_AND_ZENIN_STYLE_HEALTHSPAN_TERMINATION",
            "dataset_relation": "OVERLAPPING_UKB_COHORT;OBSERVATIONAL", "method": "MULTIVARIABLE_COX_REGRESSION",
            "rg": "NA", "se": "NA", "p": "NA", "other_numeric_result": "no self-reported snoring: fully adjusted HR=0.99, 95% CI=0.97-1.01",
            "evidence_conclusion": "The exact observational association was null after adjustment, tempering phenotypic interpretation of the positive rg.",
            "source_url": "https://pmc.ncbi.nlm.nih.gov/articles/PMC8239359/", "source_sha256": "NA",
            "verification_status": "PRIMARY_FULL_TEXT_NUMERIC_REVIEWED",
        },
    ])
    return rows


def build_summary(manifest: list[dict[str, str]], queries: list[dict[str, object]]) -> list[dict[str, object]]:
    manifest_by_id = {row["candidate_id"]: row for row in manifest}
    rows = []
    for candidate_id in EXPECTED_CANDIDATES:
        source = manifest_by_id[candidate_id]
        item = dict(SUMMARY[candidate_id])
        if item["updated_novelty_classification"] not in ALLOWED_CLASSES:
            fail(f"invalid novelty class for {candidate_id}")
        candidate_queries = [row for row in queries if row["candidate_id"] == candidate_id]
        papers = sorted({pmid for row in candidate_queries for pmid in str(row["papers_inspected"]).split(";") if pmid != "NA"})
        rows.append({
            "candidate_id": candidate_id, "sleep_trait": source["sleep_trait"], "external_trait": source["external_trait"],
            "original_novelty_classification": source["original_novelty_classification"],
            **item,
            "queries_used": ";".join(str(row["query_id"]) for row in candidate_queries),
            "databases_searched": "Europe_PMC;PubMed;PMC_full_text;checksum_pinned_source_supplements",
            "papers_inspected_pmids": ";".join(papers) if papers else "NA",
            "search_date": AS_OF_DATE,
            "claim_limit": "Literature absence is search-bounded; genetic correlation, overlap, observational association, and causality are not interchangeable.",
        })
    return rows


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", default=".")
    parser.add_argument("--validate-only", action="store_true")
    args = parser.parse_args()
    root = Path(args.root).resolve()
    manifest = verify_frozen_family(root)
    campos, jones = verify_local_supplements(root)
    queries = build_queries(root)
    evidence = evidence_rows(campos, jones)
    summary = build_summary(manifest, queries)
    if [row["candidate_id"] for row in summary] != EXPECTED_CANDIDATES:
        fail("summary candidate family/order drifted")
    if {row["candidate_id"] for row in evidence} != set(EXPECTED_CANDIDATES):
        fail("each candidate must have at least one evidence row")
    query_fields = [
        "candidate_id", "query_id", "database", "search_date", "query_scope", "query", "query_url", "hit_count",
        "returned_result_count", "papers_inspected", "accepted_evidence_pmids", "status",
    ]
    evidence_fields = [
        "candidate_id", "evidence_id", "evidence_type", "pmid", "doi", "title", "source_location", "phenotype_match",
        "dataset_relation", "method", "rg", "se", "p", "other_numeric_result", "evidence_conclusion", "source_url",
        "source_sha256", "verification_status",
    ]
    summary_fields = [
        "candidate_id", "sleep_trait", "external_trait", "original_novelty_classification", "updated_novelty_classification",
        "classification_confidence", "representative_pmids", "representative_dois", "phenotype_match", "prior_rg", "prior_se",
        "prior_p", "prior_method", "evidence_strength", "dataset_independence", "conclusion", "classification_change",
        "change_reason", "queries_used", "databases_searched", "papers_inspected_pmids", "search_date", "claim_limit",
    ]
    outputs = {
        root / "results/validation/deep_novelty_search_queries.tsv": table_text(query_fields, queries),
        root / "results/validation/deep_novelty_evidence.tsv": table_text(evidence_fields, evidence),
        root / "results/validation/deep_novelty_reaudit.tsv": table_text(summary_fields, summary),
    }
    if args.validate_only:
        for path, payload in outputs.items():
            if not path.is_file() or path.read_text(encoding="utf-8") != payload:
                fail(f"output is absent or drifted: {path.relative_to(root)}")
    else:
        for path, payload in outputs.items():
            atomic_text(path, payload)
    counts = {classification: sum(row["updated_novelty_classification"] == classification for row in summary) for classification in sorted(ALLOWED_CLASSES)}
    print(f"DEEP_NOVELTY_REAUDIT_OK candidates={len(summary)} evidence_rows={len(evidence)} query_rows={len(queries)} classes={json.dumps({k: v for k, v in counts.items() if v}, sort_keys=True)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
