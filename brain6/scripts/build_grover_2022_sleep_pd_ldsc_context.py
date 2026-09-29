#!/usr/bin/env python3
"""Record the published negative LDSC screen relevant to long sleep–PD."""
from __future__ import annotations

import csv
import hashlib
import json
import re
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
MASTER = ROOT / "brain6/manifests/gwas_master.tsv"
OUTPUT = ROOT / "brain6/results/novelty/grover_2022_sleep_pd_ldsc_context.tsv"
PROVENANCE = OUTPUT.with_suffix(".provenance.json")
ARTICLE_XML = ROOT / "brain6/qc/replication_sources/grover_2022_article_fullTextXML.xml"
SUPPLEMENTARY_FILE = "Table_1.xlsx"
FIELDS = [
    "sleep_trait", "brain_disorder", "source_pair_status", "published_rg",
    "published_se", "published_p", "study_method", "sleep_gwas_source",
    "sleep_sample_size", "disorder_gwas_source", "disorder_sample_size",
    "reported_result", "numeric_pair_estimate_reviewed", "source_overlap",
    "independent_replication", "interpretation",
]


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def read_master() -> dict[str, dict[str, str]]:
    with MASTER.open(encoding="utf-8", newline="") as stream:
        rows = list(csv.DictReader(stream, delimiter="\t"))
    return {row["trait"]: row for row in rows}


def build() -> tuple[dict[str, str], dict[str, object]]:
    sources = read_master()
    sleep = sources["longsleep"]
    pd = sources["parkinson"]
    if (sleep.get("PMID") != "30846698" or sleep.get("cases") != "34184" or
            sleep.get("controls") != "305742" or pd.get("PMID") != "31701892" or
            pd.get("cases") != "33674" or pd.get("controls") != "449056"):
        raise ValueError("Locked Brain6 long-sleep or Parkinson source metadata changed")
    article_xml = ARTICLE_XML.read_text(encoding="utf-8")
    file_record = re.search(
        r'<media\b[^>]*xlink:href="Table_1\.xlsx"[^>]*>(.*?)</media>',
        article_xml, re.DOTALL)
    file_metadata = file_record.group(1) if file_record else ""
    file_path = re.search(r"<\?cloudpmc-path\s+([^?]+Table_1\.xlsx)\?>", file_metadata)
    file_size = re.search(r"<\?size\s+(\d+)\?>", file_metadata)
    file_mime = "application/vnd.openxmlformats-officedocument.spreadsheetml.sheet"
    if (not file_path or not file_size or 'mimetype="application"' not in article_xml or
            'mime-subtype="vnd.openxmlformats-officedocument.spreadsheetml.sheet"' not in article_xml or
            Path(file_path.group(1)).name != SUPPLEMENTARY_FILE or int(file_size.group(1)) != 722516):
        raise ValueError("Official full-text XML does not bind the expected Grover supplementary workbook")
    row = {
        "sleep_trait": "longsleep",
        "brain_disorder": "parkinson",
        "source_pair_status": "EXACT_SOURCE_PAPERS_AND_MATCHING_SAMPLE_COUNTS",
        "published_rg": "NA",
        "published_se": "NA",
        "published_p": "NA",
        "study_method": "Cross-trait LDSC; pairwise screen before MR",
        "sleep_gwas_source": "Dashti et al. 2019 (PMID 30846698; DOI 10.1038/s41467-019-08917-4)",
        "sleep_sample_size": "34184 long-sleep cases / 305742 controls",
        "disorder_gwas_source": "Nalls et al. 2019 (PMID 31701892; DOI 10.1016/S1474-4422(19)30320-5)",
        "disorder_sample_size": "33674 cases / 449056 controls",
        "reported_result": "The article states that pairwise genetic-correlation analysis found no correlation of any neurodegenerative disorder with sleep or pain-related traits; it gives no pair-specific numeric LDSC estimate in the article text.",
        "numeric_pair_estimate_reviewed": "False",
        "source_overlap": "The article cites the same long-sleep and Parkinson GWAS releases and reports matching sample counts; participant-level identity is therefore shared at the summary-statistic source level.",
        "independent_replication": "False",
        "interpretation": "Prior negative direct-RG screen only. Do not infer an exact rg of zero or use this as independent replication. The Brain6 longsleep–Parkinson novelty class remains MR_ONLY because no positive direct-RG estimate was reported in the reviewed source.",
    }
    provenance = {
        "schema_version": 1,
        "status": "PASS_PRIOR_NEGATIVE_LDSC_CONTEXT",
        "scope": "One source-level literature row for the inherited significant Brain6 longsleep–Parkinson pair.",
        "builder": str(Path(__file__).resolve().relative_to(ROOT)),
        "builder_sha256": sha256(Path(__file__).resolve()),
        "input_hashes": {str(MASTER.relative_to(ROOT)): sha256(MASTER)},
        "source_file_hashes": {str(ARTICLE_XML.relative_to(ROOT)): sha256(ARTICLE_XML)},
        "source": {
            "citation": "Grover S, Sharma M. Sleep, Pain, and Neurodegeneration: A Mendelian Randomization Study. Front Neurol. 2022;13:765321.",
            "pmid": "35585838",
            "doi": "10.3389/fneur.2022.765321",
            "url": "https://pmc.ncbi.nlm.nih.gov/articles/PMC9108392/",
            "reviewed_sections": ["methods: cross-trait LDSC", "results: pairwise genetic-correlation screen", "Table 1: GWAS source/sample metadata"],
            "article_reported": "No correlation of any NDD with sleep or pain-related traits in the pairwise genetic-correlation analysis.",
            "numeric_estimate_limitation": "No pair-specific numeric rg/SE/P is stated in the article text; the linked supplementary workbook was not inspected.",
            "numeric_rg_se_p_in_article_text": False,
            "linked_supplementary_workbook_reviewed": False,
            "linked_supplementary_file": {
                "filename": SUPPLEMENTARY_FILE,
                "mime_type": file_mime,
                "size_bytes": int(file_size.group(1)),
                "cloudpmc_path": file_path.group(1),
                "source_xml": str(ARTICLE_XML.relative_to(ROOT)),
                "binary_sha256": None,
                "binary_retrieved": False,
                "retrieval_status": "OFFICIAL_FILE_METADATA_VERIFIED_BUT_BINARY_UNAVAILABLE",
                "retrieval_attempts": [
                    {"endpoint": "Europe PMC supplementaryFiles REST endpoint", "status": "TIMED_OUT_NO_BYTES"},
                    {"endpoint": "PMC linked Table_1.xlsx asset", "status": "POW_INTERSTITIAL_NOT_XLSX"},
                    {"endpoint": "Frontiers article file route", "status": "404"},
                ],
                "numeric_pair_estimate_reviewed": False,
            },
        },
        "output": {"path": str(OUTPUT.relative_to(ROOT)), "sha256": "PENDING", "rows": 1},
    }
    return row, provenance


def main() -> None:
    row, provenance = build()
    OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    with OUTPUT.open("w", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=FIELDS, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerow(row)
    provenance["output"]["sha256"] = sha256(OUTPUT)  # type: ignore[index]
    PROVENANCE.write_text(json.dumps(provenance, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(json.dumps({"output": str(OUTPUT), "sha256": sha256(OUTPUT), "status": provenance["status"]}, indent=2))


if __name__ == "__main__":
    main()
