#!/usr/bin/env python3
"""Build and verify the frozen Track B literature evidence audit.

The output is a curated primary-source ledger, not an automated claim that a
search engine proved absence.  Missing-method rows therefore say only that no
exact report was identified by the documented search protocol and date.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import io
import json
from pathlib import Path


SEARCH_DATE = "2026-09-01"
PAIR_MANIFEST = Path("results/track_b/pair_manifest.tsv")
PAIR_LOCK = Path("results/track_b/pair_manifest.lock.json")
EVIDENCE = Path("results/track_b/02_prior_evidence_map.tsv")
SEARCHES = Path("results/track_b/02_literature_search_queries.tsv")
REPORT = Path("results/track_b/02_LITERATURE_AUDIT.md")
LOCK = Path("results/track_b/02_literature_audit.lock.json")

ALLOWED_CLASSES = {
    "DIRECT_PRIOR_GLOBAL",
    "DIRECT_PRIOR_LOCAL",
    "RELATED_ONLY",
    "PRIOR_SHARED_LOCUS",
    "PRIOR_CELL_MECHANISM",
    "UNDERREPORTED",
    "UNCERTAIN",
}

EVIDENCE_FIELDS = [
    "pair_id", "sleep_trait", "external_trait", "evidence_id", "classification",
    "evidence_layer", "exact_pair_status", "sleep_phenotype_match",
    "external_phenotype_match", "cohort_relation_to_atlas", "method",
    "direct_rg_evidence", "local_evidence", "shared_locus_evidence",
    "cell_evidence", "directional_evidence", "publication_year", "title",
    "PMID", "DOI", "source_url", "primary_finding", "limitations",
    "source_verification", "audit_date",
]

SEARCH_FIELDS = [
    "pair_id", "sleep_trait", "external_trait", "method_family", "search_query",
    "phenotype_synonyms", "search_date", "search_protocol", "audit_result",
    "supporting_evidence_ids", "interpretation",
]

METHOD_FAMILIES = [
    "genetic_correlation",
    "LDSC",
    "local_genetic_correlation",
    "LAVA",
    "HESS_rho_HESS",
    "HDL_L",
    "pleiotropy",
    "PLACO",
    "conjunction_FDR",
    "fine_mapping",
    "colocalization",
    "TWAS",
    "eQTL",
    "sQTL",
    "cell_type_single_cell_single_nucleus",
    "ATAC_brain_tissue",
    "Mendelian_randomization",
    "spatial_transcriptomics",
]

PAIR_META = {
    "A": {
        "sleep": "snoring",
        "external": "parental_lifespan",
        "query_sleep": 'snoring OR habitual snoring OR sleep-disordered breathing',
        "query_external": '"parental lifespan" OR "parents age at death" OR longevity',
        "synonyms": "snoring|habitual snoring|sleep-disordered breathing;parental lifespan|parents age at death|longevity",
    },
    "B": {
        "sleep": "insomnia",
        "external": "adhd",
        "query_sleep": 'insomnia OR insomnia symptoms OR sleep initiation and maintenance disorder',
        "query_external": 'ADHD OR attention-deficit/hyperactivity disorder OR attention deficit hyperactivity disorder',
        "synonyms": "insomnia|insomnia symptoms|sleep initiation and maintenance disorder;ADHD|attention-deficit/hyperactivity disorder",
    },
    "CONTROL": {
        "sleep": "insomnia",
        "external": "frailty",
        "query_sleep": 'insomnia OR insomnia symptoms OR sleep initiation and maintenance disorder',
        "query_external": 'frailty OR frailty index OR frail',
        "synonyms": "insomnia|insomnia symptoms|sleep initiation and maintenance disorder;frailty|frailty index",
    },
}


def row(pair: str, evidence_id: str, classification: str, layer: str,
        exact: str, sleep_match: str, external_match: str, cohort: str,
        method: str, rg: str, local: str, shared: str, cell: str,
        directional: str, year: str, title: str, pmid: str, doi: str,
        url: str, finding: str, limitations: str,
        verification: str = "PRIMARY_SOURCE_ABSTRACT_OR_FULL_TEXT_REVIEWED") -> dict[str, str]:
    meta = PAIR_META[pair]
    return {
        "pair_id": pair,
        "sleep_trait": meta["sleep"],
        "external_trait": meta["external"],
        "evidence_id": evidence_id,
        "classification": classification,
        "evidence_layer": layer,
        "exact_pair_status": exact,
        "sleep_phenotype_match": sleep_match,
        "external_phenotype_match": external_match,
        "cohort_relation_to_atlas": cohort,
        "method": method,
        "direct_rg_evidence": rg,
        "local_evidence": local,
        "shared_locus_evidence": shared,
        "cell_evidence": cell,
        "directional_evidence": directional,
        "publication_year": year,
        "title": title,
        "PMID": pmid,
        "DOI": doi,
        "source_url": url,
        "primary_finding": finding,
        "limitations": limitations,
        "source_verification": verification,
        "audit_date": SEARCH_DATE,
    }


def evidence_rows() -> list[dict[str, str]]:
    rows = [
        row("A", "A_CAMPOS_2020_GLOBAL", "DIRECT_PRIOR_GLOBAL", "GLOBAL",
            "EXACT_SLEEP_RELATED_EXTERNAL_SOURCE_DEFINITION", "EXACT", "RELATED_SOURCE_DEFINITION",
            "SAME_OR_OVERLAPPING_SNORING_GWAS;EXTERNAL_DATASET_DIFFERENT_OR_UNRESOLVED",
            "LDSC", "rg=-0.1279;SE=0.0633;P=0.0432", "NONE", "NONE", "NONE", "NONE",
            "2020", "Insights into the aetiology of snoring from observational and genetic investigations in the UK Biobank",
            "32060260", "10.1038/s41467-020-14625-1", "https://pubmed.ncbi.nlm.nih.gov/32060260/",
            "Supplementary Data 2 directly reported a negative snoring correlation with parents age at death.",
            "Not independent of the atlas snoring GWAS; external source definition and cohort cannot be treated as exact independent replication; nominal P=0.0432."),
        row("A", "A_CAMPOS_2020_TRAIT_ANNOTATION", "RELATED_ONLY", "MOLECULAR_TISSUE",
            "SLEEP_TRAIT_ONLY_NOT_PAIR_MECHANISM", "EXACT", "NOT_TESTED_AS_PAIR",
            "SAME_OR_OVERLAPPING_SNORING_GWAS", "snoring GWAS functional annotation/eQTL/tissue review",
            "NONE", "NONE", "NONE", "TRAIT_LEVEL_ONLY", "NONE", "2020",
            "Insights into the aetiology of snoring from observational and genetic investigations in the UK Biobank",
            "32060260", "10.1038/s41467-020-14625-1", "https://pubmed.ncbi.nlm.nih.gov/32060260/",
            "The snoring paper prioritized loci and regulatory annotations across several tissues.",
            "Trait-only annotation cannot establish a snoring-parental-lifespan locus, shared signal, gene, tissue, or cell mechanism."),
        row("A", "A_TIMMERS_2019_LIFESPAN", "RELATED_ONLY", "EXTERNAL_GWAS",
            "EXTERNAL_TRAIT_ONLY", "NOT_TESTED", "EXACT_ATLAS_EXTERNAL_SOURCE",
            "SAME_OR_OVERLAPPING_PARENTAL_LIFESPAN_GWAS", "parental lifespan GWAS",
            "NONE", "NONE", "NONE", "NONE", "NONE", "2019",
            "Genomics of 1 million parent lifespans implicates novel pathways and common diseases and distinguishes survival chances",
            "30642433", "10.7554/eLife.39856", "https://pubmed.ncbi.nlm.nih.gov/30642433/",
            "Large parental-lifespan GWAS supplies the atlas external phenotype.",
            "Discovery source, not replication and not an exact pairwise mechanistic study."),
        row("A", "A_WRIGHT_2019_ANCESTRYDNA", "RELATED_ONLY", "REPLICATION_CANDIDATE",
            "EXTERNAL_TRAIT_ONLY", "NOT_TESTED", "CLOSELY_MATCHED_PARENTAL_LIFESPAN",
            "ANCESTRYDNA_COMPONENT_INDEPENDENT_OF_UKB;META_ANALYSIS_INCLUDES_UKB", "parental lifespan GWAS",
            "NONE", "NONE", "NONE", "NONE", "NONE", "2019",
            "A Prospective Analysis of Genetic Variants Associated with Human Lifespan",
            "31484785", "10.1534/g3.119.400448", "https://pubmed.ncbi.nlm.nih.gov/31484785/",
            "AncestryDNA discovery exceeded 300,000 participants; a later meta-analysis included an independent UK Biobank dataset.",
            "The independent AncestryDNA component is not an analysis-ready public full summary-statistics release identified by this audit; the meta-analysis is not independent of UKB."),
        row("A", "A_TANAKA_2017_HRS", "RELATED_ONLY", "REPLICATION_CANDIDATE",
            "EXTERNAL_TRAIT_ONLY", "NOT_TESTED", "RELATED_DICHOTOMOUS_PARENTAL_LONGEVITY",
            "HRS_AND_OTHER_COHORTS;NO_IDENTIFIED_UKB_DISCOVERY_OVERLAP", "parental lifespan GWAS",
            "NONE", "NONE", "NONE", "NONE", "NONE", "2017",
            "Genome-wide Association Study of Parental Life Span",
            "27816938", "10.1093/gerona/glw206", "https://pubmed.ncbi.nlm.nih.gov/27816938/",
            "HRS-based parental longevity analysis provides an independently sampled, related phenotype.",
            "Small and phenotype definition differs; the top signal did not provide a well-powered exact pair replication and no public dense signed release was identified."),
        row("A", "A_EXACT_MECHANISM_SEARCH_NEGATIVE", "UNDERREPORTED", "LOCAL_TO_CELLULAR",
            "EXACT_PAIR_SEARCH_NO_REPORT_IDENTIFIED", "EXACT", "EXACT_OR_CLOSE_SYNONYMS",
            "NOT_APPLICABLE", "LAVA/rho-HESS/HDL-L/PLACO/conjFDR/fine-mapping/coloc/TWAS/cell-QTL/scATAC/spatial search",
            "NONE_BEYOND_CAMPOS", "NO_EXACT_REPORT_IDENTIFIED", "NO_EXACT_REPORT_IDENTIFIED",
            "NO_EXACT_REPORT_IDENTIFIED", "NO_EXACT_REPORT_IDENTIFIED", "NA",
            "No exact mechanistic paper identified in targeted search", "NA", "NA", "NA",
            "No exact snoring-parental-lifespan local, signal-aware, cellular, chromatin, or spatial analysis was identified.",
            "A targeted search cannot prove absence; indexing, terminology, and unpublished work remain limitations.",
            "TARGETED_SEARCH_WITH_PRIMARY_SOURCE_FOLLOW_UP"),

        row("B", "B_CARPENA_2021_GLOBAL", "DIRECT_PRIOR_GLOBAL", "GLOBAL",
            "EXACT_PAIR", "EXACT", "EXACT", "LIKELY_SAME_OR_PREDECESSOR_PUBLIC_GWAS_FAMILIES",
            "cross-trait LDSR", "POSITIVE_SIGNIFICANT_RG_REPORTED", "NONE", "NONE", "NONE", "NONE", "2021",
            "Sleep-related traits and attention-deficit/hyperactivity disorder comorbidity: Shared genetic risk factors, molecular mechanisms, and causal effects",
            "33821771", "10.1080/15622975.2021.1907719", "https://pubmed.ncbi.nlm.nih.gov/33821771/",
            "Directly reported a positive genomic correlation between insomnia and ADHD.",
            "Uses earlier/largely overlapping public GWAS families and is prior evidence, not independent Track B replication."),
        row("B", "B_CARPENA_2021_GENES_MR", "PRIOR_SHARED_LOCUS", "GENE_DIRECTIONAL",
            "EXACT_PAIR", "EXACT", "EXACT", "LIKELY_SAME_OR_PREDECESSOR_PUBLIC_GWAS_FAMILIES",
            "gene-based cross-trait meta-analysis; enrichment; two-sample MR", "YES_IN_SAME_PAPER", "NONE",
            "SHARED_GENES_NOT_SIGNAL_AWARE_COLOCALIZATION", "NO_CELL_SPECIFIC_TEST", "BIDIRECTIONAL_COMPONENTS_REPORTED", "2021",
            "Sleep-related traits and attention-deficit/hyperactivity disorder comorbidity: Shared genetic risk factors, molecular mechanisms, and causal effects",
            "33821771", "10.1080/15622975.2021.1907719", "https://pubmed.ncbi.nlm.nih.gov/33821771/",
            "Reported shared genes/pathways and MR evidence for sleep traits and ADHD.",
            "Gene-based overlap is not fine-mapping or trait-trait colocalization; MR is vulnerable to pleiotropy and does not prove mediation."),
        row("B", "B_DEMONTIS_2023_GLOBAL", "DIRECT_PRIOR_GLOBAL", "GLOBAL",
            "EXACT_PAIR_WITHIN_ADHD_GWAS", "EXACT", "EXACT", "SAME_ADHD_GWAS_AS_ATLAS;INSOMNIA_SOURCE_SAME_OR_OVERLAPPING",
            "LDSC in ADHD GWAS follow-up", "POSITIVE_RG_REPORTED", "NONE", "NONE", "NONE", "NONE", "2023",
            "Genome-wide analyses of attention deficit hyperactivity disorder identify 27 risk loci, refine the genetic architecture, and implicate several cognitive domains",
            "36702997", "10.1038/s41588-022-01285-8", "https://pubmed.ncbi.nlm.nih.gov/36702997/",
            "The ADHD discovery paper included insomnia among genetically correlated phenotypes.",
            "It uses the exact atlas ADHD source and therefore cannot provide independent replication."),
        row("B", "B_XUE_2026_GLOBAL", "DIRECT_PRIOR_GLOBAL", "GLOBAL",
            "EXACT_PAIR_IN_MULTI_DISORDER_STUDY", "EXACT", "EXACT", "SAME_OR_OVERLAPPING_PUBLIC_GWAS_FAMILIES",
            "bivariate MiXeR/genetic correlation", "SIGNIFICANT_POSITIVE_RELATION_REPORTED", "NONE", "NONE", "NONE", "NONE", "2026",
            "Detection of pleiotropic genetic factors and critical brain cell types linking insomnia with psychiatric disorders",
            "41065713", "10.1093/sleep/zsaf317", "https://pubmed.ncbi.nlm.nih.gov/41065713/",
            "Insomnia showed significant genetic correlation with ADHD and six other psychiatric disorders.",
            "Transdiagnostic design and overlapping GWAS sources limit pair-specific novelty and independence."),
        row("B", "B_XUE_2026_SHARED_LOCI", "PRIOR_SHARED_LOCUS", "PLEIOTROPY",
            "EXACT_PAIR_IN_MULTI_DISORDER_STUDY", "EXACT", "EXACT", "SAME_OR_OVERLAPPING_PUBLIC_GWAS_FAMILIES",
            "conjFDR; ASSET; MAGMA", "YES", "NONE", "ADHD_INCLUDED;SEVEN_NOVEL_ASSOCIATIONS_REPORTED", "NONE", "NONE", "2026",
            "Detection of pleiotropic genetic factors and critical brain cell types linking insomnia with psychiatric disorders",
            "41065713", "10.1093/sleep/zsaf317", "https://pubmed.ncbi.nlm.nih.gov/41065713/",
            "Across insomnia and psychiatric disorders, 70 shared loci were reported; seven novel associations were assigned to ADHD.",
            "Shared-locus methods do not by themselves establish a single shared causal variant or mediation."),
        row("B", "B_XUE_2026_CELL", "PRIOR_CELL_MECHANISM", "CELL_ENRICHMENT",
            "EXACT_PAIR_INCLUDED_BUT_CELL_RESULTS_TRANS_DIAGNOSTIC", "EXACT", "EXACT", "SAME_OR_OVERLAPPING_PUBLIC_GWAS_FAMILIES",
            "SEISMIC on 36 brain cell types", "YES", "NONE", "YES", "EIGHT_CORTICAL_NEURON_SUBTYPES", "NONE", "2026",
            "Detection of pleiotropic genetic factors and critical brain cell types linking insomnia with psychiatric disorders",
            "41065713", "10.1093/sleep/zsaf317", "https://pubmed.ncbi.nlm.nih.gov/41065713/",
            "Four GABAergic and four glutamatergic cortical neuron subtypes were implicated across the insomnia-psychiatric analysis.",
            "Abstract-level evidence does not establish that every subtype is ADHD-specific; enrichment is not cell-specific eQTL colocalization."),
        row("B", "B_XUE_2026_ADHD_CELL_SPECIFICITY", "UNCERTAIN", "CELL_SPECIFICITY",
            "EXACT_PAIR_INCLUDED", "EXACT", "EXACT", "SAME_OR_OVERLAPPING_PUBLIC_GWAS_FAMILIES",
            "reported transdiagnostic cell enrichment", "YES", "NONE", "YES", "PAIR_SPECIFIC_SUBTYPE_ASSIGNMENT_UNRESOLVED", "NONE", "2026",
            "Detection of pleiotropic genetic factors and critical brain cell types linking insomnia with psychiatric disorders",
            "41065713", "10.1093/sleep/zsaf317", "https://pubmed.ncbi.nlm.nih.gov/41065713/",
            "The study supports cortical GABA/glutamate involvement in the broader insomnia-psychiatric architecture.",
            "A headline ADHD-specific subtype claim requires inspection/replication of pair-specific results, not extrapolation from the aggregate abstract."),
        row("B", "B_ZU_2026_GLOBAL", "DIRECT_PRIOR_GLOBAL", "GLOBAL",
            "EXACT_PAIR_AS_EXPLORATORY_ADULT_EXTENSION", "EXACT_ADULT_UKB_INSOMNIA", "EXACT", "HIGH_INSOMNIA_AND_ADHD_GWAS_OVERLAP_WITH_ATLAS",
            "LDSC", "rg=0.32;SE=0.03;P<4.8e-37", "NONE", "NONE", "NONE", "NONE", "2026",
            "Dynamic relationship and pleiotropic loci of attention deficit hyperactivity disorder with sleep traits",
            "42297780", "10.1038/s41398-026-04166-4", "https://www.nature.com/articles/s41398-026-04166-4",
            "Reported strong positive ADHD-insomnia genetic correlation in an exploratory adult UKB insomnia analysis.",
            "Adult result was an adjunct to pediatric ABCD analyses and has substantial source overlap with the atlas."),
        row("B", "B_ZU_2026_PLACO", "PRIOR_SHARED_LOCUS", "PLEIOTROPY",
            "EXACT_PAIR_AS_EXPLORATORY_ADULT_EXTENSION", "EXACT_ADULT_UKB_INSOMNIA", "EXACT", "HIGH_INSOMNIA_AND_ADHD_GWAS_OVERLAP_WITH_ATLAS",
            "PLACO after Z-score decorrelation", "YES", "NONE", "14_LEADS_12_LOCI", "NONE", "NONE", "2026",
            "Dynamic relationship and pleiotropic loci of attention deficit hyperactivity disorder with sleep traits",
            "42297780", "10.1038/s41398-026-04166-4", "https://www.nature.com/articles/s41398-026-04166-4",
            "PLACO identified 14 lead SNPs across 12 insomnia-ADHD risk loci.",
            "PLACO is statistical pleiotropy; clumping and gene mapping do not establish biological mediation."),
        row("B", "B_ZU_2026_COLOC", "PRIOR_SHARED_LOCUS", "TRAIT_COLOCALIZATION",
            "EXACT_PAIR_AS_EXPLORATORY_ADULT_EXTENSION", "EXACT_ADULT_UKB_INSOMNIA", "EXACT", "HIGH_INSOMNIA_AND_ADHD_GWAS_OVERLAP_WITH_ATLAS",
            "Bayesian coloc", "YES", "NONE", "FOUR_LOCI_PP4_OVER_0.90", "NONE", "NONE", "2026",
            "Dynamic relationship and pleiotropic loci of attention deficit hyperactivity disorder with sleep traits",
            "42297780", "10.1038/s41398-026-04166-4", "https://www.nature.com/articles/s41398-026-04166-4",
            "Four loci exceeded PP4 0.90, including signals near/in FOXP2, LSAMP, PTPRD and lincRNA RP11-6N13.1.",
            "No signal-aware fine-mapping or prior-sensitivity evidence was identified in the main text; PP4 supports a shared signal but not causality."),
        row("B", "B_ZU_2026_CELL", "PRIOR_CELL_MECHANISM", "TISSUE_CELL_ENRICHMENT",
            "EXACT_PAIR_AS_EXPLORATORY_ADULT_EXTENSION", "EXACT_ADULT_UKB_INSOMNIA", "EXACT", "HIGH_INSOMNIA_AND_ADHD_GWAS_OVERLAP_WITH_ATLAS",
            "MAGMA; WebCSEA; Lake et al single-cell enrichment", "YES", "NONE", "YES", "EXCITATORY_INHIBITORY_NEURONS_AND_ONE_MICROGLIAL_CELL_REPORTED", "NONE", "2026",
            "Dynamic relationship and pleiotropic loci of attention deficit hyperactivity disorder with sleep traits",
            "42297780", "10.1038/s41398-026-04166-4", "https://www.nature.com/articles/s41398-026-04166-4",
            "Mapped genes were enriched in brain tissue and neuronal cell populations; pooled ADHD-sleep genes also implicated one microglial cell type.",
            "Gene-set enrichment is not independent cell-atlas replication or three-way GWAS-GWAS-cell-QTL colocalization."),
        row("B", "B_ABCD_2023_CHILDHOOD", "RELATED_ONLY", "DEVELOPMENTAL_RELATED_PHENOTYPE",
            "RELATED_CHILDHOOD_DIMENSIONAL_PAIR", "CHILDHOOD_INSOMNIA_DIMENSION", "DIMENSIONAL_ADHD_SYMPTOMS",
            "ABCD_INDEPENDENT_SAMPLE_BUT_SMALL_AND_PHENOTYPICALLY_DIFFERENT", "SNP covariance/genetic correlation",
            "rG_GREATER_THAN_0.84_REPORTED", "NONE", "NONE", "NONE", "NONE", "2023",
            "Decoupling Sleep and Brain Size in Childhood: An Investigation of Genetic Covariation in the Adolescent Brain Cognitive Development Study",
            "36712562", "10.1016/j.bpsgos.2021.12.011", "https://pubmed.ncbi.nlm.nih.gov/36712562/",
            "Childhood insomnia showed significant genetic covariance with ADHD/externalizing symptoms in ABCD.",
            "N approximately 4,700, dimensional childhood traits, and different estimator/phenotype prevent use as standardized independent LDSC replication."),
        row("B", "B_EXACT_CELL_QTL_SPATIAL_SEARCH_NEGATIVE", "UNDERREPORTED", "CELL_QTL_CHROMATIN_SPATIAL",
            "EXACT_PAIR_SEARCH_NO_CONVERGENT_REPORT_IDENTIFIED", "EXACT", "EXACT", "NOT_APPLICABLE",
            "cell-specific eQTL/sQTL coloc; scATAC; spatial transcriptomics search", "NONE", "NONE", "NONE",
            "NO_THREE_WAY_GWAS_GWAS_CELL_QTL_REPORT_IDENTIFIED", "NONE", "NA",
            "No exact convergent cell-QTL/chromatin/spatial report identified", "NA", "NA", "NA",
            "No exact insomnia-ADHD report was identified that joined trait-trait fine-mapping with cell-specific QTL and same-cell chromatin plus spatial validation.",
            "Targeted search cannot prove absence; recent and unpublished studies may not be indexed.",
            "TARGETED_SEARCH_WITH_PRIMARY_SOURCE_FOLLOW_UP"),

        row("CONTROL", "C_SONG_2024_GLOBAL", "DIRECT_PRIOR_GLOBAL", "GLOBAL",
            "EXACT_PAIR", "EXACT_OR_CLOSE_UKB_INSOMNIA", "EXACT_FRAILTY_INDEX", "HIGH_UKB_OVERLAP;SAME_FRAILTY_GWAS_AS_ATLAS",
            "LDSC", "rg=1.14_unconstrained;rg=0.61_constrained;P=4.63e-189_constrained", "NONE", "NONE", "NONE", "NONE", "2024",
            "Investigating the shared genetic architecture between frailty and insomnia",
            "38425786", "10.3389/fnagi.2024.1358996", "https://pmc.ncbi.nlm.nih.gov/articles/PMC10903740/",
            "Directly reported strong positive global genetic correlation between insomnia and frailty.",
            "Both GWAS contain UKB participants; unconstrained rg above one is a warning, and constrained-intercept sensitivity is not a substitute for independent replication."),
        row("CONTROL", "C_SONG_2024_LOCAL", "DIRECT_PRIOR_LOCAL", "LOCAL",
            "EXACT_PAIR", "EXACT_OR_CLOSE_UKB_INSOMNIA", "EXACT_FRAILTY_INDEX", "HIGH_UKB_OVERLAP;SAME_FRAILTY_GWAS_AS_ATLAS",
            "rho-HESS", "YES", "3p21.31;P=9.45e-06", "NONE", "NONE", "NONE", "2024",
            "Investigating the shared genetic architecture between frailty and insomnia",
            "38425786", "10.3389/fnagi.2024.1358996", "https://pmc.ncbi.nlm.nih.gov/articles/PMC10903740/",
            "A Bonferroni-significant local correlation was reported at 3p21.31.",
            "rho-HESS/1000 Genomes result should be independently reproduced with the Track B method/reference and scrutinized for local h2 and overlap."),
        row("CONTROL", "C_SONG_2024_SHARED_LOCI", "PRIOR_SHARED_LOCUS", "PLEIOTROPY",
            "EXACT_PAIR", "EXACT_OR_CLOSE_UKB_INSOMNIA", "EXACT_FRAILTY_INDEX", "HIGH_UKB_OVERLAP;SAME_FRAILTY_GWAS_AS_ATLAS",
            "MTAG; CPASSOC", "YES", "YES", "rs34290943_and_rs10865954_highlighted", "NONE", "NONE", "2024",
            "Investigating the shared genetic architecture between frailty and insomnia",
            "38425786", "10.3389/fnagi.2024.1358996", "https://pmc.ncbi.nlm.nih.gov/articles/PMC10903740/",
            "Cross-trait analyses highlighted shared risk SNPs including rs34290943 and rs10865954.",
            "Cross-trait association is not fine-mapping or signal-aware trait colocalization."),
        row("CONTROL", "C_SONG_2024_TISSUE", "RELATED_ONLY", "TISSUE",
            "EXACT_PAIR", "EXACT_OR_CLOSE_UKB_INSOMNIA", "EXACT_FRAILTY_INDEX", "HIGH_UKB_OVERLAP;SAME_FRAILTY_GWAS_AS_ATLAS",
            "S-LDSC; MAGMA; GTEx v8", "YES", "YES", "YES", "BROAD_BRAIN_REGIONS_NOT_CELL_TYPES", "NONE", "2024",
            "Investigating the shared genetic architecture between frailty and insomnia",
            "38425786", "10.3389/fnagi.2024.1358996", "https://pmc.ncbi.nlm.nih.gov/articles/PMC10903740/",
            "Jointly associated SNPs were enriched in multiple brain regions in tissue-level analyses.",
            "Broad tissue enrichment is not a cell mechanism and may reflect generic expression/annotation architecture."),
        row("CONTROL", "C_SONG_2024_SMR", "RELATED_ONLY", "MOLECULAR",
            "EXACT_PAIR", "EXACT_OR_CLOSE_UKB_INSOMNIA", "EXACT_FRAILTY_INDEX", "HIGH_UKB_OVERLAP;SAME_FRAILTY_GWAS_AS_ATLAS",
            "SMR with GTEx eQTL and HEIDI", "YES", "YES", "GENE_LEVEL_OVERLAP", "NO_CELL_SPECIFIC_QTL", "NONE", "2024",
            "Investigating the shared genetic architecture between frailty and insomnia",
            "38425786", "10.3389/fnagi.2024.1358996", "https://pmc.ncbi.nlm.nih.gov/articles/PMC10903740/",
            "SMR/HEIDI highlighted four reported functional genes in enriched tissues.",
            "Bulk-tissue SMR is not pairwise three-way colocalization and does not establish a mediator."),
        row("CONTROL", "C_SONG_2024_MR", "RELATED_ONLY", "DIRECTIONAL",
            "EXACT_PAIR", "EXACT_OR_CLOSE_UKB_INSOMNIA", "EXACT_FRAILTY_INDEX", "HIGH_UKB_OVERLAP;SAME_FRAILTY_GWAS_AS_ATLAS",
            "bidirectional two-sample MR", "YES", "YES", "YES", "NONE", "BIDIRECTIONAL_MR_REPORTED", "2024",
            "Investigating the shared genetic architecture between frailty and insomnia",
            "38425786", "10.3389/fnagi.2024.1358996", "https://pmc.ncbi.nlm.nih.gov/articles/PMC10903740/",
            "The paper reported bidirectional MR evidence.",
            "Sample overlap and horizontal pleiotropy remain plausible; MR cannot rescue contradictory local/coloc evidence."),
        row("CONTROL", "C_EXACT_DEEP_MECHANISM_SEARCH_NEGATIVE", "UNDERREPORTED", "FINE_CELL_SPATIAL",
            "EXACT_PAIR_SEARCH_NO_DEEP_REPORT_IDENTIFIED", "EXACT", "EXACT", "NOT_APPLICABLE",
            "LAVA/HDL-L/PLACO/conjFDR/signal-aware fine-mapping/coloc/cell-QTL/scATAC/spatial search",
            "NONE_BEYOND_SONG", "RHO_HESS_ONLY", "MTAG_CPASSSOC_ONLY", "NO_CELL_SPECIFIC_REPORT_IDENTIFIED", "MR_EXISTS", "NA",
            "No exact deep cell-resolved control analysis identified", "NA", "NA", "NA",
            "The known control has global/local/shared-locus/tissue/MR evidence, but no exact signal-aware and cell-resolved chain was identified.",
            "Targeted search cannot prove absence; the control is intended for concordance, not novelty.",
            "TARGETED_SEARCH_WITH_PRIMARY_SOURCE_FOLLOW_UP"),
    ]
    return rows


SUPPORT = {
    "A": {
        "genetic_correlation": "A_CAMPOS_2020_GLOBAL", "LDSC": "A_CAMPOS_2020_GLOBAL",
        "eQTL": "A_CAMPOS_2020_TRAIT_ANNOTATION", "ATAC_brain_tissue": "A_CAMPOS_2020_TRAIT_ANNOTATION",
    },
    "B": {
        "genetic_correlation": "B_CARPENA_2021_GLOBAL;B_DEMONTIS_2023_GLOBAL;B_XUE_2026_GLOBAL;B_ZU_2026_GLOBAL",
        "LDSC": "B_CARPENA_2021_GLOBAL;B_DEMONTIS_2023_GLOBAL;B_ZU_2026_GLOBAL",
        "pleiotropy": "B_XUE_2026_SHARED_LOCI;B_ZU_2026_PLACO", "PLACO": "B_ZU_2026_PLACO",
        "conjunction_FDR": "B_XUE_2026_SHARED_LOCI", "colocalization": "B_ZU_2026_COLOC",
        "eQTL": "B_ZU_2026_CELL", "cell_type_single_cell_single_nucleus": "B_XUE_2026_CELL;B_ZU_2026_CELL",
        "Mendelian_randomization": "B_CARPENA_2021_GENES_MR",
    },
    "CONTROL": {
        "genetic_correlation": "C_SONG_2024_GLOBAL", "LDSC": "C_SONG_2024_GLOBAL",
        "local_genetic_correlation": "C_SONG_2024_LOCAL", "HESS_rho_HESS": "C_SONG_2024_LOCAL",
        "pleiotropy": "C_SONG_2024_SHARED_LOCI", "eQTL": "C_SONG_2024_SMR",
        "ATAC_brain_tissue": "C_SONG_2024_TISSUE", "Mendelian_randomization": "C_SONG_2024_MR",
    },
}


def search_rows() -> list[dict[str, str]]:
    evidence_by_id = {r["evidence_id"]: r for r in evidence_rows()}
    rows: list[dict[str, str]] = []
    for pair, meta in PAIR_META.items():
        for method in METHOD_FAMILIES:
            supporting = SUPPORT.get(pair, {}).get(method, "")
            if supporting:
                source_rows = [evidence_by_id[value] for value in supporting.split(";")]
                exact = any(r["exact_pair_status"].startswith("EXACT_PAIR") for r in source_rows)
                result = "DIRECT_OR_PAIR_INCLUDED_EVIDENCE_IDENTIFIED" if exact else "RELATED_ONLY_EVIDENCE_IDENTIFIED"
                interpretation = "See cited evidence rows; method-level strength and limitations are not collapsed across sources."
            else:
                result = "NO_EXACT_METHOD_REPORT_IDENTIFIED"
                interpretation = "No exact-pair report was identified in the targeted search; this is not proof of absence."
            query_method = method.replace("_", " ")
            rows.append({
                "pair_id": pair,
                "sleep_trait": meta["sleep"],
                "external_trait": meta["external"],
                "method_family": method,
                "search_query": f'({meta["query_sleep"]}) AND ({meta["query_external"]}) AND ("{query_method}")',
                "phenotype_synonyms": meta["synonyms"],
                "search_date": SEARCH_DATE,
                "search_protocol": "FRESH_WEB_INDEX_SEARCH_THEN_PUBMED_PMC_OR_PUBLISHER_PRIMARY_SOURCE_REVIEW",
                "audit_result": result,
                "supporting_evidence_ids": supporting or "NA",
                "interpretation": interpretation,
            })
    return rows


def tsv_text(fields: list[str], rows: list[dict[str, str]]) -> str:
    buffer = io.StringIO(newline="")
    writer = csv.DictWriter(buffer, fieldnames=fields, delimiter="\t", lineterminator="\n")
    writer.writeheader()
    writer.writerows(rows)
    return buffer.getvalue()


def sha256_text(value: str) -> str:
    return hashlib.sha256(value.encode("utf-8")).hexdigest()


def atomic_text(path: Path, value: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(value, encoding="utf-8")
    temporary.replace(path)


def validate_pair_lock() -> None:
    if not PAIR_MANIFEST.is_file() or not PAIR_LOCK.is_file():
        raise SystemExit("ERROR: frozen Track B pair manifest is missing")
    rows = list(csv.DictReader(PAIR_MANIFEST.open(newline="", encoding="utf-8"), delimiter="\t"))
    observed = {r["pair_id"]: (r["sleep_trait"], r["external_trait"]) for r in rows}
    expected = {k: (v["sleep"], v["external"]) for k, v in PAIR_META.items()}
    if observed != expected:
        raise SystemExit(f"ERROR: pair manifest drifted: {observed}")
    lock = json.loads(PAIR_LOCK.read_text(encoding="utf-8"))
    if lock.get("pair_b_identity_sha256") != sha256_text("B\tinsomnia\tadhd\n"):
        raise SystemExit("ERROR: Pair B lock identity drifted")


def report_text(evidence: list[dict[str, str]]) -> str:
    counts: dict[str, int] = {}
    for r in evidence:
        counts[r["classification"]] = counts.get(r["classification"], 0) + 1
    count_lines = "\n".join(f"- `{key}`: {counts.get(key, 0)}" for key in sorted(ALLOWED_CLASSES))
    return f"""# Track B fresh literature audit

Audit date: {SEARCH_DATE}

This audit searched each frozen pair using the exact phenotypes, documented synonyms, and all method families requested in the Track B specification. Search-engine discovery was followed by inspection of PubMed, PMC, or publisher primary-source records. A negative search row means only **no exact report was identified by this protocol on this date**; it is not proof that no paper exists.

## Pair A — snoring / parental lifespan

The global relationship is not novel: Campos et al. directly reported a negative correlation between snoring and “parents age at death” (`rg=-0.1279`, `SE=0.0633`, `P=0.0432`). It uses the same or overlapping snoring GWAS and an unresolved/different parental-lifespan source, so it is prior global evidence rather than independent replication. Independent or partly independent parental-lifespan cohorts exist, but this audit did not identify public, dense, signed, closely matched summary statistics sufficient for the requested standardized LDSC replication. No exact local-to-cellular mechanistic chain was identified.

## Pair B — insomnia / ADHD

The global relationship and generic shared-genetic claim are already well established. Carpena et al. reported direct global correlation, shared gene-level signals and MR. Demontis et al. included insomnia in the ADHD GWAS correlation analysis. Xue et al. reported conjunction-FDR/ASSET loci and cortical GABAergic/glutamatergic cell enrichments across insomnia and psychiatric disorders. Zu et al. subsequently reported `rg=0.32`, 12 PLACO loci, four coloc signals with `PP4>0.90`, and brain/cell enrichment. These studies materially narrow novelty: a defensible Track B advance would require independent replication and stronger signal-aware, cell-specific QTL/chromatin convergence, not another generic pleiotropy claim.

## Positive control — insomnia / frailty

Song et al. already reported global LDSC, a rho-HESS locus at 3p21.31, MTAG/CPASSOC shared variants, broad brain-tissue enrichment, bulk GTEx SMR and bidirectional MR. The control is therefore appropriate for pipeline concordance but not novelty. This audit did not identify an exact signal-aware, independently cell-resolved mechanism.

## Evidence classification counts

{count_lines}

## Guardrails

- Literature findings do not change any downstream threshold.
- Dataset overlap is carried forward as a limitation, not silently relabeled as replication.
- Shared-locus association is not colocalization; bulk eQTL/SMR is not a cell-specific mediator.
- Cell enrichment is not evidence that a cell type is causal.
- The 2026 literature means Pair B novelty must be judged against already published PLACO, coloc, and cell-enrichment results.
"""


def build_texts() -> tuple[str, str, str, str]:
    validate_pair_lock()
    evidence = evidence_rows()
    if len({r["evidence_id"] for r in evidence}) != len(evidence):
        raise SystemExit("ERROR: duplicate evidence_id")
    if {r["classification"] for r in evidence} - ALLOWED_CLASSES:
        raise SystemExit("ERROR: invalid evidence classification")
    searches = search_rows()
    evidence_text = tsv_text(EVIDENCE_FIELDS, evidence)
    search_text = tsv_text(SEARCH_FIELDS, searches)
    report = report_text(evidence)
    lock = {
        "schema_version": 1,
        "audit_date": SEARCH_DATE,
        "pair_manifest_sha256": hashlib.sha256(PAIR_MANIFEST.read_bytes()).hexdigest(),
        "pair_b_identity_sha256": sha256_text("B\tinsomnia\tadhd\n"),
        "evidence_row_count": len(evidence),
        "search_row_count": len(searches),
        "method_families_per_pair": len(METHOD_FAMILIES),
        "prior_evidence_map_sha256": sha256_text(evidence_text),
        "search_queries_sha256": sha256_text(search_text),
        "report_sha256": sha256_text(report),
        "negative_search_interpretation": "NO_EXACT_REPORT_IDENTIFIED_IS_NOT_PROOF_OF_ABSENCE",
        "threshold_policy": "LITERATURE_DOES_NOT_CHANGE_PREDEFINED_ANALYSIS_THRESHOLDS",
    }
    lock_text = json.dumps(lock, indent=2, sort_keys=True) + "\n"
    return evidence_text, search_text, report, lock_text


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--verify", action="store_true")
    args = parser.parse_args()
    expected = build_texts()
    paths = [EVIDENCE, SEARCHES, REPORT, LOCK]
    if args.verify:
        for path, value in zip(paths, expected):
            if not path.is_file() or path.read_text(encoding="utf-8") != value:
                raise SystemExit(f"ERROR: literature audit missing or drifted: {path}")
        print(f"verified Track B literature audit: {len(evidence_rows())} evidence rows, {len(search_rows())} searches")
        return
    for path, value in zip(paths, expected):
        atomic_text(path, value)
    print(f"wrote Track B literature audit: {len(evidence_rows())} evidence rows, {len(search_rows())} searches")


if __name__ == "__main__":
    main()
