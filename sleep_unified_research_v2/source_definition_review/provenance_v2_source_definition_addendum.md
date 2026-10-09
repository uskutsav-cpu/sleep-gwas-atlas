# MVP GIA insomnia source-definition addendum

Scientific source review; retrieval date 2026-10-09 UTC. The sealed v1 evidence at `c94ca550875b5ddeaae0b7f6eb99d861c4150d16` was not edited. This addendum contains no new GWAS body, pair result, native GWAS analysis, significance threshold, or manuscript.

**Finding:** public primary evidence resolves the objection that the documented insomnia case/control rules belong only to HARE. The final Science supplement describes GIA and the clinical coding rules; the official CIPHER insomnia record supplies the exact code list and explicitly states that PheCode exclusions were not applied. The exact GIA EUR analysis and GWAS Catalog accession share the trait identifier and counts. The source is recorded clinical insomnia; it remains a related-phenotype transport test for the core UKB frequent self-report insomnia trait.

## Exact accession and methods lineage

| Evidence | Verified identity or role |
| --- | --- |
| [Catalog GCST90475826](https://www.ebi.ac.uk/gwas/rest/api/v2/studies/GCST90475826) | Insomnia (PheCode 327.4), MVP, European ancestry, PMID 39024449; 78,566 cases and 329,572 controls; full summary statistics available. |
| [Catalog source metadata](https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90475001-GCST90476000/GCST90475826/GCST90475826.tsv.gz-meta.yaml) | N=408,138; GRCh38, one-based positions, GWAS-SSF v1.0; source MD5 `2e9c12624b653aac527444fc426e037b`; not harmonised or sorted. Metadata modified 2025-04-28. No body was downloaded in this review. |
| [dbGaP pha010373.1](https://www.ncbi.nlm.nih.gov/projects/gap/cgi-bin/analysis.cgi?study_id=phs002453.v1.p1&pha=10373) | Phe_327_4.EUR.GIA, Insomnia, European MVP participants; N=408,138, cases=78,566, controls=329,572. |
| [dbGaP dictionary phd008759.1](https://www.ncbi.nlm.nih.gov/projects/gap/cgi-bin/document.cgi?phd=8759&study_id=phs002453.v1.p1) | Its GIA Phe_327_4 row gives those exact EUR counts. The separate HARE row gives EUR 80,037 cases and 334,480 controls, so the two releases should not be conflated. |
| [Final Science supplemental methods](https://pmc-oa-opendata.s3.amazonaws.com/PMC12857194.1/NIHMS2091849-supplement-Supplemental_Text.pdf) | PDF pages 2–3 describe GIA population assignment, the September 2020 EHR horizon, CIPHER phenotype sourcing, and binary ICD ascertainment. |
| [Official CIPHER insomnia metadata](https://phenomics.va.ornl.gov/web/api/phenotype/14926) / [human-facing record](https://phenomics.va.ornl.gov/web/cipher/phenotype-viewer/details?id=14926) | ID 14926, Insomnia (gwPheWAS), keyword Phe_327_4, MVP enrollees; exact code list and algorithm rules. |

The linkage is a documented resource/trait/count chain, not a claim that CIPHER itself contains the GCST accession or participant identifiers. Catalog supplies the final-paper PMID and exact EUR counts; dbGaP supplies the exact GIA analysis name and matching counts; the final methods identify CIPHER as the initial phenotype source; CIPHER identifies its insomnia definition as used in the MVP phenome-wide GWAS and supplies Phe_327_4. No participant-level reconstruction or independent validation of the phenotyping code was performed.

The [official README](https://www.ncbi.nlm.nih.gov/projects/gap/cgi-bin/GetPdf.cgi?id=phd008760.1) does label its detailed appendix “Methods gwPheWAS [HARE]”; its EHR horizon is September 2019. That appendix alone cannot certify GIA. However, the [2023 primary preprint](https://pmc.ncbi.nlm.nih.gov/articles/PMC10327290/) already describes genetically inferred population assignment and the same binary-code rules. The recovered final supplement now provides the stronger final-publication evidence. Thus the fresh review's HARE-only limitation was a real access limitation, not evidence that no public GIA definition exists.

## Clinical ascertainment and exact code list

The final GIA supplement and CIPHER record define cases by at least two instances of mapped ICD codes and controls by zero instances. The supplement describes EHR follow-up through September 2020. CIPHER explicitly states “Phecode exclusions were not applied” and that MVP manually mapped ICD-9 to ICD-10. These facts exclude borrowing the standard PheWAS package's default related-phenotype exclusions. [Final methods, PDF pages 2–3](https://pmc-oa-opendata.s3.amazonaws.com/PMC12857194.1/NIHMS2091849-supplement-Supplemental_Text.pdf), [CIPHER insomnia algorithm](https://phenomics.va.ornl.gov/web/api/phenotype/14926).

| Vocabulary | Codes in the official CIPHER insomnia definition |
| --- | --- |
| ICD-9-CM | 307.41; 307.42; 327.0; 327.00; 327.01; 327.02; 327.09; 780.52 |
| ICD-10-CM | F51.01; F51.02; F51.03; F51.04; F51.05; F51.09; G47.0; G47.00; G47.01; G47.09 |

The CIPHER algorithm creation date is 2020-09-01. The retrieved record is revision 3, `versionInfo=1.2`, modified 2025-05-13; the revision comment adds KESER/VADC metadata and `majorRevision=false`. Do not interpret `versionInfo=1.2` as proof of the exact upstream PheWAS map software version. No executable phenotype code is attached, and the record marks the algorithm as not validated.

The recorded codes include transient, persistent, primary, psychophysiological, medical-condition-related, and mental-disorder-related insomnia categories. The public rule supplies no chronicity, symptom-frequency, daytime-impairment, or questionnaire requirement. It therefore does not certify exact equivalence to the core UKB self-report question. Codes in this PheCode group also must not be confused with the numerically identical ICD-9 code 327.4: the pinned reference PheWAS map places ICD-9 327.4 in PheCode 327.5 (parasomnia).

## Verified rules and remaining implementation uncertainty

| Field | Evidence status | Interpretation |
| --- | --- | --- |
| Clinical ascertainment | VERIFIED_PRIMARY_GIA_AND_CIPHER | Structured EHR ICD-based binary insomnia, not a survey phenotype. |
| Case occurrence threshold | VERIFIED_PRIMARY_GIA_AND_CIPHER | At least two mapped ICD instances. |
| Control occurrence threshold | VERIFIED_PRIMARY_GIA_AND_CIPHER | Zero mapped target-code instances. |
| One-instance participants | INFERRED_FROM_EXPLICIT_RULE | One instance satisfies neither stated group; exclusion follows logically. No executable implementation was available to audit. |
| Related-PheCode control exclusions | VERIFIED_CIPHER_EXPLICIT_NONE | PheCode exclusions were not applied. Do not import standard package default exclusions. |
| EHR end horizon | VERIFIED_FINAL_PUBLICATION | September 2020; the HARE-only README instead states September 2019. |
| Start date, minimum lookback or encounter density | UNRESOLVED_NOT_STATED | CIPHER's `dataUsedStart` and `dataUsedEnd` are null; no source specifies a minimum observation duration. |
| Distinct dates versus repeated codes on one date | UNRESOLVED_NOT_STATED | The source says instances/codes, not distinct days. Do not rewrite the criterion as two separate visits or two dates. |
| Minimum interval between code instances | UNRESOLVED_NOT_STATED | No stated interval; no assumed chronicity criterion. |
| Parent/child rollup implementation and historical mapping artifact | UNRESOLVED_IMPLEMENTATION | Exact included ICD list is public; executable extraction/rollup code and a release-pinned historical map were not retrieved. |
| Chart adjudication | NO_VALIDATION_REPORTED | CIPHER marks validation false and provides no adjudication result. This is an ascertainment property, not a claim that diagnoses are clinically invalid. |

The exact official code list happens to match the parent-plus-child insomnia list in the pinned [PheWAS maintainer map](https://github.com/PheWAS/PheWAS/tree/55dd1c24e228851922400cfba8d7db474565ccc7). That is a reference comparison, not the source of the MVP code list or proof that MVP used the package defaults. The package permits configurable aggregation and defaults to related-PheCode exclusions; the CIPHER algorithm explicitly says it did not apply those exclusions.

## Effect-allele clarification and a separate QC discrepancy

The [2025 author response](https://zenodo.org/records/15595560), PDF pages 1–2, says dbGaP REF/ALT denote tested alleles and may differ from genomic-reference REF/ALT; odds ratios/effects correspond to the explicit effect allele and its frequency. Some PheWeb displays orient effects to the minor allele. This supports matching the pair of tested alleles and using the supplied effect allele. It does not independently verify the EBI conversion, justify assuming ALT is the effect allele, or resolve the source `standard_error` scale. The authors' explanation is distinguished from an independent audit; this review performed no new variant-body comparisons.

The exact dbGaP analysis metadata reports filtering MAC<30; final supplement page 4 reports association inclusion MAC>40, while its initial genotype QC describes MAC<20 removal. These may concern different stages or deposited releases, but the public statements do not establish the explanation. Preserve the accession-specific metadata and report this discrepancy; do not silently normalize the thresholds. It is separate from the clinical phenotype definition.

## Decision and access receipt

**The prior source-definition stop can be resolved to the level of the documented public GIA phenotype.** Clinical ascertainment, the exact included codes, occurrence thresholds, the lack of PheCode exclusions, and the EHR end horizon are now supported by direct primary sources linked to the exact EUR accession. An independent MVP clinical-insomnia/FinnGen-outcome transport test remains a scientific choice to admit explicitly under the locked candidate family; this addendum authorizes no new estimator, threshold, candidate selection, or silent phenotype replacement.

Keep the implementation uncertainties above in any contract. If the locked gate requires distinct-day semantics, a minimum lookback, or an executable release-pinned phenotyping implementation, those additional requirements remain unresolved and must not be marked passed. Public documentation does not justify calling this exact self-report insomnia replication or certifying zero individual overlap across cohort sources.

The supplement access failure was resolved using the current [official PMC Cloud Service](https://pmc.ncbi.nlm.nih.gov/tools/pmcaws/), which provides anonymous HTTPS/S3 retrieval and replaced legacy article dataset services in August 2026. The public bucket listing identified the exact methods asset; no challenge was solved, authenticated endpoint used, investigator contacted, or restricted request made. Only the 2,587,868-byte methods PDF and small metadata were fetched from that bucket; no result spreadsheets were downloaded. Its MD5 matches official article JSON `afe2d8e36ab1a5213f3ec32a14c19a8a`; SHA-256 is `9b21e23c5a9a6e4df79df59ca3ac01ef1fdbf56c9daf868c48981e2e864edadb`. PDF pages 2–3 were rendered and visually inspected.

Earlier failures remain in the access receipt: PMC supplement returned HTTP 200 HTML rather than PDF; publisher and Europe PMC individual supplement routes returned 403; BioC list returned 429; the bounded Europe PMC supplement archive HEAD/range requests timed out with no bytes; obsolete OA endpoints returned 404; the reconstructed CDN media URL returned 404. An initial local write failed during disk exhaustion. Further evidence was routed into the newly authorized task-owned SSD namespace. These failures are not evidence that the public source is unavailable.

Companion artifacts: `provenance_v2_definition_evidence.tsv`, `provenance_v2_mvp_insomnia_icd_mapping.tsv`, `provenance_v2_definition_contract.json`, `provenance_v2_source_access_receipt.json`, and `provenance_v2_SHA256SUMS`. Source caches are confined to `/Volumes/Extreme SSD/sleep-unified-research-v1/recovery-2026-10-09/source_definition_review_v2/`; original SSD files were not modified.
