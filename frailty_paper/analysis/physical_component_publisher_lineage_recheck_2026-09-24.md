# Physical-component publisher lineage recheck

**Recheck date:** 2026-09-24 UTC. This recheck compares the live publisher
article and Zenodo record with the existing component source audit. It does
not change the analysis plan or eligibility decision.

## Primary-source observations

- The GeroScience article reports an individual-level analysis of 49,530 UK
  Biobank participants with cardiovascular disease using the five named
  physical-frailty components. Its Data availability section says that the
  physical-frailty GWAS summary statistics were obtained in a previous study,
  reference 24; reference 24 is Ye et al. (2023). The article does not say
  that the five separate files on Zenodo are the files from that cited study.
- Zenodo record 14011550 identifies Jiajin Chen as creator/contact person and
  describes five UK Biobank component GWAS. Its file list names five separate
  summary-statistic files. The record does not provide a generating-study
  crosswalk from those files to Ye et al., nor per-file build, ancestry/QC,
  effect model/coding, covariates or variant/sample processing methods.
- The existing deposit README's citation to the Chen et al. article therefore
  establishes a citation relationship, not a reproducible generation chain.
  The article's general Data availability statement points to Ye et al. but
  does not identify the Zenodo component files. The two records are not
  sufficient to transfer Ye et al.'s metadata to the five files.

## Additional author-linked source search

The publisher article links to the authors' `JiajinChen/PRS-HF` GitHub
repository for PRS coefficients. The current repository landing page lists a
`PRS Heart Failure` directory and a README citing PRS-CS and the article; its
visible file listing does not identify the five Zenodo GWAS files or a
component-file generation crosswalk. Exact-name web searches for two deposited
filenames returned the Zenodo record, without an additional author-linked
generation source. This is a scoped search result, not proof that no other
source exists; it adds no file-level methods or lineage evidence.

## Decision

No lineage or eligibility field is upgraded. Keep all five acquired component
GWAS blocked before build harmonization or signed-effect analysis. Retain
UNKNOWN for file-level build, ancestry/QC and effect model/coding. Do not infer
that Ye et al. generated the files from the publisher statement, matching
component labels, the deposit citation, or sample totals. The next sufficient
evidence remains an authoritative statement tied to Zenodo record 14011550
that supplies the per-file generation chain and missing methods metadata.
The existing author request remains a draft and was not sent.

## Sources

- Springer Nature publisher article, Data availability and reference 24:
  <https://link.springer.com/article/10.1007/s11357-025-01734-2>
- Zenodo record 14011550, creator, description and five-file listing:
  <https://zenodo.org/records/14011550>
- Author-linked PRS code repository named by the publisher article:
  <https://github.com/JiajinChen/PRS-HF>
- Ye et al. 2023 UK Biobank Fried Frailty Score GWAS:
  <https://pmc.ncbi.nlm.nih.gov/articles/PMC10651618/>
- Prior operational-definition and file-method audit:
  `physical_component_candidate_source_audit_2026-09-23.md`
