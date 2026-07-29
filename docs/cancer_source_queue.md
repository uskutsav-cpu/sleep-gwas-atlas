# Cancer-source verification queue

This discovery log keeps public cancer GWAS sources and their remaining
analysis blockers explicit. A verified download is not a licence to invent a
population prevalence, mix ancestries, infer rsIDs, or mark a trait
`CURATED` before the Phase 0 audit passes.

| Trait | Official release | Verified facts | Still required before curation |
| --- | --- | --- | --- |
| `breast_cancer` | [BCAC overall breast-cancer results](https://www.ccge.medschl.cam.ac.uk/breast-cancer-association-consortium-bcac/data-data-access/summary-results/gwas-summary-associations), [Zhang et al. 2020](https://doi.org/10.1038/s41588-020-0609-2) | The official BCAC page identifies this as the European-ancestry overall analysis with 133,384 cases and 113,789 controls. Its public 4,764,280,363-byte Google Drive archive was range-validated and has local SHA-256 `4f21dfc33d52753d496cd33756beed6e78e31fd6849fbd5755c406ad2c64e9c8`. It produced 8,623,039 autosomal explicit-rsID records (`breast_cancer.txt.gz` SHA-256 `900a1801455939e4bd170f3475f9f74a6fe78aab860d932da17da24391aeee0a`). The audit counted 669,140 coordinate-only records, 1,232,711 indel/multiallelic records, one injected header, one incomplete record, and 235,914 non-autosomal rows; the extractor fails if any autosomal allele conflict occurs. The registry uses the [NCI/SEER 13.0% US female lifetime risk](https://seer.cancer.gov/statfacts/html/breast.html) as an explicit liability-conversion approximation. | **Curated with caveat:** the prevalence is US-population, lifetime-risk evidence—not a cohort-specific European BCAC prevalence—so it must accompany every liability-scale interpretation. |

The source page asks downloaders to accept its stated terms. This repository
keeps only the reproducible acquisition procedure and checksums; downloaded
summary statistics remain ignored by Git.
