# Manuscript bibliography audit — 2026-09-23

The manuscript citation-key audit found 26 distinct cited keys and 26 bibliography entries, with no missing keys, unused entries, or duplicate BibTeX keys. The audit covers `paper/manuscript.md` and `paper/references.bib`. It is now reproducible with `make -C frailty_paper audit-bibliography`; the hash-linked machine-readable report is `analysis/manuscript_bibliography_audit_2026-09-23.json`, and focused regression coverage is in `tests/test_bibliography_audit.py`.

Foote et al. (2025), “Uncovering the multivariate genetic architecture of frailty with genomic structural equation modeling,” was verified against the journal record and indexed bibliographic metadata. Its BibTeX entry now includes *Nature Genetics* volume 57, issue 8, pages 1848–1859, and DOI `10.1038/s41588-025-02269-0`. [Nature Genetics article](https://www.nature.com/articles/s41588-025-02269-0); [PubMed record](https://pubmed.ncbi.nlm.nih.gov/40759756/).

Mak et al. (2025) title, journal, date, and DOI were checked against the UK Biobank publication record; its existing volume and page metadata were retained. [UK Biobank publication record](https://www.ukbiobank.ac.uk/publications/large-scale-genome-wide-analyses-with-proteomics-integration-reveal-novel-loci-and-biological-insights-into-frailty/).

This is a citation-integrity check only. It does not assess whether the cited literature is exhaustive or whether the manuscript is ready for submission.

Huang et al. (2026), the targeted publisher follow-up on OSA–frailty MR, was added as prior-work context; its bibliographic metadata was checked against the publisher record. This record remains unscreened and outside the frozen PRISMA counts. [BMC Pulmonary Medicine article](https://link.springer.com/article/10.1186/s12890-026-04178-2).

Che et al. (2025), PMID 39936043, was verified against its PubMed record and publisher citation metadata and added as prior MR evidence across FI, Fried frailty and sleep traits. The captured queue row remains unscreened. [Publisher record](https://www.dovepress.com/bidirectional-causal-associations-between-frailty-measures-and-sleep-d-peer-reviewed-fulltext-article-NSS); [PubMed record](https://pubmed.ncbi.nlm.nih.gov/39936043/).

The expanded focused prior-genetic audit added the 2024 Li and Lu papers, the 2025 Zhang and Gao papers, and the 2026 Zhang paper identified within the frozen genetic-search corpus. Bibliographic metadata was cross-checked against their PubMed records; all associated screening rows remain unscreened.

The regenerated report records 27 citation uses across 26 unique keys; all 26 match one bibliography entry. Missing, unused and duplicate key lists are empty. The JSON report SHA-256 is `27ec003e4a6a271332b573ac2ac55130660ba4c1361d6dc49aa5852e703dc8f6`; the checker SHA-256 is `136612bcc3dce0b016e97b45a530bfb3c52751f95b0e1601287c34efec7cb842`. Its three focused regression tests pass, and the full pinned Python 3.11.11 suite passes 77/77 (log SHA-256 `5bcc7b1c130f577583f4eec6116c8e7bf2e26afbd1074e9f2066306db6ed8c2d`), and the project Python 3.13.13 suite passes 77/77 (log SHA-256 `4d1888149e7dea6f33d1ab51fa532e30e2cedf8594a25c3958e67e16b5ce755f`).
