# FinnGen CASCADE release-metadata recheck — 2026-09-25

## Finding

FinnGen's release page and the CASCADE browser describe related but differently scoped counts. The release page reports 27,294 genes and 297,024 chromatin peaks analyzed across the released molecular resource. CASCADE v1.0.0/R12 reports 20,829 cis-eGenes and 210,584 cis-caPeaks, i.e. feature sets with cis-QTL results. Those counts should remain separate; they are not evidence that either source is missing the other features. The release page labels its 33 cell types as 9 L1 plus 24 L2, while the CASCADE browser describes L0 PBMC plus 8 L1 lineages and 24 L2 subtypes. The equal total but different level labels are another source-specific taxonomy description; no one-to-one crosswalk is exposed, so retain each source’s hierarchy labels as reported.

A version/count discrepancy remains for regulatory links. CASCADE's current browser metadata reports 593,765 peak-gene links, with 104,780 peaks linked to 12,018 genes. The 2025 preprint reports 496,488 enhancer-gene links. Both sources report 1,108 donors, 20,829 genes, 210,584 peaks, and 119,094 fine-mapped variants, but the inspected public metadata does not establish whether the link-total difference reflects an updated analysis, filtering, or terminology. Preserve both source-specific counts; do not merge or select one as the canonical count without a versioned file manifest or methods clarification.

## Source-specific record

| Source/version | Reported units and counts | Access/provenance |
|---|---|---|
| FinnGen release page, public release 2026-08-13 | 1,108 donors; 10,612,211 snRNA nuclei; 7,100,899 snATAC nuclei; 33 cell types (9 L1, 24 L2); 27,294 genes and 297,024 chromatin peaks analyzed; 17,167,198 variants | Official release/access page; browser is public, bulk summary-statistic download instructions require the official form and are emailed. |
| CASCADE browser, v1.0.0/R12 (current browser metadata) | L0 PBMC + 8 L1 lineages + 24 L2 subtypes; 20,829 cis-eGenes; 210,584 cis-caPeaks; 119,094 fine-mapped molQTL variants (PIP > 0.5); 593,765 peak-gene links (104,780 peaks × 12,018 genes) | Official interactive browser and methods page; public targeted lookup. |
| Kanai et al. medRxiv preprint v1 (2025) | 51,083 cis-eQTLs for 20,829 genes; 338,100 cis-caQTLs for 210,584 peaks; 119,094 fine-mapped variants; 496,488 enhancer-gene links | Preprint abstract/version; this older count is retained as reported and not silently replaced by the browser total. |

The FinnGen release page says the public release occurred on 2026-08-13 and links both the CASCADE browser and the Kanai et al. preprint. Its 27,294/297,024 feature totals are not the same denominator as CASCADE's cis-eGene/caPeak counts. The total of 33 cell types uses different level labels across the two pages (release page: 9 L1 + 24 L2; browser: L0 + 8 L1 + 24 L2), with no public crosswalk. Conversely, the two link totals may be more directly comparable but remain unreconciled: the public pages do not expose a crosswalk or versioned export manifest that would demonstrate equivalence.

## Relevance to this project

No eligible frailty-specific loci currently pass the locked local-sharing gate, and the running LAVA sensitivity results are incomplete. Therefore no CASCADE locus query, bulk-data request, colocalization, or mechanistic claim was made. CASCADE remains optional, post-gate functional context, not independent frailty replication. No FinnGen form was submitted and no files were downloaded.

## Sources

- [FinnGen access/results and 2026 single-nucleus multiome release](https://www.finngen.fi/en/access_results)
- [FinnGen CASCADE browser and methods](https://cascade.finngen.fi/about)
- [Kanai et al., medRxiv preprint v1, DOI 10.1101/2025.11.25.25340489](https://doi.org/10.1101/2025.11.25.25340489)
