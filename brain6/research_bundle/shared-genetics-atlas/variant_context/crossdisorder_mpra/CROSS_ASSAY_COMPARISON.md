# Comparing two public reporter assays for selected variants

The author's published Lee table has since been recovered. The [updated published cross-assay comparison](published_table/CROSS_ASSAY_RESULTS.md) replaces the provisional Lee threshold values below, and adds reference LD for rs62369151. Its published allele FDR is 0.0114, but it is not classified as an `emVar`; Siraj's row also lacks the required activity evidence. The historical reconstruction is retained for audit.

We joined the existing [Lee neural-progenitor screen](RESULTS.md) and [Siraj SK-N-SH neural-derived cell-line screen](../siraj_2026_mpra/RESULTS.md) by exact GRCh37 chromosome, position, reference and alternate allele, checking that rsIDs agree. Four of the previously ranked candidates occur in both screened tables. These are selected candidates, not a genome-wide or independent replication sample. The assays use different cell contexts and constructs; neither directly measures gene regulation at the endogenous locus.

| Variant | Lee provisional BH, two model versions | Siraj library and source allele-skew FDR | Context |
| --- | --- | --- | --- |
| rs12658451 | 0.851, 0.174 | OL33: 0.445 | Neither source value is below 0.05. |
| rs323509 | 0.492, 0.929 | OL32: 0.622 | Neither source value is below 0.05. |
| rs60271 | 0.670, 0.924 | OL32: 0.905 | Neither source value is below 0.05. |
| rs62369151 | 0.0114, 0.000344 | OL42: 0.0442 | Candidate for follow-up, with material limitations below. |

The Lee numbers are **our unvalidated reconstructions**, not published allele-effect estimates; the tested set and normalization differ across versions and publication QC has not been reproduced. They use BH across the reconstructed Lee model family. Siraj's FDR is the value supplied in the source allele-skew table and is not a correction across all of Brain6's explorations. The Lee and Siraj values are not a single combined test; no cross-study significance or replication claim follows from this table.

For **rs62369151**, the Siraj OL42 row reports log2 skew 0.933 (SE 0.323), but both allele-level Bonferroni-adjusted activity p-values are 1.0, so the row does not satisfy the stated activity criterion. It should not be labeled an author-defined expression-modulating variant on the basis of skew alone. Lee's published table reports an allele-level result for it, but mapped barcode labels did not match those found in the released raw count files in our [barcode audit](barcode_audit/RESULTS.md); the reason remains unresolved. The two studies' effect-orientation conventions have not been reconciled here, so signed effects are not presented as directional concordance.

This locus is distinct from **rs62264764 on chromosome 3**, which was not among either assay's matched design or matrix entries. The non-match does not rule out an effect. The post-selected chromosome 5 candidate is not evidence of an insomnia–ADHD causal variant, target gene or independent GWAS replication. The [joined table](cross_assay_candidate_comparison.tsv), [code](compare_public_assays.py), and [hash receipt](cross_assay_receipt.json) preserve the exact source rows. A next step is to recover the published Lee allele table or a verified construct/barcode crosswalk, establish activity and allele orientation under the authors' criteria, and test relevant variants in an appropriate endogenous context.
