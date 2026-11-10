# FOXP2-region insomnia and hospital-frailty colocalization with larger ancestry-specific references

## Result

A new **exploratory** analysis compared the MVP European-ancestry insomnia GWAS with FinnGen hospital frailty (HFRS) using public gnomAD non-Finnish European (NFE) and Finnish (FIN) signed LD references, respectively. In a narrower one-megabase GRCh37 window centered on rs12536335, 1,257 allele-matched SNPs were available in both references and both GWAS. SuSiE reported one candidate signal for each trait, both with lead rs35624276. With p1=p2=1e-4 and p12=1e-5, the signal-pair PP.H4 was **0.979** and PP.H3 was **0.0205**. With p12=1e-6 or 5e-5, PP.H4 was 0.826 or 0.996, respectively. Using total MVP N rather than the working effective N changed these results negligibly. The original rs12536335 is present in both candidate credible sets; it is not the highest-PIP variant (MVP PIP approximately 0.025; HFRS PIP approximately 0.123).

This is stronger *conditional evidence for one shared regional association signal* than the earlier reference-limited comparison, but it is **not proof** of a shared causal variant, a causal role for FOXP2, a single mechanism, or a validated individual risk score. It is a postselected locus and the priors, LD and fine-mapping model affect the posterior. HFRS is an ICD-based frailty measure, not the original FI or FFS; sleep-adjacent codes and phenotype differences matter. MVP and FinnGen are distinct named cohorts, but complete participant and relatedness non-overlap has not been verified. The analysis does not include a molecular QTL, so it does not establish a variant–gene link or cell type.

## What was checked

Both larger-reference matrices were extracted from the public gnomAD v2.1.1 common-variant BlockMatrix release; NFE was used for MVP and FIN for HFRS. They are ancestry-oriented external references, **not in-sample** LD. Per-variant gnomAD allele numbers in this window had medians of 15,394 for NFE and 3,472 for FIN; those allele counts vary by site and are not an independently established exact LD-panel sample size. Neither full multi-gigabyte source matrix nor full source index was downloaded or checksum-verified.

The full two-megabase comparison spans an unstored matrix block and the wider retrieved matrices were not positive semidefinite (minimum eigenvalues -8.48 and -3.79); missing correlations were not filled with zero. The narrowed one-megabase window avoids the unstored block, contains both lead markers, and yields positive-semidefinite matrices to numerical tolerance (FIN minimum eigenvalue about -2.4e-14; NFE about 8.2e-6). It necessarily omits the distal HFRS signal investigated in the earlier wider analysis, so this result must not be generalized to every signal across the region.

For the 1,257 SNPs, mean absolute effect-allele-frequency differences were 0.0106 between MVP and gnomAD NFE and 0.0080 between HFRS and gnomAD FIN; no difference exceeded 0.2. For 1,199 variants returned by the separate SISu 4.2 pairwise-LD service around rs35624276, gnomAD FIN versus SISu r² had mean absolute difference 0.0104 and correlation 0.997. Those variants are correlated, not independent checks, and these comparisons cannot prove source identity, signed phase or every allele match. Z/LD inconsistency estimates were 0.00167 for MVP/NFE and 0.00216 for HFRS/FIN, much smaller than in the broad, small-reference fits, but there is no universal cutoff. Both SuSiE fits converged; convergence does not validate the model.

On a further common 1,200-SNP subset compatible with 1000 Genomes FIN99, the ancestry-oriented gnomAD fit gave PP.H4=0.980 at p12=1e-5. Fits using EUR503 or FIN99 for *both* traits on that same subset gave 0.979 and 0.984, respectively. FIN99 is nested in EUR503 and its Z/LD diagnostics remained substantially worse; these are correlated sensitivity analyses, not independent confirmations. Their results do not erase the wider-window FIN99 instability, which involved a different support and additional distal signals.

## Inputs and reproduction

The input GWAS files were checksum-verified and are described in `stage6_foxp2_independent_cohort_coloc_review.md`. The MVP file has no usable SE values; its Z scores use the sign of its reported odds ratio and the reported two-sided p value, with effective case-control N=253,769 as a working approximation and total N=408,138 as sensitivity. HFRS uses its beta/SE and N=500,737. The source field named `alt` in the MVP file was not assumed to be the genomic alternate allele; both studies were aligned to the MVP effect allele and reference LD signs were aligned to that allele.

The public Hail variant-index files were decoded using `stage6_gnomad_fin_index_extract.py` and `stage6_gnomad_nfe_index_extract.py`. Every decoded partition was checked against published partition counts and sequential matrix indices, and variants were matched by uniquely lifted position and allele pair. The custom index decoder has **not** been independently checked with Hail. LD block decoding uses the public [ldcov](https://github.com/mkanai/ldcov) implementation at commit `92b8aad6b914ace1dad65eb119bae301f9448836`; block downloads and hashes are recorded in the FIN and NFE index/LD receipts. Full-broad-window FIN scratch matrices were removed after the narrow matrices were created to conserve disk space; rerun `stage6_foxp2_gnomad_signed_ld.py left` to reconstruct the required broader source slice. `stage6_foxp2_gnomad_nfe_narrow_ld.py`, `stage6_foxp2_gnomad_narrow_validation.py`, and `stage6_foxp2_gnomad_external_qc.py` generate and check the narrow matrices.

The official susieR 0.14.2 / coloc 5.2.3 analysis is `stage6_foxp2_gnomad_ancestry_matched_coloc.R`, with outputs `stage6_foxp2_gnomad_ancestry_matched_diagnostics.csv`, `stage6_foxp2_gnomad_ancestry_matched_coloc.csv`, and the per-trait `stage6_foxp2_gnomad_*_candidate_sets.csv`. Same-support sensitivity scripts/results use `stage6_foxp2_gnomad_common1200_coloc.R` and `stage6_foxp2_gnomad_narrow_1000g_sensitivity.R`. The exact variant order is `stage6_foxp2_gnomad_narrow_variants.csv`. Further work should verify the custom Hail-index decoding with the official reader, assess source overlap and wider-window signals with in-sample or better matched LD, and test molecular colocalization independently. Existing AlphaGenome, QTL and reporter evidence should remain hypotheses until linked by such analyses.

## Wider shared-prior sensitivity

Keeping the previously fitted signals, variants, LD and single-trait priors fixed, we additionally evaluated lower shared-variant priors. This is a post hoc stress test, not a determination of which prior is correct.

| Per-variant shared prior p12 | PP.H3: distinct variants | PP.H4: shared variant |
|---:|---:|---:|
| 1e-8 | 0.951 | 0.045 |
| 1e-7 | 0.675 | 0.322 |
| 1e-6 | 0.173 | 0.826 |
| 1e-5 | 0.021 | 0.979 |

The main analysis uses p1 = p2 = 1e-4; their product is 1e-8. Under this particular low shared prior, the fitted regional evidence favors H3 over H4, while the earlier higher prior favors H4. The appropriate prior is not established by this analysis; variants and signals were selected, and the coloc priors are modeling assumptions rather than observed disease probabilities. The posterior is conditional on the model, fixed precomputed SuSiE fits and one-megabase locus, and it does not incorporate uncertainty about LD, ascertainment, study overlap or locus selection. This sensitivity materially limits statements of strong shared-causal-variant evidence. It does not establish that the signals are different or that either gene or phenotype is causal.

`Rscript work/foxp2_rare_shared_prior.R` evaluates the saved fit with the original coloc version, verifies the original p12=1e-5 output and writes the [prior-sensitivity table](stage6_foxp2_low_shared_prior_sensitivity.csv). No new participant data or genetic-score outcomes were used.
