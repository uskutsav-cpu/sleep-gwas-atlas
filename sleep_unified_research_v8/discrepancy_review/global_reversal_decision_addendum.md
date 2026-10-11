# Decision on a uniform disease-direction reversal

11 October 2026. This addendum answers a finite question: may the three preserved de Lange munged inputs be globally re-expressed toward disease risk, and what can that operation establish? No fit, source edit, historical edit, or correction of released results was performed.

**Decision: YES to a separately labelled conditional reversal diagnostic; NO to an authoritative source-wide correction of all 36 released directions on the evidence currently recovered.** Keep the existing 36-row direction hold in place.

## Evidence supporting the diagnostic

The 22 eligible primary Table 1 controls span 12 chromosomes. Twelve have the published risk allele in raw `Allele1`; ten have it in raw `Allele2`. All require the same reversal after expressing the raw effect for the published risk allele. Thus the observation is not just a collection of negative raw effects or an arbitrary choice of allele order. It appears across IBD, Crohn's disease and ulcerative colitis, and source P/OR matches identify the primary comparisons. The source convention is highly plausibly uniform, but uniformity across the full files has not been established by source-owner documentation or the original conversion code.

Either of the two leading explanations—`Effect` referring to `Allele2`, or `Effect` for `Allele1` toward the opposite case/control phenotype—would imply disease-risk Z for retained `A1` equals minus the current Z. A diagnostic may test that explicit conditional statement. It must not describe the transformed body as an authoritatively corrected GWAS release.

Twenty-two nonrandom published loci cannot prove that every raw variant used the same convention. They were novel loci identified in one study, not a random sample of all variants. A mixed convention introduced during input conversion or merging remains logically possible. The single shared schema/header and uniformly reversed controls make that possibility less plausible; they do not exclude it. Reporting a probability such as 2^-22 would be invalid because independence, random sampling, and a 50% sign-error model were not established.

## Finite additional source search

The official publisher page was recovered. Its small Supplementary Table 2 workbook describes 241 known and novel loci. The actual header at `table_of_known_loci!A9:S9` contains coordinates, rsIDs, source publications, phenotype, association P values, heterogeneity P values and gene annotations. It supplies **no effect allele, signed beta, or risk OR**, so it cannot extend direction controls by itself. Supplementary Table 3's `Summary!A3:T3` provides fine-mapping/P-value/credible-set/coverage metadata, without a disease effect-direction contract. These workbooks were inspected read-only and preserved. The original source-directory listing contains the raw IBD body and a harmonised folder, with no README. This does not prove that a README never existed.

The Methods state SNPTEST v2.5 additive tests and METAL standard-error meta-analysis. SNPTEST beta uses allele_B; METAL requires explicit effect labels. These official software conventions make the source-label explanation plausible but do not reveal the original source-specific conversion command. An independent author-lab repository uses these Allele1/Allele2 columns for single-trait heritability munging; that is not direction validation because heritability is invariant to uniform sign reversal. It cannot settle this question.

The search stopped at this point. No additional broad association calculation, large download, or native fit is warranted to answer a missing source-semantic metadata question.

Sources: [de Lange primary article](https://www.nature.com/articles/ng.3760); [official Supplementary Table 2](https://media.springernature.com/original/springer-static/esm/art%3A10.1038%2Fng.3760/MediaObjects/41588_2017_BFng3760_MOESM276_ESM.xlsx); [official Supplementary Table 3](https://media.springernature.com/original/springer-static/esm/art%3A10.1038%2Fng.3760/MediaObjects/41588_2017_BFng3760_MOESM277_ESM.xlsx); [Catalog source directory](https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST004001-GCST005000/GCST004131/); [SNPTEST documentation](https://mathgen.stats.ox.ac.uk/genetics_software/snptest/snptest_v2.5.html); [METAL documentation](https://genome.sph.umich.edu/wiki/METAL_Documentation); [author-lab heritability workflow](https://github.com/andersonlab/sc_ti_atlas).

## Minimal bounded native diagnostic

If the root chooses to execute the diagnostic, use **long sleep–Crohn's disease only**. It is one of the fixed cases, has a clear historical nonzero interval and a strong opposite-sign HGA candidate. Its role is to confirm a transformation contract, not to claim an independent biological validation or estimator agreement with HGA.

1. Create a new, privately retained diagnostic copy of the **exact preserved Crohn munged body**, negating only Z while retaining SNP, A1, A2 and N, row order and precision. Preserve the original body. Explicitly label the new copy `CONDITIONAL_DISEASE_RISK_REEXPRESSION_UNCONFIRMED_SOURCE_CONVENTION`.
2. Verify all rows satisfy the transformation exactly, with unchanged keys/alleles/N and no extra/missing records. This whole-input check is necessary for this operation; it does not repeat unrelated integrity audits.
3. Run one native pair with the original long-sleep input, exact original LDSC/reference/filter/options, one worker, and existing guards. Record that this is a diagnostic; do not add it to an association testing family.
4. Compare against the captured historical full precision: consumed SNP and allele sets, N and single-trait h2 unchanged; rg, Z, genetic covariance, unrestricted cross-trait intercept and mean Z1×Z2 negated; their SE and two-sided P unchanged; interval [L,U] becomes [-U,-L]. Check available delete-block vectors with the same re-expression. Use normal numerical tolerances appropriate to the deterministic existing runner, and investigate any deviation instead of adjusting thresholds.

Passing confirms that the preserved input/runner behaves as expected under a global reversal. It cannot establish that the reversal is scientifically correct at every source variant. An all-36 rerun would repeat the same algebra and is unnecessary until a defensible source-wide contract is obtained. Once that contract is established, the remedy can document exact signed re-expression and reconcile dependent outputs; it must retain a historical/semantic distinction rather than rewrite history.

## Source-wide correction gate

A source-wide correction becomes defensible when one of the following establishes a uniform convention for each exact 20161107 body:

- Source-owner confirmation of the final effect-allele and case/control coding, including whether it is uniform across all rows; ideally with the original conversion or METAL configuration.
- An authoritative source README or source-specific preprocessing script establishing that contract, corroborated by the published controls.
- Equivalent primary source evidence that explicitly resolves the final release's coding across variants and cohorts, rather than another package copying its column names.

If such evidence confirms a mixed convention, a uniform reversal is rejected and affected rows must be repaired using the documented per-row convention or a different source. Without confirmation, hold the directions and describe the diagnostic conditionally. This is a scientific metadata gap requiring a precise clarification, not a reason to build another controller.
