# Chromosome 5: common-support comparison of three pairwise models

On the same **1,633 variants**, each pairwise single-signal model supports a regional shared signal under the default prior, but their **95% conditional SNP sets do not all overlap**. The testis expression–GWAS comparisons favor rs416223, while the insomnia–ADHD comparison favors a different group of linked markers. This makes the existing pairwise results insufficient to identify one variant shared by all three traits. It also does not rule out such a variant: these are model-dependent, correlated results with incomplete coverage.

| Pair | H4 at p12=10⁻⁵ | H4 at p12=10⁻⁶ | 95% SNP set conditional on H4 |
|---|---:|---:|---|
| Testis expression–insomnia | 0.9568 | 0.6888 | rs416223, rs325481, rs325500 |
| Testis expression–ADHD | 0.9825 | 0.8490 | rs416223, rs325500 |
| Insomnia–ADHD | 0.9911 | 0.9174 | rs77960, rs2431108, rs30266 |

Expression refers to RP11-6N13.1 in GTEx testis. The expression–insomnia conditional SNP distribution gives rs416223 probability 0.726; expression–ADHD gives it 0.873. The insomnia–ADHD conditional distribution assigns its three listed variants 0.330, 0.319 and 0.307, respectively. These values are **conditional on the pair's H4 hypothesis**, not standalone causal probabilities. The pairwise conditional distributions must not be multiplied into a three-trait posterior; the datasets and models share inputs.

## Reference LD context

All six variants in the union of the three sets matched by rsID, position and allele pair in the 1000 Genomes Phase 3 European reference. Within 503 reference individuals, the three expression-prioritized variants are strongly linked to one another (pairwise r² 0.988–0.996), as are the three trait–trait variants (r² 0.995–1.000). Across those two groups, r² is approximately **0.696–0.706**; for example, rs416223–rs77960 has r²=0.703. This shows substantial but incomplete tagging in a small reference sample. It does not establish independent signals or identify which variant is causal. The reference is not the actual GWAS or GTEx cohort; LD and frequencies can differ across cohorts.

## Interpretation and limitations

The 1,633 variants are the intersection of two prior QTL–GWAS matches, not all variants in the chromosome 5 region. Prior results on their original supports had 1,673 expression–insomnia and 1,633 expression–ADHD variants. Repeating all three pairwise models on this intersection controls this particular coverage difference but cannot restore missing causal variants or address multiple independent signals. Single-signal ABF models can give misleadingly precise conditional sets when these assumptions fail. Study populations, participant overlap, gene selection, window boundaries, incomplete source allele-frequency checks, effect priors and GTEx sample size remain limitations. Testis is not brain; frontal cortex had no usable rows for the gene in the preceding extraction. No AlphaGenome score was used, no experimental validation or independent replication was performed, and no result repairs the failed LAVA gate.

[Protocol](PROTOCOL.md), [all pairwise posteriors](pair_posteriors.tsv), [conditional SNP distributions](conditional_snps.tsv), [reference LD pairs](pairwise_ld.tsv), [reference receipt](reference_receipt.json), [model code](run.R) and [independent arithmetic and support verification](verification.json) document the calculation. The verifier reproduces all six H0–H4 vectors within 2.1×10⁻¹⁵ and confirms the conditional-set memberships; these checks do not validate the model. The LD query was chosen after inspecting the conditional sets, and the remote indexed reference file was not whole-file hashed.

A [multi-trait method readiness check](../hyprcoloc_readiness/RESULTS.md) found no completed HyPrColoc fit; source-build and usable disk-space issues prevented the planned exploratory run.
