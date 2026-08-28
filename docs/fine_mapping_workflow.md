# Fine-mapping and trait-trait colocalization workflow

The primary family contains every `PRIMARY_PHASE1` locus supported by both
PLACO+ and conjunction FDR in the same locked LAVA block. Loci are not ranked,
capped, or filtered by global genetic correlation. The exact family is frozen
before any SuSiE or colocalization result is read.

SuSiE-RSS 0.14.2 and coloc 5.2.3 are source-archive and version pinned in
`config/fine_mapping_sources.tsv`. The fine-mapping policy reuses the signed
European GRCh37 LAVA UK Biobank v1.1 correlation matrix. Both GWAS betas are
oriented to the reference A1 allele, ambiguous strand SNPs are excluded, and
the GWAS/LD rows and columns must be identical. HapMap3-prefiltered GWAS are
rejected because they are not dense fine-mapping inputs.

The reference downloader records both archive checksums and hashes for all 44
extracted `.info`/`.bcor` files. The preflight verifies those hashes, all 45
full-resolution inputs, the exact R package entrypoints, memory, disk, and the
presence of real cross-method loci. It currently reports blockers without
starting analysis:

```bash
python3 scripts/56_finemapping_preflight.py --report-only
```

After the upstream gates pass, the production sequence is:

```bash
python3 scripts/57_prepare_finemapping_loci.py
python3 scripts/58_materialize_finemapping_locus.py SHARED_LOCUS_ID --materialize
python3 scripts/59_run_finemapping_locus.py SHARED_LOCUS_ID --execute
python3 scripts/60_collate_finemapping.py
python3 scripts/60_collate_finemapping.py --validate-only
```

Materialization estimates quantitative-trait phenotype SD from MAF, SE, and N
using the coloc model when no authoritative SD is available. Binary traits use
their locked case fraction and total N. SuSiE uses at most ten signals, 95%
credible sets, `min_abs_corr=0.5`, 1,000 iterations, flat priors, and no
residual-variance estimation. Both models must converge. RSS-LD inconsistency,
kriging allele-switch outliers, and credible-set purity remain explicit.

Trait-trait coloc-SuSiE is evaluated at the complete locked p12 grid. Its
aggregate is intentionally named `trait_trait_colocalization.tsv`; it cannot
satisfy the final colocalization gate by itself. The final
`colocalization.tsv` must also contain completed eQTL, sQTL, and pQTL searches
and signal-level results or evidence-backed unavailable outcomes.

Fine-mapping PIPs and colocalization posteriors are model-based evidence, not
proof of a causal variant, gene, mechanism, mediation path, or direction.
