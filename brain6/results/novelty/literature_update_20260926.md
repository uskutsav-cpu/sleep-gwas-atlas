# Literature novelty update — 2026-09-26

## Fan et al. (2026): broad sleep–disease genetic-correlation atlas

Fan et al. reported LD score regression analyses pairing 34 European-ancestry sleep GWAS with 113 disease/trait GWAS. The paper reports genetic correlations involving psychiatric disorders, including depression, schizophrenia, and bipolar disorder, and also includes Alzheimer’s and Parkinson’s disease in its broader disease analyses. It was published in *Communications Medicine* on 20 May 2026 (version of record 11 August 2026; DOI 10.1038/s43856-026-01656-w).

This is prior context, not independent replication: the article reused previously published GWAS, including source data that overlap Brain6 inputs, while exact participant-level overlap is not established by the article-level summary.

I matched Supplementary Data 11 against the locked Brain6 map using the sleep-study metadata in Data 2 and disease-study metadata in Data 4. Twelve of 35 inherited significant Brain6 pairs have exact sleep-GWAS/disease-GWAS source matches; all 12 directions agree and all 12 pass Fan et al.’s reported FDR 0.05 threshold. Two additional chronotype pairs have four possible rows because the Jones et al. 2019a/2019b release suffixes cannot be assigned unambiguously; all four candidate directions agree and three pass the article FDR threshold. The other 21 Brain6 pairs have no exact source-pair row in this table. This does not establish that related published correlations do not exist.

The corrected audit normalizes repeated whitespace in publisher labels. The earlier 11-pair output is preserved as historical v1; v2 adds the long-sleep–Parkinson exact match omitted because Data 11 contains a repeated space in the disease label. No Brain6 FDR family or inherited novelty class changed, and no first-ever claim is supported. See `fan_2026_source_match_context_v2.tsv` and its provenance for matched rows and checksums. The builder is `../../scripts/audit_fan_2026_rg_context.py`. Publisher workbooks were read unmodified from the source URLs and checksums recorded in provenance; they were not copied into the repository.

**Disposition:** use as prior-context evidence in the next pair-level novelty refresh; preserve existing novelty classes pending a full current literature search.

## Grover et al. (2022): prior negative long-sleep–Parkinson LDSC context

Grover et al. report a pairwise cross-trait LDSC analysis across seven sleep/circadian traits and five neurodegenerative conditions. Their article states that the pairwise analysis found no genetic correlation between any neurodegenerative disorder and sleep or pain-related traits. Table 1 cites the same Dashti et al. 2019 long-sleep GWAS and Nalls et al. 2019 Parkinson GWAS used by the corresponding Brain6 pair, with matching sample counts (34,184/305,742 and 33,674/449,056). The article text does not state the long-sleep–Parkinson pair's numeric `rg`, `SE`, or `P`. Official Europe PMC full-text XML identifies the linked `Table_1.xlsx` supplement as an OpenXML workbook of 722,516 bytes at PMC storage path `8fc4/9108392/54eee41c5296/Table_1.xlsx`; the workbook bytes could not be retrieved, so its contents remain unreviewed and no exact estimate or zero effect is inferred.

This is prior negative context from the same source GWAS, not independent replication. It does not change the inherited Brain6 estimate or FDR. The Brain6 novelty crosswalk's `MR_ONLY` class remains scoped to positive published direct-rg findings; Supplementary Table S35 now records this non-significant direct-LDSC context separately.

Source: Grover S, Sharma M. *Sleep, Pain, and Neurodegeneration: A Mendelian Randomization Study*. Front Neurol. 2022;13:765321. [PMC full text](https://pmc.ncbi.nlm.nih.gov/articles/PMC9108392/); DOI 10.3389/fneur.2022.765321.

## The MULTI Consortium / SleepChart (2026): long-sleep LDSC context

The SleepChart paper reports rounded positive LDSC estimates for long sleep with ADHD (`rg=0.28`, `SE=0.04`, `P=2.24e-12`), MDD (`0.29`, `0.04`, `2.57e-11`), schizophrenia (`0.28`, `0.03`, `3.47e-16`), and bipolar disorder (`0.21`, `0.03`, `1.09e-7`). All four directions agree with the locked Brain6 global estimates. These are not exact source-pair or phenotype matches: the paper's exposure is UK Biobank field 1160, long sleep >8 h versus normal sleep 6–8 h (25,049 cases; 300,420 controls), while Brain6 uses the same UK Biobank source population with a ≥9 h versus 7–8 h contrast (34,184 cases; 305,742 controls). Because the exposure cohorts overlap and the threshold differs, this is prior context rather than independent two-trait replication. The SleepChart estimate also has fewer long-sleep cases and is not a power-upgrade candidate. Exact PGC outcome-release identity is unresolved. The values in Table S40 are rounded article-text numbers; the linked supplementary genetic-correlation workbook and summary-statistics files were not downloaded or reviewed. A 30 July 2026 author correction adds two omitted references; it does not report estimate changes.

Source: The MULTI Consortium, O'Toole CK, Song Z, et al. *Sleep chart of biological ageing clocks in middle and late life*. Nature (2026). [Publisher article](https://www.nature.com/articles/s41586-026-10524-5); DOI 10.1038/s41586-026-10524-5. [Author correction](https://www.nature.com/articles/s41586-026-10920-x); DOI 10.1038/s41586-026-10920-x.

## Source

Fan Z, Yang Y, Guo Y, et al. Multi-organ imaging and genetics show the impact of sleep patterns on the human brain and body. *Communications Medicine*. 2026;6:435. Published 20 May 2026. [Publisher article](https://www.nature.com/articles/s43856-026-01656-w).

The article states that 34 sleep GWAS and 113 disease/trait GWAS were used for LDSC and directs readers to Supplementary Data 11 for the sleep–disease genetic-correlation results.


## Xue et al. (2026): insomnia–psychiatric shared architecture

Xue et al. report significant genetic correlations between insomnia and ADHD, bipolar disorder, MDD, and schizophrenia. Their abstract reports insomnia Neff=314,149 and, across the study, 70 shared loci containing 97 candidate SNPs, with per-disorder novel-association counts of 7, 5, 15, and 19 for ADHD, bipolar disorder, MDD, and schizophrenia, respectively. These source-reported counts were not compared against Brain6 loci. The abstract does not provide pair-specific rg/SE/P or the exact insomnia GWAS release. Brain6's locked insomnia source-card counts (109,402 cases; 277,131 controls) imply Neff=313,750 under 4/(1/cases+1/controls), close to the reported value but insufficient to establish identical source data or participants. The three insomnia pairs among Brain6's 35 inherited significant rows (ADHD, bipolar disorder, MDD) already have prior direct-rg classifications; insomnia–SCZ is not significant under the inherited 396-pair family. This article adds source-level context only. It is not independent replication, a locus-level validation, or a reason to change the locked novelty or FDR labels.

Source: [publisher abstract](https://academic.oup.com/sleep/article-abstract/49/1/zsaf317/8279895); [PubMed record](https://pubmed.ncbi.nlm.nih.gov/41065713/). The source-specific supplementary coefficients and GWAS release identities were not verified in this refresh.
