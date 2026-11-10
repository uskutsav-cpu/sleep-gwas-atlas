# SMAD5 brain-eQTL frequency discrepancies: reference check

Four matched variants in the original cis-window Mak FI comparison had absolute QTL–GWAS allele-frequency differences above 0.1. The colocalization analysis already excluded these four under its frequency-difference filter. We checked their reported GRCh38 alternate allele against exact two-allele records at the lifted GRCh37 positions in 503 European-labeled 1000 Genomes phase 3 participants:

| Variant | GRCh38 alternate allele | Brain QTL non-reference AF | Mak FI aligned AF | Reference-panel AF |
|---|---|---:|---:|---:|
| rs4976271 | T | 0.238 | 0.602 | 0.585 |
| rs56345521 | C | 0.535 | 0.822 | 0.832 |
| rs116687508 | A | 0.278 | 0.794 | 0.786 |
| rs113343755 | T | 0.016 | 0.183 | 0.223 |

The external reference is closer to the Mak frequency than to the reported QTL frequency for all four. This makes source-level investigation worthwhile but **does not establish which dataset is wrong**, whether the QTL field is correctly interpreted, or whether the reference matches the QTL population. Different ancestry, sample composition, imputation and representation can affect frequencies; some QTL values are not simply the complement of the reference frequency. Exact position and allele-pair matches are useful checks, not proof of strand or effect-direction harmonization. These variants were excluded by the existing >0.1 filter, so this check does not change the reported filtered coloc result or establish a shared causal signal.

The comparison used a locally saved GWAS/QTL alignment table and the public indexed reference VCF; the latter was accessed by region and not fully checksum verified. The [receipt](stage6_niagads_smad5_frequency_discrepancy.json) records matching alleles, sample counts, inputs and hashes. Reproduce with `python3 work/niagads_smad5_frequency_discrepancy_check.py`. Investigating the QTL provider's exact `non_ref_af` definition, allele representation, ancestry and per-variant QC remains necessary before reintroducing these records. This check does not evaluate a genetic score or participant outcomes.

## Source-record representation follow-up

We inspected the QTL archive's `original_variant_id` for all **4,617** matched cis-window Mak FI rows. Simple common-suffix and common-prefix trimming produces the reported single-base contrast for every row. Twenty-five original records use longer allele strings, including two of the four frequency-discrepant records. This check supports internal string consistency; it does not independently verify the reference genome sequence, source genotype calls, multi-allelic context, or effect orientation. Reproduce with `python3 work/niagads_qtl_original_variant_audit.py`; the [receipt](stage6_niagads_qtl_original_variant_audit.json) records input hashes and counts.

The QTL record at GRCh38 `chr5:136850352:G:T` names rs10037001 while the older reference VCF names rs113343755. On 29 October 2026, Ensembl returned both identifiers at the same GRCh38 position and G/T contrast ([rs10037001](https://rest.ensembl.org/variation/human/rs10037001?content-type=application/json), [rs113343755](https://rest.ensembl.org/variation/human/rs113343755?content-type=application/json)); the identifier difference alone is therefore not evidence of a different base contrast. It also does not explain the observed frequency difference. The QTL record's reported minor-allele count and sample count are consistent with its listed frequency for this row, but that arithmetic does not verify the underlying calls. The four records remain excluded in the filtered analysis.
