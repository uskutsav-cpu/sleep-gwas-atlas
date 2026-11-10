# Transcript annotations for five shared-risk candidates

Five of the 22 post-selected markers had positive estimates and nominal p < 0.05 in all eight inspected GWAS results. I queried the GRCh37 Ensembl Variant Effect Predictor (VEP) for those exact positions and the alternate allele of each score contrast. This is **annotation of variants**, not evidence that their nearest genes mediate insomnia or frailty.

| Variant | Shared-count allele / other allele | VEP alternate allele | Selected returned transcript annotations |
|---|---|---|---|
| rs12536335 | A / G | G | Intronic and noncoding transcript annotations in FOXP2; also an NMD-transcript annotation. |
| rs817509 | G / A | G | Intronic and noncoding transcript annotations in RP11-384F7.1. |
| rs2838788 | T / C | C | Intronic and noncoding transcript annotations in ADARB1; also an NMD-transcript annotation. |
| rs34290943 | T / C | C | Synonymous annotations in LAMB2 transcripts; noncoding-transcript and upstream/downstream annotations also occur, including upstream annotations for USP19. |
| rs12671330 | C / T | C | Downstream-gene annotation for FOXP2 transcripts. |

For rs34290943, VEP reports a synonymous consequence in two LAMB2 transcripts: this describes a coding substitution without an amino-acid change for those transcripts, not proof of altered function or proof that LAMB2 is the effector gene. The shared-count allele T is the reference, whereas the annotated alternate is C; the annotation must not be mistaken for the direction of the GWAS effect. NMD-transcript annotation identifies a transcript category, not evidence that the variant causes nonsense-mediated decay. Transcript annotations can differ by isoform, and a blank annotation would not establish lack of regulatory activity.

The five candidates were selected after examination of prior cross-trait and follow-up data. The eight association results include overlapping or potentially overlapping samples and unadjusted tests; they are not eight independent replications. Two markers are in the FOXP2 region and show nonzero reference-panel LD. The Song source genome build remains assumed, and multiallelic representations and source-specific effect-allele coding require verification. VEP position and allele matches confirm only the queried reference records, not causal sharing across traits. The FI includes sleep-related items; hospital frailty is a distinct diagnosis-based measure.

## What would strengthen the gene link

First, resolve variant identity and source alleles and use dense, ancestry-appropriate LD for conditional analysis and fine-mapping in each trait. Then assess colocalization with adequately powered molecular QTL data in relevant tissues, with sensitivity to multiple signals and overlapping samples. Only supported loci should proceed to cell-type-specific functional follow-up. An annotation alone cannot identify a causal variant, tissue, cell type, effect on gene expression, or usable personal risk score.

## Reproducibility

The [extracted table](shared5_allele_specific_vep.csv) includes positions, exact alleles, all returned transcript terms and gene symbols, URLs and response SHA-256 hashes. We required one matching coordinate record and a reference/alternate pair matching the score contrast. Results were retrieved from Ensembl's GRCh37 VEP ID endpoint for [rs12536335](https://grch37.rest.ensembl.org/vep/human/id/rs12536335?content-type=application/json), [rs817509](https://grch37.rest.ensembl.org/vep/human/id/rs817509?content-type=application/json), [rs2838788](https://grch37.rest.ensembl.org/vep/human/id/rs2838788?content-type=application/json), [rs34290943](https://grch37.rest.ensembl.org/vep/human/id/rs34290943?content-type=application/json) and [rs12671330](https://grch37.rest.ensembl.org/vep/human/id/rs12671330?content-type=application/json). Retrieval and extraction script: `work/vep_shared_five.py`; CSV SHA-256 `cc74b6e9eda9045699d5e345df559128b8c4f1887f6d359cac5318bff0e35a34`. See the [association review](shared22_eight_source_evidence_review.md) and [regional LD review](shared22_reference_ld_interpretation.md).
