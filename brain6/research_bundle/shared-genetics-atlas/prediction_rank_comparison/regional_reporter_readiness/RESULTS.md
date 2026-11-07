# Expanded regional reporter-set feasibility

The public Siraj OL33 SK-N-SH file contains **47 biallelic SNV rows** in the two specified Brain6 GRCh37 windows: 12 on chromosome 3 and 35 on chromosome 5. Under the prespecified feasibility rules, **one** row is a candidate allele-skew positive and **three** have an unadjusted 95% skew interval within a selected ±0.25 log₂ margin. All four have at least one allele meeting the working positive-activity screen.

| GRCh37 assay variant | Feasibility category | Reported log₂ skew (SE) | Source skew FDR |
| --- | --- | ---: | ---: |
| chr5:103904399 G>A | Candidate positive | +0.717 (0.156) | 0.000098 |
| chr3:117780054 T>C | Near-zero comparison | −0.070 (0.089) | 0.730 |
| chr5:103902161 G>C | Near-zero comparison | −0.048 (0.081) | 0.817 |
| chr5:103981726 T>G (rs40465) | Near-zero comparison | −0.015 (0.066) | 0.944 |

These categories are **analyst-defined proxies**, not verified author emVar labels or biological ground truth. The near-zero category uses a chosen practical margin and unadjusted confidence intervals; it does not prove absence of an effect. The candidate positive is outside the original selected 425 variants, and its LD and trait associations have not been verified here. None of these rows is a direct functional test of rs62264764 or proof of a shared causal disease variant. The earlier [rs40465 study](../../variant_context/rs40465_mpra/RESULTS.md) describes other contexts and assay limitations.

A later [genetic context check](chr5_positive_followup/RESULTS.md) identifies the candidate as rs12659965, with moderate reference LD to the chromosome 5 leads and marginal associations in FinnGen insomnia and UKB FFS. The lookup remains post-selected and is not a corrected shared-locus discovery; the candidate was absent from the archived Brain6 ADHD regional input.

One candidate positive and three comparisons from two post-selected, correlated regions are far too few for a reliable AlphaGenome performance estimate or a generalizable train/test exercise. We did not fit a model, calculate AUROC or make an independent-validation claim. Adult-frontal DNase predictions and neural-derived short-construct expression are different readouts, source study training overlap is unknown, and the exact complete author QC and label definitions remain unverified. A useful benchmark needs a larger independent set, prespecified validated labels and negatives, allele/build harmonization, LD-block-aware evaluation and documented training-data overlap. The [AVI mirror](https://huggingface.co/datasets/just-dna-seq/alphagenome_avi) is an authorized research source under the user's accepted terms; it is not used here to train a variant-effect model.

The [pre-result protocol](PROTOCOL.md), [scan](scan.py), [47 source rows](regional_rows.tsv), [source receipt](receipt.json), [independent proxy check](verification.json), and [run log](run.log) preserve the analysis. The single archive member passed ZIP CRC and gzip completion; the complete archive was not downloaded or whole-file hashed. Coordinates follow the source's working GRCh37 mapping, which is not a full independent assembly audit.
