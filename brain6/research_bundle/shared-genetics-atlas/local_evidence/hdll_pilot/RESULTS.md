# HDL-L feasibility pilot: insomnia and frailty

This post-selected pilot ran the authors' [HDL-L method](https://www.nature.com/articles/s41588-025-02123-3) on one chromosome 3 region. It returned positive local correlation estimates under three hypothetical sample-overlap settings. **These are exploratory, uncorrected estimates from the same data, not independent replication or a rescue of the failed LAVA quality gates.** The frailty phenotype includes a sleep-related component, which complicates interpretation.

| Assumed shared participants | Local genetic correlation | Method-reported interval | Uncorrected p-value |
| --- | ---: | ---: | ---: |
| 0 | 0.8710 | 0.5042–0.9756 | 0.000647 |
| 50,000 | 0.9002 | 0.6026–0.9864 | 0.000694 |
| 100,000 | 0.9002 | 0.6026–0.9864 | 0.000694 |

Actual participant overlap is unknown; these values are assumptions, not measurements. The agreement between the two positive-overlap runs does not establish robustness to overlap or residual confounding. Reversing the order of the traits at 50,000 produced the same correlation, interval and p-value to displayed precision, and swapped the univariate results as expected; this software check does not validate the model. All runs returned without recorded R warnings. The FI univariate p-value was printed as zero; its numerical cause is unverified and it should not be interpreted as literally zero.

The region was chosen after earlier analyses, and neither the p-values nor intervals account for selection or multiple testing. As context only, 0.05/2,468 is approximately 0.0000203, below all three reported p-values. This is not an actual family-wide result: only one HDL-L block was tested, and an appropriate testing family was not established for this pilot. A [subsequent preprint](https://www.preprints.org/manuscript/202507.0966) disputes aspects of the HDL-L paper's method comparison and reports reference-panel sensitivity; it has not been presented here as a settled resolution.

## Inputs and reproducibility

The selected UK Biobank reference block is chr3.44, GRCh37 positions 47,598,663–50,380,661, which overlaps but is not identical to the previously examined LAVA locus 464 (47,588,462–50,387,742). HDL-L and LAVA have different local random- and fixed-effect targets, respectively, and different LD partitions; the results cannot be directly equated. The reference contains 681 markers, of which 678 matched FI and 677 matched insomnia by identifier and allele pair. The median input N is 175,226 for FI and 313,750 for insomnia; the latter is an effective binary-trait sample size, not a count of distinct participants. The analysis uses the authors' default reference sample size 335,272 and 99% eigenvalue setting. Reference and trait cohorts may overlap; this was not established.

The full small BIM archive was verified against its publisher MD5. The selected LD member was fetched from the [published reference archive](https://zenodo.org/records/14825987) by byte range and checked using its ZIP CRC and SHA-256; the 1.83 GB archive was not downloaded or checked as a whole. The authors' HDL-L source was pinned to commit `cadb27c9aaa0e64a06a480cf2b0390a38f4f1d85`. The reference and source hashes were rechecked locally. See the [protocol](PROTOCOL.md), [download receipt](fetch_receipt.json), [analysis script](run.R), [raw results](run_results.json), [trait-order script](swap.R), and [trait-order results](swap_results.json). The small reference members and source are on the SSD in `/Volumes/Brain6Work/hdll-chr3-44-20261028`.

A defensible follow-up needs a frozen testing family and calibration, appropriate LD and participant-overlap information, a phenotype without the sleep component where available, and independent data. A missing or failed local test should not be treated as evidence of no overlap.
