# GENOA public methylation-QTL follow-up

The [GENOA meQTL study](https://doi.org/10.1038/s41467-023-37961-4) provides public cis-QTL summary statistics from whole blood in 961 African American participants; the study analysed variants within 50 kb of each CpG. Its [Zenodo archive](https://doi.org/10.5281/zenodo.7697509) supplies chromosome-level files on GRCh37. We downloaded the chromosome 15 file through HTTPS, verified its MD5 against the repository (`679b3497ae07e5694c41c9e6aaadaeb8`), and scanned all 7,750,441 rows. There were 369 rows for cg12594244.

The project candidate rs11636381 is 83,871 bp from the site, outside GENOA's ±50-kb cis window. The full chromosome file has no row for that exact candidate–CpG pair, as expected from its design. This is **not a null association result**. The study is not an individual-level score validation and does not measure sleep or frailty in this analysis.

Two previously identified GoDMC site-associated variants were present in the GENOA rows:

| Variant | GRCh37 position | Reported alleles (minor/major) | GENOA beta (SE) | GENOA p |
|---|---:|---|---:|---:|
| rs4077284 | 74228391 | A/G | −0.0574 (0.0602) | 0.340 |
| rs2165241 | 74222202 | T/C | 0.0140 (0.0568) | 0.806 |

The source README labels allele1 as minor and allele0 as major, but does not explicitly say which allele its beta counts; precise direction harmonization remains unverified. These are marginal estimates from an additional population, not a formal cross-population replication test. Lack of a significant association does not prove no effect, especially with the differing sample sizes, genetic backgrounds and LD. The analysis does not independently authenticate the previous GoDMC response, which came over HTTP. The smallest p among all 369 rows is post hoc and not a prespecified test; it is not evidence that our candidate variant shares a causal signal.

Project implication: the nearby methylation result remains a hypothesis, and current data do not link the candidate to this CpG. A full mQTL resource covering the candidate–site distance, a verified effect allele and suitable ancestry-matched regional data would be needed to test a link and then consider colocalization.

Artifacts: [site-level extracted rows](stage6_genoa_cg12594244.csv), [receipt](stage6_genoa_cg12594244_receipt.json), `work/stage6_genoa_cpg_extract.py`, and the checksum-verified source file `work/genoa_mqtl/meQTL_summarystat_chr15.txt.gz`.

A separate [variant-centered scan](stage6_genoa_rs11636381_mqtl_review.md) identified 38 tested CpGs within the GENOA cis window, including a post hoc small-p association at cg03789871. That site differs from cg12594244 and has not been established as frailty-associated.
