# Brain6 × AlphaGenome: six-region AVI pilot

**Status:** 425 of 425 prepared candidate–variant rows have an overall AlphaGenome Variant Impact (AVI) score from a pinned third-party re-encoding of Google DeepMind’s static SNV dataset. Subsequently, frontal-cortex reference and alternate predictions were obtained for one variant per region; see the separate [prediction report](../jbrowse_predictions/FRONTAL_CORTEX_REPORT.md).

AVI is a composite prioritization score that includes AlphaGenome regulatory predictions and other inputs. Its PHRED scale expresses a genome-wide rank, not a probability of causality or an estimate of disease risk. A higher score does not identify a tissue, gene, direction of effect or shared causal mechanism. Results concern the union of the previously reported single-signal ABF credible sets, not every regional GWAS variant.

## Chromosome 5 case study

The chr5:103,447,968–104,447,968 GRCh37 candidate contains 71 selected variants. The highest mirrored AVI score in this selected set belongs to rs2431108 (GRCh38 chr5:104612267:T>C). Its PHRED score is 19.99941; the prior Brain6 conditional ABF PIPs are 0.131 for insomnia and 0.244 for ADHD. Both PIPs depend on the single-causal-variant model; none of these values establishes that the variant is causal. The exploratory trait–trait PP.H4 was 0.991 at p12=10⁻⁵ and 0.917 at p12=10⁻⁶.

| Variant | GRCh38 variant | AVI PHRED interval | Insomnia ABF PIP | ADHD ABF PIP | Credible-set membership |
|---|---|---:|---:|---:|---|
| [rs2431108](https://deepmind.google.com/science/alphagenome/atlas?q=chr5%3A104612267%3AT%3EC&m=variant) | chr5:104612267:T>C | 19.999 | 0.131 | 0.244 | insomnia, adhd |
| [rs62362443](https://deepmind.google.com/science/alphagenome/atlas?q=chr5%3A104576413%3AA%3EC&m=variant) | chr5:104576413:A>C | 17.070 | 0.003 | 0.000 | insomnia |
| [rs40465](https://deepmind.google.com/science/alphagenome/atlas?q=chr5%3A104646025%3AT%3EG&m=variant) | chr5:104646025:T>G | 14.418 | 0.040 | 0.009 | insomnia, adhd |
| [rs185260](https://deepmind.google.com/science/alphagenome/atlas?q=chr5%3A104629522%3AA%3EC&m=variant) | chr5:104629522:A>C | 13.146 | 0.004 | 0.001 | insomnia |
| [rs254024](https://deepmind.google.com/science/alphagenome/atlas?q=chr5%3A104608319%3AG%3ET&m=variant) | chr5:104608319:G>T | 10.630 | 0.004 | 0.001 | insomnia |
| [rs10052804](https://deepmind.google.com/science/alphagenome/atlas?q=chr5%3A104584121%3AA%3EC&m=variant) | chr5:104584121:A>C | 9.760 | 0.002 | 0.000 | insomnia |
| [rs30266](https://deepmind.google.com/science/alphagenome/atlas?q=chr5%3A104636656%3AG%3EA&m=variant) | chr5:104636656:G>A | 9.585 | 0.117 | 0.261 | insomnia, adhd |
| [rs325485](https://deepmind.google.com/science/alphagenome/atlas?q=chr5%3A104659667%3AA%3EG&m=variant) | chr5:104659667:A>G | 9.487 | 0.002 | 0.026 | insomnia, adhd |
| [rs416223](https://deepmind.google.com/science/alphagenome/atlas?q=chr5%3A104655775%3AC%3EA&m=variant) | chr5:104655775:C>A | 9.083 | 0.022 | 0.067 | insomnia, adhd |
| [rs6869862](https://deepmind.google.com/science/alphagenome/atlas?q=chr5%3A104550512%3AT%3EG&m=variant) | chr5:104550512:T>G | 7.630 | 0.004 | 0.000 | insomnia |

Ordering is by the mirror’s exact integer raw AVI score, with no threshold. PIP and AVI are different types of evidence and are not combined into a causal probability. PHRED intervals reflect ties introduced by the published rounded raw score.

## All six regions

| Candidate region (GRCh37) | Selected variants | Highest AVI in selected set | AVI PHRED interval | PP.H4, p12=10⁻⁵ | PP.H4, p12=10⁻⁶ |
|---|---:|---|---:|---:|---:|
| chr5:103447968–104447968 | 71 | [rs2431108](https://deepmind.google.com/science/alphagenome/atlas?q=chr5%3A104612267%3AT%3EC&m=variant) | 19.999 | 0.991 | 0.917 |
| chr11:28076261–29076261 | 88 | [rs4923550](https://deepmind.google.com/science/alphagenome/atlas?q=chr11%3A28615816%3AT%3EC&m=variant) | 14.661 | 0.003 | 0.000 |
| chr14:98137479–99137479 | 69 | [rs76247422](https://deepmind.google.com/science/alphagenome/atlas?q=chr14%3A98162916%3AC%3ET&m=variant) | 17.978 | 0.545 | 0.107 |
| chr3:117142005–118142005 | 128 | [rs62264764](https://deepmind.google.com/science/alphagenome/atlas?q=chr3%3A117920728%3AG%3EA&m=variant) | 20.486 | 0.921 | 0.537 |
| chr5:87185500–88185500 | 3 | [rs4916723](https://deepmind.google.com/science/alphagenome/atlas?q=chr5%3A88558577%3AA%3EC&m=variant) | 15.992 | 0.138 | 0.016 |
| chr7:113571035–114571035 | 66 | [rs1456021](https://deepmind.google.com/science/alphagenome/atlas?q=chr7%3A114688105%3AT%3EG&m=variant) | 17.163 | 0.843 | 0.350 |

The chr11, chr14 and lower-coordinate chr5 insomnia units did not pass the original association-posterior gate for an insomnia credible set. Their selected variants come from the ADHD credible sets; reporting an insomnia PIP in the full table does not change that gate. The original LAVA family failed QC, no candidate met the project’s positive evidence-class gates, and none has compatible independent two-trait locus replication.

## Provenance and checks

The mirror is pinned at revision `04820dafc5af625d1d941e080b5ae7d2b17dae87`; its claimed original source checksum and input/output hashes are recorded in `avi_provenance.json`. Its README describes a long table, but the delivered Parquet uses three score columns per position. The assumed alphabetical alternative-allele order was spot-checked against published original-file examples for reference A and T and against Broad SpliceAI Lookup for all three alternatives at one C-reference and one G-reference site. An independent Broad lookup also returned PHRED 20.00 for rs2431108, consistent with the mirror’s 19.99941 after rounding. These spot checks do not prove every mirror record; the original 88.5 GB file was not independently hashed.

AlphaGenome and AVI are research tools, not clinical decision tools. The static AVI download and the API have different terms; see `mirror_LICENSE.txt` and the [AlphaGenome output terms](https://deepmind.google.com/science/alphagenome/output-terms). This pilot does not train a variant-effect model. Review the separate tissue-specific predictions against independent functional evidence; further Atlas track scoring can use an API key.

Sources: [AlphaGenome Atlas research report](https://www.medrxiv.org/content/10.64898/2026.09.16.26363192v1.full), [AlphaGenome model paper](https://www.nature.com/articles/s41586-025-10014-0), [pinned AVI mirror](https://huggingface.co/datasets/just-dna-seq/alphagenome_avi/tree/04820dafc5af625d1d941e080b5ae7d2b17dae87), [Broad score lookup](https://github.com/broadinstitute/SpliceAI-lookup/tree/main/google_cloud_run_services/score_lookup).
