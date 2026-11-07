# AVI cross-check across six regions

This check compares the highest mirrored AVI variant from each of the six candidate regions with the Broad Institute SpliceAI Lookup service. The selection was fixed from the mirror before querying, with ties resolved by the GRCh38 variant string. All six agree within the stated rounding tolerances. Both services use the same underlying AlphaGenome scores, so this checks data access and allele mapping—not biological replication.

| Variant | GRCh38 variant | Mirror PHRED interval | Broad PHRED | Result |
|---|---|---:|---:|---|
| [rs4923550](https://deepmind.google.com/science/alphagenome/atlas?q=chr11:28615816:T%3EC&m=variant) | chr11:28615816:T>C | 14.66099 | 14.66 | AGREES_TO_ROUNDING |
| [rs76247422](https://deepmind.google.com/science/alphagenome/atlas?q=chr14:98162916:C%3ET&m=variant) | chr14:98162916:C>T | 17.97764 | 17.98 | AGREES_TO_ROUNDING |
| [rs62264764](https://deepmind.google.com/science/alphagenome/atlas?q=chr3:117920728:G%3EA&m=variant) | chr3:117920728:G>A | 20.48566 | 20.49 | AGREES_TO_ROUNDING |
| [rs2431108](https://deepmind.google.com/science/alphagenome/atlas?q=chr5:104612267:T%3EC&m=variant) | chr5:104612267:T>C | 19.99941 | 20.0 | AGREES_TO_ROUNDING |
| [rs4916723](https://deepmind.google.com/science/alphagenome/atlas?q=chr5:88558577:A%3EC&m=variant) | chr5:88558577:A>C | 15.99192 | 15.99 | AGREES_TO_ROUNDING |
| [rs1456021](https://deepmind.google.com/science/alphagenome/atlas?q=chr7:114688105:T%3EG&m=variant) | chr7:114688105:T>G | 17.16267 | 17.16 | AGREES_TO_ROUNDING |

The comparison allows 0.000501 for raw score and 0.00501 around the mirror PHRED interval, reflecting the precision exposed by the Broad service. The original source genome-wide data were not downloaded or hashed locally. Six checks cannot verify all 425 variant rows.

The hash-bound cross-check table and lookup details are in `broad_crosscheck.tsv`; selection, hashes and timestamp are in `broad_crosscheck_provenance.json`. The [Broad lookup source](https://github.com/broadinstitute/SpliceAI-lookup/tree/main/google_cloud_run_services/score_lookup) describes its upstream tables.
