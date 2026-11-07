# Neural-progenitor reporter candidates near rs62264764

Four previously selected chromosome 3 variants occur in the public neural-progenitor reporter matrices examined for Brain6. All four have **low pairwise LD with rs62264764** in the 503-person 1000 Genomes European reference: r² is approximately 0.009–0.010. Being in the same broad region is therefore not sufficient to treat their assay rows as direct tests of the prioritized marker.

| Matrix candidate | r² with rs62264764 | 95% reference bootstrap interval |
| --- | ---: | ---: |
| rs16827974 | 0.0105 | 0.0025–0.0237 |
| rs56376353 | 0.0088 | 0.0017–0.0212 |
| rs56385956 | 0.0099 | 0.0021–0.0228 |
| rs73187934 | 0.0096 | 0.0020–0.0224 |

The earlier [Lee assay analysis](../RESULTS.md) found rs16827974 and rs73187934 threshold crossings in at least one **provisional** reconstruction, but the authors' complete QC, activity criterion and allele-effect results were not reproduced. Those provisional calculations must not be presented as validated published effects. A regulatory result at any of these markers could still be biologically interesting, but the low pairwise LD does not support using it as direct functional validation of rs62264764. This does not exclude a different shared signal elsewhere in the region or establish which variant or gene affects sleep or ADHD.

The r² values describe phased haplotypes in one finite, European-labeled reference panel; the discovery cohorts' LD can differ. Bootstrap intervals resample reference individuals, not independent cohorts, and the query is limited to four post-selected variants. The original LAVA gates remain failed. Direct allele-specific testing in a relevant cellular context, full author-QC results and well-calibrated, study-matched fine-mapping are still needed.

The [pre-result plan](PROTOCOL.md), [analysis](analyze.py), [LD counts](ld.tsv), [receipt](receipt.json) and [independent count-based verification](verification.json) preserve the work. All five exact rsID/GRCh37 position/allele pairs matched the [UCSC-hosted 1000 Genomes reference](https://hgdownload.soe.ucsc.edu/gbdb/hg19/1000Genomes/phase3/). The indexed VCF was not downloaded or whole-file hashed. The first script attempt stopped on a coordinate-parsing error ([log](run.log)); the corrected [retry](retry.log) completed and saved no individual genotypes. Verification checks the aggregate haplotype arithmetic, not phasing, whole-file provenance or bootstrap regeneration.
