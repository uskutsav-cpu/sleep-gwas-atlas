# Same-marker FI sensitivity with an independent LD reference

We repeated the selected locus 464 **FI univariate calculation** with a 1000 Genomes Phase 3 European-labeled reference and compared it with the original UKB LAVA reference on the same matched markers. Of the original 2,741 markers, **2,738** matched by GRCh37 position, rsID and allele pair in the public reference; three were not matched. Both calculations use the same FI statistics, marker order, sample size and 99% LD-variance retention rule.

| LD reference | Markers | Retained components | FI genetic-variance diagonal |
|---|---:|---:|---:|
| Original UKB LAVA reference, matched support | 2,738 | 216 | −1.822×10⁻⁹ |
| 1000 Genomes European panel | 2,738 | 183 | +1.881×10⁻⁶ |

With UKB LD, FI remained below the original failure boundary. The 1000 Genomes calculation produced a positive variance on the matched support. Thus a change of reference representation is sufficient to change this **selected univariate** result, even after matching marker support. It does not establish which estimate is correct, why they differ, or whether the reference approximation caused the original failure. The two panels differ in size, ancestry composition, genotypes and estimation; the public panel contains only 503 people and has rank at most 502. The 99% thresholds therefore retain different numbers of components. No bivariate local correlation was computed, and the frozen LAVA family gates remain failed.

A [30-resample diagnostic](bootstrap/RESULTS.md) shows sizeable changes in the estimate and component count within this small reference. Bootstrap draws repeat people and are not independent replication or calibrated intervals.

A [five-population holdout](population_holdout/RESULTS.md) leaves a positive FI variance in each overlapping reference subset, with smaller sample sizes and fewer retained components. It does not identify an ancestry effect or validate the local result.

A [direct LAVA runtime check](lava_runtime/RESULTS.md) reproduces the FI-only SVD result and reports uncorrected univariate values for both traits, using the public dosage LD but the original UKB allele frequencies. It does not establish a bivariate finding or rescue a QC family.

## Methods and limitations

The [export](export.R) uses LAVA's processed FI Z statistics and independently checks them against raw statistics multiplied by the package's allele-alignment factor. LAVA retains original allele labels after flipping the statistics, so this check is necessary before using the UKB reference A1 orientation. We streamed the indexed [UCSC-hosted VCF](https://hgdownload.soe.ucsc.edu/gbdb/hg19/1000Genomes/phase3/) and calculated A1 diploid dosages for 503 European-labeled panel participants. Missing calls, if any, were mean-imputed per marker; nonmatching or ambiguous records were excluded. The [Python calculation](eur.py) uses singular-value decomposition of standardized dosages and the saved LAVA univariate formulas; it is a numerical reproduction of the calculation, not a run of LAVA with a fully validated alternative reference. The UKB comparison was performed with the pinned LAVA runtime and exactly the same matched SNP list.

The [protocol](PROTOCOL.md), [aligned FI data](aligned_fi.tsv), [matched support](support.tsv), [excluded markers](excluded.tsv), [public-panel result](eur_result.json), [UKB replay code](matched_ukb.R), [UKB result](matched_ukb.json), and [verification](verify.py) preserve the analysis. The VCF and index were not downloaded or hashed in full; no individual genotypes were saved. Verification checks marker sets, input order and saved arithmetic, not independent genotype retrieval, study overlap, source provenance or false-positive calibration. Both inputs are post-selected and neither outcome should be chosen based on its sign.
