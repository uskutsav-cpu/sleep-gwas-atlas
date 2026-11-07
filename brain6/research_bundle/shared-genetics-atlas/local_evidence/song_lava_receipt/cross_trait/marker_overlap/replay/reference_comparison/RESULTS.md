# Independent LD reference comparison at locus 464

We selected 60 evenly spaced markers from the existing 2,741-marker FI–insomnia replay and compared the saved UK Biobank LAVA reference with the [1000 Genomes Phase 3 GRCh37 reference](https://hgdownload.soe.ucsc.edu/gbdb/hg19/1000Genomes/phase3/) in 503 European-labeled individuals. All 60 markers matched the selected rsID, position and allele pair, and all 1,770 pairwise dosage comparisons were available. The two sources differ in participant composition, sample size, LD computation and storage; no claim of independence from every GWAS participant is made.

Across these selected pairs, the median absolute difference in squared LD values was **0.0044**. For the 156 pairs with UKB squared value at least 0.2, the median difference was **0.021**, with three differences above 0.1. Most pairs across a 2.8 Mb region have low LD, so the all-pairs median can look small without demonstrating close agreement where fine-mapping is most sensitive. These paired observations are correlated and this post-selected sample of markers is not representative of the full genome or necessarily of all 2,741 locus markers.

## Reference-matrix quality finding

The LAVA reader returned **1.001693** for the off-diagonal entry between rs62262671 and rs55754265. An independent two-marker extraction reproduced the same value. Since a true correlation cannot exceed one, this is a numerical approximation or other reference/reader issue requiring investigation; we did not silently clip it. The 60×60 symmetric matrix has unit diagonal but a minimum eigenvalue **−0.001695**, so this sampled submatrix is not positive semidefinite. The 1000 Genomes dosage correlation for the pair was approximately 0.9958. The UKB reference documentation describes an eigen-compressed, truncated representation, but these checks do not establish exactly why the value exceeds one or whether it influenced the FI negative variance. Because this is a principal submatrix of the extracted locus LD, its negative eigenvalue also shows that an exact symmetric full matrix containing the same entries cannot be positive semidefinite.

**This is not a reference repair or a local genetic-correlation result.** The original local gates remain failed. The small independent panel cannot establish which LD source is correct or exclude sampling differences, and agreement for sampled pairs cannot validate the full LD decomposition or address sample overlap. An ancestry- and study-matched LD resource, appropriate handling of numerical approximations, and a frozen family-wide calibration are still needed.

A [one-marker-at-a-time replay](pair_removal/RESULTS.md) found that the FI variance remained negative after removing either member of the identified pair. This does not establish the cause of the failure or validate the remaining matrix.

A [whole-locus spectrum check](full_spectrum/RESULTS.md) finds 2,629 out-of-range entries and 484 eigenvalues below −10⁻⁸ in the 2,741-marker extracted reference matrix. The 216 components retained by the original 99% rule have positive eigenvalues; the cause and impact of the approximation remain unverified.

A [two-adjustment PSD sensitivity](psd_sensitivity/RESULTS.md) shows that eigenvalue clipping and a diagonal ridge give different FI variance and component counts. Neither is a validated correction or a rescue of the original family gates.

A [same-marker 1000 Genomes sensitivity](eur_lava_sensitivity/RESULTS.md) shows positive FI variance with 503 European-labeled reference participants while the matched-marker UKB reference remains negative. The references have different precision and composition, and neither result validates a local correlation or repairs the family gates.

## Audit trail

The [protocol](PROTOCOL.md) precedes the comparison. [UKB extraction](ukb.R), [selected markers](markers.tsv), [reference matrix](ukb_ld.tsv), [1000 Genomes comparison code](compare.py), [pair results](pairwise.tsv), [marker status](marker_status.tsv) and [receipt](receipt.json) document the work; no individual genotypes were saved. [Matrix QC](matrix_qc.py) records eigenvalues and an input hash in [its receipt](matrix_qc.json), while the separate [two-marker replay](pair_check.R) and [result](pair_check.json) reproduce the out-of-range entry. [Verification](verify.py) checks pair counts and summary arithmetic. The remote VCF and index were not fully downloaded or hashed; the full original LAVA reference provenance and source-GWAS sample matching remain unverified.
