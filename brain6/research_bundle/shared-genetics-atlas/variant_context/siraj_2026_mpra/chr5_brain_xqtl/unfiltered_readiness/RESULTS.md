# Unfiltered brain methylation–GWAS pilot at chromosome 5

The previously identified parietal-cortex methylation signal at CpG `cg24905254` was followed up using **unfiltered** Knight-ADRC SNP results from the open [NIAGADS NG00184.v1 atlas](https://dss.niagads.org/datasets/NG00184/). Indexed byte-range extraction returned 3,335 rows for this CpG in the selected GRCh38 1-Mb Brain6 window. Among the earlier significant-only portal rows, all **121 SNP rows** in the query interval matched the unfiltered extract's p values and betas; the other ten were indels, outside the SNP-only source file. This cross-check supports the lookup but does not authenticate every extracted record or prove complete regional coverage.

A post-selected, **single-signal approximate Bayes-factor** comparison used 2,360 allele-compatible SNPs for insomnia and 2,277 for ADHD. With per-SNP priors p1=p2=10⁻⁴, shared prior p12=10⁻⁵ or 10⁻⁶, methylation effect prior SD from 0.001 to 0.02, and case–control GWAS effect prior SD 0.2, the shared-variant posterior H4 ranges were:

| Pair | H4 across tested priors | Dominant model under these settings |
|---|---:|---|
| CpG methylation–insomnia | 0.000075–0.00160 | Both associated, different single variants (H3) |
| CpG methylation–ADHD | 0.0000020–0.000418 | Both associated, different single variants (H3) |

The models favor distinct dominant signals **conditional on their assumptions**, not a proof that secondary signals cannot be shared. The CpG and region were chosen after looking at earlier results. The 1-Mb window truncates the source's wider cis window; GWAS frequencies were unavailable for independent harmonization checks, and allele matching excludes palindromic or ambiguous SNPs. Possible multiple causal variants, incomplete coverage, source differences, unknown study overlap, modest QTL sample size and selected effect priors limit inference. The tested SD range is illustrative, not calibrated for the methylation phenotype. We did not establish valid LD–Z diagnostics for multi-signal fine-mapping. The analysis does not identify a causal variant or gene, demonstrate mediation or independent locus replication, correct for the full project search, or repair the failed LAVA gates.

## Source p-value diagnostic

The source Z agrees with beta/SE. Across 1,020 selected CpG rows with 2<|Z|<9, the reported p-values closely match a **two-sided t distribution with 370 degrees of freedom**: the largest absolute discrepancy in log p is about 2.2×10⁻⁶. For rs12658032, t(370) gives 1.0131×10⁻⁶, matching the source-reported 1.0131×10⁻⁶; a standard-normal tail gives 6.61×10⁻⁷. The degrees of freedom were inferred from these rows, not confirmed from the original model documentation. The coloc calculation treats beta/SE as approximately Gaussian, an additional limitation; the source p-value was not silently replaced. [Diagnostic code](diagnose.py) and [receipt](diagnosis.json) preserve the comparison.

## Provenance and verification

The 1.47-TB methylation archive and its 25.3-GB chromosome 5 member were not downloaded in full. The [publisher manifest](mqtl_manifest.tsv) matches its published MD5. Its 170,910-byte chromosome index matches publisher MD5 `57fe8e4233f4a4ecea73bb321a20bf2e`; indexed byte-range retrieval read about 30.4 MB of compressed blocks and retained only the selected CpG rows. The [locator](locate.py), [index receipt](locate_receipt.json), [extractor](extract.py), [range receipt](extract_receipt.json), [unfiltered selected rows](selected_cpg.tsv), and [matching receipt](matching_receipt.json) preserve positions and hashes. The full source member and archive were not hashed, and the range has no independent publisher checksum. Regional boundary extraction uses an index approach previously checked against a fully downloaded atlas file, not a complete checksum verification of this member.

[Posterior table](coloc_posteriors.tsv), [top conditional variants](conditional_top.tsv), [matching code](match.py), [model code](coloc_pilot.R), and [independent arithmetic check](verification.json) document the pilot. The 16 posterior rows were recomputed independently to within 4.5×10⁻¹⁶; this verifies arithmetic, not statistical calibration. Per-pair matched GWAS files are on the SSD-backed work volume, not bundled here.

Data were prepared, archived and distributed by the National Institute on Aging Alzheimer’s Disease Data Storage Site (NIAGADS) at the University of Pennsylvania (U24-AG041689); the atlas accession is NG00184. The xQTL atlas data depend on the participation of research volunteers and the work of the contributing investigators and staff.

A [tail-sensitivity diagnostic](t_tail_sensitivity/RESULTS.md) replaces the methylation statistic with an equal-p normal score while retaining the original SE and prior. H3 remains dominant across tested settings, but this heuristic is not a calibrated t-likelihood.
