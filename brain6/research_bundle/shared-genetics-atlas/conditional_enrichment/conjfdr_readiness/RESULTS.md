# Conjunctional FDR: readiness and precision audit

**Status: no real-data conjunction-FDR discovery run.** Synthetic software and calibration stress tests were subsequently run; see the [R implementation assessment](r_implementation/RESULTS.md) and [mixed-signal simulation](r_implementation/synthetic_ld_overlap/mixed_signal/RESULTS.md). The existing max-p BH/BY analysis and the conditional-enrichment diagnostic are not the published pleioFDR conjunctional FDR method. This audit identifies concrete requirements and input limitations; it does not estimate conjunctional FDR or produce validated discoveries.

## Current resources and official requirements

The [official pleioFDR repository](https://github.com/precimed/pleiofdr) specifies MATLAB (tested versions ≥2015) and a workstation with at least **16 GB RAM**. It documents an experimental, unsupported Octave path. At this check, neither `matlab` nor `octave` resolved on this computer's command path; the computer reports **8 GiB** of physical memory. MATLAB might exist outside the command path; a license and executable have not been verified. The official genome-wide European LD reference URL reports **2,383,912,974 bytes** (about 2.22 GiB); the system volume had approximately 835 MiB available and the external volume approximately 9.2 GiB at the time of the check. These values can change. The reference alone exceeds the free system volume space and is not sufficient by itself to run the complete analysis. It was not downloaded or validated, and the user subsequently designated the external SSD for heavy work. A scan of the external volume did not locate the named official reference. No files were deleted.

The repository's `master` was observed at commit `0da963cac22fe8de9030166d7aea974bcbb0a367`; future runs should pin and archive their actual version. The README asks users to increase random pruning from its default 20 to 500 iterations for its example. The [default configuration](https://github.com/precimed/pleiofdr/blob/0da963cac22fe8de9030166d7aea974bcbb0a367/config_default.txt) contains `stattype=conjfdr`, actual `fdrthresh=0.01` and an explanatory recommendation of 0.05 for conjFDR; the analyst must freeze the intended threshold rather than assume those two statements are identical. It also specifies an hg19 MHC fitting exclusion and distinguishes exclusion from model fitting from exclusion from discovery. Random pruning, genome-build alignment, LD-reference fit and exclusion choices require explicit review and must not be changed after seeing results merely to obtain hits.

## Precision of available common-support files

A fresh scan of the already-used aligned HapMap3 inputs shows the following (each has 1,138,262 rows):

| Input | Distinct displayed absolute Z values | Rows tied with at least one other | Rows whose p=0.01 classification could change solely within a ±0.0005 Z rounding interval |
|---|---:|---:|---:|
| FinnGen insomnia | 4,483 | 1,137,872 | 60 |
| Atkins FI | 5,034 | 1,137,790 | 53 |
| Ye FFS | 856 | 1,138,161 | 1,387 |

The two-sided normal p-value calculations assume Z was rounded to three decimals. For the conventional p=5×10⁻⁸ cutoff, five FFS rows and no insomnia or FI rows in these files straddle the cutoff under that rounding interval. These counts are sensitivity checks on the **displayed Z values**, not proofs of errors in the original GWAS or estimates of FDR bias. Ties alone do not invalidate an empirical FDR procedure. The especially small number of FFS values reflects the available input representation; source-level p-value precision and the BOLT association convention require separate review. Rounded derived Z values cannot restore information absent from a source p-value. Obtain the most precise defensible source statistics and use the exact official converter with recorded alleles, build, sample sizes, missingness, and overlap checks.

The [script](precision.py), [table](precision_summary.tsv) and [receipt](receipt.json) include source hashes and precise definitions; the fixed ±0.0005 interval does not account for model differences, rounding before conversion, or all sources of numerical error.

## Conditions for a credible run

1. Use a suitable machine with verified MATLAB and sufficient memory and disk, obtain and verify the official reference and its build, and perform the repository's documented example as a reproducibility smoke test. Document the actual peak resources instead of assuming 16 GB guarantees success.
2. Reconstruct dense per-study GWAS inputs with source provenance, p-value precision, alleles and cohort-specific sample-size conventions; use the official reference converter and audit all joins, duplicate rsIDs, palindromic variants, absent variants, filtering and fit/discovery masks. The current HapMap3 intersection is not a substitute for proving adequate genome-wide coverage.
3. Document participant overlap and relatedness for each pair, ancestry compatibility and LD limitations. Sample overlap can bias conditional FDR procedures; [LeBlanc et al.](https://pmc.ncbi.nlm.nih.gov/articles/PMC6019513/) studied this issue and corrections. Review whether the chosen implementation and any correction actually address the estimated dependence; unverified overlap cannot be assumed absent.
4. Freeze pair family, directions, random-pruning seed and iterations, exclusion regions, genomic-control settings, candidate threshold and evaluation plan before inspecting outputs. Validate under realistic null simulations or negative controls that preserve LD, source variance and overlap; report all variants, LD-clumped loci, missingness, and sensitivity results. Component-excluded frailty analyses and independent cohorts remain separate requirements.

This plan is not a completed method validation or an instruction to run on the current constrained machine. The existing [conditional-enrichment assessment](../RESULTS.md) remains descriptive.

A later [R-implementation assessment](r_implementation/RESULTS.md) identifies a MATLAB-free candidate and checks the current SSD and reference requirements. Its documented sample-overlap correction is not implemented, so it does not remove the scientific blocker.
