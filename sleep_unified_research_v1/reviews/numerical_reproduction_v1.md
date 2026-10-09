# Independent numerical-reproduction review v1

Client audit date: 2026-10-08. Separate numerical audit of the completed insomnia–BMI native pilot; this report does not expand the earlier statistical-validity audit into a claim of complete native reproduction.

**Verdict: the processed-input native pilot is numerically consistent with stock LDSC, its stock jackknife arrays and the frozen printed result. The return-value capture wrapper has no detected estimator mutation. This is one natively estimated pair, not reproduction of the 396-pair source-QC chain or of the extension/validation families.**

## What was independently executed

`reviews/numerical_reproduction_v1.py` imports no LDSC or repository statistical implementation. Standard-library arithmetic reconstructs ratio estimates, pseudovalues, sample variances and two-sided normal tails from the existing pilot receipt and three stock delete arrays. It writes only `numerical_reproduction_v1.json` and `.tsv`. The JSON hashes all inspected files, including the capture wrapper, pilot files, frozen table/log, ten estimator source files in both checkouts, two pilot processed inputs and 44 reference members. No further native estimator run or dense-source reprocessing was launched.

Pilot: `native/core_pilot_v1/rg_insomnia__bmi.full_precision.json`. Frozen pair: `sources/recovered/results/tables/rg_matrix.tsv:14`. Original log: `/Volumes/Extreme SSD/Utsav-Research-Archive/Sleep-GWAS/FAILED-PARTIAL-2026-09-01-sleep_gwas_atlas/results/logs/rg_insomnia.log`, BMI section beginning at line 345 and summary row recorded in the review JSON.

## Capture-wrapper audit

At `scripts/native_ldsc_capture.py:33`, the original `estimate_rg` and `estimate_h2` functions are saved. The wrappers call each original once (`:36`, `:53`), then read return-object attributes, serialize a receipt and return the same object (`:50`, `:56`). There is no assignment to estimator attributes, input data or estimator arguments within these capture wrappers. The driver forwards remaining arguments and runs stock `ldsc.py` (`:68–69`).

Ten numerical source files were independently compared byte-for-byte with the archival `ldsc/` source and all match. The local estimator checkout identifies pinned commit `6c673952cee74bd5c57aef1555a03b1c015399a0`. This agrees with the parent comparison manifest; the reviewer did not rely solely on its `same:true` fields.

The additional `--print-delete-vals` operation is stock output instrumentation: `ldscore/sumstats.py:405–412` fits the estimate before writing delete arrays, and `:550–554` delegates to stock `np.savetxt`. The wrapper adds receipt writing after estimation. Static inspection supports numerical noninterference. An instrumented-versus-uninstrumented second native run was not performed, so observational equivalence is supported by source inspection and the independent numerical identities below rather than by a second complete fit.

`NATIVE_ESTIMATE_RETURNED` records that an object was returned; downstream scientific success must additionally check complete pair cardinality and finite, valid estimate/SE/P fields. The current one-pair receipt passes those checks. The wrapper's JSON does not itself enforce source-admission hashes; that responsibility belongs to its runner and input ledger.

## Independent ratio-jackknife verification

Each stock `.hsq1.delete`, `.hsq2.delete` and `.gencov.delete` file has 200 finite numbers; both h2 delete arrays are strictly positive. They are printed at sufficient decimal precision to preserve binary64 values when read back.

For block b, compute `r_delete[b]=cov_delete[b]/sqrt(h2_1_delete[b]*h2_2_delete[b])`, and pseudovalue `u[b]=200*r_full−199*r_delete[b]`. The jackknife estimate is the pseudovalue mean; SE is `sqrt(sample_variance(u)/200)`. This directly reproduces the pinned `RatioJackknife` calculation, including uncertainty in both h2 denominators.

| Quantity | Captured native value | Independently reconstructed value | Absolute difference |
|---|---:|---:|---:|
| rg ratio | 0.17229925855319342 | 0.17229925855319342 | 0 |
| Bias-corrected jackknife estimate | 0.1721569059188145 | 0.17215690591881455 | 5.55×10⁻¹⁷ |
| rg SE | 0.023041746382148513 | 0.023041746382148558 | 4.51×10⁻¹⁷ |
| Z | 7.477699636807108 | 7.477699636807093 | 1.51×10⁻¹⁴ |
| Two-sided P | 7.563481978016354×10⁻¹⁴ | 7.563481978017247×10⁻¹⁴ | 8.92×10⁻²⁷ |

P relative difference is 1.18×10⁻¹³. Captured P comes from SciPy chi-square survival probability at Z²; the independent calculation uses `math.erfc`. Both include the small reconstructed-SE arithmetic difference. The reviewer declared relative tolerance 10⁻¹², absolute estimate tolerance 10⁻¹⁵, and absolute P tolerance 10⁻³⁰⁰ before executing the arithmetic. All quantities pass. Independently reconstructed hsq1, hsq2 and gencov SEs also match, with maximum absolute difference 8.68×10⁻¹⁹.

Stock LDSC reports/tests `rg_ratio`, while `rg_jknife` is a separately recorded bias-corrected estimate. They differ here by approximately 0.000142353. Replacing the reported ratio with `rg_jknife` would change the estimand/output convention and is not warranted by this audit.

## Frozen-result comparison and precision limits

All 11 numeric summary/scalar-P fields in the new pilot log equal the archived BMI row. The original log scalar output is rg=0.1723, SE=0.023, Z=7.4777, P=7.5635×10⁻¹⁴. The frozen TSV further serializes the already printed values with `%.4g`, yielding rg=0.1723, SE=0.023, Z=7.478 and P=7.563×10⁻¹⁴. The independent checker reproduces this staged serialization exactly.

The native-versus-frozen differences are −7.41447×10⁻⁷ for rg, +4.17464×10⁻⁵ for SE, −0.000300363 for Z and +4.81978×10⁻¹⁸ for P. These differences are consistent with the known display/serialization precision; they are not evidence of a changed numerical estimate. The old logs and TSV do not retain the old unrounded result, so **identical at full precision cannot be established**. Use a comparison class such as `NATIVE_AGREES_WITH_HISTORICAL_PRINTED_PRECISION`, not `IDENTICAL_FULL_PRECISION`.

The new receipt retains the current binary64 output; “full precision” means the returned floating-point values from this run. It does not mean exact real arithmetic or recovered historical unrounded estimates. No frozen threshold or FDR family was changed.

## Input and scope limitations

The pilot uses recovered munged insomnia/BMI inputs and 22 chromosomes of European LD-score and M_5_50 reference members. The reviewer rehashed both processed files and all 44 used reference members: current processed hashes equal the ledger's current hashes; all reference hashes equal their expected ledger hashes. The processed files have **no historical expected SHA-256**, which prevents a byte-identical archival processed-input claim.

The parent raw-hash ledger records insomnia and BMI as exact expected-source SHA-256 matches. This reviewer checked ledger consistency, but did not repeat the large raw-source hashing or reconstruct raw-to-harmonized-to-munged transformations. The whole-core raw ledger has 42 matches and three missing dense raw files (LDL, HDL, triglycerides); inherited prefilter copies are a separate recovery stage and do not close the missing raw-QC boundary.

Native runtime receipt: Python 3.9.23, NumPy 1.21.5, pandas 1.3.3 and SciPy 1.7.3. The checkpoint's NumPy/pandas entries describe pipeline numerical/table tooling; they do not independently establish the historical dedicated LDSC runtime library versions. The native runtime is now recorded, but exact historical library-environment equivalence remains unverified.

The pilot and old log each retain 1,006,820 valid-allele SNPs. This supports the processed-input comparison, while not validating all prior harmonization decisions or source phenotype/effect coding. No extension or replication native estimate is certified by this report. The larger run's free-space/swap guard is not a scientific failure and should not be bypassed to manufacture a completion claim.

The stock arrays retain delete values but no SNP identifiers or genomic deletion coordinates. They suffice for this within-pair SE check. They are insufficient to establish cross-pair block alignment for discovery/validation covariance; capture or otherwise verify the actual SNP order and deletion boundaries before reusing arrays for such inference. The valid pilot cannot repair the archived multivariate fixed-scale rg-uncertainty issue identified separately in `statistical_validity_v1.md`.
