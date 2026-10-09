# Independent numerical reproduction review, v1

Verdict: **one processed-input insomnia–BMI pilot passes independent arithmetic reconstruction and current byte-provenance checks. Full source-chain and 396-/1,200-family native reproduction remain unverified.**

This is a new reviewer implementation, distinct from the existing statistical/numerical reviewer. It does not read or import that reviewer's code, reports, or computed results. Formulas were checked against the actual pinned LDSC source. Computation uses Python standard-library `math.fsum`, direct delete-vector variance, and `math.erfc`; neither LDSC nor the capture wrapper is imported or executed. No native GWAS analysis was rerun and original SSD inputs/frozen outputs were not modified.

The capture contains one estimate, and each of the three stock delete files contains 200 finite scalar values. All delete h2 products are positive. Independently recomputed values are:

| Statistic | Independently calculated | Captured reference | Absolute difference |
|---|---:|---:|---:|
| rg ratio | 0.17229925855319342 | 0.17229925855319342 | 0 |
| Bias-corrected jackknife rg | 0.17215690591881838 | 0.1721569059188145 | 3.886e-15 |
| rg SE | 0.023041746382148544 | 0.023041746382148513 | 3.123e-17 |
| z | 7.4776996368070971 | 7.4776996368071078 | 1.066e-14 |
| two-sided P | 7.5634819780170295e-14 | 7.5634819780163542e-14 | 6.753e-27 |

The displayed LDSC point estimate is `rg_ratio=cov/sqrt(h2_1*h2_2)`. Its SE comes from the stock ratio jackknife. The bias-corrected jackknife mean differs by 0.00014235263437503698; it is not the displayed point estimate. With B=200 and delete ratio d_i, the independent SE calculation is `sqrt((B-1)/B * sum((d_i-mean(d))^2))`. An independently evaluated pseudovalue expression agrees within 1.388e-17. Heritability and genetic-covariance SEs reconstructed separately from their stock deletes also pass.

Capture comparisons use relative tolerance 1e-12 and absolute tolerance 1e-14; P uses relative tolerance 1e-12 and zero absolute tolerance. These limits and the historical tolerances were declared in the review code before reading pilot estimates. There are 26 passing comparisons in the detailed TSV.

The immutable 396-row `rg_matrix.tsv` hashes to `161756ac61775ad3393572cf5fbdf51690eb47ff010b8a50dbcc5fabf79bc077`, matching the frozen checkpoint. Its insomnia–BMI row records rg=0.1723, SE=0.023, z=7.478, P=7.563e-14; all tested fields agree within declared print-precision intervals. For rg/SE the interval is ±5e-5; z and h2 intercept use ±5e-4; P uses ±5e-18. The historical and new native scalar logs both retain z=7.4777 and P=7.5635e-14, with all scalar comparisons passing half-unit intervals at the final printed place. Full-precision historical identity cannot be tested from rounded archived output. The fixed-width summary P=0.0000 is a display truncation and must not replace its scalar scientific-notation P. FDR was not recomputed by this review.

Actual inputs were independently SHA-256 hashed again rather than accepting ledger assertions:

| Input stage | File | Actual SHA-256 | Evidence |
|---|---|---|---|
| core_raw | insomnia.txt.gz | `32848cd92a6324c9048cad4a804cf1546299af1d0281981c6ca43c2046a9065b` | historical expected hash matched |
| core_munged | insomnia.sumstats.gz | `b4f66dde69d0c793f77b674dca45526659e69e6fbb31a5a4f3e57870077e6fc4` | current ledger hash matched; historical hash unavailable |
| core_harmonized | insomnia.harmonized.tsv.gz | `23e06f4a682477064e73551eeab78c292af4d2058fe636ddf328c4bccf560ea4` | current ledger hash matched; historical hash unavailable |
| core_raw | bmi.txt.gz | `0d6ed0ea97870916b830ccae349df94ca6f3cc68c025c79e784729af7f7136a4` | historical expected hash matched |
| core_munged | bmi.sumstats.gz | `481bcfaf59fc6459853c6319414507d67cbf102dc77d93e527705ed23429b078` | current ledger hash matched; historical hash unavailable |
| core_harmonized | bmi.harmonized.tsv.gz | `c9b0828c40deb29b5b34f4a1f99fedb510ffc17552ec0086698b7832d2e09965` | current ledger hash matched; historical hash unavailable |
| core_source_archive | Meta-analysis_Locke_et_al+UKBiobank_2018_UPDATED.txt.gz | `0d6ed0ea97870916b830ccae349df94ca6f3cc68c025c79e784729af7f7136a4` | historical expected hash matched |
| core_source_archive | Insomnia_sumstats_Jansenetal.txt.gz | `32848cd92a6324c9048cad4a804cf1546299af1d0281981c6ca43c2046a9065b` | historical expected hash matched |

The official 33,357,890-byte reference archive hashes to `9537f00eb0d163a935aaa2cf04b358b7cf21852279b9c7925802526f6060b069`. All 44 consumed chromosome LD-score/default M_5_50 members match their archive-member expected hashes and the input ledger. The reference and weight prefixes point to the same checked directory. Active LDSC HEAD is `6c673952cee74bd5c57aef1555a03b1c015399a0`, its worktree is clean, and ten actual code files each match both the pinned Git payload and archived SSD code bytes. The JSON receipt records every hash. Byte matching is distinct from validating reference ancestry/design assumptions.

Static wrapper review found no estimator modification: the wrapper invokes the saved original estimator, copies return attributes, writes JSON, and returns the original estimator object. The pilot lacks a single atomic receipt that binds at-execution source/code hashes to output hashes, so this review verifies current byte provenance rather than asserting that post-hoc hashes establish every runtime property.

The pilot log reports 1,133,335 insomnia SNPs, 1,217,311 BMI SNPs, 1,130,007 after reference/SNP merging, and 1,006,820 with valid alleles. Those counts were observed in the log; filtering, allele orientation, genome mapping, sample definition, effect coding and overlap assumptions were not independently rerun. Raw-source/archive hashes match historical expectations; harmonized/munged hashes match current receipts but have no historical expected hashes. This supports the pilot's processed-input numerical reproduction and does not complete the original QC chain.

The 200 scalar output rows have no SNP identities, genomic coordinates or separator map. Their numerical agreement cannot establish genomic alignment between discovery/validation runs, justify cross-run paired jackknife covariance, or establish an independent biological replication. Such claims require matched SNP/block evidence and source/cohort review.

The full-family plan's latest inspected resource receipt records 258,363,392 internal bytes free versus a 3,221,225,472-byte (3 GiB) floor, so its resource gate fails. No full family run is admitted by that receipt. This review performed one numerical pilot reconstruction; a 396-row frozen table and a job manifest are not completion of 396 native estimates or 1,200 extension estimates. Raw-chain rerun, all-family execution, genomic alignment/covariance calibration and biological interpretation remain unresolved human/scientific review requirements.

Reproduce this audit with `python3 sleep_unified_research_v1/reviews/independent_numerical_v1.py`. The script, report, JSON receipt and TSV are sealed by `independent_numerical_v1.sha256`. Scope is research evidence only; no manuscript text is generated.
