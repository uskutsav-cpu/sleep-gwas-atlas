# Sample-size conventions in the main LDSC inputs

The four previously rerun input files matched their recorded SHA-256 hashes. A complete scan of their HapMap3-shaped rows found:

| Input | Nonmissing N rows | Distinct N values | N representation observed |
|---|---:|---:|---|
| Brain6 insomnia | 1,133,335 | 1 | Constant 313,750.036. |
| Brain6 ADHD | 1,087,785 | 1,241 | 87,373.146–128,213.795; most common 128,213.795. |
| Frailty index | 1,169,496 | 1 | Constant 175,226. |
| Fried frailty score | 1,165,370 | 1 | Constant 386,565. |

Each file has 1,217,311 total rows and N is missing on the remaining rows; no nonmissing N was nonfinite or nonpositive. The first two counts agree with the archived munging logs' nonmissing-beta counts; the fresh LDSC run itself uses fewer markers after cross-trait and LD matching. Missing N in these munged tables is not evidence that the original study had a missing sample size at those variants.

For five deterministic ADHD markers representing its five most frequent N values, each uniquely matched by rsID and allele pair in the archived raw file, munged N agrees within 0.0005 with **4×Ncase×Ncontrol/(Ncase+Ncontrol)** computed from that row's case and control counts. This strongly supports an effective-N convention for those rows, not a proof of all 1,087,785 values or authentication of the raw source. For the first valid insomnia marker, the munged N is 313,750.036 whereas its matching archived raw row reports N=385,506; the preparation recipe and whether a constant effective-N transformation is appropriate across the source study remain unverified. The FI and FFS constants numerically match the project’s reported source-study sample sizes, but this does not prove per-variant counts or the full source transformation.

Therefore the LDSC log's automatic **“observed scale”** label is insufficient to interpret the Brain6 insomnia and ADHD h² outputs as validated observed-scale heritability. They depend on the supplied N values and binary-phenotype conventions; no liability-scale conversion is established. The reproduced genetic correlation remains an LDSC estimate conditional on these inputs, ancestry, alleles and model; the separate constant-rescaling experiment on the archived FinnGen-derived file does not validate variable-N or cross-study handling. This audit does not change the 396-test correction, establish replication, or rescue LAVA QC.

The [protocol](PROTOCOL.md), [complete N scan](scan.py), [N receipt](receipt.json), [raw-row lookup](raw_sample.py), [selected comparisons](raw_sample.json), and [scan log](scan.log) preserve the assessment. The selected raw ADHD and insomnia files match previously recorded project hashes; the exact five-marker arithmetic is a targeted check, not a whole-file verification of every N transformation.

A [same-marker N sensitivity](insomnia_n_sensitivity/RESULTS.md) compares the constant insomnia N with its raw per-variant values; r_g changes from 0.4033 to 0.4034 at displayed precision, without establishing which N convention is correct.
