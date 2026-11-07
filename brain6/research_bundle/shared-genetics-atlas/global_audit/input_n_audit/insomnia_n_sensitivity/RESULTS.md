# Insomnia sample-size sensitivity for the Brain6 genetic correlation

The original Brain6 insomnia–ADHD genetic correlation is **numerically stable in this specific same-marker N sensitivity**. We matched all 1,133,335 valid archived munged insomnia rows uniquely to raw rows by rsID and nonpalindromic allele pair. Both copies have identical SNPs, alleles and Z values; only the insomnia N column differs. The raw per-variant N divided by the archived constant N ranges from 1.169 to 1.232, with median 1.228. Both LDSC runs ultimately use 1,054,978 allele-compatible SNPs.

| Insomnia N used | r_g | SE | p | Insomnia h² output |
|---|---:|---:|---:|---:|
| Archived constant 313,750.036 | 0.4033 | 0.0318 | 5.836×10⁻³⁷ | 0.0571 (0.0025) |
| Archived raw per-variant N | 0.4034 | 0.0318 | 5.9754×10⁻³⁷ | 0.0466 (0.0020) |

The matched constant-N run reproduces the full earlier reexecution's r_g, SE, p and final SNP count at displayed precision. The r_g changes by 0.0001 at displayed precision; the two estimates share data and this is not an independent replication or a statistical test of their difference. The h² values are LDSC's printed, N-dependent outputs, not validated observed- or liability-scale heritabilities. The corresponding cross-trait intercepts round to −0.0035 and −0.0034, each SE 0.0072; these do not establish participant nonoverlap.

The raw N may be total rather than effective sample size; the correct input for the original meta-analysis is not established. This analysis does **not** endorse either convention or establish source authenticity, variable-N treatment, ancestry fit, cohort independence, local colocalization, or failed LAVA quality gates. A constant scaling experiment on a different source cannot justify conclusions about this per-variant sensitivity. The original 396-test BH correction applies to the original analysis family, not automatically to this post-hoc alternate run.

The [protocol](PROTOCOL.md), [matching code](build.py), [matching receipt](build_receipt.json), [run script](run.sh), [constant-N log](rg_constantN.log), [raw-N log](rg_rawN.log) and [verification](verification.json) retain the comparison. Both derived inputs are stored alongside these files for reproduction; original source hashes matched the prior receipts. The original ADHD input, LD scores, implementation and LDSC flags are unchanged. This is a stability observation, not proof of a valid genetic correlation under all assumptions.
