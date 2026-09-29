# Canonical LAVA binary-N crosswalk — read-only audit

The seven canonical trait inputs came from `scripts/01_harmonize.py`, whose binary-trait branches write conventional effective `N` for LDSC. The pinned LAVA 0.1.5 binary path instead consumes that per-SNP `N` with a trait-wide observed case fraction from `input_info.tsv`, setting reconstructed cases to `N * case.prop`. The first locus's exact files and all seven source-card hashes are bound in `provenance.json`; the locked dense-input audit certifies the full harmonized file hashes.

| Trait | NOT_RUN | source fields relevant to N | study cases / controls | first locus N | N / total | modeled cases |
|---|---:|---|---:|---:|---:|---:|
| adhd | 223 | Nca;Nco | 38691 / 186843 | 93,086.6 | 0.413 | 15,969.3 |
| bipolar | 168 | NCAS;NCON;NEFFDIV2 | 41917 / 371549 | 95,582.2 | 0.231 | 9,690.1 |
| insomnia | 671 | N | 109402 / 277131 | 313,750.0 | 0.812 | 88,801.9 |
| longsleep | 1291 | NONE | 34184 / 305742 | 122,985.4 | 0.362 | 12,367.8 |
| mdd | 589 | NONE | 170756 / 329443 | 449,855.9 | 0.899 | 153,570.1 |
| parkinson | 636 | N_cases;N_controls | 33674 / 449056 | 98,344.3 | 0.204 | 6,860.2 |
| scz | 142 | NCAS;NCON;NEFF | 53386 / 77258 | 114,216.8 | 0.874 | 46,673.2 |

Five source cards label the harmonized `N` as `total`; inspection of their inputs and the harmonizer shows effective `N` instead. The two cards labeled `effective` are internally consistent but still coupled to an observed aggregate case fraction in LAVA. Source archives expose variant-level case/control counts for ADHD, bipolar disorder, Parkinson disease, and schizophrenia, and variant-level total `N` for insomnia. The Dashti long-sleep and Howard MDD releases inspected here lack variant-level sample counts.

This is a **definite input-model semantic mismatch**, not a proven explanation for any particular low-local-h² `NOT_RUN` cell. Association Z, LD, and true local signal may dominate testability. No canonical output, input, or frozen threshold has been altered. A separate source-specific method decision and outcome-blinded pilot lock must precede any repaired run; do not substitute total `N` for a heterogeneous meta-analysis solely to increase testability.
