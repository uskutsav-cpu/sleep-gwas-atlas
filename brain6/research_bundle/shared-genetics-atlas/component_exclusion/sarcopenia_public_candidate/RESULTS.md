# Public sarcopenia-factor GWAS: source audit and selected-marker screen

The public [Zenodo sarcopenia-factor release](https://doi.org/10.5281/zenodo.17069328) is a useful candidate for a multiscale physical-function comparison. The depositor describes a genetic factor based on muscle mass, strength and function, estimated with genomic SEM. **It is not a frailty GWAS with sleep and exhaustion removed**: the release does not identify its individual input GWAS or fully document their phenotypes, factor direction, model, overlap or the meaning of its constant N=522,632. Consequently, the file alone does not establish that every underlying input excludes sleep or fatigue, that N counts unique people, or that a score was measured in individuals.

We streamed the full 269,614,908-byte file twice without storing it, and both passes matched the publisher MD5; the first also recorded SHA-256. Its 7,319,145 rows included 1,135,028 rows with rsIDs found in the project's GRCh37 HapMap reference. **All those reference-matched rows** had matching chromosome, position and unordered allele pair. These are row counts, not necessarily independent variants; this restricted reference check does not validate the full genome or the factor model. A second scan found 7,396 rows with missing effect-allele frequency (`NA`); the other checked numeric fields were valid under the simple finite/range checks. In the first 10,000 valid numeric rows, reported p-values agreed with normal Wald p-values computed from beta/SE to numerical precision. This checks arithmetic on an ordered prefix, not study calibration.

## Selected-marker comparison

The six previously selected Brain6 insomnia–ADHD leads each had an exact GRCh37 coordinate and allele-pair match in the factor file and an existing MVP clinical-insomnia record. For two, both marginal p-values are small enough that their larger p-value multiplied by six is below 0.05:

| Lead | MVP insomnia p | Sarcopenia-factor p | max-p × 6 |
| --- | ---: | ---: | ---: |
| rs2582894, chr11 | 0.000933 | 0.0000390 | 0.00560 |
| rs1476535, chr7 | 1.20×10⁻⁷ | 0.00620 | 0.0372 |

These are **post-selected, source-dependent signals for follow-up**, not confirmed joint discoveries. The original candidate search was larger than six, and the later cross-source selection and other comparisons are not captured by a sixfold multiplier. The MVP insomnia and sarcopenia factor studies' participant or relative overlap is unverified. The [complete six-lead table](selected_comparison.tsv) includes the four non-crossing markers. For example, rs2431108 has factor p=8.92×10⁻⁷ but MVP insomnia p=0.0895; rs62264767 has factor p=0.367. The separately prioritized, strongly linked rs62264764 has factor p=0.369; it is not a seventh independent locus or part of this benchmark. A weak p-value does not prove absence of an effect.

The factor's sign and measurement are not sufficiently documented to assign a protective direction, a physical-frailty mechanism, a shared causal allele or a gene. A shared-marker p-value is not colocalization or two-trait replication. Further work needs the factor's underlying GWAS definitions, loading/sign convention, sample overlap and correct N semantics, followed by appropriately controlled genome-wide analysis and independent data. In parallel, the exact requested matched FI/FFS with sleep and exhaustion removed remains unavailable in the checked public sources.

[Pre-run protocol](PROTOCOL.md), [full-stream audit](audit.py), [first receipt](receipt.json), [numeric-failure audit](invalid_audit.py) and [second receipt](invalid_receipt.json), [raw selected rows](selected_source_rows.json), and [comparison code](compare.py) document the source, checks and results. The project GRCh37 reference file's SHA-256 was independently recorded in the first receipt after the stream. The entire source was not retained locally. The six-lead benchmark assumes valid marginal tests; the source's missing model metadata prevents accepting it as a calibrated discovery claim.

A subsequent [four-trait LDSC feasibility check](ldsc_feasibility/RESULTS.md) returned numerical results but cannot establish calibrated sleep–frailty genetic correlations because the factor inputs, direction and N semantics remain unverified.
