# Exhaustion GWAS: repeated-ID and allele-mismatch audit

The full exhaustion archive was streamed a second time and matched both its publisher MD5 and the previous run's SHA-256. The independently recomputed chromosome counts exactly match the [first whole-file audit](../zenodo_full_exhaustion_audit/RESULTS.md).

Of 873 source rows whose rsID and GRCh37 position matched the HapMap reference but whose unordered allele pair differed:

- **862 rows**, representing **861 rsIDs**, have another source row with the same rsID and the reference allele pair.
- **11 rows**, representing **11 rsIDs**, have no source row with that same rsID and reference allele pair.
- Every mismatched row shares exactly one allele with the HapMap pair. This pattern is consistent with alternate alleles at multi-allelic sites, but a reference panel and full variant normalization are needed to characterize each site. It does not establish effect direction or decide which allele is biologically relevant.

All 861 repeated reference-matched IDs had at least one mismatched row and a corresponding exact-pair row. A join by rsID alone can therefore duplicate or misassign associations. A defensible future pipeline must match build, position, and allele representation, retain per-allele identity, explicitly handle multiple records and palindromic variants, and report exclusions. The 11 IDs lacking a matching source allele pair should be held for reference/normalization review rather than guessed or strand-flipped automatically. An exact pair is an eligible *candidate for further QC*, not evidence that the study's effect convention or phenotype has been validated.

The [mismatch table](allele_mismatches.tsv) preserves all 873 rows and flags reference-pair companions. [Code](audit.py), [protocol](PROTOCOL.md), [chromosome counts](by_chromosome.tsv) and [receipt](receipt.json) record the source, build, reference and hashes. The 861 and 11 figures count rsIDs; 862 and 11 count source rows. The underlying reference is a filtered HapMap map, not an exhaustive representation of all alleles. The full archive was not saved.

This is a source-data harmonization diagnostic. Ancestry, genotype model, phenotype coding, sample overlap and per-variant N remain unresolved, and no sleep–exhaustion association or modified frailty GWAS has been estimated. The five component GWAS files do not themselves provide FFS with exhaustion removed.

A separate [prefix-only model and frequency diagnostic](../zenodo_frequency_model_audit/RESULTS.md) records compatible simple case–control standard-error scale and selected European-reference frequencies; neither verifies the source metadata.
