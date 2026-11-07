# Original-versus-harmonized source cross-check

A partial byte-range query of the original GRCh37 file found exactly one row at each previously specified coordinate and allele pair: chr5:103947968 T/C and chr5:103964585 G/A. This original file's header contains no rsID, so the row-to-rsID association uses the independently recorded project coordinate and allele pairs, together with the harmonized query; it is not a direct rsID match in the original row.

For both selected variants, the harmonized GRCh38 release uses the opposite effect allele. The corresponding beta changes sign and the effect-allele frequency is complementary; SE, p-value and N are identical. Thus the negative harmonized betas in the [selected lookup](../RESULTS.md) are consistent with the positive original betas on their opposite alleles. The [recorded rows](receipt.json) and [arithmetic verification](verification.json) document this check.

Both queries are partial and do not independently verify whole-file content. A separate earlier full-stream read of the original file matched its publisher MD5, but that prior receipt does not verify the full harmonized file or every row at the time of this query. This is a source-orientation check, not independent replication, colocalization, or evidence of a causal variant.
