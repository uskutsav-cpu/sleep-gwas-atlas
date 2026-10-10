# Result-independent FinnGen preprocessing detail

This supplements the existing R13 insomnia feasibility protocol before full-body decompression or any new fit. It changes no historical analysis, biological testing family or QC threshold. The earlier bounded header inspection and its sequencing deviation remain in their original receipt.

## Literal input and build checks

Require the exact acquired R13 F5_INSOMNIA body, original size/MD5 and sealed SHA-256, and the literal 13-column header recorded in the bounded inspection. Require the provenance-bound historical nonambiguous HM3 allele list and GRCh37 chromosome/position reference, together with the original UCSC hg38ToHg19 chain and unchanged fail-closed point-liftover implementation. Each candidate variant must have exactly one HM3 rsID, a unique autosomal coordinate mapping, and the identical GRCh37 chromosome/position. Complement REF/ALT on a reverse-strand mapping before the original nonambiguous allele-orientation check. Reject candidate rows that do not meet these requirements and retain rejection counts; no source-wide coordinate concordance percentage is invented.

Source order is preserved. Every already retained rsID is excluded as a duplicate. Do not invent rsIDs from position or use a newer SNP map. Retain only finite supplied beta/sebeta, positive sebeta, finite ALT frequency with MAF strictly above 0.01, nonambiguous compatible alleles and finite Z=beta/sebeta. Exclude chromosome 6 positions 25,000,000 through 34,000,000 inclusive in GRCh37. Every source row receives exactly one terminal disposition; diagnostic counts are additional and cannot be summed as rejection dispositions.

Use the historical FinnGen direct-munged schema SNP/A1/A2/Z/N and its 12-significant-digit Z/N formatting. N is the explicitly assumed constant effective size 4/(1/51643+1/446273). No INFO field or per-variant sample-size verification is manufactured. Write a deterministic gzip derivative with mtime zero and empty stored filename. Require the historical minimum 700,000 retained variants as a necessary diagnostic condition, without treating that condition as scientific source admission.

## Diagnostics and preservation

For retained rows, record supplied-P versus two-sided normal beta/sebeta-P diagnostic bins: missing/invalid P, literal P=0, absolute relative difference above 0.1, and maximum finite relative difference using a denominator floor of 1e-300. These bins never filter or select rows. A zero supplied P is kept as a distinct diagnostic, not proof of consistency. Record oriented Z and finite N streams, ordered SNP/allele identity, full decompressed stream SHA, gzip CRC/EOF, full source hashes before/after and every filter count. Missing annotations remain missing.

One supervised worker, BLAS threads one, the shared heavy-worker mutex, 2 GiB aggregate owned RSS, 3 GiB internal free, 5 GiB SSD free and the previously reserved 2 GiB feasibility namespace remain mandatory. Decompression and the single observed-scale h2 fit follow completion and numerical review of the original 190 commands. Freeze code, assets, exact commands, output namespace and independent admission before execution. Preserve partial files and immutable failure receipts; no automatic retry, overwrite, alternate release or new rg is authorized.

The output is a qualified source-feasibility diagnostic. INFO/N omissions, clinical construct differences, overlap, observation horizon and human clinical interpretation remain unresolved even if all computational checks pass. No population prevalence, liability conversion, independent replication or biological-null inference is admitted.
