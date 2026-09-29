# HFRS check against the official FinnGen R12 manifest

**Checked:** 2026-09-24 (UTC).
**Purpose:** Determine whether the public standard R12 summary-statistics manifest identifies the custom HFRS or HFRS-without-dementia GWAS reported by Mak et al.

## Reproducible check

Streamed the official manifest URL into an in-memory TSV scan; no summary-statistics files were downloaded and no access form was submitted:

```sh
curl -fsSL --max-time 60 'https://storage.googleapis.com/finngen-public-data-r12/summary_stats/finngen_R12_manifest.tsv' \
  | python3 -c 'import sys,csv,re,hashlib; data=sys.stdin.buffer.read(); print("manifest_sha256="+hashlib.sha256(data).hexdigest()); rows=list(csv.DictReader(data.decode("utf-8").splitlines(),delimiter="\t")); print("manifest_rows="+str(len(rows))); hits=[r for r in rows if re.search(r"frail|HFRS|hospital frailty", " ".join(r.values()), re.I)]; print("frailty_matches="+str(len(hits))); [print("\t".join(r.get(k,"") for k in ("phenocode","phenotype","category","num_cases","num_controls","path_https"))) for r in hits]'
```

Observed output:

```text
manifest_sha256=730fc3f54f6a73d37dc6fc9c6a0fd13159f4134bbc274e93c4db980f93120772
manifest_rows=2469
frailty_matches=0
```

The scan searched all values in every manifest row for `frail`, `HFRS`, or `hospital frailty`, case-insensitively. It is a metadata-name check, not a search of every R12 result file.

## Interpretation

The primary article reports a custom continuous score formed from 109 weighted ICD-10 codes, analyzes an HFRS-without-dementia sensitivity phenotype, and reports FinnGen N=500,737. Its methods say the custom GWAS used SAIGE v0.35.8.8. FinnGen's public R12 methods documentation describes the standard core release as 2,502 endpoints (2,499 binary and three quantitative) in 500,348 participants, analyzed with Regenie. These are distinct analysis descriptions; a similarly named standard endpoint or an R13/R14 release result cannot be assumed equivalent to the custom HFRS.

The public standard R12 manifest contains no entry whose available metadata labels match HFRS/frailty. This narrows the standard-manifest route but does not prove that no separately delivered, opaque-named, or author-held custom result exists. The exact full HFRS and no-dementia summary-statistics files remain unverified. The archived publisher tables remain thresholded reported-hit tables only; they are not genome-wide inputs for LDSC or LAVA.

## Sources

- Mak et al., *Nature Aging* (published 2025-08-05), [primary article](https://doi.org/10.1038/s43587-025-00925-y), including Methods and Data Availability.
- FinnGen, [R12 data-download documentation](https://finngen.gitbook.io/documentation/data-download), which links the public manifest and describes the form/email route for summary-statistic downloads.
- FinnGen, [R12 GWAS documentation](https://finngen.gitbook.io/documentation/methods/phewas), describing the standard release analysis and endpoint/sample counts.
- Public manifest: <https://storage.googleapis.com/finngen-public-data-r12/summary_stats/finngen_R12_manifest.tsv>.

No web form, email, login, or protected route was used. Preserve the existing request draft as unsent; use only the official provider route if the exact custom statistics are still required.
