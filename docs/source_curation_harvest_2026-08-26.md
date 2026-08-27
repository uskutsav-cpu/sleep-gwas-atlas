# Selective source-provenance harvest — 2026-08-26

The remote `claude/neuro-immune-cardio` branch expands the project to 149
traits, so it cannot be merged into the locked atlas-v1.0 scope. Commit
`8023d143f2bde431ce0f54fa827302cd9d37adcd` nevertheless contains useful
file-level provenance from complete downloads that were later evicted to save
disk. This review considered only the 28 source-pending traits in the locked
45-trait manifest.

## Promoted exact matches

Ten selected releases matched the locked publication, phenotype stratum,
sample-size intent, European ancestry, and GRCh37 requirement. Their official
GWAS Catalog landing/download URLs, exact byte counts, SHA-256 checksums, and
branch-recorded build-anchor validation are now in
`config/public_gwas_sources.tsv`:

- de Lange 2017 combined IBD, Crohn's disease, and ulcerative colitis;
- Ishigaki 2022 European-only rheumatoid arthritis;
- Yengo 2018 GIANT+UKB BMI;
- Graham 2021 European HDL, LDL, and log-triglycerides strata;
- Mishra 2022 GIGASTROKE European any stroke; and
- Nielsen 2018 atrial fibrillation.

This branch harvest advances source verification from 17/45 to 27/45. A
separate direct verification of Campos 2020 snoring brings the current total to
28/45. Subsequent direct verification of the selected colorectal- and
lung-cancer releases brings the current total to 30/45. These counts do not
claim that the archives are currently materialized. Direct verification of the
FinnGen R9 sleep-apnoea endpoint subsequently brings the total to 31/45; its
GRCh38 build remains a separate harmonization blocker. Official file-level
verification of the Jansen 2019 UK Biobank-only insomnia release and Jones 2019
UK Biobank morning-person release brings the current total to 33/45. Complete
verification of the Bellenguez 2022 Alzheimer disease Stage I file then brings
the current total to 34/45. The public genome-wide file has 85,934 clinically
diagnosed or proxy cases and 401,577 controls; the paper's larger headline
total includes Stage II follow-up and must not be attached to every genome-wide
row. Direct verification of the Nalls 2019 public European Parkinson disease
and UKB proxy-case release brings the total to 35/45. Its exact N is 482,730,
not the inconsistent 1,474,097 previously carried in the manifest, and its
literal raw schema requires a GRCh37 coordinate-and-allele-to-rsID mapping.
Direct verification of OpenGWAS `ieu-b-38`, the exact Evangelou 2018 European
UKB+ICBP systolic-blood-pressure meta-analysis, then brings the total to 36/45.
Its 216,210,855-byte GWAS-VCF passed complete SHA-256 and multi-member gzip
checks and declares 7,088,083 harmonized GRCh37 variants. A source-specific
downloader resolves the public endpoint's short-lived URL from the stable
dataset ID, and a strict materializer maps ALT-relative ES/SE/LP/AF/SS fields
to the atlas TSV contract. The primary phenotype is untransformed SBP in mmHg;
UK Biobank averaged available readings and added 15 mmHg for participants
taking antihypertensive medication.
Verification of Aragam 2022 `GCST90132314`, the exact primary CAD discovery
meta-analysis rather than its Biobank-Japan extension, then brings the total to
37/45. The official README and ranged file inspection establish the GRCh37
schema, while an earlier complete branch download supplies the independently
recorded 3,250,443,774-byte SHA-256 and nine-anchor build check. This promotion
does not erase the ancestry limitation: the paper reports more than 95%
European participants, but the Catalog explicitly labels the source European
plus ancestry not reported. The manifest preserves that label, so CAD remains
blocked from the strictly EUR harmonization path. Its coordinate-only
`markername` also requires an audited rsID mapping.
The exact Han 2020 main-UKB asthma file then brings the total to 38/45. Its
complete prior download, current EBI byte/header checks, upstream MD5, primary
paper, and archive README agree on a 64,538-case/329,321-control hg19 GWAS with
OR, confidence interval, P, EAF, INFO, and N fields. Source verification is not
an ancestry waiver: the paper says the main analysis included participants
regardless of ancestry and gives different 53,924/276,523 counts for its
European sensitivity analysis. The manifest therefore remains `MIXED`, and the
schema records the required log-OR-SE derivation from the confidence interval.
These counts also do not claim that every source is on GRCh37, binary
prevalence is resolved, or LDSC is ready. Shared samples are retained as
explicit analysis concerns: IBD/CD/UC overlap, the three lipid traits share
cohorts, and BMI overlaps UKB-derived sleep phenotypes.

## Deferred branch candidates

These records were not promoted:

| Trait | Reason |
|---|---|
| major depression | Branch substituted Wray 2018 for the selected, larger Howard 2019 release because the latter is controlled. Ease of access is not the selection criterion. |
| schizophrenia | Accession/publication/phenotype notes conflict and require primary-source reconciliation. |
| bipolar disorder | Branch substituted a roughly five-fold smaller 2016 GWAS for Mullins 2021. |
| ADHD | Branch substituted a later accession whose phenotype and relationship to the locked Demontis 2023 source require review. |
| multiple sclerosis | Branch substituted a roughly ten-fold smaller 2016 German GWAS for IMSGC 2019. |
| type 2 diabetes | Branch substituted Scott 2017; although genome-wide, its bare coordinate marker lacks the rsID mapping needed by the current LDSC path. |

No branch candidate record existed for snoring, sleep apnea, colorectal cancer,
lung cancer, or melanoma. Snoring, colorectal cancer, and lung cancer were
resolved independently from official GWAS Catalog sources, and sleep apnea was
resolved from FinnGen's R9 manifest, endpoint definition, and data dictionary.
The selected Landi melanoma accession `GCST010304` has `fullPvalueSet=false` in the official API
and no directory in the Catalog summary-statistics FTP tree, so it remains
pending rather than being represented by an invented download. Sleep apnea
requires a validated build-conversion path; melanoma still requires new
file-level primary-source work.

The rule for the next pass is unchanged: prefer the largest appropriate and
scientifically compatible release, even when it requires controlled access,
documented liftover, or explicit variant-ID mapping. A smaller public file is a
candidate or sensitivity dataset, not an automatic replacement.
