# Remaining source-access blockers

All 45 locked traits now have a selected-source record and a schema record.
This ledger describes why three records remain `SOURCE_PENDING`. A paper,
accession, request page, or top-hit extract is not sufficient evidence for
`SOURCE_VERIFIED`: the exact genome-wide file must be obtainable, identified,
integrity-checked, and mapped to its literal schema.

## Multiple sclerosis (`ms`)

Selected scientific dataset: the 15-cohort genome-wide discovery meta-analysis
from IMSGC 2019, represented by [GWAS Catalog
GCST009597](https://www.ebi.ac.uk/gwas/studies/GCST009597) and the [primary
paper](https://pmc.ncbi.nlm.nih.gov/articles/PMC7241648/).

- The genome-wide discovery analysis has 14,802 cases, 26,703 controls, and
  8,278,136 analyzed SNPs.
- The paper-wide 47,429-case/68,374-control total adds MS-Chip and ImmunoChip
  replication at 4,842 prioritized effects. It is not a genome-wide sample size.
- The Catalog reports `fullPvalueSet=false` and European plus
  ancestry-not-reported samples.
- The paper identifies controlled cohort data and offers only ANZGENE,
  Rotterdam, and Berkeley cohort summary statistics upon request. Those three
  components are not the complete discovery meta-analysis.
- The [current IMSGC request
  form](https://imsgc.net/request-summary-statistics/) lists “MS Chip (Science
  2019),” which is the targeted follow-up rather than a clearly identified full
  discovery GWAS.

Required resolution: obtain a consortium-confirmed complete discovery
meta-analysis file or make an explicitly reviewed dataset substitution. Then
record the immutable file identity, byte count, SHA-256, build, ancestry,
sample-size behavior, and literal allele/effect schema. No request has been
submitted.

## Type 2 diabetes (`t2d`)

Selected scientific dataset: the full BMI-unadjusted European DIAMANTE
meta-analysis from Mahajan 2018b, available through the [official DIAGRAM
download page](https://www.diagram-consortium.org/downloads.html), form route
`check_check.Mahajan2018b.T2D.3.php`.

- The selected analysis has 74,124 cases and 824,006 controls.
- The official README names
  `Mahajan.NatGenet2018b.T2D.European.gz`, states GRCh37, and documents
  forward-strand EA/NEA, EAF, log-odds Beta, SE, Pvalue, and Neff.
- DIAGRAM publishes SHA-1
  `112e0e789179e15049649151163f5a517ff3ed7e` and MD5
  `53855134fd0862c402cf49ab6233ed6d` for the download.
- The form requires the downloader to affirm that they will not attempt
  re-identification. The README also prohibits reposting.

Required resolution: the user must decide whether to accept the official terms.
After an authorized download, record the actual archive name, bytes, SHA-256,
container integrity, member header, and complete row/schema audit. The
repository must not redistribute the data. No checkbox has been affirmed.

## Cutaneous melanoma (`melanoma`)

Selected scientific dataset: the histopathologically confirmed-only GWAS at
[dbGaP phs001868 / pha004971](https://www.ncbi.nlm.nih.gov/projects/gap/cgi-bin/analysis.cgi?study_id=phs001868.v1.p1&pha=4971),
described by [Landi et al.
2020](https://pmc.ncbi.nlm.nih.gov/articles/PMC7255059/).

- The selected genome-wide analysis has 30,134 confirmed cases and 81,415
  controls. The headline total of 36,760/375,188 additionally includes
  self-reported UK Biobank and 23andMe strata and has no full public release.
- dbGaP release notes name the full archive
  `phs001868.MelanomaSusceptibility.analysis-PI.MULTI.tar.gz`, but the study
  routes it through Authorized Access.
- The public FTP file `phs001868.pha004971.txt` is a truncated top-hit
  companion: 1,200,378 bytes, 16,898 data rows, and SHA-256
  `2074f42190ec3f3b004e54823054512e245ae4c59352eefd0ea1c3f02ab75073`.
  It cannot support LDSC or other genome-wide analyses.
- dbGaP states GRCh38 and labels the analysis Multi Ancestry, whereas the
  primary methods describe European-PCA sample QC. The controlled archive must
  resolve that metadata conflict and will require a validated GRCh38-to-GRCh37
  path if retained.

Required resolution: an eligible institutional requester must obtain access to
the full archive. Afterward, verify archive/member hashes, the literal schema,
per-variant N behavior, ancestry, and build before harmonization. No data-access
request has been submitted.
