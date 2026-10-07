# Independent genotype LD reconstruction: locked research protocol

Locked 2026-10-07 before any new reference genotype/LD outcomes. The candidate
regions were chosen retrospectively from the inherited research program. This
is reference reconstruction and an explanatory sensitivity experiment, not
new GWAS discovery or independent two-trait replication.

## Sources and resource stop rules

Use the official EMBL-EBI 1000 Genomes Phase 3 GRCh37 release 20130502,
autosomal v5b VCFs, matching TBI indexes, sample panel and current pedigree
correction/known-issues READMEs. Public availability without an embargo is
documented by IGSR and the NIH-managed AWS registry. Never mix AWS v5a
compressed indexes with EBI v5b objects. Record exact source URL, ETag,
Last-Modified, bytes and SHA256 for each retrieved body/range.

Regions: chr5:103447968–104447968; chr6:100630147–102636772;
chr11:112459489–114257728, all 1-based inclusive GRCh37. These are the full
inherited C, B and A windows; no outcome-selected subwindows. Three reference
regions support four pair-region hypotheses; this protocol performs no
four-hypothesis validation test and does not change its 0.05/4 denominator.

Initially permit <=30,000,000 total compressed genotype/header bytes; metadata
<=1,000,000 bytes. Require >=1 GiB disk free immediately before acquisition,
and bound matrices to <512 MiB each with streamed genotype parsing. HTTP
range responses must be 206 with matching Content-Range; fail rather than
silently accepting a whole chromosome. Existing retained TBIs are used to
plan exact ranges, but current remote indexes must match SHA256 before use.

## Scientific/source admission

VCF header and source README must establish GRCh37, release, sample IDs,
phased genotype format, and processing history. Require exact unique VCF
sample IDs matching all 2504 panel entries. Use panel EUR populations CEU,
FIN, GBR, IBS and TSI. Audit the corrected pedigree; if included EUR samples
have confirmed first-degree relationships or unresolved identity problems,
do not label the resulting reference an unrelated EUR reference. Retain the
full sample provenance and report unresolved source caveats.

Use diploid biallelic A/C/G/T SNPs with PASS filter and complete GT in every
selected EUR sample. Count ALT alleles as 0,1,2; report every rejected record
with its reason. Require unique chr:POS:REF:ALT identities and unique rsID
where an rsID is provided. Require observed EUR MAF>=0.01. Exclude A/T and
C/G variants from the primary harmonization-ready set regardless of frequency;
retain a separate available-record audit. Order by POS,REF,ALT,rsID and save
that ordering and ALT frequency. No genotype imputation, LD clipping,
diagonal changes or PSD projection.

## Numerical diagnostics and interpretation

Construct signed Pearson LD directly from centered ALT dosage columns in
float64. Numerical acceptance requires finite entries, diagonal error<=1e-12,
max asymmetry<=1e-12, |r|<=1+1e-12, and a Gram spectrum consistent with PSD
within 1e-8 times the largest eigenvalue. A p>n reference is necessarily
rank-deficient: PSD is adequate for a reference diagnostic, but no claim of
positive definiteness or original SuSiE-gate repair is permitted. Check an
independent scalar Pearson calculation against the matrix implementation on
100 deterministically seeded SNP pairs plus the two chr5 leads; tolerance
1e-12. Save full genotype factors enabling exact LD regeneration. Check
subpopulation allele frequencies and rs2431108–rs77960 signed r/r² separately
in EUR, CEU, FIN, GBR, IBS, TSI. These population comparisons are exploratory
reference sensitivity with no disease-association P value.

## GWAS alignment and fine-mapping gates

If retained variant-level signed Z/alleles are present, audit exact rsID,
coordinate and allele concordance to REF/ALT; orient Z to ALT using exact,
swap, complement and complement-swap relationships, excluding palindromes.
Never treat an unsigned or allele-free PLACO output as aligned GWAS input.
Record missing frequency/INFO, trait-specific effect/SE/N and original source
receipts. Lack of original dense trait contracts and valid analyzed/effective
N semantics blocks SuSiE/coloc-SuSiE even if a new matrix is numerical PSD.
Full-window SNP coverage and GWAS/reference frequency/summary-statistic
consistency need their own prospective method-specific protocol before a
native fine-mapping experiment. Do not replace frozen historical 0/32 LD
failure status. No posterior analysis is authorized by numeric LD passage.

The rs2431108–rs77960 comparison can determine whether different lead IDs
strongly tag one another in this reference. Even high r² does not prove one
causal variant; low r² does not prove a new signal without conditional analysis.

Sources: https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ ;
https://www.internationalgenome.org/category/data-reuse/ ;
https://registry.opendata.aws/1000-genomes/ .
