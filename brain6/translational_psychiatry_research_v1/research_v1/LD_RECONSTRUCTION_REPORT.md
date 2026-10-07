# Actual independently sourced LD reconstruction

Three new full-window signed LD references were reconstructed from actual
public 1000 Genomes Phase3 genotypes. This is new empirical reference work,
not a rescaling of the failed UKB/LAVA factors and not a synthetic result.
The historical32/32 failed LD-matrix status remains immutable.

|GRCh37 region|EUR samples|admitted nonpalindromic MAF≥.01 SNPs|effective rank|numeric status|
|---|---:|---:|---:|---|
|chr5:103447968–104447968|503|2694|502|PASS_REFERENCE_NUMERIC_QC|
|chr6:100630147–102636772|503|4878|502|PASS_REFERENCE_NUMERIC_QC|
|chr11:112459489–114257728|503|4627|502|PASS_REFERENCE_NUMERIC_QC|

All matrices are finite symmetric signed ALT-dosage correlations. Maximum
unit-diagonal error is1.9984e-14; maximum absolute r exceeds1 by at most
1.88e-14 (floating-point roundoff within the frozen1e-12 tolerance).
No reconstructed matrix values or diagonals were modified to repair the
reference. Native packages handle roundoff-level negative eigenvalues at
their documented numerical tolerance; the genotype-derived source LD
matrix remains unchanged. An independent scalar
Pearson implementation and SciPy Pearson recomputation on100 deterministically
selected actual SNP pairs per region agree within4.72e-15. The p>503 matrices
are necessarily rank-deficient; their Gram spectra are PSD and rank502.
Numeric PSD cannot prove GWAS/LD compatibility or establish fine-mapping
calibration. Source-qualified trait contracts, frequency/INFO/N checks and
pre-fit summary-statistic consistency remain separate admission gates.

## A decisive chr5 identity sensitivity

Brain6 lead rs2431108 corresponds to chr5:103947968:T:C; prior Zu2026 lead
rs77960 is chr5:103964585:G:A. Both exact allele identities were checked
against the primary EnsemblGRCh37 variation endpoint and the actual reference
VCF. rs2431108's current dbSNP mapping also includes T/A; that alternate was
not collapsed into the T/C dosage used here.

Across503EUR samples, signed ALT correlation is0.9977377424798654 and
r²=0.9954806027688181. CEU99, FIN99, GBR91 andIBS107 each have r²=1.0;
TSI107 has r²=0.9787615361749753. Different lead IDs therefore provide no
credible basis for a distinct chr5 signal in this reference. This does not
prove a single causal variant, shared causal effects or independent disease
replication; those require conditional inference and new trait data.
Exact frequencies/alleles/population estimates are retained in
`ld/chr5_prior_vs_brain6_lead_LD.tsv`.

## Source, acquisition and inclusion audit

The current official EBI release20130502 autosomal v5b was used. Metadata,
including the2025 corrected pedigree and2020 known-issues list, total754117
bytes. Fifteen exact HTTP206 byte ranges/header segments total26443591
compressed bytes; no full chromosome was downloaded. Every response range,
source URL, date, ETag, Last-Modified, SHA256 and byte count is recorded.
Current remote TBI SHA256 agrees with the retained original indexes for all
three chromosomes. VCF header specifies hs37d5 reference and assembly=b37.

All2504 unique VCF IDs match the official sample panel. Five503EUR population
sample counts and ordering are saved. Current pedigree records contain no
selected-EUR parent/child/sibling/second/third-order edges; every selected
sample carries phase3-genotype flag1. Known2020ASW and2025PEL corrections do
not affect the selectedEUR set. This is a source/pedigree audit, not a newly
computed genome-wide kinship analysis.

Only PASS diploid biallelic A/C/G/T SNPs, complete EUR genotypes, observed
MAF≥.01 and unambiguous primary nonpalindromic alleles were retained. Every
excluded record has an explicit reason. Frequencies computed directly from
genotypes agree with source-rounded EUR_AF to≤4.991e-5. Genotype counts are
ALT0/1/2; all signed matrices share exact saved variant order.

## Preservation and reproducibility

`LD_SOURCE_MANIFEST.json` and `LD_VALIDATION.tsv` link exact evidence. The
three compressed genotype-factor files are small self-contained research
artifacts, permitting exact full signed-LD regeneration. Full float64
matrices are retained locally under ignored`work/ld_genotypes_research_v1/`;
checksums are preserved. A factor-based rebuild avoids reacquisition.

The original UKB100000/LAVA0.1.5 BCOR/INFO matrices and original normalized
dense trait SQLite snapshots still point to the unmounted ExtremeSSD. No
factor-based historical reconstruction has been promoted. New503EUR source
is a much smaller external reference with its own sampling uncertainty; it
cannot silently replace the original scientific experiment.

New native fine-mapping was governed separately by
`ld/FINEMAP_PRE_FIT_PROTOCOL.md`; no result was admitted from numeric LD passage
alone. Exact original source identities were restored, 2185 shared SNPs
passed allele/AF/INFO/N gates, and both native RSS diagnostics passed the
locked operational consistency gate. Eight fixed native SuSiE fits, twelve
coloc-SuSiE comparisons and three same-universe ABF sensitivities were
executed. Primary H4 is0.991568229456347 and lower-prior H4 is
0.921629512291253; each trait has one qualifying credible set, containing
70insomnia and8ADHD variants. This is exploratory source-qualified support
for a known region, not independent GWAS replication or a causal variant.
`ld/FINEMAP_RESULTS_AND_LIMITATIONS.md` and complete numerical tables
retain the binary-model, effective-N, finite-reference and source-overlap
limits. No historical gate or result was changed.
