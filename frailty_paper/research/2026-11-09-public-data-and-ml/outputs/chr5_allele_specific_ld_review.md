# Chromosome 5: allele-specific reference LD for insomnia and FFS candidates

I retrieved the three exact GRCh37 sites from the publicly available [1000 Genomes phase 3 phased chromosome 5 VCF](https://hgdownload.soe.ucsc.edu/gbdb/hg19/1000Genomes/phase3/), using indexed region reads, and selected the 503 individuals labeled EUR in the [phase 3 sample panel](https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/integrated_call_samples_v3.20130502.ALL.panel). All three reference records were biallelic in this VCF; all 1,006 sampled haplotypes had phased calls. This does not eliminate the broader multiallelic rsID ambiguity in other releases or in the GWAS sources.

I calculated the Pearson correlation of phased haplotype allele indicators for the specific alleles shown below (r² is its square):

| Pair of counted alleles | Signed r | r² | Both alleles on same reference haplotype |
|---|---:|---:|---:|
| rs323509 A and rs7703746 G | +0.692 | 0.480 | 306 / 1,006 |
| rs323509 A and rs1592757 C | +0.844 | 0.712 | 301 / 1,006 |
| rs7703746 G and rs1592757 C | +0.803 | 0.644 | 358 / 1,006 |

The signed correlations show that the listed alleles tend to occur together in this reference panel; the final column is an observed haplotype count, not a case risk. The allele-specific r² for rs323509–rs7703746 (0.480) is close to the earlier Ensembl rsID-level lookup (0.484); rs7703746–rs1592757 agrees closely with 0.643. Differences can reflect reference subset or LD calculation details, which were not established here.

The project FI and FFS sources associate rs323509 A with higher measured scores; the original UKB insomnia source associates A with higher insomnia liability. The published FFS locus table reports a positive coefficient for rs7703746 with `refA=G`; the interpretation of that field as the effect allele is still a source-format inference. The Ye supplement records a separate insomnia study associating rs1592757 C with higher insomnia odds (OR 1.02, p=2×10⁻¹⁰). These patterns and the signed reference LD are coherent as a **regional hypothesis**; they are not evidence that a single allele causes both traits or that these individual estimates replicate independently. The FinnGen and UKB insomnia/frailty data may overlap; the FI UKB and TwinGene component effects for rs323509 disagree in sign under the project's inferred adapter.

This reference check is distinct from a study-specific LD matrix, fine-mapping or a colocalization test. It cannot estimate individual risk, establish a causal gene, or justify combining these markers into a score. The summary counts, all allele frequencies, direct VCF URL and panel checksum are in `chr5_allele_specific_reference_ld.json`; no per-person genotypes were saved. Reproducible code is `work/chr5_reference_haplotypes.py`.
