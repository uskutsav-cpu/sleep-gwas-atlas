# Candidate markers near published FFS loci

I compared the 15 candidate positions with the 37 conditionally identified FFS lead variants in the Ye et al. supplement. For each candidate, `candidate_nearest_ffs_lead_ld.csv` lists the closest reported lead on the same chromosome; for the four pairs within 500 kb I queried the Ensembl GRCh37 1000 Genomes phase 3 European pairwise LD endpoint. This is a proximity screen against **those 37 reported leads**, not a genome-wide search for every FFS signal or proof that more distant variants are independent.

| Candidate | Nearest reported FFS lead | Reported position difference | EUR reference r² | Interpretation |
|---|---|---:|---:|---|
| rs323509 | rs7703746 | 174,842 bp | 0.484 | Moderate reference LD; candidate may tag part of the previously identified FFS region |
| rs10865954 | rs71080535 | 404,064 bp | 0.321 | Moderate reference LD; exact indel representation needs review |
| rs11152363 | rs660010 | 253,781 bp | No result returned | Cannot infer r²=0 or independence |
| rs1942262 | rs660010 | 437,652 bp | No result returned | Cannot infer r²=0 or independence |

The source’s lead positions for rs7703746 and rs660010 match the Ensembl GRCh37 coordinates. For rs71080535, the Ye clumping table represents a G/GA insertion at chr3:49,616,053, whereas Ensembl represents an insertion beginning at chr3:49,616,054 with alternate strings including A; this is consistent with an anchored versus unanchored insertion representation, but the precise study/reference allele mapping is not independently confirmed. Both rs323509 and rs7703746 have multiallelic Ensembl records, as does rs10865954; the LD service returns rsID-level summaries without sufficient detail here to verify the exact study allele contrasts. Thus the r² values are approximate reference-panel context, not allele-specific evidence of same signed effects.

Because Ye's original FFS GWAS is also the FFS source used for the candidate associations, this comparison is **not independent replication**. LD does not establish a shared causal variant or that FFS and insomnia colocalize. The chromosome 5 pair is worth analyzing alongside the previously noted insomnia catalog variant rs1592757, but pairwise r² values do not establish a three-variant haplotype or colocalization. Score development should use the actual target ancestry's genotype reference, normalized variants and a prespecified LD-aware procedure, rather than simply keeping or dropping these markers on this exploratory lookup.

The 15 candidates were post-selected, and no individual score validation or genetic risk probability is produced. The original FFS locus list and its source workbook are documented in `fi_ffs_risk_followup.md`; API URLs and response hashes are included in the CSV. No records were imputed for empty reference queries.
