# Developmental chromatin lookup for the shared chromosome 5 region

A targeted search of the previously used [human cortical neuron and progenitor dataset](https://link.springer.com/article/10.1186/s13059-023-02959-0) found 26 exact position-and-allele variant–peak records for three post-selected candidates. **None passed a descriptive Bonferroni correction across those 26 returned rows.** The smallest nominal p was 0.0448 for rs416223 in neurons; multiplying by 26 gives a value capped at 1. This post-hoc correction does not cover the broader project search, and the result is not evidence for or against a shared causal mechanism.

| Cell type | Candidate | Matching variant–peak rows | Smallest nominal p | Interpretation |
|---|---|---:|---:|---|
| Neuron | rs2431108 | 0 | Not available | No exact row in this selected file. |
| Neuron | rs77960 | 0 | Not available | No exact row in this selected file. |
| Neuron | rs416223 | 6 | 0.0448 | Nominal only; no returned row passes the 26-row benchmark. |
| Progenitor | rs2431108 | 7 | 0.253 | No returned row passes the benchmark. |
| Progenitor | rs77960 | 7 | 0.253 | No returned row passes the benchmark. |
| Progenitor | rs416223 | 6 | 0.521 | No returned row passes the benchmark. |

The two insomnia–ADHD candidates rs2431108 and rs77960 have identical reported p and beta for their seven progenitor peaks; they are strongly linked in the reference panel and are not independent confirmations. The progenitor peak `chr5_104610871_104612270` contains rs2431108 under the assumed zero-based, half-open peak convention, but its nominal p is 0.253. The source effect allele at that row is T; the reported alternate-C-aligned beta is −0.1006. The peak is measured developmental accessibility, but this association is not strong evidence of an allele effect and does not assign an effector gene. Coordinate and peak conventions should be confirmed before interpreting individual bases.

The downloaded files are a selected **Model1A** candidate set, not all possible brain caQTL tests. Missing neuron rows do not mean the variants were never tested elsewhere or have no effect; non-significant results do not demonstrate absence of effects. The [study](https://doi.org/10.1186/s13059-023-02959-0) used cortical progenitors and differentiated neurons, with 75 and 57 donors in the paired-modality subset. Those cultures differ from adult frontal tissue and the GTEx testis QTL context. The lookup is not a GWAS–QTL colocalization, a measured validation of AlphaGenome, independent locus replication, or a rescue of failed LAVA QC.

The neuron file contained 4,359,743 rows and progenitor file 4,645,460; both complete gzip reads passed source-member CRC and SHA-256 checks. The [protocol](PROTOCOL.md), [candidate rows](candidate_rows.tsv), [retrieval receipt](retrieval_receipt.json), [scan receipt](scan_receipt.json), [code](scan.py) and [independent count and arithmetic verification](verification.json) preserve the analysis. The selected Zenodo ZIP members were range-retrieved to the SSD; the 6.3 GB archive's full publisher MD5 was not verified. See the existing [testis expression follow-up](../chr5_dense_qtl/RESULTS.md) and [trait–trait locus atlas](../../locus_atlas/LOCUS_ATLAS.md).

A [regional-readiness check](regional_readiness/RESULTS.md) finds incomplete overlap with the previously matched GWAS files and no reported QTL standard errors or allele frequencies. A new colocalization is not admitted on these inputs.
