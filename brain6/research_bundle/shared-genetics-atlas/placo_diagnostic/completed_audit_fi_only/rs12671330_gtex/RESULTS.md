# GTEx v8 context for rs12671330

For exact GRCh38 variant **chr7:114693988 T>C**, the GTEx v8 significant single-tissue eQTL and sQTL endpoints each returned zero records across the tissues they expose. Separately, dynamic FOXP2 eQTL tests returned the following nominal results:

| Tissue | Sample count | ALT C normalized expression slope | Nominal p |
|---|---:|---:|---:|
| Cortex | 205 | +0.0840 | 0.0562 |
| Frontal cortex BA9 | 175 | +0.0449 | 0.192 |
| Adrenal gland | 233 | −0.0000654 | 0.999 |

These tests were chosen after examining a previously published candidate and are not a prespecified multiple-testing family. The dynamic endpoints return nominal tests, including nonsignificant ones; the cortex result is not a significant association. The normalized effect is an ALT-versus-REF slope on normalized expression, not a fold change. A missing significant-pair record or a nonsignificant test does not establish absence of an effect, especially in unmeasured cell types. Conversely, neither a nominal result nor the marker's position near FOXP2 establishes an effector gene, shared causal variant, or independent replication. Dense tissue-appropriate GWAS–QTL colocalization and experimental follow-up remain necessary.

`query_status.tsv` records all five successful requests, exact URLs and response hashes; `Eqtl_response.json` and `Sqtl_response.json` retain the static responses. Dynamic summaries are retained, but per-donor arrays in those API responses were not saved, so their hashes cannot by themselves reproduce the full response bytes. The [pre-query protocol](PROTOCOL.md) and `query.py` document scope and limitations. This is measured QTL context, distinct from the AlphaGenome model prediction and the calibration-limited FI diagnostic.
