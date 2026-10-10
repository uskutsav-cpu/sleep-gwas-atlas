# Candidate release schema

All tabular files are UTF-8, tab-delimited, with one header row. Empty cells mean the field is not applicable or was unavailable in the frozen source; they are not zero. Numeric estimates are serialized from full-precision V4 result tables. Confidence intervals are marginal 95% intervals.

| File | Rows | Contents |
|---|---:|---|
| `core_rg.tsv` | 396 | Frozen core correlation family; q-values use the complete 396-pair BH family. |
| `extension_rg.tsv` | 1,200 | Frozen extension correlation family; q-values use the complete 1,200-pair BH family. |
| `external_validation_rg.tsv` | 41 | Estimated outcome-side external checks; threshold labels refer to the original 217-pair Bonferroni family. |
| `heritability.tsv` | 158 | All standalone core, extension and validation h² estimates, with observed/reported scale fields kept distinct. |
| `sensitivities.tsv` | 62 | Prespecified fitted-specification sensitivities; deltas are descriptive and difference uncertainty is not calibrated. |
| `sleep_constructs.tsv` | 12 | Sleep-source construct/coding metadata, with unknown reliability and overlap caveats retained. |
| `validation_217_class_counts.tsv` | 4 | Candidate counts by frozen/current validation class. |
| `outcome_phenotypes.tsv` | 158 | Core, extension and estimated validation-outcome labels and definitions for deterministic table joins. |

The RG field `discovery_outcome_trait_id` links estimated external validation rows back to the original extension family. `outcome_source_id` links each row to the corresponding phenotype definition/source record. See `schema.json` for the field allowlists and `release_manifest.json` for row counts and SHA-256 digests. No file contains local paths or source-body hashes.

## Record keys and joins

Use `(analysis_family, trait_id)` to join `heritability.tsv` to `outcome_phenotypes.tsv`; require `source_id` to agree across those two records. For correlation tables, use `(analysis_family, sleep_trait, outcome_trait)` as the pair key and verify the sleep and outcome source IDs against the metadata for that row's family. To link an external estimate to its discovery pair, match `(sleep_trait, discovery_outcome_trait_id)` to `(sleep_trait, outcome_trait)` in `extension_rg.tsv`. The validation and discovery outcome source IDs can differ; check each against its own family-specific metadata rather than requiring equality across sources. Source IDs alone are not unique phenotype identifiers. Sensitivity rows use `(job_id, kind, sleep_trait, outcome_trait)`; `job_id` alone is not unique.
