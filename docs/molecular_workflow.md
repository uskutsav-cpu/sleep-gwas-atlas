# Molecular-QTL and TWAS workflow

This layer starts only after the complete primary fine-mapping locus family is
checksum locked. Its source, query, feature-selection, model, and variance-control
rules are frozen in `config/molecular_analysis_policy.json`. A molecular
association is evidence for a model-dependent gene hypothesis, not proof of
mediation, a causal gene, or a mechanism.

## Exact source family

- eQTL Catalogue release 7: every dataset using the locked `ge`, `leafcutter`,
  or `aptamer` quantification method, with no result-ranked dataset cap.
- The 49 GTEx v8 eQTL datasets imported by the eQTL Catalogue.
- Fenland cis-pQTL v1 as a plasma-protein replication source.
- PsychENCODE release 3 as supplementary brain evidence; a dated,
  checksum-bound `ACCESS_BLOCKED` outcome is valid when an exact file cannot be
  accessed without user-controlled credentials or terms acceptance.
- The 49-tissue GTEx v8 elastic-net PredictDB release containing gene-specific
  phi values. The older MASHR eQTL and sQTL archives are excluded because they
  lack the phi values required by the locked variance-control method.

The current GTEx portal API identifies the collapsed GENCODE v26 GRCh38 gene
model as exactly 134,502,408 bytes. The workflow pins its SHA-256 and uses it
only for interval overlap after mapping each locked GRCh37 locus to GRCh38 with
the pinned UCSC chain.

## Result-free locks and execution

Run the metadata audit first:

```bash
bash scripts/61_fetch_molecular_metadata.sh
python3 scripts/61_molecular_preflight.py
python3 scripts/62_prepare_molecular_search_plan.py
```

Each tabix task is executed by `scripts/63_run_molecular_search.py`. Curator-
mediated sources must instead be recorded with
`scripts/63_record_molecular_search.py`, together with the exact source snapshot
and either normalized data or a permitted terminal outcome. Only after every
planned source unit is recorded does `scripts/64_lock_molecular_features.py`
freeze the analyzable feature family. Scripts 65 and 66 then materialize and run
every trait-feature SuSiE-RSS/coloc-SuSiE comparison.

The result-free PredictDB inventory contains 98 files, 49 model databases and 49
covariance files, totaling 3,135,665,776 bytes (2.920 GiB). The downloader
refuses this large transfer unless the user has first been warned and the command
includes `--acknowledge-large-download`:

```bash
bash scripts/67_fetch_twas_resources.sh --models --acknowledge-large-download
bash scripts/67_fetch_twas_resources.sh --runtime
python3 scripts/68_index_twas_models.py --materialize
python3 scripts/69_prepare_twas_manifest.py --preflight results/tables/twas_preflight.json
```

The model index excludes each model-feature missing phi before association
results are read. S-PrediXcan receives the exact locked GWAS N and LDSC h2.
Traits with h2 outside `0 < h2 <= 1` are `NOT_APPLICABLE`; estimates are never
capped or substituted. Scripts 70–72 map the dense GWAS, run all eligible
trait-by-context tests, retain raw and corrected statistics, and apply BH FDR
within each trait/model family.

## Integration and absence handling

`scripts/73_collate_molecular.py` validates all upstream hashes and recomputes
the integrated tables. A QTL signal is supported only when the same signal pair
passes `PP.H4 >= 0.8` and `PP.H4/PP.H3 >= 5` at every locked p12 value. TWAS
evidence enters a locus only through exact GENCODE interval overlap and uses the
variance-controlled statistic.

The canonical outputs are:

- `results/tables/molecular_qtl_colocalization.tsv`
- `results/tables/colocalization.tsv`
- `results/tables/twas.tsv`
- `results/tables/molecular_locus_coverage.tsv`
- `results/tables/molecular_evidence.tsv`
- `results/atlas/genes.tsv`
- `results/atlas/molecular.provenance.json`

`genes.tsv` contains only genes with at least one supported molecular stream.
Every primary locus still has complete four-modality coverage. A locus with no
supported gene remains an explicit evidenced absence in
`molecular_locus_coverage.tsv`; the workflow never fabricates a nearest-gene or
placeholder row to satisfy the integrated-atlas schema.

Validate an existing complete layer without changing it:

```bash
python3 scripts/73_collate_molecular.py --validate-only
```
