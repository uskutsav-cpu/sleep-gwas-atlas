# Brain6 × AlphaGenome pilot

Predictions and derivatives are distributed under the [AlphaGenome Output Terms and binding use notice](TERMS_OF_USE.md). They are research predictions, not measured biological effects or clinical risk estimates.

**State: AVI scoring and frontal-cortex reference/alternate AlphaGenome predictions completed for all 425 selected variants across six regions.** This extension does not alter Brain6's frozen claims, evidence classes, inputs or final tables.

## Results

- `results/jbrowse_predictions/FRONTAL_CORTEX_REPORT.md` and `six_region_predictions.tsv`: actual AlphaGenome predictions through the public JBrowse service, for one AVI-selected variant in each region.
- `results/avi_mirror/AVI_PILOT_REPORT.md` and `variant_evidence_avi.tsv`: six-region overall AVI results and joined evidence, from a pinned third-party mirror. See its provenance and spot-check limitations below.
- `results/CHR5_CASE_STUDY.md` and `results/chr5_pip.svg`: exploratory chromosome 5 case study and plot.
- `results/regions.tsv`: all six insomnia–ADHD intervals, per-trait set statuses, variant coverage, and colocalization prior sensitivity.
- `results/variant_evidence.tsv`: union of the published per-trait 95% single-signal ABF credible sets, with source beta/SE, conditional PIP, membership, genome build mapping and Atlas links. The unit is candidate × SNP. A beta is for source A1; mapped REF/ALT does not change its sign. Empty PIPs are missing; `NO` membership alone is not evidence against association. The region table retains the per-trait gate status.
- `results/brain6_final_*.tsv`: unchanged, pair-specific context rows from the final package for candidate status, regulatory annotation and component QTL tests. Coordinate overlap and nearest genes are descriptive. QTL tests where present are component-GWAS tests, not pair-level colocalization.
- `results/provenance.json`: source and output hashes; `sequence_cache/` records UCSC hg19 and hg38 sequence responses used to check reference bases.

Each candidate is analyzed separately even if coordinates or variants overlap. The 425 rows (71 at the chromosome 5 case-study locus) are **selected credible-set variants**, not all GWAS variants. All 425 passed the mapping gates in this run. Genomic mapping does not independently prove the rsID or a causal effect. The underlying posteriors assume a single causal variant per trait, are conditional on the specified model, and are not SuSiE PIPs. Brain6's primary local genetic correlation family failed its QC gate; no protected candidate has a positive validation class or independent two-trait locus replication. Results are research hypotheses, not clinical guidance.

## Reproduce preparation

Run from a checkout at source commit `303c905c229104b7bf623df6de7431abb87e0727`, with the `brain6/alphagenome_pilot` directory present. The original SSD inputs are not included in the GitHub repository; the script checks each input against the frozen manifest and checks the full PIP file against the final provenance. The chain and UCSC reference sequence are also hashed. Install `pyliftover==0.4.1` in a suitable Python environment. If Python's CA setup is missing on this Mac, set `SSL_CERT_FILE=/etc/ssl/cert.pem`.

```bash
python3 brain6/alphagenome_pilot/prepare.py \
  --external '/Volumes/Extreme SSD/brain6-work/exploratory-finemap-coloc-v2' \
  --chain '/Volumes/Extreme SSD/brain6-work/independent-ld-reference-v1/hg19ToHg38.over.chain.gz' \
  --out brain6/alphagenome_pilot/results
Rscript brain6/alphagenome_pilot/plot.R brain6/alphagenome_pilot/results
python3 brain6/alphagenome_pilot/report.py
python3 brain6/alphagenome_pilot/test_pilot.py
```

The prepare script caches UCSC sequence responses. It uses unique same-chromosome chain mappings, resolves source alleles against the GRCh37 reference, handles chain strand, and requires the mapped GRCh38 reference to agree. Palindromic, non-SNV, nonunique, mismatched or changed-reference variants are held. Cached sequences should be preserved along with the provenance receipt. Check input and output hashes whenever the pipeline is rerun. `test_pilot.py` verifies the included run, not the remote Atlas service.

## Obtain and run Atlas scores

1. Visit [AlphaGenome](https://deepmind.google.com/science/alphagenome/) to complete the account profile, review and accept the applicable terms yourself, and request an API key. Do not put the key in source control or chat.
2. In a Python environment with enough disk space install the [official SDK](https://github.com/google-deepmind/alphagenome), for example `python -m pip install alphagenome==0.10.0`. Set `ALPHAGENOME_API_KEY` securely in the environment of the process running the commands. The pilot's SDK integration has been checked against SDK source and mocked locally; it has **not** been verified with a live key.
3. Retrieve track metadata using `python brain6/alphagenome_pilot/atlas_score.py metadata --out /path/to/new-metadata-dir`. Select the available brain-related ontology IDs and relevant scorer names from `scorer_metadata.json`. Tissue track availability is a property of the service; avoid assuming that a proxy tissue establishes brain-cell specificity.
4. Start with chromosome 5 using a fresh output directory. Replace both placeholders with IDs present in metadata and add repeatable options as needed:

```bash
python brain6/alphagenome_pilot/atlas_score.py score \
  --evidence brain6/alphagenome_pilot/results/variant_evidence.tsv \
  --out /path/to/new-chr5-run --scorer SCORER_NAME --ontology ONTOLOGY_ID
```

Use `--limit 2` for an initial service smoke test. To include all six regions use `--region all` after reviewing the smoke test. Each distinct ready variant is queried once per run. The tool writes the configuration, SDK version, all returned per-score track and observation metadata, raw and calibrated scores when present, and a status for each query. It does not silently convert a failed or empty query into zero. Inspect `query_status.tsv` and review errors before interpretation; a new directory is required for a rerun. The service's model release is not exposed by this client and is recorded as unknown. Review scorer-specific definitions before comparing signs or magnitudes; expression, chromatin and splicing scores use different semantics. This script does not recalculate colocalization or turn an AI score into evidence of disease causation.

**Terms:** As of the review, the [official AlphaGenome repository](https://github.com/google-deepmind/alphagenome#terms) says that, except as expressly provided for specified downloadable artifacts, API outputs and other Atlas information are for non-commercial use and should not be used to train other machine-learning models. Check the current terms and applicable access route before use. This pilot does not train or fit any ML model on AlphaGenome outputs. Its Atlas links enable manual browsing but account onboarding may still be required.

After scoring, create a long-form join while retaining failed, held and unqueried rows:

```bash
python brain6/alphagenome_pilot/integrate.py \
  --evidence brain6/alphagenome_pilot/results/variant_evidence.tsv \
  --run /path/to/new-chr5-run --out /path/to/new-integrated-output
```

`atlas_query_status` describes the actual run; the original `alphagenome_status` column records the earlier preparation snapshot. Track and observation metadata stay attached to each score. The integration receipt hashes the input and returned score files. This join is not a composite causal score.

## Overall AVI mirror query

With the user's authorization to use the mirrored data under its terms, `fetch_avi.py` queries the publicly hosted, pinned re-encoding of the permissive static AVI SNV download using DuckDB HTTP range requests. The output covers 425/425 selected candidate–variant rows. It contains overall AVI raw scores and PHRED **intervals**, not tissue- or gene-specific AlphaGenome scores. A PHRED value is a genomic rank and is not a probability of pathogenicity or shared causality.

The mirror's README describes a long layout, whereas its pinned files have `alt0`, `alt1`, `alt2` columns. The script assumes alphabetical alternative alleles excluding REF and fails on four reference-base spot-check anchors against published original-file examples and independent Broad results. The C and G anchor checks from Broad are rounded to three decimal places; matching mirror values were preserved at source precision. An additional independent Broad lookup for rs2431108 agrees after rounding. These checks do not establish correctness of every record, and the original whole-genome source file was not downloaded or hashed locally. PHRED is reconstructed from the mirror's published knot intervals; no midpoint or cutoff is imposed.

Run `fetch_avi.py` with DuckDB and its `httpfs` extension installed in an explicit scratch directory; `--out` must be a new directory. Then run `avi_report.py` if the output is at `results/avi_mirror`. Record the mirror revision, source checksum claimed by the mirror, license text, output hashes, and reproduction environment. Do not train genetic variant-effect models on the scores. Atlas tissue-specific scores remain separate and require the API workflow above.

`python3 brain6/alphagenome_pilot/compare_rankings.py` builds `RANK_COMPARISON.md`, `variant_ranks.tsv` and `ranking_overlap.tsv`. These compare AVI and existing single-trait ABF PIP ranks within each selected credible-set union, with original association-posterior gates retained. They do not test significance or train a model; smaller sets and ties have explicitly reported denominators.

`results/avi_mirror/BROAD_CROSSCHECK.md` documents an additional one-variant-per-region check against Broad SpliceAI Lookup. It verifies access and allele agreement for those six variants only; both lookups derive from the same AlphaGenome data, so this is not independent biological replication. Its raw results and source-binding hashes are saved alongside the report. `validate_avi_regions.py` respects the service's published rate guidance; reruns require network access.

`python3 brain6/alphagenome_pilot/molecular_context.py` binds the existing GTEx v8 exact-lead lookup to the six AVI regions. See `results/avi_mirror/MOLECULAR_CONTEXT.md` and its two tables. QTL associations are measured source data and are not tissue-specific AlphaGenome predictions; nominal lead-only lookups do not constitute QTL colocalization or causal-gene evidence.

The automated tests now exercise a simulated Atlas score run end-to-end: a successful variant, a failed query, a mapping hold, metadata retention, integration, and preventing an error message from leaking a test key. This is an offline mock; it does not verify service authentication, quotas, track availability, or real prediction values.

For a local alternative to the API, [AlphaGenome Research](https://github.com/google-deepmind/alphagenome_research#model-requirements) says its model weights require acceptance of separate non-commercial model terms and recommends at least an NVIDIA H100 GPU. The inspected computer has an Apple M1 GPU, and its main disk had about 181 MiB free during this check, so local inference has not been attempted. The API remains the prepared route when the user has obtained access and sufficient space for the SDK.

## Frontal-cortex predictions

With explicit user permission for non-commercial research, the public JBrowse service was used to predict one highest-AVI variant per region (six total). `results/jbrowse_predictions/FRONTAL_CORTEX_REPORT.md` describes the real reference and alternate DNase and RNA tracks, their 501 bp local summaries, exact track metadata, and limitations. That earlier report is a six-variant feasibility pilot; the later full expansion covers all 425 selected variants. Predictions are not measured associations or causal evidence. The JBrowse service did not expose a model release identifier. Signed service URLs are deliberately excluded; only 2,001 bp around each variant and source-byte-range hashes were retained. Reuse of AlphaGenome API outputs is subject to its non-commercial and other output terms; this pilot does not train an ML model on them.

`extract_jbrowse_window.py` extracts a window from a **ready** JBrowse response with live signed URLs, `run_jbrowse_five.py` retains request handles and avoids duplicate submission after an unknown outcome, `plot_jbrowse.R` plots the case-study tracks, and `build_jbrowse_report.py` builds the six-region report. The signed response is scratch data with expiring URLs; a new approved service query is needed to reproduce extraction after expiration.

## Full tissue expansion

The full-set analysis is stored under `results/jbrowse_expanded/`. Its live or final report and coverage table explicitly distinguish completed variants from ones not yet extracted. The evidence table retains all 425 input variants and their original evidence, and adds four local track summaries. The 2,001-bp float32 windows are compressed per variant; `decode_expanded_window.py` converts an individual folder to a TSV. The JBrowse model release is not exposed, and the local RNA window sum is not a gene-level score or a calibrated association test.

`expand_tissues.py` checkpoints a submitted job before polling and retains an unknown-submission marker if a network error could have occurred after acceptance. Resume with the same scratch and output paths; do not delete submission markers or automatically resubmit an unknown request. Only run concurrent workers for **disjoint** region sets, with no worker permitted to advance into another worker's region; there is no cross-process lock. Honor service errors and quota responses. On this computer it stops if free disk space falls below 80 MiB. `summarize_expanded.py` audits completed compressed windows and rebuilds the coverage/evidence tables; `test_expanded.py` independently verifies values, checksums and input identity.
