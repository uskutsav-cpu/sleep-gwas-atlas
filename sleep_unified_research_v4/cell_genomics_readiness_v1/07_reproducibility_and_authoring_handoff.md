# Reproducibility and authoring handoff

## Frozen inputs and outputs

The candidate-release builder reads the completed, versioned inputs listed in `release_candidate/release_manifest.json`. Each input is identified by repository-relative name and SHA-256. It extracts an explicit field allowlist; it does not rerun analyses or alter the original tables. The full-precision source results remain in `sleep_unified_research_v4/tables/`. Raw GWAS source bodies remain on the external SSD and are not needed to rebuild this compact candidate if the listed result tables and source metadata are available.

The original result families remain separate: 396 core, 1,200 extension, and 217 validation candidates (41 estimated). The 62 sensitivities remain descriptive; neither their differences nor the Cov0 displays are calibrated cross-estimate tests. The 23-page V4 evidence figure set and its source index already exist. `figure_selection.tsv` identifies a compact main-figure selection and its source tables without rerendering completed figures.

## Rebuild and check

From the repository root:

```bash
python3 sleep_unified_research_v4/cell_genomics_readiness_v1/scripts/build_candidate_release.py
python3 sleep_unified_research_v4/cell_genomics_readiness_v1/tests/test_candidate_release.py
```

The builder uses Python 3 standard-library modules only. A clean-checkout run requires the listed committed result tables, `config/analysis_panel.tsv`, `discovery_extension/config/candidate_traits.tsv`, and the 12-row sleep construct dictionary. It does not require GWAS downloads, access credentials, R, PLINK, or an LDSC rerun. The candidate manifest gives input/output digests and builder digest. A different output location can be supplied with `--out`.

The pushed commit `a5c955ae` was also tested from a fresh source archive on the external SSD containing only committed builder inputs and package files. The builder and integrity test passed there, and all 11 generated candidate files matched the committed release byte-for-byte. No ignored raw-source data or workspace scratch files were needed.

## Scientific interpretation

The current evidence is a qualified association resource, not a new primary discovery or validated inference method. Independent result-level reviews found no newly generated association, locus, mechanism, calibrated contrast, or fully independent two-trait replication. Prior-art benchmarking bounds the novelty claim. Source-specific definition and harmonization limitations accompany the results. All correlation intervals are marginal within-fit intervals. The estimated discovery/validation pairs share the sleep GWAS and lack calibrated sampling covariance.

## Critical path to human manuscript authoring

The package can support a manuscript-writing decision after the independent benchmark, source harmonization, uncertainty review, and independent usability review are finalized and incorporated. Investigators then need to approve the central contribution and outlet; source owners/investigators must confirm rights, phenotype metadata, and participant overlap; the institution must determine applicable ethics requirements; and authors must provide authorship, funding, conflicts, AI-use, and deposition declarations. The candidate release must not be publicly deposited until source terms are cleared. The remaining core/extension raw-to-processed replay tails are optional provenance enhancements for the completed estimates, not a prerequisite to writing a carefully qualified resource paper.

No claim of journal competitiveness, new scientific discovery, causal mechanism, measurement reliability calibration, or independent replication is supported by this handoff. The benchmark supports only a qualified resource contribution for investigator evaluation.

Actual result-level reviews are preserved at `../../sleep_integrated_discovery_methods_v1/reviews/completed_results_scientific_interpretation_v1.md` and `../../sleep_integrated_discovery_methods_v1/reviews/genetic_epidemiologist_result_level_statistical_review_v1.md`. These are independent automated readings of the actual frozen outputs, not human approval. The source-specific review is `02_source_harmonization.md`; the prior-resource comparison is `01_resource_benchmark.md`.
