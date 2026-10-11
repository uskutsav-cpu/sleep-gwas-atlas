# Exact source-linked candidate lookup

This small standard-library Python reader exposes existing frozen evidence for researchers. It performs no GWAS/LDSC fit and changes no estimate, QC rule, threshold or source record. The candidate has not been cleared for public deposit.

Use Python 3 and supply the existing candidate directory:

```sh
python3 query_atlas.py --release /path/to/release_candidate search 'abdominal pain'
python3 query_atlas.py --release /path/to/release_candidate estimate --family validation --sleep insomnia --outcome gwas_catalog_GCST90479148
python3 query_atlas.py --release /path/to/release_candidate estimate --family core --sleep accel_sleep_duration --outcome hdl
```

Search returns metadata candidates, with family and exact trait IDs. A text match does not mean phenotype equivalence. Estimate retrieval requires an exact family/sleep/outcome key; source agreement is enforced when joining. Output includes the source construct and coding, frozen family/QC classifications, marginal uncertainty and relevant sensitivity estimates. External outcome estimates explicitly link to discovery and disclose sleep-source reuse. Absent rows produce an abstention rather than a biological-null classification.

The six scientific examples demonstrate evidence retrieval, including a source-reused validation estimate, an HDL sensitivity, an ineligible melanoma h² record and a held bowel-source record with blank signed fields. They are execution examples, not evidence of human usability or newly discovered associations. `scientific_query_examples.json` records the returned fields. V8 queries read the required hold from the release directory; a missing required hold fails retrieval. The old-release fallback is the adjacent tables/source_direction_hold.json. Use the V8 interpretation view for current work.

`build_interpretive_release.py` builds the separate V8 candidate view from the preserved historical candidate and reviewed V8 tables. It withholds 36 signed core results without altering the original complete testing families or historical files. `plot_core_interpretive.R` shows the new masked core view; `plot_transport.R` shows three source-defined estimates. Both use base R, require no added packages and perform no statistical estimation.
