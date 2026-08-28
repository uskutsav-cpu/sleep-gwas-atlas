# Integrated atlas and immutable release workflow

The atlas release is deliberately impossible to create until every scientific
gate other than the release itself passes. File presence alone is insufficient.
`config/atlas_table_schema.json` freezes the exact ten-table field order,
primary keys, row-count minima, probability constraints, and the non-result
markers that are forbidden in scientific output.

The integrated validator checks the exact 45-trait and 396-pair families,
autosomal locus and variant coordinates, allele validity, primary-key
uniqueness, probability bounds, and foreign-key links from variants through
regulatory elements, genes, cells, pathways, causal tests, and graph edges.
Passing this validator establishes structural traceability; it does not turn an
association or annotation into causal evidence.

The robustness table has one row per major conclusion and each of the nine
locked robustness families. Every applicable row must point to a real local
evidence file with its exact SHA-256. Non-applicability requires an explanation,
and a material contradiction cannot remain unresolved.

Validation commands are safe and read-only:

```bash
python3 scripts/52_validate_integrated_atlas.py
python3 scripts/53_validate_robustness.py
python3 scripts/99_atlas_acceptance.py
```

The release builder first requires a clean tracked worktree and 22 passing
non-release gates. It has no overwrite mode. A real release is created only by
an explicit execution flag:

```bash
python3 scripts/54_build_release.py
python3 scripts/54_build_release.py --execute
python3 scripts/55_validate_release.py
```

`releases/atlas-v1.0` contains the exact panel, 45-row source provenance,
pinned tool registry, configurations and analysis code, QC logs and final
tables, the producing Git commit, a byte-and-SHA-256 manifest, and a canonical
`checksums.sha256`. Creation happens through a staging directory followed by
an atomic rename. An existing release is never replaced in place.

Raw and harmonized GWAS files are not redistributed in the release. Their
registered source identities and checksums, canonical-data checksums in
`traits.tsv`, and generating code are retained instead.
