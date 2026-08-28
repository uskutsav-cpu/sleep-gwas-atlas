# Integrated atlas and immutable release workflow

The atlas release is deliberately impossible to create until every scientific
gate other than the release itself passes. File presence alone is insufficient.
`config/atlas_table_schema.json` freezes the exact ten-table field order,
primary keys, row-count rules, probability constraints, and the non-result
markers that are forbidden in scientific output. Positive-evidence tables may
be header-only after a genuine complete null analysis, but only when their
separate checksum-bound task coverage is exact and contains no unresolved
access blocker.

Validated null families are first-class release evidence. When no Genomic SEM
candidate passes the frozen held-out criteria, the release contains header-only
factor-GWAS and Q_SNP tables plus terminal provenance, not invented association
rows. A zero primary shared-locus family is handled analogously by the
fine-mapping and downstream locus tables.

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

## Complete production target

`atlas_v1_release` is the single full-DAG entry point. It explicitly requires
the terminal Genomic SEM, LAVA, MiXeR, integrated-atlas, and robustness outputs;
the latter two transitively require the complete pleiotropy, fine-mapping,
molecular/TWAS, and interpretation layers. From the current checkpoint its dry
run resolves 3,042 jobs before dynamic bivariate MiXeR and later checkpoint
families are expanded:

```bash
snakemake --cores 1 -n atlas_v1_release
snakemake --cores 16 atlas_v1_release
```

The second command belongs on the documented production host, not the current
laptop. Large acquisitions remain disabled by false-by-default acknowledgement
keys in `config/workflow.yaml`; each key authorizes only its named transfer.
The MiXeR 64-file reference is manual-stage-only and is checksum sealed before
any task manifest can be written. The release rule runs last and still requires
a clean tracked worktree and all 22 non-release acceptance gates.
