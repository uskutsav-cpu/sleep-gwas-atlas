# Brain6 v0.2: executable handoff

This extends a recovered v0.1 overlay. It has **not** been merged into a verified live repository. The original checkout, inputs, code, A/B/CONTROL, results, and locked testing families remain authoritative. The new six-disorder panel was chosen after the global screen; describe it as follow-up, not a retroactively preregistered discovery study.

## 1. Protect and identify the actual project

Apply the outer package from your real repository root. Default is read-only. `--upgrade` permits only exact v0.1 payload replacement. Any locally edited extension file conflicts and stops the installer. Core files and data are outside the install allowlist. Never force-push.

From `extensions/brain6`, install Python and test:

```bash
python3 -m venv .venv
source .venv/bin/activate
python -m pip install -e '.[test]'
python -m pytest -q
python -m brain6 audit --repo ../.. --out work/local-audit-v02.json
```

The executable package needs Python 3.11 or newer. Native statistical programs and reference assets are separately required. Python dependency ranges describe compatibility, not a validated native production lock. Four native smoke checks are deliberately skipped without their runtimes.

## 2. Rank measured correlations; do not invent sleep partners

Use README commands to map the **actual** final matrix columns, rank all 72 sleep–brain rows retaining original family FDR, and create the reviewed decision lock. Insomnia–ADHD is a protected legacy pair. The prior proposed SCZ/BIP/AD/PD partners were examples, not verified top rows. A disorder can be `NO_ELIGIBLE_PAIR`; six null decisions are now valid and remain visible. Do not force a positive AD/PD connection.

Copy templates into ignored `work/`, never replace existing core manifests. Run each version into a fresh output directory: code-version changes intentionally invalidate resume fingerprints. Legacy evidence must be imported with its original method/input identity, not silently recomputed or renamed PLACO+.

## 3. Shared data preparation and native discovery

`assemble-preparation` builds a shared normalization/pair/index/chunk DAG. `pipeline` runs synchronously; `workflow/Snakefile` can schedule the same tasks. `plan-placo` + `execute-placo` estimate genome-wide PLACO+ parameters once and execute bounded chunks. `deep-followup` binds native PLINK clumping, signed LD, allele matching, SuSiE and coloc. These native paths need real native executables and approved data, and were **not** executed in this delivery.

LDSC uses its native munge/h2/rg CLI. Local analysis uses the provided LAVA per-locus adapter: low local h² is not a program failure. The prior rejected LAVA analysis is never reused. HDL-L/SUPERGNOVA are **not implemented here**; no blanket switch to them is scientifically justified by a prior QC failure. GenomicSEM constructs actual S/V/I then fits a reviewed model. Its adapter is not a complete factor-GWAS campaign. The MATLAB adapter requires the official pleioFDR programs, valid aligned MAT inputs, reviewed settings and a license.

## 4. MR: exposure-only candidate selection, both directions

Run candidate selection on the full normalized exposure, **not PLACO hits**:

```bash
python -m brain6 mr-candidates --exposure work/preparation-v02/normalized_sleep \
  --root work/mr-v02 --name forward_candidates
```

Create a native clump job with `associations=.../forward_candidates/candidates.tsv`, `p_column=P`, and `p1` equal to the candidate threshold. Use `configs/mr_clump.template.json` (default r² ≤ 0.001, window ≥ 10,000 kb). The harmonizer enforces those separate instrument defaults unless explicitly reviewed overrides are supplied. Choose and document an **MR-specific LD r² and window**; do not blindly reuse the pleiotropy clumping threshold. Pass its accepted native artifact to:

```bash
python -m brain6 mr-harmonize --candidates work/mr-v02/forward_candidates \
  --clump work/mr-v02/forward_clump --outcome work/preparation-v02/normalized_disease \
  --root work/mr-v02 --name forward_harmonized
```

The output has exposure/outcome tables and attrition reasons; outcome P is not a filter. Bind these tables to `configs/mr.template.json` and execute the native MR adapter through `make-job` / `run-job`. Reverse direction starts again with the disease as exposure, using disease-selected instruments.

Copy `configs/mr_family.template.json`; bind pair ID, direction, native artifact fingerprint and harmonization receipt. `mr-family` validates trait/direction identity and native settings for empirical inputs. It preserves 2K family slots and applies BH to the primary IVW family without treating absent jobs as tested nulls. Method-specific sensitivity outputs remain separate. MR-PRESSO, CAUSE and Steiger were discussed previously but are **not automated in this package**. Agreement of estimators does not prove IV assumptions.

## 5. Replication: source and cohort independence

Copy `replication_family.template.json`. Every actual estimate needs a one-row `rg,se,p` result file and SHA, a source, a hash-bound cohort audit, power/QC and phenotype-match review. `DISEASE_ONLY` means the sleep GWAS is reused. The `discovery_cohorts` / `replication_cohorts` lists refer to the side(s) asserted independent; for BOTH_TRAITS include both sets. Disjoint labels alone do not prove individuals are independent. Unavailable, weak and opposite results remain distinct.

```bash
python -m brain6 replication-family --manifest work/replication.reviewed.json \
  --root work/replication-v02 --name summary
python -m brain6 mr-family --manifest work/mr-family.reviewed.json \
  --root work/mr-v02 --name family_summary
```

## 6. Full-cis molecular QTLs through native colocalization

Start with authorized, complete cis-summary exports, not significant-eQTL-only tables. Normalize source column names explicitly via the QTL source card. This indexer is for one tissue/study/molecular phenotype per source card. Multi-allelic variants, palindromes and duplicated feature-coordinate rows are conservatively excluded; review attrition and resolve source aliases in a documented upstream step if necessary. It does not silently lift genome builds or assume all QTL betas have the same scale.

```bash
python -m brain6 qtl-index --source work/qtl-source.reviewed.json \
  --root work/qtl-v02 --name cortex_expression
python -m brain6 qtl-locus --gwas work/preparation-v02/normalized_disease \
  --qtl work/qtl-v02/cortex_expression --feature ACTUAL_FEATURE_ID \
  --region work/actual-region.json --root work/qtl-v02 --name actual_locus
```

A cis-only pair is explicitly forbidden from genome-wide PLACO parameter estimation. GWAS N and QTL N are retained separately. A local pair can enter `index-pair` and `prepare-locus` using an actual signed LD manifest.

For all preselected feature/tissue/locus queries, configure `molecular_followup.template.json` and run:

```bash
python -m brain6 molecular-followup --settings work/molecular.reviewed.json
```

This executes QTL locus preparation, signed-LD alignment, trait-specific native SuSiE and `coloc.susie`, with resumable units and every failed/blocked query in the ledger. Supply actual QTL expression SD and study N; do not invent `sdY=1`. H4 maxima are descriptive; inspect all signal comparisons, H3 and sensitivity. QTL colocalization is not proof of mediation.

## 7. Cell annotation and cross-disorder comparison

Provide BED0 half-open regulatory intervals, source SHA, genome build and citation per cell class. The new scorer converts one-based SNP positions correctly, merges overlaps and reports omitted non-autosomal intervals.

Input variants require `pair_id,locus_id,SNP,CHR,BP,PIP,selected,stratum`. Include appropriately matched independent control loci, not just discoveries. A score is the fraction of total PIP mass overlapping the annotation, **not** a cell-specific causal probability. Conditioning on PIPs and selecting controls must be scientifically reviewed.

```bash
python -m brain6 annotation-scores --variants work/locus-variants.tsv \
  --annotations work/cell-annotations.reviewed.json --root work/cells-v02 --name scored
python -m brain6 cell-enrichment --table work/cells-v02/scored/scores.tsv \
  --out work/cell-enrichment.tsv --permutations 10000
```

`integrate` and `cross-disorder` operate on normalized sourced evidence tables. No raw scRNA/scATAC preprocessing, new spatial model, differential-expression experiment or wet-lab validation was performed.

## 8. Completion cannot be inferred from successful commands

```bash
python -m brain6 campaign-template --pair-lock work/brain6-pairs.lock.json \
  --out work/campaign.reviewed.json
python -m brain6 campaign-audit --manifest work/campaign.reviewed.json \
  --root work/campaign-v02 --name coverage
```

Six selected pairs produce 72 pair-stage slots plus GenomicSEM and conjFDR globals (74 slots). Secondary pairs increase the count. Slots stay explicit when data, software or compute are absent. Bind a summary of **all** relevant loci and assert reviewed scope, not one passing locus. This coverage checklist never certifies publication readiness or mechanisms, even with all receipts present.

## Native test boundary

`python -m pytest -q tests/test_native_runtime.py` invokes actual programs if installed. It parses R scripts, executes synthetic SuSiE/coloc, runs actual PLINK clumping/signed LD, and executes synthetic TwoSampleMR. It does not fake their results. Missing programs/packages skip explicitly. The supplied CI native CRAN/PLINK job has **not run** in this session. LDSC/LAVA/PLACO/GenomicSEM/MATLAB still require separate production smoke/reproduction runs on your host.
