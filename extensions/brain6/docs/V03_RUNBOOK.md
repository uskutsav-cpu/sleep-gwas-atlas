# Brain6 v0.3 computational execution runbook

This release extends the exact v0.2 payload; it does not replace the main atlas,
root scripts, core panel, A/B/CONTROL contracts, or any measured GWAS result.
The task is computation, not figures or manuscript preparation. The code has
been exercised here on artificial fixtures and available archived numerical
outputs. Full production execution is NOT complete.

## 1. Install on a separate review branch

The outer `apply_brain6.py` defaults to checksum/conflict checks only. `--upgrade`
permits replacement only when existing bytes match the known v0.2 delivery.
Unrecognized local edits, wrong repositories, symlink output targets, unfinished
merges and unrelated staged files stop application/commit. No force pushes,
branch resets, root-file replacements, automatic commits, or main merges occur.

```bash
REPO="/Users/swethasunilkumar/Documents/Codex/2026-08-26/go/work/sleep-gwas-atlas"
DELIVERY="$HOME/Downloads/brain6_v03_delivery"
python3 "$DELIVERY/apply_brain6.py" "$REPO" --upgrade &&
python3 "$DELIVERY/apply_brain6.py" "$REPO" --upgrade --apply
```

After a successful application, use a subshell so a failed command stops later
commands rather than running in another project directory:

```bash
(
  set -e
  cd "$REPO/extensions/brain6"
  .venv/bin/python -m pip install -e '.[test]'
  .venv/bin/python -m pytest -q
  bash scripts/audit_on_host.sh "$REPO"
)
```

The existing extension virtual environment is reused. A missing interpreter
stops the sequence; create a Python >=3.11 environment at that exact location
before retrying. A host audit writes a new timestamped folder in Downloads;
it does not download data, execute production genetics, edit decisions, commit,
or push. `HOST_AUDIT.json`, `steps_1_24.tsv`, and `atlas_review/brain72.tsv`
are the first results to inspect. Run `host-audit24 --hash-dense --scan-shards`
explicitly for full byte hashing and the available 176-shard scan.

## 2. Important correction from measured results

The recovered September 4 archive contains 396 global rows (372 primary and
24 QC-failed sensitivity). Its original 396-family FDR reproduces exactly;
153 primary rows pass (98 positive,55 negative). The 72 brain rows contain 35
significant relationships. The draft selection is:

| Disorder | Draft partner | Rationale |
|---|---|---|
| ADHD | Insomnia | Preserve previously selected Track-B Pair B; short sleep has the smallest FDR, but does not retroactively replace B |
| MDD | Insomnia | Smallest original family FDR among eligible measured pairs |
| Schizophrenia | Long sleep | Smallest original family FDR, not chronotype by default |
| Bipolar disorder | Long sleep | Smallest original family FDR, not chronotype by default |
| Parkinson's | Long sleep | Only significant sleep partner in this archived matrix |
| Alzheimer's | NO_ELIGIBLE_PAIR | None of its 12 global tests pass original-family FDR |

These are proposals, NOT a new frozen scientific choice. Re-evaluate the current
host matrix rather than trusting an older archive. Alzheimer's remains in the
six-disorder reporting universe. A globally null result does not exclude local
sharing. To investigate Alzheimer's anyway, define a separate local/pleiotropic
exploratory family and its multiplicity/power requirements; do not falsely call
insomnia–AD a significant pre-existing global discovery or weaken the core FDR.
The current primary-pair freezer intentionally refuses globally ineligible pairs.
A separately approved null-track selection interface is not implemented here.

After review, use the existing `freeze` command with the actual host matrix,
`atlas_profile.json`, and the reviewed decision table. Do not use a fabricated
reviewer's name. `plan24 --pair-lock ... --out ...` constructs all 23 prerequisites
plus the final release gate. Missing data cannot be bypassed by deleting a pair.

## 3. Production runtime and data

`doctor24 --out ...` reports actual executable hashes, Python/R versions and R
imports. Found/importable is distinct from statistically validated. In this
container R, PLINK, MATLAB and Snakemake were not available, and network/DNS
failure prevented installation. These limitations describe this environment,
not the user's Mac or a lab server.

`bootstrap_native.R NEW_EMPTY_LIBRARY` is supplied for explicit host use. It
installs CRAN dependencies and three commit-pinned upstream packages into a NEW
isolated R library. It refuses a nonempty target rather than updating a running
project's packages. It records package versions and session information. CRAN
versions must subsequently be frozen with the host's environment tooling; this
is NOT an already reproduced cross-platform environment lock. System R/compiler,
PLINK, LDSC, its reference files, MATLAB licensing, and genomic data are separate
requirements. The bootstrap script was not parsed or executed by R here.

Use `validate24` in `native_environment` mode against
`extensions/brain6/tests/test_native_runtime.py` for actual native functional
checks on artificial data. The mode requires that exact shipped test file;
skips/missing imports cannot yield PASS. These checks cover syntax, SuSiE/coloc,
PLINK and TwoSampleMR artificial examples, not full LAVA/GenomicSEM/PLACO
production calibration. Production stages retain their own semantic gates.

`source-inventory` reads the ORIGINAL panel, schema and source registry; it
reports actual bytes, source URLs/checksums, coordinate, ancestry and sample-size
conventions. It does not guess columns, invent INFO, or turn HapMap3 shards into
dense whole-genome data. Supply reviewed source cards before existing
`assemble-preparation`/`pipeline` can run normalization and pair preparation.
The original signed effect and standard error must use compatible scales.

### Source-scale gates

The archive specifies short/long-sleep BOLT binary-LMM coefficients, not ordinary
logistic log odds; bipolar `2*NEFFDIV2` and schizophrenia `2*NEFF` require
source-specific N handling. Coordinate-poor MDD inputs require an audited map.
Alzheimer's direct/proxy cases and GRCh38 source require an explicit definition
and build reconciliation. Do not use cohort headline N when the registered
per-variant sample size is required.

Production SuSiE jobs now require `scale_reviewed`, `effect_scale`,
`n_semantics`, and a sample-size justification. Case-control data must be
log-odds scaled; quantitative data need a sourced sdY. Treating binary-LMM
coefficients as quantitative is an approximation and requires an explicit
justification; it is not enabled silently. MR jobs similarly declare exposure
and outcome units. An odds ratio cannot be requested for a quantitative outcome.
Existing v0.2 job settings must be upgraded, not reused under a different scale.

## 4. Full-family analysis commands

The existing native adapters still run the actual named programs. None of the
new Python collectors substitutes a home-made test for LDSC/LAVA/PLACO/SuSiE,
coloc, GenomicSEM or TwoSampleMR.

| Command | Contract/output |
|---|---|
| `joint-placo-family` | Full selected-pair x variant joint BH with disk-backed sorting; failed variant slots retained; extreme-Z exclusions explicit; exact native/legacy output hashes; mixed PLACO/PLACO+ requires compatibility review |
| `local-family` | Complete planned pair x LD-block denominator, explicit missing/error versus low-local-heritability outcomes |
| `pathways24` | Competitive gene-set overlap at independent LD-block level, exact stratified hypergeometric convolution; not MAGMA/S-LDSC |
| `cell-family24` | Complete expected pair x cell family of matched-locus permutation tests, including unavailable cells; finite permutation resolution reported |
| `compare24` | Common-testable-block overlap and Jaccard; descriptive-only by default because shared sleep GWAS/cohorts induce dependence |
| `robustness24` | Full frozen specification grid; missing runs retained; consistent negative results not reported as positive evidence |
| `provenance24` | Transitive input/output receipt verification, including underlying input-byte drift |
| `validate24` | Executes declared actual pytest suites; saves raw logs and JUnit counts; skips cannot make a complete validation PASS |
| `audit24` | Full 24-step/unit coverage with method identities, input fingerprints, explicit locus/feature universes and no single-locus-for-whole-family shortcut |
| `seal24` | Refuses incomplete/synthetic contamination; exports only explicitly redistribution-reviewed aggregate summaries |

Each collector has a JSON template under `configs/`. `reviewed:false`, empty
families and UNRESOLVED paths are intentionally not runnable defaults.
Specify the entire planned family before data interpretation. All output
folders are immutable: use a new run identifier when inputs, code, or policies
change. Transaction publication now rechecks input and implementation hashes
at the end to catch concurrent modification.

The pathway null fixes selected-locus counts within each matched stratum. It
requires control exchangeability, appropriate locus-size/LD/gene-density
matching, and reviewed variant-to-gene evidence; gene overlap alone does not
prove mechanism. Cross-disorder inferential overlap remains disabled unless a
specific shared-GWAS dependence/exchangeability justification is recorded.
Do not use nominal independent-set Fisher tests on reused GWAS results.

## 5. Coverage and scientific boundaries

A 24-step draft expands global, trait, pair, direction, locus and QTL-feature
requirements. Scope manifests pin the actual tested locus/feature universe.
Empty scopes need verified upstream NO_SIGNAL evidence. Missing, partial,
underpowered, unavailable, discordant and failed analyses are not converted into
successful biological results. An explicit prespecified scope omission is
reported as such, never as an analysis that ran. Accepted artifact classes and
native-method identities are restricted per step; production units additionally
pin expected input bytes.

The existing GenomicSEM model adapter returns INSUFFICIENT_EVIDENCE pending
model-fit review. A source-faithful model-acceptance adapter is still required
before it can satisfy final computational completion. Do not edit receipt flags
to bypass that boundary. Optional HDL-L/SUPERGNOVA, CAUSE/MR-PRESSO/Steiger,
full factor-GWAS, raw single-cell and spatial preprocessing are not newly
implemented/validated by this upgrade. Keep existing root methods where
appropriate rather than claiming this extension contains them.

## 6. Computational finishing sequence

Preserve existing compatible global/Track-B artifacts rather than rerunning for
appearance. Bind and validate native tools and dense inputs; review/freeze the
measured family; prepare shared trait files once; finish source-aware replication;
run native local and genome-wide pleiotropy analyses; apply full-family accounting;
map independent loci; align signed LD; fine-map BOTH traits; run signal-wise and
prior-sensitive coloc; run full-cis QTL fine-mapping/coloc; perform regulatory,
cell/pathway and cross-disorder comparison; exposure-only instrument MR in both
directions; validated GenomicSEM and declared complementary conjFDR; complete
robustness and transitive provenance; execute all declared tests; seal aggregate
outputs only after every required unit resolves under the frozen protocol.

The archive scan in this delivery verifies data already available here. It does
not finish those missing native computations. The final real-data study remains
blocked on the actual host data and native execution.
