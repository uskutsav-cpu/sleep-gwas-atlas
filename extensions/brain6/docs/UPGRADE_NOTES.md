# v0.1 → v0.2: changes and compatibility

## Actual provenance

The current remote repository returned 404. The baseline was the user's previously saved `brain6-extension.patch`, recovered into an isolated work directory. No current `main` files were read or overwritten. All work in this delivery is an extension of that baseline. The upgrade patch is against that recovered v0.1 payload, not against an asserted current remote commit.

## New implementation

`mr_workflow.py`: exposure-only candidates; native clump settings must bind the exact candidate table. MR instrument r²/window defaults are separate from pleiotropic-locus defaults. Outcome P is not a filter. Actual result-family imports require native MR and harmonization receipts, source identities, direction and input settings. Both directions of all pairs remain in the family.

`replication.py`: a new release is not automatically independent. Keep per-side cohort scope, exact/comparable phenotype matching, heritability/QC status, hash-bound evidence and the complete testing denominator. Underpowered pairs stay unpromoted. Estimates in empirical manifests must match a hashed one-row output file.

`qtl.py`, `molecular.py`: build a bounded-memory full-cis SQLite index, conservatively exclude duplicated feature coordinates, align QTL/GWAS effects and sample sizes, create native SuSiE/coloc jobs, and account for all frozen feature/locus queries. Local molecular pairs cannot enter genome-wide PLACO.

`annotations.py`: correct BED0/GWAS1 coordinates, merged regulatory intervals, excluded contig counts, PIP-weighted locus scores. This is a prepared-annotation analysis, not raw single-cell sequencing processing or cell-mechanism proof.

`campaign.py`: explicit stage coverage. One passing locus is not a completed multi-locus stage. Null or data-blocked disorders are preserved. No checklist issues a scientific/publication seal.

## Safety fixes

Strict JSON rejects duplicate keys and non-finite constants. Receipt verification checks nonempty/nonduplicated output inventories and file byte counts/hashes. R and Python native source changes participate in fingerprints. Ambiguous-position deletion uses SQLite rather than loading the entire conflicting-position list into RAM. All-null pair decisions can be frozen without division by zero. LAVA custom LD references can be bound using an explicit file/SHA inventory; old rejected analyses are untouched. Native MR and coloc adapters have stronger input/output checks.

## Applying the update

The outer `apply_brain6.py` defaults to read-only. Fresh installs add only `extensions/brain6/` and the separate CI file. With `--upgrade`, replacements are allowed only if each existing file matches the exact recovered v0.1 digest. Local edits conflict and stop before file installation. `--apply`, `--commit` and `--push` are explicit, use existing identity/credentials and do not force-push. Review current branch/history before pushing: a feature branch also contains its existing parent history.

The standalone delta patch may be inspected with `git apply --check`; it also must not be forced over modified local files. Keep backup copies of your data-host checkout outside the installation process. An interrupted multi-file installation can be rerun; identical installed files are idempotent and modified files remain conflicts.

## Run versioning and native validation

Use a **new v0.2 output root**. Broader code fingerprints intentionally invalidate v0.1 resume compatibility. Do not delete old results to make a new run pass; import compatible legacy artifacts with the explicit import pathway. Nothing in this patch asks you to rerun Pair A or Pair B, restart CONTROL or reuse failed LAVA science.

The four native smoke tests use real binaries only when available: R syntax, synthetic SuSiE/coloc, synthetic PLINK clumping/LD and synthetic TwoSampleMR. They skipped in this environment. Native CRAN/PLINK CI is provided but has not run on GitHub. Its environment resolves dependencies at test time and is not a production lock.

Original inferential decisions need review. The 72-row brain subset retains the original 396-test FDR. Six-disorder selection follows an observed screen; additional pair thresholds do not remove selective-inference issues. Do not replace old thresholds, PLACO versions, data sources or method names silently.
