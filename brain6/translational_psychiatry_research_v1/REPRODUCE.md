# Reproduction and evidence boundaries

Use a checkout of `brain6/translational-psychiatry-research-v1` from the project repository. The handoff directory contains review artifacts; executable scripts resolve historical inputs from the repository root. It is not a replacement for authorized access to the original GWAS sources. Never put new raw ADHD source files, derived beta/SE/Z vectors, native RDS or filtered molecular association slices in a public branch or ZIP.

Python 3.13.13; NumPy 2.4.3; SciPy 1.18.1; pandas 3.0.5; pytest 9.1.1; matplotlib 3.11.2. R 4.6.1, isolated susieR 0.14.2 and coloc 5.2.3. Full package locks, exact source tarball hashes and R session records are in `research_v1/ld/`. Thirty source tarballs total 34,293,244 bytes; the isolated library is 64 MiB. No global R library was altered. TeX uses the built-in editor and an existing Tectonic export runtime.

## Source-free checks and retained-table replay

Run from the repository root. These checks require no new external acquisition and do not certify biological inference. The molecular replay uses the already inherited committed association input slices; no new restricted input is copied into the handoff.

```sh
python -m pytest -q -c brain6/translational_psychiatry_research_v1/research_v1/statistics/pytest.ini -m 'not integration' brain6/translational_psychiatry_research_v1/tests brain6/translational_psychiatry_research_v1/research_v1/statistics/test_statistics.py brain6/translational_psychiatry_research_v1/research_v1/statistics/test_parent_review.py brain6/translational_psychiatry_research_v1/research_v1/ld/test_ld_research.py
python brain6/translational_psychiatry_research_v1/scripts/verify_historical_integrity.py
python brain6/translational_psychiatry_research_v1/scripts/global_contrasts.py run
python brain6/translational_psychiatry_research_v1/scripts/molecular_sensitivity.py run
python brain6/translational_psychiatry_research_v1/scripts/validate_research_package.py
```

Expected: 31 source-free tests pass, two local integration tests explicitly deselected; 5,666 historical files have zero drift; 72 contrasts/648 grid rows and 60 molecular posterior rows replay; 23 artifact/manuscript validation checks pass. Exact hosted results and commit are in the GitHub receipt. Test results are separate from scientific status.

## Actual data-dependent integration and native runs

The separate local integration invocation requires the actual local input objects and fails if absent. It passed two tests on the present authorized local data; a public clean checkout cannot silently certify these inputs.

```sh
python -m pytest -q -c brain6/translational_psychiatry_research_v1/research_v1/statistics/pytest.ini -m integration brain6/translational_psychiatry_research_v1/research_v1/statistics/test_parent_review.py
```

`research_v1/ld/REPRODUCE.md` contains exact acquisition, factor reconstruction, source harmonization, native RSS, eight SuSiE fits, 12 coloc-SuSiE conditions, three ABF conditions and synthetic API commands. Its native artifact/source manifests bind every input and fitted object. The three small public genotype factors rebuild full matrices sequentially. No reference is clipped or PSD-projected. Native synthetic API tests (three) remain distinct from eight real GWAS fits and their independently recomputed RDS/PIP/CS/ELBO/posteriors.

`research_v1/novelty/REPRODUCE_SOURCE_RESEARCH.md` describes primary-source acquisition, all five paginated literature queries, exact supplement extraction and the four-slot FinnGen protocol. `research_v1/statistics/README.md` describes the pinned-source calibration, 24 × 10,000 simulations, 12 unchanged-source fixtures, historical numeric audit and independent review. The Gaussian variance cross-check script has its own 2 × 100,000-draw settings/provenance. These simulations are truth-known software/method diagnostics, not estimated empirical false-positive rates.

## Publication material generation

```sh
python brain6/translational_psychiatry_research_v1/scripts/build_publication_materials.py
python brain6/translational_psychiatry_research_v1/scripts/build_handoff_ledgers.py
tectonic --keep-logs brain6/translational_psychiatry_research_v1/MANUSCRIPT/brain6_reassessment.tex
```

Use the built-in LaTeX editor for source edits and compilation diagnostics. Vector/PNG outputs use deterministic metadata and original saved numerical results. Figures do not run a scientific model. Rebuilding ledgers binds existing results; it does not prospectively freeze or select hypotheses. Regenerate the package checksums only after authorized new work, preserving original protocol freezes and amendments.

The original source/overlap contracts, full canonical worker path, participant independence, human review and author/ethics/permissions approvals remain outside this reproduction claim. Historical packages and unrelated projects stay immutable.
