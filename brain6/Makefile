# Reproducible Brain6 entry points. Run from any directory with `make -C brain6 <target>`.
# Expensive production targets require explicit input/output paths and do not form an implicit `all` target.
REPO_ROOT := $(abspath ..)
BRAIN6_PYTHON ?= $(REPO_ROOT)/extensions/brain6/.venv/bin/python
PYTHON ?= $(BRAIN6_PYTHON)
R_RUNTIME_DIR ?= /Volumes/Extreme SSD/brain6-work/r-env-4.3.3-short/bin
RUNTIME_RSCRIPT ?= $(R_RUNTIME_DIR)/Rscript
PREPARATION_ROOT ?= /Volumes/Extreme SSD/brain6-work/preparation-v1
LAVA_INPUT_ROOT ?= $(PREPARATION_ROOT)/lava-inputs-v1
LAVA_OUTPUT_ROOT ?= /Volumes/Extreme SSD/brain6-work/lava-results-v1
LAVA_RETRY_LOCUS ?= 740
LAVA_RETRY_PAIR ?= longsleep__parkinson
LAVA_RETRY_TMPDIR ?= /private/tmp/brain6-lava-roundoff-tmp
LAVA_ROUNDOFF_AUDIT ?= brain6/qc/lava_production_checkpoint_audit_20260924T2139Z.json
MPLCONFIGDIR ?= $(shell $(PYTHON) -c 'import tempfile; print(tempfile.gettempdir())')/brain6-matplotlib

.PHONY: help analysis-deps global refresh-figures sensitivity profile-sensitivity supplement-software supplement-context supplement-sensitivity supplement-dense-qc validate validate-manuscript validate-final-report audit-checkpoints audit-source-archives script-tests extension-tests test check \
		placo-publish placo-candidates placo-candidate-figure placo-directions placo-signed-ld placo-valid-ld placo-valid-signed-effects placo-validate lava-run lava-roundoff-run lava-roundoff-retry lava-details lava-roundoff-details lava-canonical-finalize \
		lava-not-run-causes lava-not-run-summary lava-roundoff-compare lava-pair-context-audit cross-layer-evidence \
		power-screen-aggregate validate-finngen-r9 novelty-crosswalk published-rg-context

help:
	@printf '%s\n' \
	  'Brain6 workflow (run: make -C brain6 TARGET [VARIABLE=value])' \
	  '  analysis-deps     Check NumPy, pandas, SciPy and Matplotlib in the default Brain6 venv' \
	  '  global             Rebuild the locked 12x6 map without a new FDR family' \
	  '  refresh-figures    Re-render figures after validating the locked tables' \
	  '  profile-sensitivity Rebuild the frozen leave-one-sleep-trait-out profile table' \
	  '  sensitivity        Rebuild the interim source-hashed claim sensitivity ledger' \
	  '  supplement-software Build Supplementary Table S17 from the software manifest' \
	  '  supplement-context Build Supplementary Table S18 with overlap-limited published rg context' \
	  '  supplement-sensitivity Build interim Supplementary Table S16 from the claim sensitivity ledger' \
	  '  supplement-dense-qc Rebuild dense QC summaries and provenance-bound Supplementary Table S1' \
	  '  lava-not-run-causes Build the checksum-bound canonical v3 NOT_RUN cause table' \
	  '  lava-not-run-summary Build trait and locus summaries of canonical v3 NOT_RUN causes' \
	  '  lava-pair-context-audit Compare prior pair-context univariates with canonical trait-only LAVA v3' \
	  '  cross-layer-evidence Build the descriptive five-pair evidence table' \
	  '  power-screen-aggregate Combine and validate the completed Brain6-only sleep-duration screen batches' \
	  '  validate           Verify current Brain6 outputs and evidence hashes' \
	  '  validate-manuscript Check key global-map prose claims against machine-readable tables' \
	  '  validate-final-report Check the 35-row provisional final report and source hashes' \
	  '  audit-checkpoints  Read-only verification of available PLACO/LAVA production receipts' \
	  '  audit-source-archives Compare cached raw GWAS files to frozen source sizes and SHA256 values' \
	  '  script-tests       Run Brain6 orchestration/output tests' \
	  '  extension-tests    Run extension suite with pinned R 4.3.3 first on PATH' \
	  '  test               Run both Brain6 test suites' \
	  '  check              Rebuild sensitivity, run tests, validate outputs and receipts' \
	  '  placo-publish      Require PAIR=<frozen pair id>; validate and publish a completed pair' \
	  '  placo-candidates   Regroup currently published PLACO candidates with frozen LD rules' \
	  '  placo-candidate-figure Render/verify the explicitly partial four-pair PLACO candidate overview' \
	  '  placo-directions    Build the checksum-bound partial candidate direction table' \
	  '  placo-signed-ld     Export signed LAVA-reference LD for partial candidate-to-lead pairs' \
	  '  placo-valid-ld      Rebuild partial clumps excluding invalid raw LD edges' \
	  '  placo-valid-signed-effects Align partial PLACO effect signs to valid signed LD' \
	  '  lava-run           Run/resume frozen LAVA family; override LAVA_INPUT_ROOT/LAVA_OUTPUT_ROOT as needed' \
	  '  lava-roundoff-run  Run/resume the separately locked LAVA 0.1.5 roundoff-stability family' \
	  '  lava-roundoff-retry Retry one verified environmental failure after its coordinator stops; set LAVA_RETRY_ID' \
	  '  lava-roundoff-compare Compare complete, receipt-verified baseline and roundoff pair-locus outputs descriptively' \
	  '  lava-roundoff-details Materialize details from a completed corrected run; set LAVA_ROUNDOFF_RUN_ID' \
	  '  lava-details       After a complete LAVA run, restore local-h2/CI detail outputs from verified receipts' \
	  '  lava-canonical-finalize Create the immutable v3 family decision after all 17,465 cells are verified' \
	  '  lava-not-run-causes Build a source-bound table for all frozen canonical v3 NOT_RUN cells' \
	  '  cross-layer-evidence Build the source-bound five-pair Brain6 evidence integration' \
	  '  novelty-crosswalk  Bind the 35 inherited Brain6 hits to the dated parent-atlas literature audit' \
	  '  published-rg-context Compare the 21 published sleep–psychiatric rg rows with locked Brain6 estimates descriptively'

analysis-deps:
	@$(PYTHON) -c 'import importlib.util,sys; required=("numpy","pandas","scipy","matplotlib"); missing=[p for p in required if importlib.util.find_spec(p) is None]; print("Missing analysis dependencies: "+", ".join(missing) if missing else "Brain6 analysis dependencies are available"); sys.exit(bool(missing))' || { printf '%s\n' "Install the analysis extra with: $(BRAIN6_PYTHON) -m pip install -e '$(REPO_ROOT)/extensions/brain6[analysis]'"; exit 2; }

global: analysis-deps
	mkdir -p "$(MPLCONFIGDIR)"
	cd $(REPO_ROOT) && MPLCONFIGDIR="$(MPLCONFIGDIR)" $(PYTHON) brain6/scripts/build_locked_global.py

refresh-figures: analysis-deps
	mkdir -p "$(MPLCONFIGDIR)"
	cd $(REPO_ROOT) && MPLCONFIGDIR="$(MPLCONFIGDIR)" $(PYTHON) brain6/scripts/build_locked_global.py --refresh-figures

profile-sensitivity:
	cd $(REPO_ROOT) && $(PYTHON) brain6/scripts/build_profile_leave_one_out.py

sensitivity: supplement-context profile-sensitivity
	cd $(REPO_ROOT) && $(PYTHON) brain6/scripts/analyze_placo_candidate_mhc_exclusion.py
	cd $(REPO_ROOT) && $(PYTHON) brain6/scripts/build_headline_claim_sensitivity.py
	cd $(REPO_ROOT) && $(PYTHON) brain6/scripts/build_supplement_sensitivity.py

supplement-software:
	cd $(REPO_ROOT) && $(PYTHON) brain6/scripts/build_software_supplement.py

supplement-context:
	cd $(REPO_ROOT) && $(PYTHON) brain6/scripts/build_published_context_table.py

supplement-sensitivity:
	cd $(REPO_ROOT) && $(PYTHON) brain6/scripts/build_supplement_sensitivity.py

lava-not-run-causes:
	cd $(REPO_ROOT) && $(PYTHON) brain6/scripts/build_lava_canonical_v3_not_run_causes.py

lava-not-run-summary: lava-not-run-causes
	cd $(REPO_ROOT) && $(PYTHON) brain6/scripts/build_lava_canonical_v3_cause_summaries.py

cross-layer-evidence:
	cd $(REPO_ROOT) && $(PYTHON) brain6/scripts/build_cross_layer_evidence.py

novelty-crosswalk:
	cd $(REPO_ROOT) && $(PYTHON) brain6/scripts/build_brain6_novelty_crosswalk.py

published-rg-context:
	cd $(REPO_ROOT) && $(PYTHON) brain6/scripts/build_brain6_published_rg_context.py

power-screen-aggregate:
	cd $(REPO_ROOT) && $(PYTHON) brain6/scripts/aggregate_power_optimized_sleep_screen_v1.py --write

validate: lava-not-run-causes lava-not-run-summary cross-layer-evidence placo-candidate-figure placo-directions placo-signed-ld placo-valid-ld placo-valid-signed-effects validate-finngen-r9 validate-pgc-mdd2025-rg
	cd $(REPO_ROOT) && $(PYTHON) brain6/scripts/validate_current_outputs.py
	cd $(REPO_ROOT) && $(PYTHON) brain6/scripts/build_brain6_novelty_crosswalk.py --validate
	$(MAKE) validate-manuscript
	$(MAKE) validate-final-report

validate-manuscript:
	cd $(REPO_ROOT) && $(PYTHON) brain6/scripts/validate_manuscript_claims.py

validate-pgc-mdd2025-rg:
	cd $(REPO_ROOT) && $(PYTHON) brain6/scripts/validate_pgc_mdd2025_rg_sensitivity.py

validate-final-report:
	cd $(REPO_ROOT) && $(PYTHON) brain6/scripts/validate_final_brain6_report.py

validate-finngen-r9:
	cd $(REPO_ROOT) && $(PYTHON) brain6/scripts/validate_finngen_r9_replication_candidates.py

audit-checkpoints:
	cd $(REPO_ROOT) && $(PYTHON) brain6/scripts/audit_production_checkpoints.py

audit-source-archives:
	cd $(REPO_ROOT) && $(PYTHON) brain6/scripts/audit_external_source_archives.py

supplement-dense-qc: audit-source-archives
	cd $(REPO_ROOT) && $(PYTHON) brain6/scripts/build_dense_qc_report.py --preparation-root "$(PREPARATION_ROOT)"

script-tests:
	mkdir -p "$(MPLCONFIGDIR)"
	cd $(REPO_ROOT) && MPLCONFIGDIR="$(MPLCONFIGDIR)" $(PYTHON) -m pytest -q brain6/scripts/tests

extension-tests:
	cd $(REPO_ROOT)/extensions/brain6 && PATH="$(R_RUNTIME_DIR):$$PATH" $(BRAIN6_PYTHON) -m pytest -q -p no:cacheprovider tests

test: script-tests extension-tests

check:
	$(MAKE) sensitivity
	$(MAKE) supplement-software
	$(MAKE) test
	$(MAKE) validate
	$(MAKE) audit-checkpoints

placo-publish:
	@test -n "$(PAIR)" || { echo 'Set PAIR to a frozen PLACO pair id'; exit 2; }
	cd $(REPO_ROOT) && $(PYTHON) brain6/scripts/publish_placo_results.py --pair "$(PAIR)"

placo-candidates:
	cd $(REPO_ROOT) && $(PYTHON) brain6/scripts/clump_placo_candidates.py --rscript "$(RUNTIME_RSCRIPT)"

placo-candidate-figure:
	cd $(REPO_ROOT) && MPLCONFIGDIR="$(MPLCONFIGDIR)" $(PYTHON) brain6/scripts/plot_partial_placo_candidate_loci.py

placo-directions:
	cd $(REPO_ROOT) && $(PYTHON) brain6/scripts/build_partial_placo_direction_table.py

placo-signed-ld:
	cd $(REPO_ROOT) && $(PYTHON) brain6/scripts/build_partial_placo_signed_ld.py

placo-valid-ld:
	cd $(REPO_ROOT) && $(PYTHON) brain6/scripts/rebuild_partial_placo_clumps_valid_ld.py

placo-valid-signed-effects:
	cd $(REPO_ROOT) && $(PYTHON) brain6/scripts/build_partial_placo_signed_effect_ld.py

placo-validate:
	cd $(REPO_ROOT) && $(PYTHON) brain6/scripts/audit_placo_family_outputs.py

lava-run:
	cd $(REPO_ROOT) && $(PYTHON) brain6/scripts/run_lava_family.py \
	  --input-root "$(LAVA_INPUT_ROOT)" --output-root "$(LAVA_OUTPUT_ROOT)" \
	  --rscript "$(RUNTIME_RSCRIPT)"

lava-roundoff-run:
	cd $(REPO_ROOT) && $(PYTHON) brain6/scripts/run_lava_family_roundoff_v1.py \
	  --input-root "$(LAVA_INPUT_ROOT)" --output-root "$(LAVA_OUTPUT_ROOT)" \
	  --rscript "$(RUNTIME_RSCRIPT)"

lava-roundoff-retry:
	@test -n "$(LAVA_RETRY_ID)" || { printf '%s\n' 'Set LAVA_RETRY_ID to a unique immutable retry label'; exit 2; }
	cd $(REPO_ROOT) && $(PYTHON) brain6/scripts/retry_lava_roundoff_environmental_locus.py \
	  --retry-id "$(LAVA_RETRY_ID)" --locus-id "$(LAVA_RETRY_LOCUS)" \
	  --pair-id "$(LAVA_RETRY_PAIR)" --temp-dir "$(LAVA_RETRY_TMPDIR)"

lava-roundoff-compare:
	cd $(REPO_ROOT) && $(PYTHON) brain6/scripts/compare_lava_baseline_roundoff.py \
	  --audit "$(LAVA_ROUNDOFF_AUDIT)" \
	  --output brain6/results/lava/lava_baseline_roundoff_pair_locus_comparison_v1.tsv \
	  --summary-output brain6/results/supplement/table_S24_lava_roundoff_comparison.tsv \
	  --provenance brain6/results/lava/lava_baseline_roundoff_pair_locus_comparison_v1.provenance.json

lava-pair-context-audit:
	cd $(REPO_ROOT) && $(PYTHON) brain6/scripts/audit_lava_pair_context_univariates.py \
	  --audit "$(LAVA_ROUNDOFF_AUDIT)" \
	  --canonical work/lava-canonical-v3-production/d730debf45266d298401564f3260bdecb14739d1c1f1835a5aebd615c83fa60b/results/canonical_family_results.tsv \
	  --output brain6/results/lava/lava_pair_context_univariates_v1.tsv \
	  --summary-output brain6/results/supplement/table_S27_lava_pair_context_univariates.tsv \
	  --provenance brain6/results/lava/lava_pair_context_univariates_v1.provenance.json

lava-roundoff-details:
	@test -n "$(LAVA_ROUNDOFF_RUN_ID)" || { printf '%s\n' 'Set LAVA_ROUNDOFF_RUN_ID from brain6/manifests/lava_roundoff_v1_run.json'; exit 2; }
	cd $(REPO_ROOT) && $(PYTHON) brain6/scripts/materialize_lava_family_details.py \
	  --input-root "$(LAVA_INPUT_ROOT)" --output-root "$(LAVA_OUTPUT_ROOT)" \
	  --reference "$(REPO_ROOT)/ref/lava/ukb_v1.1" --rscript "$(RUNTIME_RSCRIPT)" \
	  --run-id "$(LAVA_ROUNDOFF_RUN_ID)"

lava-details:
	cd $(REPO_ROOT) && $(PYTHON) brain6/scripts/materialize_lava_family_details.py \
	  --input-root "$(LAVA_INPUT_ROOT)" --output-root "$(LAVA_OUTPUT_ROOT)" \
	  --reference "$(REPO_ROOT)/ref/lava/ukb_v1.1" --rscript "$(RUNTIME_RSCRIPT)"

lava-canonical-finalize:
	@test -n "$(LAVA_CANONICAL_RUN_DIR)" || { printf '%s\n' 'Set LAVA_CANONICAL_RUN_DIR to the completed canonical v3 run directory'; exit 2; }
	@test -n "$(LAVA_REFERENCE_PROVENANCE)" || { printf '%s\n' 'Set LAVA_REFERENCE_PROVENANCE to the sealed official reference.provenance.json'; exit 2; }
	cd $(REPO_ROOT) && $(PYTHON) brain6/scripts/finalize_lava_canonical_v3.py \
	  --run-dir "$(LAVA_CANONICAL_RUN_DIR)" --input-root "$(LAVA_INPUT_ROOT)" \
	  --reference-provenance "$(LAVA_REFERENCE_PROVENANCE)"
