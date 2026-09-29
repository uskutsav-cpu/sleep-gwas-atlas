# Run from the repository root with `make -C frailty_paper <target>`.
# These targets reuse current frozen inputs; none starts a network download.
PYTHON ?= .venv/bin/python
ROOT := ..

.PHONY: help setup-env audit-frozen manifest verify-manifest audit-manifest-metadata audit-hlma-portal-crosswalk storage-configure storage-check validate-plan validate-overlap scan-gwas-sources plot-frozen-fi plot-frozen-fi-preview physical-component-qc fi-correction-sensitivity latent-correction-sensitivity phase16-sensitivity-summary supplementary-tables checklists import-review-exports validate-review audit-bibliography audit-missing-abstract-sources audit-manuscript-claims audit-latent-qsnp audit-reporting-locators audit-hfrs-published-leads test preflight smoke

help:
	@printf '%s\n' \
	  'Frailty package targets (run from repository root):' \
	  '  Set PYTHON=path/to/python to select a validated interpreter; default: .venv/bin/python' \
	  '  make -C frailty_paper setup-env FRAILTY_PYTHON=... [FRAILTY_VENV_DIR=...] Set up pinned Python 3.11.11 without overlaying a different existing venv' \
	  '  make -C frailty_paper audit-frozen  Verify and extract frozen atlas and extension results' \
	  '  make -C frailty_paper manifest      Rebuild the checksum inventory of verified acquired files' \
	  '  make -C frailty_paper verify-manifest Verify every listed file size and SHA-256 without writing' \
	  '  make -C frailty_paper audit-manifest-metadata Require complete 22-column acquisition metadata' \
	  '  make -C frailty_paper audit-hlma-portal-crosswalk Rebuild link-only HLMA portal/OMIX metadata crosswalks from the saved listing snapshot' \
	  '  make -C frailty_paper storage-configure FRAILTY_EXTERNAL_STORAGE_ROOT=... FRAILTY_ANALYSIS_WORKSPACE=... Link empty derived outputs to an external volume' \
	  '  make -C frailty_paper storage-check FRAILTY_ANALYSIS_WORKSPACE=... Check the actual derived-output volume has at least 20 GiB free' \
	  '  make -C frailty_paper validate-plan Check frozen plan YAML, semantics, and input hashes' \
	  '  make -C frailty_paper validate-overlap Validate source-level cohort ledger against frozen panel' \
	  '  make -C frailty_paper scan-gwas-sources Stream-check acquired GWAS gzip files, headers, and basic numeric domains without filtering' \
	  '  make -C frailty_paper plot-frozen-fi   Plot frozen FI rg extract after checking pinned Python/Matplotlib versions' \
	  '  make -C frailty_paper plot-frozen-fi-preview Allow version mismatch for a provenance-labeled preview only' \
	  '  make -C frailty_paper physical-component-qc Summarize structural QC while keeping provenance-blocked components ineligible' \
	  '  make -C frailty_paper fi-correction-sensitivity Compare inherited all-396 BH with same-family Bonferroni' \
	  '  make -C frailty_paper latent-correction-sensitivity Compare secondary all-84 BH with same-family Bonferroni' \
	  '  make -C frailty_paper phase16-sensitivity-summary Build the headline conclusion-by-sensitivity status matrix' \
	  '  make -C frailty_paper supplementary-tables Rebuild the current source-linked supplementary package' \
	  '  make -C frailty_paper checklists     Regenerate item-level reporting checklists' \
	  '  make -C frailty_paper import-review-exports Normalize untouched licensed RIS/CSV exports' \
	  '  make -C frailty_paper validate-review Validate dual-screening decisions and exclusions' \
	  '  make -C frailty_paper audit-bibliography Check Pandoc citation keys against references.bib' \
	  '  make -C frailty_paper audit-missing-abstract-sources Confirm empty abstracts are absent from all retained PubMed XML occurrences' \
	  '  make -C frailty_paper audit-manuscript-claims Check selected quantitative manuscript claims against source tables' \
	  '  make -C frailty_paper audit-latent-qsnp Audit source-recommended Q_SNP pruning for latent factors' \
	  '  make -C frailty_paper audit-reporting-locators Validate checklist line references against source files' \
	  '  make -C frailty_paper audit-hfrs-published-leads Reproduce the publisher-table UKB p-value counts for reported HFRS lead variants' \
	  '  make -C frailty_paper test           Run package unit/integration checks' \
	  '  make -C frailty_paper preflight FRAILTY_EXTERNAL_STORAGE_ROOT=... Verify resources, frozen contracts, review ledger, and package checks' \
	  '  make -C frailty_paper smoke          Run read-only frozen audits and plan validation'

setup-env:
	FRAILTY_PYTHON="$(FRAILTY_PYTHON)" FRAILTY_VENV_DIR="$(FRAILTY_VENV_DIR)" bash scripts/00_setup.sh $(ROOT)

audit-frozen:
	$(PYTHON) scripts/11_audit_frozen_atlas_frailty.py
	$(PYTHON) scripts/12_audit_frozen_extension.py --repo $(ROOT)

manifest:
	$(PYTHON) scripts/10_build_acquired_resource_manifest.py --repo $(ROOT)

verify-manifest:
	$(PYTHON) scripts/16_verify_acquired_resource_manifest.py --repo $(ROOT) $(if $(FRAILTY_EXTERNAL_STORAGE_ROOT),--external-storage-root "$(FRAILTY_EXTERNAL_STORAGE_ROOT)",)

audit-manifest-metadata:
	$(PYTHON) scripts/45_audit_acquisition_manifest_metadata.py --manifest manifests/all_acquired_resources.tsv

audit-hlma-portal-crosswalk:
	$(PYTHON) scripts/46_build_hlma_portal_archive_crosswalk.py --repo $(ROOT)

storage-configure:
	@test -n "$(FRAILTY_EXTERNAL_STORAGE_ROOT)" -a -n "$(FRAILTY_ANALYSIS_WORKSPACE)" || { echo 'set FRAILTY_EXTERNAL_STORAGE_ROOT and FRAILTY_ANALYSIS_WORKSPACE'; exit 2; }
	mkdir -p "$(FRAILTY_ANALYSIS_WORKSPACE)"
	$(PYTHON) scripts/21_configure_external_storage.py --repo $(ROOT) --storage-root "$(FRAILTY_EXTERNAL_STORAGE_ROOT)" --derived-output-root "$(FRAILTY_ANALYSIS_WORKSPACE)" --apply

storage-check:
	@test -n "$(FRAILTY_ANALYSIS_WORKSPACE)" || { echo 'set FRAILTY_ANALYSIS_WORKSPACE'; exit 2; }
	bash scripts/09_require_storage.sh $(ROOT) "$(FRAILTY_ANALYSIS_WORKSPACE)"

validate-plan:
	$(PYTHON) scripts/13_validate_analysis_plan.py --repo $(ROOT)

validate-overlap:
	$(PYTHON) scripts/20_validate_cohort_overlap.py --repo $(ROOT)

scan-gwas-sources:
	$(PYTHON) scripts/19_scan_acquired_gwas.py --repo $(ROOT)

plot-frozen-fi:
	$(PYTHON) scripts/22_plot_frozen_fi_rg.py

plot-frozen-fi-preview:
	$(PYTHON) scripts/22_plot_frozen_fi_rg.py --allow-unpinned

physical-component-qc:
	$(PYTHON) scripts/34_build_physical_component_qc_summary.py --repo $(ROOT)

fi-correction-sensitivity:
	$(PYTHON) scripts/35_build_fi_correction_sensitivity.py --repo $(ROOT)

latent-correction-sensitivity:
	$(PYTHON) scripts/40_build_latent_correction_sensitivity.py --repo $(ROOT)

phase16-sensitivity-summary:
	$(PYTHON) scripts/41_build_phase16_sensitivity_summary.py --repo $(ROOT)

supplementary-tables: latent-correction-sensitivity phase16-sensitivity-summary
	$(PYTHON) scripts/39_build_supplementary_tables.py --repo $(ROOT)

checklists:
	$(PYTHON) scripts/14_build_reporting_checklists.py

import-review-exports:
	$(PYTHON) scripts/15_import_manual_review_exports.py --repo $(ROOT)

validate-review:
	$(PYTHON) scripts/17_validate_review_screening.py --review-dir review

audit-bibliography:
	$(PYTHON) scripts/42_audit_bibliography.py --repo $(ROOT) --output frailty_paper/analysis/manuscript_bibliography_audit_2026-09-23.json

audit-missing-abstract-sources:
	$(PYTHON) scripts/47_audit_missing_abstract_sources.py --repo $(ROOT) $(if $(FRAILTY_EXTERNAL_STORAGE_ROOT),--external-storage-root "$(FRAILTY_EXTERNAL_STORAGE_ROOT)",)

audit-manuscript-claims:
	$(PYTHON) scripts/48_audit_manuscript_quantitative_claims.py --repo $(ROOT)

audit-latent-qsnp:
	$(PYTHON) scripts/43_audit_latent_qsnp_pruning.py --repo $(ROOT) --output frailty_paper/analysis/latent_factor_qsnp_pruning_audit.tsv --provenance frailty_paper/analysis/latent_factor_qsnp_pruning_audit.provenance.json

audit-reporting-locators:
	$(PYTHON) scripts/44_audit_reporting_locators.py --repo $(ROOT)

audit-hfrs-published-leads:
	$(PYTHON) scripts/49_audit_hfrs_published_lead_replication.py --repo $(ROOT)

test:
	$(PYTHON) -m unittest discover -s tests -v

preflight: verify-manifest audit-manifest-metadata validate-plan validate-overlap validate-review audit-bibliography audit-missing-abstract-sources audit-manuscript-claims audit-latent-qsnp audit-reporting-locators test

smoke: audit-frozen validate-plan
