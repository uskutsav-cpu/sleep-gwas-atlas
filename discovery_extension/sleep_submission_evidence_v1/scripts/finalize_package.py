#!/usr/bin/env python3
"""Inventory the authorized release; retain scientific failures as explicit gates."""
import argparse
import collections
import csv
import hashlib
import json
import subprocess
from pathlib import Path

P = Path(__file__).resolve().parents[1]
R = P.parents[1]
BASE = '303c905c229104b7bf623df6de7431abb87e0727'


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def rows(path):
    with path.open(newline='') as handle:
        return list(csv.DictReader(handle, delimiter='\t'))


def read(name):
    return json.loads((P / 'logs' / name).read_text())


def git(*args):
    return subprocess.check_output(['git', '-C', str(R), *args], text=True).strip()


def release_files():
    output = subprocess.check_output(['git', '-C', str(R), 'ls-files', '-z', '--cached', '--others', '--exclude-standard', '--', str(P.relative_to(R))])
    return sorted({R / x.decode() for x in output.split(b'\0') if x and (R / x.decode()).is_file()})


def verify():
    errors = []
    manifest = P / 'hashes/PACKAGE_SHA256SUMS'
    for line in manifest.read_text().splitlines():
        expected, rel = line.split('  ', 1)
        path = R / rel
        if not path.is_file() or sha(path) != expected:
            errors.append(rel)
    if errors:
        raise SystemExit('FAIL_RELEASE_HASHES: ' + ', '.join(errors))
    print(json.dumps({'status': 'PASS_RELEASE_HASHES', 'files': len(manifest.read_text().splitlines())}))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--verify', action='store_true')
    parser.add_argument('--ci-json', type=Path)
    args = parser.parse_args()
    if args.verify:
        verify()
        return
    numerical = read('numerical_execution.json')
    independent = read('numerical_independent.json')
    source = read('provenance_execution_receipt.json')
    literature = read('literature_audit_receipt.json')
    workbook = read('workbook_independent_validation.json')
    assert independent['status'] == 'PASS' and workbook['status'] == 'PASS'
    assert read('package_consistency.json')['status'] == 'PASS_PACKAGE_CONSISTENCY'
    ci = json.loads(args.ci_json.read_text()) if args.ci_json else None
    if ci:
        assert ci['conclusion'] == 'success' and ci['headBranch'] == 'publication/sleep-phenome-evidence-v1'
        (P / 'logs/github_evidence_ci.json').write_text(json.dumps(ci, indent=2) + '\n')
    elif (P / 'logs/github_evidence_ci.json').exists():
        ci = read('github_evidence_ci.json')
    gates = []

    def gate(name, status, boundary, evidence):
        assert status in ['PASS', 'FAIL', 'BLOCKED', 'NOT_APPLICABLE']
        gates.append(dict(gate=name, status=status, scope_or_boundary=boundary, evidence=evidence))

    gate('all_applicable_existing_tests', 'FAIL', 'Full checkout:459 tests;12 failures,37 errors,3 skips; missing historical artifacts/R. Five extension scripts pass;two R-dependent scripts blocked.', 'logs/existing_tests_complete_checkout.log;logs/existing_extension_tests_final.json')
    gate('all_new_tests', 'PASS', '25 source-free numerical/evidence tests, including failure modes;not native science', 'logs/final_tests.log')
    gate('independent_numerical_recalculation', 'PASS', '1200 discovery rows,217 validation rows,41 printed-log estimates;archived precision', 'logs/numerical_execution.json;logs/numerical_independent.json')
    gate('frozen_extension_output_hashes', 'PASS', '9/9 original frozen extension outputs exact', 'logs/provenance_frozen_nine_recheck.json')
    gate('complete_original_core_checkpoint', 'BLOCKED', 'One of ten formerly missing artifacts recovered;nine original files still missing;environment order policy fails', '02_SOURCE_AND_PROVENANCE_AUDIT.md')
    gate('original_locked_outputs_unchanged', 'PASS', 'Only new evidence namespace,new workflow,and isolated synthetic test change from authoritative main', 'logs/repository_release_scope.json')
    gate('missing_values_not_zero', 'PASS', 'TSV original fields exact;all117507 workbook cells independently compared;typed blanks explicitly documented', 'logs/package_consistency.json;logs/workbook_independent_validation.json')
    gate('complete_test_universes', 'PASS', 'Exact1200 Cartesian and217 frozen candidate identities retained', 'logs/package_consistency.json')
    gate('novelty_classifications_consistent', 'PASS', '1200 classifications reproduced;zero first-ever novelty promotions;this gate is consistency,not priority', 'logs/final_literature_canonical.log;06_NOVELTY_AND_PRIOR_ART.md')
    gate('exhaustive_current_pair_priority', 'BLOCKED', 'Current primary full tables reviewed;complete per-pair current database searches not executed;521 pairs unresolved', 'logs/literature_audit_receipt.json')
    gate('figure_source_mappings', 'PASS', 'Six source exports;1200 atlas,23 forest,7 heterogeneity,18 domain rows', 'logs/numerical_figure_source_audit.json')
    gate('visual_and_editable_exports', 'PASS', 'Seven figure pages and six one-page numerical Word tables inspected;9-sheetXLSX verified', 'logs/visual_qa.json;logs/workbook_independent_validation.json')
    gate('public_branch_release_scope', 'PASS', 'Derived evidence only;no new raw GWAS,reference bodies,publisher cache,credentials or manuscript staged', '16_PUBLIC_RELEASE_AND_LICENSE_AUDIT.md;logs/repository_release_scope.json')
    gate('license_and_archive_approval', 'BLOCKED', 'Investigator source-specific rights,project-code license and any archive/DOI require human decisions', '16_PUBLIC_RELEASE_AND_LICENSE_AUDIT.md')
    gate('lint_and_compilation', 'PASS', 'Ruff correctness rules F,E9;Python compileall;Node syntax;git diff whitespace', 'logs/lint_final.log;logs/compilation.json')
    clean = read('clean_checkout_receipt.json') if (P / 'logs/clean_checkout_receipt.json').exists() else None
    gate('clean_checkout_reproduction', 'PASS' if clean and clean['status'] == 'PASS' else 'BLOCKED', 'Clean archived checkout,25 tests,independent arithmetic,figure regeneration,canonical evidence and synthetic e2e', 'logs/clean_checkout_receipt.json')
    gate('actual_branch_github_actions', 'PASS' if ci else 'BLOCKED', 'Dedicated evidence workflow at recorded actual commit;legacy global workflow evaluated separately', 'logs/github_evidence_ci.json' if ci else 'Await first authorized push and actual run')
    gate('main_branch_unchanged', 'PASS', BASE, 'logs/repository_release_scope.json')
    gate('no_manuscript_generated_or_modified', 'PASS', 'Six Word files are standalone editable numerical tables;no manuscript sections,significance or cover letter', 'logs/repository_release_scope.json')
    gate('current_remote_source_bodies', 'PASS', 'Two MVP complete-body SHA256/MD5 matches;raw retention0;other111 retain historical receipts', 'logs/provenance_remote_hash_gwas_catalog_GCST90479148.json;logs/provenance_remote_hash_gwas_catalog_GCST90479330.json')
    gate('native_gwas_ldsc_replay', 'BLOCKED', 'Native reruns0;original sleep/munged outcomes and EUR LD/weights missing', '02_SOURCE_AND_PROVENANCE_AUDIT.md')
    gate('fully_independent_two_trait_replication', 'FAIL', 'All23 positives reuse sleep GWAS;qualified outcome-side validation survives', '04_REPLICATION_INDEPENDENCE_AUDIT.md')
    gate('heterogeneity_calibrated_covariance', 'BLOCKED', 'Seven nominal positive-pair flags assume zero covariance;shared sleep estimator covariance unknown', '05_HETEROGENEITY_VALIDATION.md')
    gate('smoking_adiposity_conditional_analysis', 'BLOCKED', 'Valid joint covariance/source inputs unavailable;no synthetic observed estimates', '09_NEGATIVE_AND_BLOCKED_RESULTS.md')
    gate('local_architecture_placo_finemapping_coloc', 'BLOCKED', 'Dense source-verified stats,matched LD,QTL/tooling absent;not automatic prerequisite to global paper', '09_NEGATIVE_AND_BLOCKED_RESULTS.md')
    gate('mr_reporting_checklist', 'NOT_APPLICABLE', 'No MR conducted;literature MR comparisons do not constitute new MR', '12_STROBE_STREGA_CHECKLIST.md')
    gate('investigator_human_approvals', 'BLOCKED', '15 unsupplied author/ethics/rights/review decisions;no invented approvals', '15_ETHICS_AUTHOR_APPROVAL_CHECKLIST.md')
    gate('sleep_editorial_novelty', 'FAIL', 'No adequately supported first-ever result;prior screens materially overlap;NO_GO as currently supported', '17_SLEEP_EDITORIAL_READINESS.md')
    path = P / 'tables/final_acceptance_gates.tsv'
    with path.open('w', newline='') as handle:
        writer = csv.DictWriter(handle, fieldnames=list(gates[0]), delimiter='\t', lineterminator='\n')
        writer.writeheader(); writer.writerows(gates)
    table_inventory = {x.name: len(rows(x)) for x in sorted((P / 'tables').glob('*.tsv'))}
    receipt = dict(schema_version='1.0', execution_date='2026-10-08', authoritative_main_sha=BASE,
                   branch='publication/sleep-phenome-evidence-v1', overall_readiness='NO_GO_FOR_SLEEP_AS_CURRENTLY_SUPPORTED',
                   human_manuscript_handoff='CONDITIONAL', package_completion='FEASIBLE_NON_MANUSCRIPT_WORK_COMPLETE_WITH_EXPLICIT_EXTERNAL_AND_HUMAN_BLOCKERS',
                   numbered_audits=18, sensitivity_analysis_categories=len({r['analysis_id'] for r in rows(P/'tables/sensitivity_results.tsv')}),
                   sensitivity_rows=numerical['sensitivity_rows'], numerical=numerical, independent=independent,
                   source_integrity_rows=source['source_integrity_rows'], source_metadata_rows=source['source_metadata_rows'],
                   recovered_previously_missing_checkpoint_artifacts=source['exact_original_missing_checkpoint_artifacts_recovered'],
                   recovered_project_original_files=3, current_full_raw_bodies_verified=2, raw_bodies_retained=0,
                   native_scientific_reruns=0, source_body_bytes=source['current_full_body_bytes'], literature=literature,
                   figure_types=6, figure_pages=7, editable_word_tables=6, workbook_sheets=9,
                   workbook_cells_verified=workbook['total_checked_table_cells'], new_tests=25,
                   existing_suite=dict(tests=459, failures=12, errors=37, skipped=3),
                   canonical_tsv_tables=table_inventory, gate_status_counts=dict(collections.Counter(x['status'] for x in gates)),
                   branch_ci=ci, clean_checkout=clean,
                   manuscript_written=False, journal_contacted=False, charges_incurred=False,
                   restricted_inputs_released=False, main_modified=False,
                   hash_inventory_excludes=['FINAL_EXECUTION_RECEIPT.json', 'hashes/PACKAGE_SHA256SUMS'])
    (P / 'FINAL_EXECUTION_RECEIPT.json').write_text(json.dumps(receipt, indent=2) + '\n')
    (P / 'README.md').write_text('''# Sleep GWAS research and submission support

Overall SLEEP assessment: **NO_GO_FOR_SLEEP_AS_CURRENTLY_SUPPORTED**. Human manuscript handoff: **CONDITIONAL**. This is a completed non-manuscript evidence package with explicit external-data, scientific and human blockers; no manuscript or submission was created.

Start with [18 Manuscript author handoff](18_MANUSCRIPT_AUTHOR_HANDOFF.md), [17 Readiness](17_SLEEP_EDITORIAL_READINESS.md), [Final blockers](FINAL_BLOCKERS_AND_RESTART_ACTIONS.md), and [Final execution receipt](FINAL_EXECUTION_RECEIPT.json). The 18 numbered audits cover science, prior art, sources, ethics, reporting and journal rules. The main-source baseline is recorded in the receipt. Historical outputs remain unchanged.

All 1,200 discovery comparisons and all 217 validation candidates were checked numerically at archived precision. The 603 BH discoveries remain. All 23 original positives remain Bonferroni positive but are **external-outcome-side replication**, with the same sleep GWAS; none is fully independent two-trait replication. Seven nominal heterogeneity flags assume zero estimator covariance. Nine original core checkpoint files remain missing. Native GWAS/LDSC reruns: zero.

The `tables/` directory contains canonical TSVs, a nine-sheet editable `supplementary_data.xlsx`, and six editable numerical tables under `tables/word/`. The `figures/` directory contains five research figure types and one preliminary graphical abstract, seven pages total, with PDF/vector/300-dpi PNG, source data and alt-text metadata. These are separate support assets, not a manuscript. Final figure selection, print scaling and journal portal supplement rules need human review.

Independent receipts: `logs/numerical_independent.json`, `logs/numerical_figure_source_audit.json`, `logs/provenance_independent_review.json`, `logs/workbook_independent_validation.json`, `logs/visual_qa.json`, `logs/clean_checkout_receipt.json`. All 25 new tests pass. The historical full test suite has 12 failures and 37 errors, retained explicitly; passing source-free tests cannot certify native reproduction.

Use [REPRODUCE](REPRODUCE.md) and `scripts/finalize_package.py --verify` to check the released hash inventory from the repository root. A standalone output copy is a review packet; reproduction uses the dedicated Git branch because canonical historical inputs remain in their original repository paths. New full publisher/supplement/search caches and raw GWAS are excluded from public release; exact acquisition URLs and hashes remain in `sources/`. No new project-code license, public archive or DOI has been asserted.

Machine gates: `tables/final_acceptance_gates.tsv`. PASS concerns the named evidence level only. FAIL, BLOCKED and NOT_APPLICABLE rows are deliberate scientific boundaries and cannot be promoted by software test success.
''')
    # Inventory does not hash itself or the final receipt, avoiding a circular claim.
    excluded = {P/'FINAL_EXECUTION_RECEIPT.json', P/'hashes/PACKAGE_SHA256SUMS', P/'logs/release_inventory.json'}
    files = [x for x in release_files() if x not in excluded]
    inventory = [dict(path=str(x.relative_to(R)), bytes=x.stat().st_size, sha256=sha(x)) for x in files]
    inv_path = P / 'logs/release_inventory.json'
    inv_path.write_text(json.dumps(dict(files=len(inventory), total_bytes=sum(x['bytes'] for x in inventory), artifacts=inventory), indent=2) + '\n')
    files.append(inv_path)
    (P / 'hashes/PACKAGE_SHA256SUMS').write_text(''.join(sha(x)+'  '+str(x.relative_to(R))+'\n' for x in sorted(files)))
    print(json.dumps(dict(status='RELEASE_INVENTORIED_WITH_QUALIFICATIONS', files=len(files), gates=len(gates), tsv_tables=len(table_inventory), gate_status_counts=receipt['gate_status_counts'])))


if __name__ == '__main__':
    main()
