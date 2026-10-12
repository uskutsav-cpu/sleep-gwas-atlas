#!/usr/bin/env python3
"""Publish aggregate real feasibility outcomes without gold or model answers.

Post-execution reporting only. Does not modify the frozen instrument/scorer.
Requires every registered cell to have an observed, preserved definitive record.
"""
import argparse
import csv
import hashlib
import json
from pathlib import Path


def write_tsv(path, columns, rows):
    with path.open('w', newline='') as stream:
        writer = csv.DictWriter(stream, fieldnames=columns, delimiter='\t', lineterminator='\n')
        writer.writeheader()
        writer.writerows(rows)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--scores', type=Path, required=True)
    parser.add_argument('--registration', type=Path, required=True)
    parser.add_argument('--exclusion-ledger', type=Path, required=True)
    parser.add_argument('--out', type=Path, required=True)
    args = parser.parse_args()
    report = json.loads((args.scores / 'score_report.json').read_text())
    cells = json.loads((args.scores / 'scored_cells.json').read_text())
    registration = json.loads(args.registration.read_text())
    ledger = json.loads(args.exclusion_ledger.read_text())
    if len(registration['selected_task_ids']) != 10 or report['planned_cells'] != 90:
        raise ValueError('This reporter is specific to the registered ten-task feasibility.')
    if ledger['original_definitive_cells'] != 90 or ledger['eligible_cells'] != report['observed_cell_records']:
        raise ValueError('Complete the registered definitive inventory and reconcile exclusions.')
    if ledger['excluded_cells'] != report['unobserved_planned_cells']:
        raise ValueError('Scoring missingness differs from the preserved exclusion ledger.')
    if report['dependence_sensitivity'] is not None:
        raise ValueError('Ten-item feasibility must not have a confidence interval.')
    args.out.mkdir(parents=True, exist_ok=True)
    (args.out / 'aggregate_results.json').write_text(json.dumps(report, indent=2) + '\n')
    (args.out / 'provenance_exclusion_ledger.json').write_text(json.dumps(ledger, indent=2) + '\n')
    rows = []
    for summary in report['summaries']:
        relevant = [r for r in cells if (r['model_alias'], r['condition']) ==
                    (summary['model_alias'], summary['condition'])]
        rows.append(dict(
            model_alias=summary['model_alias'], provider='OpenAI',
            condition=summary['condition'], planned=summary['expected_cells'],
            completed=len(relevant), strict_correct=sum(r['strict_correct'] for r in relevant),
            schema_valid=sum(r['schema_valid'] for r in relevant),
            status_correct=sum(r['status_correct'] for r in relevant),
            mean_field_fraction=summary['mean_field_fraction'],
            numeric_correct=sum(r['numeric_correct'] for r in relevant),
            numeric_applicable=sum(r['numeric_fields'] for r in relevant),
            unsafe_tasks=sum(bool(r['unsafe_indicators']) for r in relevant)))
    columns = ['model_alias', 'provider', 'condition', 'planned', 'completed',
               'strict_correct', 'schema_valid', 'status_correct', 'mean_field_fraction',
               'numeric_correct', 'numeric_applicable', 'unsafe_tasks']
    write_tsv(args.out / 'model_condition_outcomes.tsv', columns, rows)
    paired = [{k: x[k] for k in ['model_alias', 'paired_tasks',
               'strict_task_rate_difference_C_minus_B', 'field_fraction_difference_C_minus_B']}
              for x in report['paired_B_C']]
    write_tsv(args.out / 'paired_B_C_descriptive.tsv',
              ['model_alias', 'paired_tasks', 'strict_task_rate_difference_C_minus_B',
               'field_fraction_difference_C_minus_B'], paired)
    cell_rows = [{
        'model_alias': x['model_alias'], 'condition': x['condition'], 'task_id': x['task_id'],
        'stratum': x['stratum'], 'strict_correct': int(x['strict_correct']),
        'schema_valid': int(x['schema_valid']), 'status_correct': int(x['status_correct']),
        'field_fraction': x['field_fraction'], 'numeric_correct': x['numeric_correct'],
        'numeric_applicable': x['numeric_fields'],
        'unsafe_indicators': '|'.join(x['unsafe_indicators']),
        'response_digest_as_scored': x['response_sha256']}
        for x in cells]
    write_tsv(args.out / 'cell_outcomes_no_answers.tsv',
              list(cell_rows[0]) if cell_rows else ['model_alias', 'condition', 'task_id'], cell_rows)
    lines = [
        'SLEEP GWAS ATLAS — ACTUAL CODEX FEASIBILITY RESULTS',
        'Registered ten-task subset, 12 October 2026; automated scores provisional.',
        '',
        f"Preserved {ledger['original_definitive_cells']}/{report['planned_cells']} actual completed cell envelopes; "
        f"{report['completed_cells']} provenance-eligible responses scored; "
        f"{ledger['excluded_cells']} excluded from primary scoring for returned packet-hash mismatch.",
        'Excluded original envelopes/answers remain unchanged; no fabricated hash correction or retry.',
        'The frozen scorer calls these two absent from its eligible projection; they were actually executed.',
        'Three requested model aliases, all OpenAI: gpt-6.1-sol, gpt-6-sol, gpt-6-luna.',
        'Backend version, effective temperature/seed, latency, tokens and cost are not authenticated.',
        'Primary outcome is strict all-required-fields agreement with the frozen structured key.',
        'This is finite exact retrieval/admissibility correctness, not validated broad interpretation accuracy.',
        '', 'Model / condition / strict correct / schema valid / unsafe requested-field tasks']
    lines += [f"{r['model_alias']} / {r['condition']} / {r['strict_correct']}/{r['completed']} / "
              f"{r['schema_valid']}/{r['completed']} / {r['unsafe_tasks']}/{r['completed']}" for r in rows]
    lines += ['', 'Matched-evidence paired comparisons, C minus B (descriptive):']
    lines += [f"{r['model_alias']}: {r['strict_task_rate_difference_C_minus_B']:+.3f} strict rate; "
              f"{r['field_fraction_difference_C_minus_B']:+.3f} mean field fraction; "
              f"{r['paired_tasks']} paired tasks." for r in paired]
    if all(x['strict_task_rate_difference_C_minus_B'] == 0 for x in paired):
        lines += ['No observed paired strict-rate difference; zero difference does not establish equivalence.']
    a = report['A_answerability_and_abstention']
    lines += ['', 'Conventional-resource condition A is a distinct practical comparison.',
              f"Exact target-source support coverage: {a['exact_target_source_coverage']}; "
              f"required abstention fields: {a['required_target_abstention_fields']}; "
              f"valid abstention-field rate: {a['valid_abstention_field_rate']}.",
              f"Supported conventional factual fields: {a['supported_conventional_factual_fields']}; "
              f"correct factual yield: {a['supported_conventional_factual_yield']}; "
              f"unsupported target-substitution task rate: {a['unsupported_target_substitution_task_rate']}.",
              'A correct abstention does not establish equivalent coverage. A/C differences cannot isolate presentation.',
              '', 'Limits and claim boundaries',
              'Ten final tasks were selected before responses, one per failure-mode stratum.',
              'The complete 60-item inventory is prepared; this run is not 60-task execution.',
              'Related source facts recur. Models are not independent providers, and answers are not biological replicates.',
              'No confidence interval, significance test or broad superiority claim is made for this feasibility subset.',
              'Instructions restricted solver access, but technical filesystem/tool blinding and actual traces were unavailable.',
              'Closed-field safety scores do not certify the absence of unsafe claims in all free prose.',
              'Separate execution logs preserve all launch attempts and any permitted transport retry.',
              'Independent AI scientific review and accountable human semantic scoring remain distinct requirements.',
              'No human participants, consent, researcher responses or human utility outcomes are supplied by this experiment.',
              'No GWAS estimate, testing family or biological claim was altered by these evaluations.',
              '', 'Reproduce from preserved actual response envelopes using the frozen score_benchmark.py,',
              'registration, task inventory and private condition-specific key; then run this reporter.',
              f"Registration SHA256: {hashlib.sha256(args.registration.read_bytes()).hexdigest()}"]
    (args.out / 'ACTUAL_FEASIBILITY_REPORT.txt').write_text('\n'.join(lines) + '\n')
    print(json.dumps({'out': str(args.out), 'completed': report['completed_cells']}))


if __name__ == '__main__':
    main()
