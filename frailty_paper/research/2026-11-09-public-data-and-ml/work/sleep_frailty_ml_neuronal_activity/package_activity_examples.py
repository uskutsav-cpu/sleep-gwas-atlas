"""Package audited measured-activity examples without inventing allele targets."""
import collections
import csv
import gzip
import hashlib
import json
import math
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / 'outputs/sleep_frailty_ml_neuronal_activity'


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def rows(path, delimiter=','):
    with gzip.open(path, 'rt') as handle:
        yield from csv.DictReader(handle, delimiter=delimiter)


def main():
    audit_path = OUT / 'additional_activity_independent_no_prefix_audit.json'
    audit = json.loads(audit_path.read_text())
    assert audit['passed']
    screen = OUT / 'additional_activity_closed_core_candidate_screen.csv.gz'
    assert sha(screen) == audit['candidate_file_sha256']
    kept = {r['name']: r for r in rows(screen) if r['passes_shared_core_and_duplicate_screen'] == 'True'}
    original = json.loads((OUT / 'additional_neuronal_activity_inventory_receipt.json').read_text())
    correction = OUT / 'activity_target_interpretation_correction.json'
    assert sha(correction) == original['source_target_description_correction']['sha256']
    metadata = OUT / Path(original['sources'][0]['path']).name
    fitted = OUT / Path(original['sources'][1]['path']).name
    assert sha(metadata) == original['sources'][0]['sha256']
    assert sha(fitted) == original['sources'][1]['sha256']
    # Independent file join repeats exact original coefficient identity for
    # this packaged subset, retaining source string precision.
    coefficients = {r['name']: float(r['logFC']) for r in rows(fitted, '\t')}
    examples = []
    for source in rows(metadata, '\t'):
        name = source['name']
        if name not in kept:
            continue
        prior = kept[name]
        value = float(source['logFC'])
        assert math.isfinite(value) and abs(value - coefficients[name]) < 1e-10
        assert source['logFC'] == prior['source_logFC']
        sequence = source['sequence'].upper()
        reverse = sequence.translate(str.maketrans('ACGT', 'TGCA'))[::-1]
        assert sha_bytes(min(sequence, reverse).encode()) == prior['physical_sequence_sha256_both_orientations']
        assert source['ref'] == 'GRCh38' and source['category'] == 'element'
        examples.append({'source_construct_name': name, 'sequence_in_original_assay_orientation': sequence,
            'original_BCalm_activity_log2_RNA_over_DNA': source['logFC'],
            'activity_target_semantics': 'Unshifted fitted construct activity; not a difference from scrambled controls',
            'allelic_effect_target': '', 'insomnia_or_frailty_risk_target': '',
            'source_study_role': 'Salomon', 'cell_context': 'NGN2_induced_excitatory_neuron',
            'assay': 'lentiMPRA', 'declared_genome_build': source['ref'],
            'chromosome': source['chr'], 'interval_start_0_based': str(int(float(source['start']))),
            'interval_end_exclusive': str(int(float(source['end']))), 'source_declared_strand': source['strand'],
            'full_sequence_matches_declared_reference_in_one_orientation': True,
            'component_id_for_source_exclusion': prior['source_component_after_shared_145bp_closure'],
            'canonical_full_sequence_sha256': prior['physical_sequence_sha256_both_orientations'],
            'target_kind': 'auxiliary_construct_activity_only',
            'source_metadata_sha256': original['sources'][0]['sha256'],
            'source_fitted_results_sha256': original['sources'][1]['sha256']})
    assert len(examples) == len(kept) == audit['candidate_rows']
    assert len({r['source_construct_name'] for r in examples}) == len(examples)
    assert len({r['sequence_in_original_assay_orientation'] for r in examples}) == len(examples)
    target = OUT / 'audited_additional_neuron_activity_examples.csv.gz'
    with gzip.open(target, 'wt') as handle:
        writer = csv.DictWriter(handle, fieldnames=list(examples[0]))
        writer.writeheader()
        writer.writerows(examples)
    # Re-read the actual saved artifact, including intentionally absent targets.
    saved = list(rows(target))
    assert len(saved) == len(examples)
    assert all(r['allelic_effect_target'] == r['insomnia_or_frailty_risk_target'] == '' for r in saved)
    assert all(r['original_BCalm_activity_log2_RNA_over_DNA'] == kept[r['source_construct_name']]['source_logFC'] for r in saved)
    receipt = {'prepared': True, 'example_rows': len(examples),
        'distinct_canonical_full_sequences': len({r['canonical_full_sequence_sha256'] for r in examples}),
        'genomic_components': len({r['component_id_for_source_exclusion'] for r in examples}),
        'sequence_lengths': dict(collections.Counter(len(r['sequence_in_original_assay_orientation']) for r in examples)),
        'interpretation': original['source_target'],
        'permitted_role': 'Candidate additional construct-activity supervision in a separately specified comparison. Entire source omitted whenever Salomon is held out. Never counted as a new independent study.',
        'missing_allelic_and_disease_targets_preserved_as_empty': True, 'model_predictions_used_as_labels': False,
        'actual_DNA_or_genomic_context_invented': False, 'significance_or_effect_size_filter_applied': False,
        'opposite_assay_orientations_averaged': False, 'coefficient_values_changed': False,
        'current_comparison_changed': False, 'new_models_trained_with_this_artifact': False,
        'remaining': 'Fix the additional-data comparison, activity centering and sampling rules before fitting. Evaluate actual allelic effect size and cell/nearby differences on unchanged source-held-out development data; preserve the reserved study for a justified frozen candidate.',
        'independent_no_prefix_audit_sha256': sha(audit_path), 'interpretation_correction_sha256': sha(correction),
        'artifact_sha256': sha(target), 'code_sha256': sha(Path(__file__))}
    (OUT / 'audited_additional_neuron_activity_package_receipt.json').write_text(json.dumps(receipt, indent=2) + '\n')
    (OUT / 'README.txt').write_text(
        'Additional measured neuron activity for the Sleep–Frailty molecular model\n\n'
        f'{len(examples):,} original source constructs; 3,329 canonical DNA sequences; 2,535 connected region groups.\n'
        'The examples describe whole-sequence activity in induced excitatory neurons.\n'
        'The original BCalm coefficients are unshifted fitted log2 RNA/DNA activity.\n'
        'Scrambled controls determine statistical testing; their mean is not subtracted from these released coefficients.\n'
        'The original ambiguous description is archived; activity_target_interpretation_correction.json documents the correction.\n\n'
        'Eligibility uses original DNA/build checks, entire-source connected groups, a conservative 2kb history screen,\n'
        'all historical DNA and the entire reserved study, exact shared 145bp cores in both orientations, and duplicate-ID checks.\n'
        'An independent scanner repeated every window without the prefix shortcut or sequence deduplication.\n'
        'No significance or effect-size threshold selected examples. No reserved experimental outcomes were read.\n\n'
        'The allelic-effect and insomnia/frailty-risk target fields are intentionally empty, never zero.\n'
        'Do not count opposite assay orientations as independent biological regions.\n'
        'The source is already represented in development data. It is not a new independent validation study.\n'
        'Exclude the whole artifact when the Salomon study is held out.\n'
        'The original raw counts and historical BCalm runtime have not been refitted.\n\n'
        'These examples are packaged for a future prespecified activity-supervision comparison.\n'
        'They have not been inserted into the currently running fixed comparison.\n')
    print(json.dumps(receipt), flush=True)


def sha_bytes(data):
    return hashlib.sha256(data).hexdigest()


if __name__ == '__main__':
    main()
