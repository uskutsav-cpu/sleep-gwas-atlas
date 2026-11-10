"""Lock the fresh-study benefit test without reading experimental outcomes."""
import csv
import datetime
import hashlib
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT/'outputs/sleep_frailty_ml_resumable'
AVI = OUT/'AVI_preparation'


def sha(path):
    digest = hashlib.sha256()
    with Path(path).open('rb') as handle:
        for chunk in iter(lambda:handle.read(1024**2),b''):
            digest.update(chunk)
    return digest.hexdigest()


def main():
    assert not list((OUT/'AVI_development').glob('*_record.json')), 'Must lock this protocol before new scientific AVI fits.'
    assert not list((OUT/'AVI_development').glob('*_started.json'))
    guard_path = ROOT/'outputs/sleep_frailty_ml_bounded/AVI_exposure/AVI_development_chromosome_exposure_receipt.json'
    guard = json.loads(guard_path.read_text())
    assert guard['passed_publication_and_identity_screen'] and not guard['new_study_outcomes_read']
    candidates = guard_path.with_name('Ulirsch_candidates_outside_AVI_development_chromosomes.csv')
    assert sha(candidates) == guard['candidate_file_sha256']
    with candidates.open() as handle:
        rows = list(csv.DictReader(handle))
    assert len(rows) == guard['candidate_tests_after'] == 183
    reconstruction = ROOT/'outputs/sleep_frailty_ml_measurement/Ulirsch_count_reanalysis_audit.json'
    assert json.loads(reconstruction.read_text())['passed']
    count_protocol = reconstruction.with_name('Ulirsch_count_reanalysis_protocol.json')
    sealed_receipt = reconstruction.with_name('Ulirsch_sealed_reconstruction_receipt.json')
    assert sha(count_protocol) == json.loads(reconstruction.read_text())['protocol_sha256']
    result = {
        'time_UTC':datetime.datetime.now(datetime.timezone.utc).isoformat(),
        'scope':'Additional-input benefit on molecular reporter-effect magnitude. No patient risk, disease causation, genetic-score or cross-cell validation claim.',
        'primary_study':'Ulirsch 2016, basal K562, GEO GSE87711; GATA1 overexpression remains outside this test.',
        'candidate_rule':'Use the already declared native-gate candidate rule: adapted reader/common epoch if its native gate passes; otherwise the original frozen reader at exactly three epochs for this separate input question.',
        'training_and_development_requirement':'Complete the unchanged finite native comparison and the separate matched AVI development gate. A failed AVI development gate stops this test before reserved outcomes are joined.',
        'matched_arms':'Same expanded model, initialization, supervised observations, epochs, optimizer, updates and availability flag. Only permitted composite AVI raw/PHRED inputs differ.',
        'identity_pool_before_sealed_count_QC':{'constructs':183,'variants':64,'connected_genomic_components':47,
                                             'candidate_file':str(candidates.relative_to(ROOT)),
                                             'candidate_sha256':sha(candidates)},
        'eligibility':'Retain the existing complete-source DNA/history/reader-benchmark and AVI train/validation chromosome exclusions, followed only by the already frozen count-availability rule. Do not select by significance, activity, measured size, SE, prediction or error.',
        'primary_target':'Absolute regularized basal RNA/DNA allelic contrast from the frozen reanalysis; not the original published fitted coefficient or denoised truth.',
        'primary_unit':'One unique physical allele-pair construct. Different source tiles are distinct constructs, not independent people. Report variant-balanced and component-balanced sensitivity estimates separately.',
        'primary_metric':'Unweighted mean absolute error in observed log2-effect magnitude on exactly the same complete set of eligible constructs for both arms.',
        'baseline':'Fixed training-only exact K562 line median when supported, otherwise the fixed training-only K562 cell median. No baseline is estimated from reserved outcomes.',
        'minimum_evaluable':{'constructs':100,'distinct_variants':40,'connected_components':20},
        'paired_prediction_completeness':'Every count-QC-eligible identity requires finite predictions from both arms and the fixed baseline. Missing predictions make the primary test inconclusive; they cannot be silently dropped. Missing AVI is a neutral input plus the identical availability flag in both arms, not a reason to remove cases.',
        'primary_retention_criteria':{'score_arm_relative_MAE_reduction_vs_scoreless_at_least':0.01,
                                      'score_arm_relative_MAE_reduction_vs_fixed_baseline_at_least':0.01,
                                      'paired_component_bootstrap_relative_gain_CI95_lower_strictly_above':0.0},
        'uncertainty':{'method':'Paired percentile bootstrap of whole existing connected genomic components, carrying all their constructs in each sampled component.',
                       'replicates':5000,'seed':20261109,'interval_quantiles':[0.025,0.975],
                       'interpretation':'Within-study clustered sampling uncertainty; components are not proven ancestry-specific LD blocks, and bootstrap intervals do not establish transport to other assays.'},
        'required_sensitivity_reports':['MAE with equal weight per distinct variant and per connected component.','Per-chromosome results, without model or threshold changes.','Descriptive nearby distinct-variant ordering where available; do not treat multiple tiles of the same variant as distinct DNA changes.','Source count availability, eligible and excluded identities, and complete prediction coverage.'],
        'confidence_boundary':'Freeze any development-qualified confidence policy and its source hashes before external prediction/outcome joining. If no method qualifies, explicitly lock no confidence release: the paired point-prediction test may proceed, but no trustworthy-interval, useful-abstention or project-ranking success follows. All cases withheld never counts as success.',
        'required_before_outcome_join':['Locked native assessment, AVI plan and successful matched development decision.','Both fixed fitted checkpoints and every source/input/normalization hash.','Actual 145-bp native-construct inference verified without invented flanks.','Exact matched input identities, permitted score values/missingness and fixed baseline.','Both complete prediction files, source-mapped confidence policy or explicit no-confidence-release, and this protocol hash.'],
        'stop_rule':'Reject the added input for the primary magnitude endpoint if any retention criterion fails. Insufficient eligible groups or incomplete predictions means inconclusive, not success. Preserve both arms and all failed criteria; never tune on this study or replace it with a favorable secondary source.',
        'retention_scope':'Passing retains an experimental magnitude-input candidate only. Useful confidence, nearby/cross-cell differences and sleep–frailty relevance require their own evidence before model promotion.',
        'additional_neural_source':'Lee2025/GSE244011 is a newly catalogued candidate requiring full source/overlap/QC audits. It does not replace or rescue this primary study.',
        'reserved_outcome_values_read_by_this_script':False,
        'new_model_predictions_or_fits_made':False,
        'evidence_sha256':{str(p.relative_to(ROOT)):sha(p) for p in [guard_path,candidates,reconstruction,count_protocol,sealed_receipt]},
        'code_sha256':sha(Path(__file__))}
    destination = AVI/'fresh_study_benefit_protocol.json'
    assert not destination.exists(), 'An existing external protocol must not be overwritten.'
    destination.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({'external_protocol_frozen':True,'candidate_constructs_before_count_QC':183,
                      'outcome_values_read':False,'path':str(destination)}))


if __name__ == '__main__':
    main()
