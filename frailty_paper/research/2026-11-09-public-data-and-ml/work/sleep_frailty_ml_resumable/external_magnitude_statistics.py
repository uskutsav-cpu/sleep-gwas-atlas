"""Low-memory paired magnitude statistics for an already locked external test.

This module never opens a dataset or authorizes outcome joining. A future
evaluation runner must enforce the independently frozen model/input/outcome
locks. Self-checks use only synthetic measurements, not reserved outcomes.
"""
import collections
import hashlib
import json
import math
from pathlib import Path
import random


def quantile(values, q):
    values = sorted(values)
    position = q*(len(values)-1)
    lower = math.floor(position)
    fraction = position-lower
    return values[lower]*(1-fraction)+values[min(lower+1,len(values)-1)]*fraction


def evaluate(rows, protocol):
    identities = [str(row['construct']) for row in rows]
    if len(set(identities)) != len(identities):
        raise ValueError('Duplicate physical construct identities must not become extra observations.')
    variant_components = {}
    for row in rows:
        if row['variant'] in variant_components and variant_components[row['variant']] != row['group']:
            raise ValueError('One variant cannot be split across supposedly independent genomic components.')
        variant_components[row['variant']] = row['group']
    counts = {'constructs':len(rows),'distinct_variants':len({r['variant'] for r in rows}),
              'connected_components':len({r['group'] for r in rows})}
    result = {'counts':counts,'primary_retention_passed':False,'model_promoted':False,
              'confidence_or_abstention_success_established':False,
              'scope':'Low-level paired point-error statistics only; not an outcome-lock or biological provenance audit.'}
    if not rows:
        return {**result,'status':'inconclusive','reason':'No eligible constructs.'}
    invalid = []
    for row in rows:
        for key in ['observed','scoreless','scores','baseline']:
            value = float(row[key])
            if not math.isfinite(value) or value < 0:
                invalid.append({'construct':str(row['construct']),'field':key})
    if invalid:
        return {**result,'status':'inconclusive','reason':'Incomplete or invalid paired predictions/observations; no rows were dropped.',
                'invalid_values':invalid}
    errors = [{key:abs(float(r['observed'])-float(r[key])) for key in ['scoreless','scores','baseline']} for r in rows]
    mae = {key:math.fsum(r[key] for r in errors)/len(rows) for key in errors[0]}
    variant_groups,components = collections.defaultdict(list),collections.defaultdict(list)
    for row,error in zip(rows,errors):
        variant_groups[row['variant']].append(error)
        components[row['group']].append(error)
    balanced = {}
    for label,groups in [('variant_balanced',variant_groups),('component_balanced',components)]:
        balanced[label] = {key:math.fsum(math.fsum(row[key] for row in group)/len(group)
                                       for group in groups.values())/len(groups) for key in mae}
    gain = 1-mae['scores']/mae['scoreless'] if mae['scoreless'] > 0 else None
    skill = 1-mae['scores']/mae['baseline'] if mae['baseline'] > 0 else None
    result.update(MAE=mae,sensitivity_MAE=balanced,relative_gain_vs_scoreless=gain,
                  relative_gain_vs_fixed_baseline=skill)
    enough = all(counts[key] >= threshold for key,threshold in protocol['minimum_evaluable'].items())
    if not enough:
        return {**result,'status':'inconclusive','reason':'Insufficient eligible constructs, variants or genomic components.'}
    settings = protocol['uncertainty']
    # Resample complete components. Preaggregating sums avoids treating the
    # repeated constructs inside each component as independently resampled.
    blocks = [(len(group),math.fsum(row['scoreless'] for row in group),
               math.fsum(row['scores'] for row in group)) for key,group in sorted(components.items(),key=lambda item:str(item[0]))]
    rng = random.Random(settings['seed'])
    gains = []
    invalid_resamples = 0
    for _ in range(settings['replicates']):
        chosen = [blocks[rng.randrange(len(blocks))] for _ in blocks]
        total_control = math.fsum(block[1] for block in chosen)
        total_score = math.fsum(block[2] for block in chosen)
        if total_control <= 0:
            invalid_resamples += 1
        else:
            gains.append(1-total_score/total_control)
    interval = [quantile(gains,q) for q in settings['interval_quantiles']] if gains and not invalid_resamples else None
    criteria = protocol['primary_retention_criteria']
    checks = {'relative_MAE_gain':gain is not None and gain >= criteria['score_arm_relative_MAE_reduction_vs_scoreless_at_least'],
              'beat_fixed_baseline':skill is not None and skill >= criteria['score_arm_relative_MAE_reduction_vs_fixed_baseline_at_least'],
              'clustered_CI_supports_gain':interval is not None and interval[0] > criteria['paired_component_bootstrap_relative_gain_CI95_lower_strictly_above']}
    result.update(status='evaluated' if not invalid_resamples else 'inconclusive',
                  paired_component_bootstrap_gain_CI95=interval,
                  bootstrap_replicates=settings['replicates'],undefined_zero_control_error_resamples=invalid_resamples,
                  primary_retention_criteria=checks,primary_retention_passed=all(checks.values()))
    return result


def self_check():
    root = Path(__file__).resolve().parents[2]
    protocol_path = root/'outputs/sleep_frailty_ml_resumable/AVI_preparation/fresh_study_benefit_protocol.json'
    protocol = json.loads(protocol_path.read_text())
    rows = [{'construct':f'c{i}','variant':f'v{i//3}','group':f'g{i//6}',
             'observed':2.,'scoreless':1.,'scores':1.2,'baseline':1.} for i in range(120)]
    good = evaluate(rows,protocol)
    assert good['primary_retention_passed'] and math.isclose(good['relative_gain_vs_scoreless'],.2,abs_tol=1e-14)
    assert all(math.isclose(q,.2,abs_tol=1e-14) for q in good['paired_component_bootstrap_gain_CI95'])
    equal = evaluate([{**r,'scores':r['scoreless']} for r in rows],protocol)
    assert not equal['primary_retention_passed'] and equal['paired_component_bootstrap_gain_CI95'] == [0.,0.]
    missing = evaluate([{**r,'scores':float('nan')} if i == 0 else r for i,r in enumerate(rows)],protocol)
    assert missing['status'] == 'inconclusive' and missing['counts']['constructs'] == 120
    correlated = evaluate([{**r,'group':'one_group'} for r in rows],protocol)
    assert correlated['status'] == 'inconclusive' and not correlated['primary_retention_passed']
    clustered_rows = [{**r,'scores':1.4 if int(r['group'][1:]) < 10 else .6} for r in rows]
    clustered = evaluate(clustered_rows,protocol)
    # Twenty perfectly correlated clusters produce a much wider interval than
    # 120 independent rows. The binomial cluster reference has limits +/- .16;
    # a generous bound allows finite resampling while rejecting row bootstrap.
    assert clustered['paired_component_bootstrap_gain_CI95'][0] < -.12
    assert clustered['paired_component_bootstrap_gain_CI95'][1] > .12
    zero = evaluate([{**r,'observed':1.,'scoreless':1.,'scores':1.,'baseline':1.} for r in rows],protocol)
    assert zero['status'] == 'inconclusive' and not zero['primary_retention_passed']
    assert not good['confidence_or_abstention_success_established'] and not good['model_promoted']
    try:
        evaluate(rows+[rows[0]],protocol)
    except ValueError:
        pass
    else:
        raise AssertionError('Duplicate constructs accepted.')
    try:
        evaluate([{**r,'group':'wrong_component'} if i == 0 else r for i,r in enumerate(rows)],protocol)
    except ValueError:
        pass
    else:
        raise AssertionError('One variant split across components.')
    assert math.isclose(quantile([0.,1.,2.,3.],.25),.75)
    result = {'passed':True,'synthetic_only':True,'real_experimental_outcomes_read':False,
              'cases':['known 20% paired improvement','identical arms rejected','missing predictions inconclusive without dropping rows',
                       'one correlated component cannot validate','zero-control-error resamples inconclusive','duplicate construct rejected',
                       'correlated clusters yield wider uncertainty than row resampling','one variant split across components rejected',
                       'linear percentile interpolation','point benefit does not imply confidence or promotion'],
              'code_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
              'external_protocol_sha256':hashlib.sha256(protocol_path.read_bytes()).hexdigest()}
    target = protocol_path.with_name('external_magnitude_statistics_synthetic_checks.json')
    target.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result),flush=True)


if __name__ == '__main__':
    self_check()
