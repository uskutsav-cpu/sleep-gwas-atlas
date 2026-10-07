#!/usr/bin/env python3
"""Read-only historical evidence audit; independent arithmetic, no promotion."""
from __future__ import annotations
import csv
import hashlib
import json
import math
from collections import Counter
from pathlib import Path
import subprocess

ROOT = Path(__file__).resolve().parents[3]
AREA = ROOT / 'brain6/confirmatory_v2'
PACKAGE = ROOT / 'brain6/paper/final_package_v1'

def sha(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as f:
        for b in iter(lambda: f.read(1 << 20), b''): h.update(b)
    return h.hexdigest()

def rows(path):
    with Path(path).open(newline='') as f: return list(csv.DictReader(f, delimiter='\t'))

def write_tsv(path, records, fields=None):
    fields = fields or list(records[0])
    with Path(path).open('w', newline='') as f:
        w = csv.DictWriter(f, fields, delimiter='\t'); w.writeheader(); w.writerows(records)

def bh(values):
    n = len(values); order = sorted(range(n), key=lambda i: values[i]); result = [0.] * n; carry = 1.
    for rank in range(n, 0, -1):
        idx = order[rank-1]; carry = min(carry, values[idx]*n/rank); result[idx] = carry
    return result

def reweight(vector, factor):
    weights = list(vector); weights[4] *= factor; total = sum(weights)
    if total <= 0 or any(v < 0 or not math.isfinite(v) for v in weights): raise ValueError('invalid posterior')
    return [v / total for v in weights]

def overlap(c, block):
    region = c['candidate_locus_id'].split('_chr')[-1].split('_')
    chromosome, start, end = map(int, region)
    return c['pair_id'] == block['pair_id'] and chromosome == int(block['chr']) and start <= int(block['end']) and end >= int(block['start'])

def main():
    qc = AREA / 'qc'; qc.mkdir(parents=True, exist_ok=True)
    provenance = json.loads((PACKAGE/'BRAIN6_FINAL_PROVENANCE.json').read_text())
    checks = []
    for family, mapping, base in [('source',provenance['source_sha256'],ROOT),('output',provenance['output_sha256'],PACKAGE)]:
        for name, expected in mapping.items():
            p = Path(name) if name.startswith('/') else base/name
            observed = sha(p) if p.is_file() else ''
            checks.append({'family':family,'path':name,'expected_sha256':expected,'observed_sha256':observed,
                'status':'MATCH' if observed == expected else ('MISSING' if not observed else 'MISMATCH')})
    for directory, receipt_name, key in [
        ('brain6_alternative_local_validation_v1','run_provenance.json','output_sha256'),
        ('brain6_exploratory_finemap_coloc_v2','integrity_provenance.json','result_sha256'),
        ('brain6_exploratory_functional_v1','provenance.json','output_sha256'),
        ('brain6_exploratory_enrichment_v3','provenance.json','outputs_sha256')]:
        base = ROOT/'brain6/results'/directory
        for name, expected in json.loads((base/receipt_name).read_text())[key].items():
            p=base/name; observed=sha(p) if p.is_file() else ''
            checks.append({'family':'upstream_output','path':str(p.relative_to(ROOT)), 'expected_sha256':expected,
                'observed_sha256':observed,'status':'MATCH' if observed == expected else ('MISSING' if not observed else 'MISMATCH')})
    write_tsv(qc/'hash_checks.tsv',checks)
    assert not any(r['status']=='MISMATCH' for r in checks), 'Changed historical data'
    data = {n:rows(PACKAGE/f'BRAIN6_FINAL_{n}.tsv') for n in ['GLOBAL','CANDIDATES','REGIONS','LAVA','ALT_LOCAL_BLOCKS','TRAIT_COLOC','TRAIT_COLOC_PRIOR_SENSITIVITY','TISSUE_CELLTYPE','PATHWAYS']}
    global_rows=data['GLOBAL']; lava=data['LAVA']; blocks=data['ALT_LOCAL_BLOCKS']; candidates=data['CANDIDATES']
    significant_global=[r for r in global_rows if float(r['original_BH_FDR_q']) < .05]
    assert all((r['significance_under_original_396_family']=='True') == (float(r['original_BH_FDR_q']) < .05) for r in global_rows)
    assert len(global_rows)==72 and len({(r['sleep_trait'],r['brain_disorder']) for r in global_rows})==72
    assert len(candidates)==25 and len(data['REGIONS'])==20
    lava_counts=Counter(r['status'] for r in lava)
    assert lava_counts=={'TESTED':13745,'NOT_RUN':3720}
    assert all(r['family_decision']=='FAILED_QC_NOT_PROMOTED' for r in lava)
    alpha=.05/len(blocks); significant_blocks=[]; max_fwer_error=0.
    for b in blocks:
        if b['analysis_status']=='ESTIMATED':
            p=float(b['p']); expected=min(1,p*len(blocks))
            max_fwer_error=max(max_fwer_error,abs(float(b['p_fwer'])-expected))
            if p < alpha:
                entry=dict(b); entry['overlapping_candidates']=';'.join(c['candidate_locus_id'] for c in candidates if overlap(c,b)); significant_blocks.append(entry)
    assert len(blocks)==8465 and len(significant_blocks)==3 and max_fwer_error < 1e-12
    assert all(not b['overlapping_candidates'] for b in significant_blocks)
    write_tsv(qc/'significant_secondary_blocks.tsv',significant_blocks)
    prior_by={(r['candidate_locus_id'],float(r['p12'])):r for r in data['TRAIT_COLOC_PRIOR_SENSITIVITY']}
    replay=[]
    for r in data['TRAIT_COLOC']:
        if r['status'] != 'EXPLORATORY_TRAIT_COLOC_ABF': continue
        vector=[float(r[f'PP.H{k}']) for k in range(5)]; expected=reweight(vector,.1)
        observed=prior_by[(r['candidate_locus_id'],1e-6)]
        error=max(abs(expected[k]-float(observed[f'PP.H{k}'])) for k in range(5))
        assert abs(sum(vector)-1)<1e-10 and error < 1e-10
        replay.append({'candidate_locus_id':r['candidate_locus_id'],'default_H4':vector[4], 'low_prior_H4_replayed':expected[4], 'max_absolute_error':error, 'interpretation':'POSTERIOR_ARITHMETIC_ONLY_SINGLE_SIGNAL'})
    write_tsv(qc/'coloc_prior_replay.tsv',replay)
    enrichment={}
    for name,record_type,qname,n in [('TISSUE_CELLTYPE','GTEX_BULK_TISSUE_18_REGION_TEST','q_bh_54',54),('PATHWAYS','REACTOME_POSITIONAL_18_REGION_TEST','q_bh_1680',1680)]:
        rr=[r for r in data[name] if r['record_type']==record_type]
        assert len(rr)==n
        ps=[float(r['p_empirical_high']) for r in rr]; computed=bh(ps)
        errors=[abs(q-float(r[qname])) for q,r in zip(computed,rr)]
        empirical=max(abs(float(r['p_empirical_high'])-(int(r['null_extreme_count'])+1)/10001) for r in rr)
        assert max(errors)<1e-10 and empirical<1e-10
        enrichment[name]={'n_tests':n,'minimum_q':min(computed),'n_significant':sum(q < .05 for q in computed),'max_bh_error':max(errors),'max_empirical_p_error':empirical}
    summary={'run_identity':'brain6-continuation-20261007-v1','status':'PASS_AVAILABLE_TABLES_ORIGINAL_INPUT_REPLAY_INCOMPLETE',
        'starting_commit':subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip(),
        'hash_counts':dict(Counter(r['status'] for r in checks)), 'global_rows':72, 'inherited_significant_rows':len(significant_global),
        'significant_by_disorder':dict(Counter(r['brain_disorder'] for r in significant_global)),
        'candidate_count':25,'region_count':20,'candidate_classes':dict(Counter(r['final_evidence_class'] for r in candidates)),
        'lava_counts':dict(lava_counts),'lava_not_run_by_trait':dict(Counter(r['phen'] for r in lava if r['status']=='NOT_RUN')),
        'lava_not_run_by_reason':dict(Counter(r['reason'] for r in lava if r['status']=='NOT_RUN')),
        'lava_max_not_run':873,'lava_recovery_required':3720-873, 'secondary_counts':dict(Counter(r['analysis_status'] for r in blocks)),
        'secondary_alpha':alpha,'secondary_significant_blocks':3,'candidate_overlaps':0,'max_fwer_error':max_fwer_error,
        'coloc_estimated':len(replay),'default_H4_at_least_08':sum(r['default_H4']>=.8 for r in replay),
        'lower_prior_H4_at_least_08':sum(r['low_prior_H4_replayed']>=.8 for r in replay),'enrichment':enrichment,
        'package_provenance_sha256':sha(PACKAGE/'BRAIN6_FINAL_PROVENANCE.json'),'protocol_sha256':sha(AREA/'NEW_ANALYSIS_PROTOCOL.md'),
        'audit_script_sha256':sha(__file__), 'limitations':['No raw-GWAS or native-LD replay','Original canonical input absent in this checkout','No independent replication','Posterior replay validates arithmetic only']}
    (qc/'baseline_verification.json').write_text(json.dumps(summary,indent=2,sort_keys=True)+'\n')
    print(json.dumps(summary,sort_keys=True))

if __name__=='__main__': main()
