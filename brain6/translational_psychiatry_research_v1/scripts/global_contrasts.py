"""Retrospective source-defined rg contrasts, never independent replication."""
from pathlib import Path
import csv, hashlib, itertools, json, math, sys
import numpy as np
from scipy.stats import norm

ROOT = Path(__file__).resolve().parents[3]
OUT = ROOT / 'brain6/translational_psychiatry_research_v1/research_v1/global_contrasts'
INPUT = ROOT / 'brain6/paper/final_package_v1/BRAIN6_FINAL_GLOBAL.tsv'

def sha(path): return hashlib.sha256(path.read_bytes()).hexdigest()

def marginal_eligible(row):
    fields = ['rg','se','sleep_h2','sleep_h2_se','h2_disorder','h2_disorder_se']
    if any(not math.isfinite(float(row[x])) for x in fields): return False
    return (abs(float(row['rg'])) <= 1 and float(row['se']) > 0
            and float(row['sleep_h2']) > 0 and float(row['h2_disorder']) > 0
            and float(row['sleep_h2_se']) > 0 and float(row['h2_disorder_se']) > 0
            and float(row['sleep_h2'])/float(row['sleep_h2_se']) >= 4
            and float(row['h2_disorder'])/float(row['h2_disorder_se']) >= 4)

def bounded_contrast(delta, se1, se2):
    se_max = se1 + se2
    return se_max, 2*norm.sf(abs(delta)/se_max)

def write(name, rows):
    with (OUT/name).open('w',newline='') as f:
        writer=csv.DictWriter(f,fieldnames=list(rows[0]),delimiter='\t')
        writer.writeheader(); writer.writerows(rows)

def main():
    OUT.mkdir(parents=True,exist_ok=True)
    if sys.argv[1:] == ['freeze']:
        protocol={'id':'RETROSPECTIVE_GLOBAL_RG_CONTRASTS_V1','classification':'RETROSPECTIVE_SAME_SOURCE_DIAGNOSTIC','input_sha256':sha(INPUT),'code_sha256':sha(Path(__file__)), 'family_size':72,'alpha':.05,'disorders':['adhd','mdd','scz','bipolar'],'sleep_traits':'All12 inherited sleep traits; all6 psychiatric pairs each; no outcome selection','marginal_qc':'Finite bounded rg, positive SE/h2/h2SE and both h2 Z>=4; failed contrasts retained in72 denominator','primary_method':'Maximum two-sided asymptotic normal-Wald P over all admissible sampling covariances: Var(delta)<= (se1+se2)^2. Bonferroni72','sensitivity_correlation_grid':[-1,-.75,-.5,-.25,0,.25,.5,.75,1],'limitations':['Unknown joint estimator covariance','Marginal asymptotic normality and SE calibration assumed','Historical outcomes already known; not prospective discovery','Differences are source-defined genetic-correlation parameters, not psychiatric mechanisms','Ascertainment/cohort/phenotype/ancestry differences unresolved','Power bound assumes calibrated normal-Wald model, not actual cohort prediction power']}
        p=OUT/'protocol.json'
        if p.exists(): raise ValueError('Freeze already exists')
        p.write_text(json.dumps(protocol,indent=2)+'\n'); return
    if sys.argv[1:] != ['run']: raise SystemExit('freeze|run')
    cfg=json.loads((OUT/'protocol.json').read_text())
    assert cfg['input_sha256']==sha(INPUT) and cfg['code_sha256']==sha(Path(__file__))
    with INPUT.open() as f: data=list(csv.DictReader(f,delimiter='\t'))
    table={(x['sleep_trait'],x['brain_disorder']):x for x in data}
    sleep=list(dict.fromkeys(x['sleep_trait'] for x in data)); assert len(sleep)==12
    results=[]; sensitivity=[]; simultaneous_z=norm.ppf(1-.05/(2*72))
    for trait in sleep:
        for a,b in itertools.combinations(cfg['disorders'],2):
            x,y=table[trait,a],table[trait,b]; delta=float(x['rg'])-float(y['rg']); s1=float(x['se']);s2=float(y['se'])
            admitted=marginal_eligible(x) and marginal_eligible(y); semax,pmax=bounded_contrast(delta,s1,s2)
            results.append({'sleep_trait':trait,'disorder1':a,'disorder2':b,'rg1':x['rg'],'se1':s1,'rg2':y['rg'],'se2':s2,'difference':delta,'maximum_sampling_SE':semax,'maximum_normal_model_P':pmax if admitted else 'NA','bonferroni72_P':min(1,pmax*72) if admitted else 'NA','simultaneous_lower':delta-simultaneous_z*semax if admitted else 'NA','simultaneous_upper':delta+simultaneous_z*semax if admitted else 'NA','sufficient_abs_difference_for80pct_normal_power':(simultaneous_z+norm.ppf(.8))*semax if admitted else 'NA','status':'RETROSPECTIVE_MODEL_CONDITIONAL_DIFFERENCE' if admitted and pmax<=.05/72 else 'NO_BOUNDED_DIFFERENCE' if admitted else 'NOT_ADMITTED_MARGINAL_H2_QC','input_sha256':cfg['input_sha256'],'independence':'SAME_DISCOVERY_SOURCES'})
            for corr in cfg['sensitivity_correlation_grid']:
                se=math.sqrt(max(0,s1*s1+s2*s2-2*corr*s1*s2))
                p=2*norm.sf(abs(delta)/se) if se else (0 if delta else 1)
                sensitivity.append({'sleep_trait':trait,'disorder1':a,'disorder2':b,'assumed_sampling_error_correlation':corr,'difference_SE':se,'normal_model_P':p if admitted else 'NA','status':results[-1]['status']})
    assert len(results)==72
    write('contrasts.tsv',results);write('covariance_sensitivity.tsv',sensitivity)
    summary={'planned_contrasts':len(results),'admitted_contrasts':sum(x['status']!='NOT_ADMITTED_MARGINAL_H2_QC' for x in results),'bounded_differences':sum(x['status']=='RETROSPECTIVE_MODEL_CONDITIONAL_DIFFERENCE' for x in results),'classification':cfg['classification'],'input_sha256':cfg['input_sha256']}
    (OUT/'summary.json').write_text(json.dumps(summary,indent=2)+'\n');print(summary)

if __name__=='__main__':main()
