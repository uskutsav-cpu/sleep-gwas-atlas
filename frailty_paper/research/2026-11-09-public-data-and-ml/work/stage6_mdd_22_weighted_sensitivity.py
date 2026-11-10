import csv,json,math,hashlib
import numpy as np
from pathlib import Path
p=Path('outputs/unvalidated_22_shared_direction_score.csv');score=list(csv.DictReader(p.open()));r=json.loads(Path('outputs/shared22_reference_diagnostic.json').read_text());assert r['score_sha256'] if 'score_sha256' in r else True
assert hashlib.sha256(p.read_bytes()).hexdigest()==r['shared_score_sha256'];den=float(r['shared_allele_count_sd']);assert r['n_reference_eur']==503 and r['counted_variants']==len(score)==22
freq={x['variant']:x for x in csv.DictReader(open('outputs/22_score_reference_group_allele_frequencies.csv'))};hits={(x['chr'],x['pos'],frozenset(x['pair'])):x for x in json.loads(Path('work/mdd_symptoms/candidate_hits_private.json').read_text())}
res=[]
for trait in ['insomnia','fatigue']:
 a=[];z=[];se=[]
 for x in score:
  v=freq[x['variant']];assert v['counted_allele'] in (x['counted_allele'],x['other_allele']);p0=float(v['reference_EUR_frequency']);p0=p0 if v['counted_allele']==x['counted_allele'] else 1-p0;assert 0<p0<1
  k=(x['chr_GRCh37_assumed'],int(x['position_GRCh37_assumed']),frozenset((x['counted_allele'],x['other_allele'])));h=hits[k]['hits'][trait][0];s=float(h['se']);assert s>0;se.append(s);a.append(math.sqrt(2*p0*(1-p0)));z.append(float(h['beta'])/s*(1 if h['a1']==x['counted_allele'] else -1))
 a=np.array(a);z=np.array(z);q=np.array(se)*a;Z=float(a@z/den)
 res.append({'trait':trait,'markers':22,'weighted_z':Z,'approximate_nominal_p':math.erfc(abs(Z)/math.sqrt(2)),'se_times_reference_sd_hwe_min':float(q.min()),'se_times_reference_sd_hwe_median':float(np.median(q)),'se_times_reference_sd_hwe_max':float(q.max())})
Path('outputs/stage6_mdd_22_weighted_sensitivity.json').write_text(json.dumps({'results':res,'reference_score_sd':den,'reference':'1000 Genomes EUR n=503 cached summary; score SD is empirical with cross-chromosome correlations','assumptions':'score weights equal one; HWE SD from reference allele frequency for each numerator; empirical score SD denominator; source sample overlap and genotype variance/LD unverified','limitations':'Approximate postselected marginal association, case-enriched symptoms, not individual-level score prediction or prespecified replication.'},indent=2)+'\n');print(json.dumps(res,indent=2))
