import csv,json,numpy as np
from pathlib import Path
score15=list(csv.DictReader(open('outputs/unvalidated_15_variant_score_prototype.csv')));score22=list(csv.DictReader(open('outputs/unvalidated_22_shared_direction_score.csv')));x=np.load('work/stage6_hfrs_sparse_score_EUR503_genotypes.npz',allow_pickle=True);assert list(x['variants'])==[r['variant'] for r in score15];freq={}
for j,r in enumerate(score15):freq[r['variant']]={r['counted_allele']:float(x['genotypes'][:,j].mean()/2),r['other_allele']:float(1-x['genotypes'][:,j].mean()/2)}
for r in csv.DictReader(open('outputs/22_score_reference_group_allele_frequencies.csv')):
 target=next(y for y in score22 if y['variant']==r['variant']); a=r['counted_allele'];b=next(v for v in [target['counted_allele'],target['other_allele']] if v!=a);p=float(r['reference_EUR_frequency']);
 if r['variant'] in freq:assert abs(freq[r['variant']][a]-p)<1e-10
 else:freq[r['variant']]={a:p,b:1-p}
hits=json.loads(Path('work/mdd_symptoms/candidate_hits_private.json').read_text());results={}
for trait in ['insomnia','fatigue']:
 dif=[];info=[];neff=[]
 for r in hits:
  h=r['hits'][trait][0];dif.append((abs(float(h['fcon'])-freq[r['id']][h['a1']]),r['id']));info.append(float(h['info']));neff.append(float(h['neff']))
 results[trait]={'marker_count':len(dif),'info_range':[min(info),max(info)],'neff_range':[min(neff),max(neff)],'max_absolute_control_vs_eur_reference_frequency_difference':max(dif)[0],'markers_frequency_difference_gt_0_1':[rs for d,rs in dif if d>.1]}
Path('outputs/stage6_mdd_candidate_qc.json').write_text(json.dumps({'summary':results,'frequency_caveat':'Control-group frequency compared on reported A1 to 1000 Genomes EUR sample; informative QC only, not proof of strand, build, target population or allele identity.'},indent=2)+'\n');print(json.dumps(results,indent=2))
