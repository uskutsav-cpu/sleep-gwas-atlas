import csv,json
from pathlib import Path
rows=list(csv.DictReader(open('outputs/stage6_gsa_reference_tag_screen.csv')));index={}
for line in open('work/gepitel/GEPITEL.bim'):
 x=line.split();index.setdefault(x[1],[]).append(x)
out=[]
for r in rows:
 if float(r['best_r2'])<0.8:continue
 hit=[x for x in index.get(r['best_tag_id'],[]) if x[0]==r['chr_GRCh37'] and x[3]==r['best_tag_pos']]
 out.append({'candidate':r['variant'],'score_lists':r['score_lists'],'tag':r['best_tag_id'],'reference_r2':float(r['best_r2']),'bim_match_count':len(hit),'bim_alleles':hit[0][4:6] if len(hit)==1 else None})
summary={s:sum(s in x['score_lists'].split(';') and x['bim_match_count']==1 for x in out) for s in ['score15','score22']}
Path('outputs/stage6_gepitel_existing_tag_screen.json').write_text(json.dumps({'summary':summary,'rows':out,'interpretation':'Reference r2 only; no phase/direction, no phenotype, no marker QC or GEPITEL imputation accuracy established.'},indent=2)+'\n');print(summary)
