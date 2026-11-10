import csv,collections,json
import numpy as np
from pathlib import Path
base=Path('work/stage4_magma_reference/g1000_eur')
fam=[line.split()[1] for line in open(str(base)+'.fam')]
panel={r['sample']:r['pop'] for r in csv.DictReader(open('work/1kg_phase3_panel.tsv'),delimiter='\t')}
groups={k:np.array([i for i,s in enumerate(fam) if panel.get(s)==k]) for k in ('GBR','FIN')}
if tuple(len(groups[k]) for k in ('GBR','FIN'))!=(91,99):raise ValueError('Unexpected sample groups')
targets=list(csv.DictReader(open('outputs/stage6_nshap_gsa_v1_reference_tags.csv')))
need={(r['chr_GRCh37'],int(r[k])) for r in targets for k in ('position_GRCh37','best_tag_pos')}
bim=collections.defaultdict(list)
for i,line in enumerate(open(str(base)+'.bim')):
 c,rs,_,pos,a,b=line.split();key=(c,int(pos))
 if key in need:bim[key].append((i,rs,{a,b}))
ni=len(fam);step=(ni+3)//4;bed=np.memmap(str(base)+'.bed',mode='r',dtype='uint8')
if bytes(bed[:3])!=b'\x6c\x1b\x01':raise ValueError('BED format')
lut=np.array([2,np.nan,1,0],dtype=float)
def decode(i):
 x=bed[3+i*step:3+(i+1)*step];return lut[((x[:,None]>>np.array([0,2,4,6],dtype=np.uint8))&3).reshape(-1)[:ni]]
score={}
for fn in ('outputs/unvalidated_15_variant_score_prototype.csv','outputs/unvalidated_22_shared_direction_score.csv'):
 for r in csv.DictReader(open(fn)):score[r['variant']]={r['counted_allele'],r['other_allele']}
rows=[]
for r in targets:
 c=r['chr_GRCh37'];p=int(r['position_GRCh37']);t=int(r['best_tag_pos']);
 anchor=[v for v in bim[(c,p)] if v[2]==score[r['variant']]];tag=[v for v in bim[(c,t)] if v[1]==r['best_tag_id']]
 if len(anchor)!=1 or len(tag)!=1:raise ValueError('ambiguous marker '+r['variant'])
 a,b=decode(anchor[0][0]),decode(tag[0][0]);row={'variant':r['variant'],'score_lists':r['score_lists'],'tag':r['best_tag_id'],'eur_selected_r2':r['best_r2']}
 for name,idx in groups.items():
  x,y=a[idx],b[idx];m=np.isfinite(x)&np.isfinite(y);x,y=x[m],y[m]; row[name+'_n']=len(x)
  row[name+'_r2']=round(float(np.corrcoef(x,y)[0,1]**2),6) if len(x)>=80 and np.std(x)>0 and np.std(y)>0 else ''
 rows.append(row)
with open('outputs/stage6_nshap_selected_tag_eur_subgroups.csv','w',newline='') as f:
 w=csv.DictWriter(f,fieldnames=rows[0]);w.writeheader();w.writerows(rows)
print('rows',len(rows))
for name in groups:
 x=[r for r in rows if r[name+'_r2']!=''];print(name,'evaluable',len(x),'r2 >= 0.8',sum(r[name+'_r2']>=.8 for r in x))
