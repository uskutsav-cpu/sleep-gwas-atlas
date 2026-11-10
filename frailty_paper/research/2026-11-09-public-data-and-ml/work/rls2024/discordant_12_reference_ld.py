import csv,json,hashlib,itertools
from pathlib import Path
import numpy as np,pysam
from pyliftover import LiftOver
p=Path('outputs/stage6_rls_three_study_region.csv');rows=[r for r in csv.DictReader(p.open()) if r['first_two_same_sign']=='0'];assert len(rows)==12
panel=list(csv.DictReader(open('work/1kg_phase3_panel.tsv'),delimiter='\t'));pops={'EUR':sorted(r['sample'] for r in panel if r['super_pop']=='EUR'),'FIN':sorted(r['sample'] for r in panel if r['pop']=='FIN')};assert len(pops['EUR'])==503 and len(pops['FIN'])==99
lo=LiftOver('work/hg38ToHg19.over.chain.gz');url='https://hgdownload.soe.ucsc.edu/gbdb/hg19/1000Genomes/phase3/ALL.chr2.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz';vcf=pysam.VariantFile(url);store={};inventory=[]
for row in rows:
 rs=row['akcimen_rsid'];pos=int(row['pos']);results=lo.convert_coordinate('chr2',pos-1);same=[v for v in results if v[0]=='chr2' and v[2]=='+'];item={'rsid':rs,'pos38':pos,'ref38':row['ref'],'alt38':row['alt'],'liftover_hits':len(results),'status':''}
 if len(results)!=1 or len(same)!=1:item['status']='ambiguous liftover';inventory.append(item);continue
 pos37=same[0][1]+1;item['pos37']=pos37
 matches=[v for v in vcf.fetch('2',pos37-1,pos37) if v.pos==pos37 and v.id==rs and len(v.alleles)==2 and set(v.alleles)=={row['ref'],row['alt']}]
 if len(matches)!=1:item['status']='no unique biallelic reference match';inventory.append(item);continue
 v=matches[0];a=v.alleles.index(row['alt']);by={}
 for pop,ids in pops.items():
  d=[]
  for i in ids:
   g=v.samples[i].get('GT');d.append(sum(q==a for q in g) if g and len(g)==2 and all(q is not None for q in g) else np.nan)
  by[pop]=np.array(d,float);called=np.isfinite(by[pop]);item[pop+'_n_called']=int(called.sum());item[pop+'_alt_count']=int(np.nansum(by[pop]))
 store[rs]=by;item['status']='matched';inventory.append(item)
pairs=[]
for pop in pops:
 for a,b in itertools.combinations(store,2):
  x=store[a][pop];y=store[b][pop];keep=np.isfinite(x)&np.isfinite(y)
  if keep.sum()<80 or np.std(x[keep])==0 or np.std(y[keep])==0:continue
  corr=float(np.corrcoef(x[keep],y[keep])[0,1]);pairs.append({'reference':pop,'variant_1':a,'variant_2':b,'n':int(keep.sum()),'r':corr,'r2':corr*corr})
summary={}
for pop in pops:
 pair=[x for x in pairs if x['reference']==pop];m=[x for x in inventory if x['status']=='matched'];ids=[x['rsid'] for x in m if x.get(pop+'_alt_count',0)>0 and x.get(pop+'_alt_count',0)<2*x.get(pop+'_n_called',0)];adj={i:set() for i in ids}
 for x in pair:
  if x['r2']>=.8 and x['variant_1'] in adj and x['variant_2'] in adj:adj[x['variant_1']].add(x['variant_2']);adj[x['variant_2']].add(x['variant_1'])
 groups=[];unseen=set(ids)
 while unseen:
  todo=[unseen.pop()];g=[]
  while todo:
   v=todo.pop();g.append(v)
   for u in adj[v]&unseen:unseen.remove(u);todo.append(u)
  groups.append(sorted(g))
 summary[pop]={'matched_variable':len(ids),'pairs_evaluable':len(pair),'pairs_r2_ge_0_8':sum(x['r2']>=.8 for x in pair),'maximum_r2':max((x['r2'] for x in pair),default=None),'connected_groups_r2_ge_0_8':sorted(groups)}
out=Path('outputs/stage6_rls_discordant12_reference_ld.csv')
with out.open('w',newline='') as f:w=csv.DictWriter(f,fieldnames=['reference','variant_1','variant_2','n','r','r2']);w.writeheader();w.writerows(pairs)
receipt={'input_sha256':hashlib.sha256(p.read_bytes()).hexdigest(),'chain_sha256':hashlib.sha256(Path('work/hg38ToHg19.over.chain.gz').read_bytes()).hexdigest(),'panel_sha256':hashlib.sha256(Path('work/1kg_phase3_panel.tsv').read_bytes()).hexdigest(),'source':url,'variant_inventory':inventory,'summary':summary,'pairs_output_sha256':hashlib.sha256(out.read_bytes()).hexdigest(),'limitation':'Postselected variants; EUR503 includes FIN99, potentially related and structured samples; external LD not study data; chain mapping and exact alleles do not resolve all multiallelic study semantics; rare variant estimates can be unstable; r2 threshold grouping is descriptive and not conditional-analysis independence.'};Path('outputs/stage6_rls_discordant12_reference_ld_receipt.json').write_text(json.dumps(receipt,indent=2)+'\n');print(json.dumps({'inventory':inventory,'summary':summary},indent=2))
