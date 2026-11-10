"""Screen selected Estonian insomnia leads within 1 Mb of mapped project markers."""
import csv,json
from pathlib import Path
import numpy as np, openpyxl, pysam
root=Path(__file__).resolve().parents[1]
def read(name):
 with (root/name).open() as h:return list(csv.DictReader(h))
a={r['marker']:(int(r['chr38']),int(r['pos38'])) for r in read('outputs/stage5_marker_grch37_grch38_liftover_audit.csv') if r['status']=='unique_expected_positive_strand'}
for r in read('outputs/stage6_tpmi_sleep_candidate_inventory.csv'):
 if r['ensembl_liftover_check']=='unique_positive_match':a[r['variant']]=(int(float(r['chr38'])),int(float(r['pos38'])))
for rs in ['rs113851554','rs78206187']:
 m=json.loads((root/f'work/stage6_ensembl_{rs}_mapping.json').read_text())['mappings'][0]['mapped'];a[rs]=(int(m['seq_region_name']),m['start'])
candidates={r['variant']:r for f in ['unvalidated_15_variant_score_prototype.csv','unvalidated_22_shared_direction_score.csv'] for r in read('outputs/'+f)}
w=openpyxl.load_workbook(root/'work/sleep_preprint_supplement.xlsx',read_only=True,data_only=True);rows=[r for r in list(w['ST5'].values)[2:] if len(r)>=15 and r[0]=='Insomnia']
manual={'rs113851554':66750564,'rs1989903':114137940,'rs1093463':117744762,'rs745631':112945892}
panel={r['sample']:r for r in csv.DictReader((root/'work/1kg_phase3_panel.tsv').open(),delimiter='\t')}; cache={};result=[]
for rs,v in candidates.items():
 chr38,pos38=a[rs]
 for row in rows:
  if row[1]!=chr38 or abs(row[2]-pos38)>1e6 or row[3]==rs:continue
  lead=row[3]; lead38=int(row[2]);f=root/f'work/stage6_ensembl_grch38_{chr38}_{lead38}.json'
  if lead in manual:lead37=manual[lead]
  elif f.exists():
   maps=json.loads(f.read_text())['mappings'];assert len(maps)==1 and maps[0]['mapped']['strand']==1;lead37=maps[0]['mapped']['start']
  else:result.append({'candidate':rs,'lead':lead,'error':'lead mapping unavailable'});continue
  cpos=int(v.get('position_GRCh37_inferred') or v.get('position_GRCh37_assumed')); chrom=str(chr38)
  if chrom not in cache:
   url=f'https://hgdownload.soe.ucsc.edu/gbdb/hg19/1000Genomes/phase3/ALL.chr{chrom}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz'; t=pysam.TabixFile(url,index=str(root/f'work/stage6_1000g_chr{chrom}_phase3.tbi'));samples=next(x for x in t.header if x.startswith('#CHROM')).split('\t')[9:];idx=[i for i,s in enumerate(samples) if s in panel and panel[s]['super_pop']=='EUR'];assert len(idx)==503;cache[chrom]=(t,idx)
  t,idx=cache[chrom]; vectors=[]; alleles=[];error=None
  for name,pos in [(rs,cpos),(lead,lead37)]:
   hits=[x.split('\t') for x in t.fetch(chrom,pos-1,pos) if x.split('\t')[2]==name]
   if len(hits)!=1 or ',' in hits[0][4]:error=f'{name}: no unique biallelic VCF record';break
   h=hits[0]; alleles.append([h[3],h[4]]); vector=[]
   for i in idx:
    gt=h[9+i].split(':')[0].replace('|','/').split('/');vector.append(sum(int(z)==1 for z in gt) if len(gt)==2 and all(z!='.' for z in gt) else np.nan)
   vectors.append(np.array(vector))
  item={'candidate':rs,'lead':lead,'distance_bp_grch38':int(row[2])-pos38,'lead_p':row[14],'alleles':alleles,'error':error}
  if not error:
   mask=np.isfinite(vectors[0])&np.isfinite(vectors[1]);r=float(np.corrcoef(vectors[0][mask],vectors[1][mask])[0,1]);item.update(n=int(mask.sum()),r=r,r_squared=r*r)
  result.append(item)
(root/'outputs/stage6_estonian_all_nearby_ld.json').write_text(json.dumps(result,indent=2)+'\n')
for x in sorted(result,key=lambda x:x.get('r_squared',-1),reverse=True):print(x['candidate'],x['lead'],'r2',round(x.get('r_squared',-1),3),'p',x.get('lead_p'),'error',x['error'])
