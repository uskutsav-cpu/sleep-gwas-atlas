import csv,zipfile,io,json,collections,hashlib
import numpy as np
from pathlib import Path
P=Path('outputs'); groups={}
for label,fn in [('score15','unvalidated_15_variant_score_prototype.csv'),('score22','unvalidated_22_shared_direction_score.csv')]:
 for r in csv.DictReader((P/fn).open()):
  c=r.get('chr',r.get('chr_GRCh37_assumed')); pos=int(r.get('position_GRCh37_inferred',r.get('position_GRCh37_assumed'))); alleles=frozenset([r['counted_allele'],r['other_allele']]); k=(r['variant'],c,pos,alleles)
  groups.setdefault(k,[]).append(label)
bychr=collections.defaultdict(list)
for (rs,c,pos,pair),labels in groups.items():bychr[c].append({'id':rs,'chr':c,'pos':pos,'pair':pair,'lists':labels,'anchors':[],'tags':[]})
comp=str.maketrans('ACGT','TGCA'); manifest=collections.defaultdict(list)
with open('work/gepitel/GEPITEL.bim') as h:
 for line in h:
  c,rs,_,bp,a,b=line.split()
  if c not in bychr or not bp.isdigit() or not any(abs(int(bp)-x['pos'])<=100000 for x in bychr[c]): continue
  if len(a)!=1 or len(b)!=1 or a not in 'ACGT' or b not in 'ACGT' or a==b: continue
  manifest[(c,int(bp),frozenset((a,b)))].append(rs)
base=Path('work/stage4_magma_reference/g1000_eur');n=0
with open(str(base)+'.bim') as f:
 for n,line in enumerate(f,1):
  c,rs,_,bp,a,b=line.split();bp=int(bp);pair=frozenset((a,b)); candidates=bychr.get(c,[])
  if not candidates or len(a)!=1 or len(b)!=1 or a not in 'ACGT' or b not in 'ACGT' or len(pair)!=2:continue
  names=manifest.get((c,bp,pair))
  for x in candidates:
   if bp==x['pos'] and pair==x['pair']:x['anchors'].append((n-1,rs))
   if names and abs(bp-x['pos'])<=100000 and bp!=x['pos']:x['tags'].append((n-1,rs,bp,a,b,names))
ni=sum(1 for _ in open(str(base)+'.fam'));step=(ni+3)//4; bed=np.memmap(str(base)+'.bed',mode='r',dtype='uint8');assert bytes(bed[:3])==b'\x6c\x1b\x01' and len(bed)==3+n*step
lut=np.array([2,np.nan,1,0],dtype=np.float32)
def decode(ids):
 offsets=3+np.asarray(ids,dtype=np.int64)[:,None]*step+np.arange(step)[None,:]; b=bed[offsets]; codes=((b[:,:,None]>>np.array([0,2,4,6],dtype=np.uint8))&3).reshape(len(ids),-1)[:,:ni];return lut[codes]
rows=[]
for entries in bychr.values():
 for x in entries:
  out={'variant':x['id'],'score_lists':';'.join(x['lists']),'chr_GRCh37':x['chr'],'position_GRCh37':x['pos'],'anchor_matches':len(x['anchors']),'eligible_tag_records':len(x['tags']),'best_tag_id':'','best_manifest_probe':'','best_tag_pos':'','best_r2':'','pairwise_n':'','tags_r2_ge_0_8':'','status':'anchor missing or ambiguous'}
  if len(x['anchors'])==1 and x['tags']:
   y=decode([x['anchors'][0][0]])[0];X=decode([v[0] for v in x['tags']]); ok=np.isfinite(X)&np.isfinite(y)[None,:];N=ok.sum(1);xx=np.where(ok,X,0);yy=np.where(ok,y[None,:],0); sx=xx.sum(1);sy=yy.sum(1); cov=(xx*yy).sum(1)-sx*sy/np.maximum(N,1);vx=(xx*xx).sum(1)-sx*sx/np.maximum(N,1);vy=(yy*yy).sum(1)-sy*sy/np.maximum(N,1); den=vx*vy;r2=np.divide(cov*cov,den,out=np.full_like(cov,np.nan),where=(den>0)&(N>=400));
   if np.isfinite(r2).any():
    j=int(np.nanargmax(r2));v=x['tags'][j];out.update(best_tag_id=v[1],best_manifest_probe=';'.join(v[5]),best_tag_pos=v[2],best_r2=round(float(r2[j]),6),pairwise_n=int(N[j]),tags_r2_ge_0_8=int(np.sum(r2>=.8)),status='evaluated')
   else:out['status']='no tag with variable genotypes and at least 400 paired calls'
  elif len(x['anchors'])==1:out['status']='no matching manifest/reference tag in window'
  rows.append(out)
with (P/'stage6_gepitel_full_reference_tags.csv').open('w',newline='') as h:w=csv.DictWriter(h,fieldnames=rows[0]);w.writeheader();w.writerows(rows)
receipt={'reference':'MAGMA 1000 Genomes EUR GRCh37, 503 samples','manifest':'GEPITEL released BIM, build and strand unverified; matching assumes GRCh37' ,'window_bp':100000,'matching':'exact candidate assumed GRCh37 position and unordered allele pair; GEPITEL build and strand unverified','method':'pairwise-complete unphased dosage Pearson r squared, >=400 pairs','n_unique_candidates':len(rows),'status_counts':dict(collections.Counter(r['status'] for r in rows)),'n_r2_ge_0_8':sum(r['best_r2']!='' and r['best_r2']>=.8 for r in rows),'manifest_sha256':hashlib.sha256(Path('work/gepitel/GEPITEL.bim').read_bytes()).hexdigest(),'caveat':'Reference LD is neither cohort imputation quality nor evidence of phenotype association or independent validation; genotype phenotype linkage and array coding unknown.'}
(P/'stage6_gepitel_full_reference_tags.json').write_text(json.dumps(receipt,indent=2)+'\n');print(json.dumps(receipt,indent=2))
