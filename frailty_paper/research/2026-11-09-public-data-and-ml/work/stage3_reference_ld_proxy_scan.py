import csv,collections,json,os
import numpy as np
from pathlib import Path
m=list(csv.DictReader(open('outputs/shared22_eight_source_evidence_matrix.csv')))
c=[{'id':x['variant'],'chr':x['chr_GRCh37_assumed'],'pos':int(x['position_GRCh37_assumed']),'a':x['score_counted_allele'],'b':x['other_allele']} for x in m]
z=next(x for x in csv.DictReader(open('outputs/integrated_15_candidate_review.csv')) if x['variant']=='rs323509')
c.append({'id':'rs323509','chr':z['GRCh37_chr'],'pos':int(z['GRCh37_pos']),'a':z['counted_allele'],'b':z['other_allele']})
bychr=collections.defaultdict(list)
for x in c: x['window']=[];x['anchors']=[];bychr[x['chr']].append(x)
base=Path('work/stage4_magma_reference/g1000_eur')
with open(str(base)+'.bim') as f:
 for i,line in enumerate(f):
  v=line.split(); chrom,rs,_,pos,a,b=v;pos=int(pos)
  for x in bychr.get(chrom,[]):
   if abs(pos-x['pos'])<=250000:
    entry=(i,rs,pos,a,b)
    if a in 'ACGT' and b in 'ACGT' and len(a)==len(b)==1: x['window'].append(entry)
    if pos==x['pos'] and {a,b}=={x['a'],x['b']}:x['anchors'].append(entry)
n=i+1; nind=503; step=(nind+3)//4
bed=np.memmap(str(base)+'.bed',mode='r',dtype='uint8');assert bytes(bed[:3])==bytes([0x6c,0x1b,0x01]);assert len(bed)==3+n*step
# PLINK SNP-major codes 00=2 copies of BIM A1, 10=1, 11=0, 01=missing.
lookup=np.array([2,np.nan,1,0],dtype=np.float32)
def decode(indices):
 ix=np.asarray(indices,dtype=np.int64); offset=3+ix[:,None]*step+np.arange(step)[None,:]
 b=bed[offset];code=((b[:,:,None] >> np.array([0,2,4,6],dtype=np.uint8)) & 3).reshape(len(ix),-1)[:,:nind];return lookup[code]
rows=[];audit=[]
for x in c:
 item={'candidate':x['id'],'chr_GRCh37':x['chr'],'position_GRCh37':x['pos'],'alleles':x['a']+'/'+x['b'],'anchor_matches':len(x['anchors']),'window_biallelic_records':len(x['window']),'proxy_count_r2_ge_0_8':None,'status':''}
 if len(x['anchors'])!=1:
  item['status']='anchor unavailable or ambiguous by exact coordinate and allele pair';audit.append(item);continue
 anchor=x['anchors'][0];y=decode([anchor[0]])[0];count=0
 for j in range(0,len(x['window']),2048):
  chunk=x['window'][j:j+2048];X=decode([v[0] for v in chunk]);valid=np.isfinite(X)&np.isfinite(y)[None,:];N=valid.sum(axis=1);xf=np.where(valid,X,0);yf=np.where(valid,y[None,:],0);sx=xf.sum(1);sy=yf.sum(1);cov=(xf*yf).sum(1)-sx*sy/np.maximum(N,1);vx=(xf*xf).sum(1)-sx*sx/np.maximum(N,1);vy=(yf*yf).sum(1)-sy*sy/np.maximum(N,1);den=vx*vy;r2=np.divide(cov*cov,den,out=np.zeros_like(cov),where=den>0)
  for k in np.flatnonzero((r2>=.8)&(N>=400)):
   idx,rs,pos,a,b=chunk[int(k)];count+=1
   rows.append({'candidate':x['id'],'anchor_reference_id':anchor[1],'candidate_chr_GRCh37':x['chr'],'candidate_pos_GRCh37':x['pos'],'candidate_bim_A1_A2':anchor[3]+'/'+anchor[4],'proxy_reference_id':rs,'proxy_pos_GRCh37':pos,'proxy_bim_A1_A2':a+'/'+b,'distance_bp':pos-x['pos'],'n_pairwise_nonmissing':int(N[k]),'r2_genotype_dosage':float(r2[k]),'interpretation':'reference LD only; no association, allele orientation or causal claim'})
 item['status']='evaluated exact position and allele pair';item['proxy_count_r2_ge_0_8']=count;audit.append(item);print(x['id'],count,flush=True)
rows.sort(key=lambda r:(r['candidate'],-r['r2_genotype_dosage'],r['proxy_pos_GRCh37']))
with open('outputs/stage3_candidate_reference_ld_proxies.csv','w',newline='') as f:
 w=csv.DictWriter(f,fieldnames=rows[0].keys());w.writeheader();w.writerows(rows)
with open('outputs/stage3_candidate_reference_ld_proxy_audit.json','w') as f:json.dump({'reference':'MAGMA 1000 Genomes EUR GRCh37, 503 samples','window_bp':250000,'r2_threshold':.8,'min_nonmissing':400,'n_bim':n,'method':'pairwise-complete allele dosage Pearson r squared','candidates':audit},f,indent=2)
print('anchors',sum(v['anchor_matches']==1 for v in audit),'rows',len(rows))
