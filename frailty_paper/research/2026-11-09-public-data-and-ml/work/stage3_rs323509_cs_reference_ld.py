import csv,json,numpy as np
from bed_reader import open_bed
z=list(csv.DictReader(open('outputs/stage3_rs323509_testis_splice_credible_set.csv'))); assert len(z)==4
matches={}
for s in open('work/stage3_rs323509_cs_bim_matches.txt'):
 i,line=s.split(':',1);v=line.split();matches.setdefault(v[1],[]).append((int(i)-1,v))
idx=[];meta=[]
for x in z:
 q=matches[x['rsid']];assert len(q)==1;i,v=q[0];_,pos,ref,alt=x['variant'].split('_');assert v[0]=='5' and {v[4],v[5]}=={ref,alt};assert int(pos)-int(v[3])==664299
 idx.append(i);meta.append({'rsid':x['rsid'],'GRCh38_pos':pos,'GRCh37_pos':v[3],'alt':alt,'pip':x['pip'],'bim_A1':v[4]})
a=json.load(open('outputs/stage3_candidate_reference_ld_proxy_audit.json'));b=open_bed('work/stage4_magma_reference/g1000_eur.bed',iid_count=503,sid_count=a['n_bim'],count_A1=True,num_threads=1);G=b.read(index=np.s_[:,idx],dtype='float64',num_threads=1)
for i,x in enumerate(meta):
 if x['alt']!=x['bim_A1']:G[:,i]=2-G[:,i]
rows=[]
for i in range(4):
 for j in range(i+1,4):
  k=np.isfinite(G[:,i])&np.isfinite(G[:,j]);r=float(np.corrcoef(G[k,i],G[k,j])[0,1]);rows.append({'variant_1':meta[i]['rsid'],'variant_2':meta[j]['rsid'],'r_alt_dosage':r,'r2_alt_dosage':r*r,'n_pairwise':int(k.sum()),'interpretation':'reference LD, not disease fine mapping or shared-causality proof'})
with open('outputs/stage3_rs323509_molecular_credible_set_reference_ld.csv','w',newline='') as f:w=csv.DictWriter(f,fieldnames=rows[0]);w.writeheader();w.writerows(rows)
print(rows)
