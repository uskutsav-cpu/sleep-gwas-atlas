import csv,json,numpy as np
from bed_reader import open_bed
matches={}
for s in open('work/stage3_ld_validation_bim_matches.txt'):
 i,line=s.split(':',1);v=line.split();matches[v[1]]=int(i)-1
pairs=[('rs12536335','rs2030915'),('rs323509','rs60271'),('rs7734045','rs58125298'),('rs34290943','rs78807522')]
a=json.load(open('outputs/stage3_candidate_reference_ld_proxy_audit.json'))
b=open_bed('work/stage4_magma_reference/g1000_eur.bed',iid_count=503,sid_count=a['n_bim'],count_A1=True,num_threads=1)
rows=list(csv.DictReader(open('outputs/stage3_candidate_reference_ld_proxies.csv')));out=[]
for x,y in pairs:
 g=b.read(index=np.s_[:,[matches[x],matches[y]]],dtype='float64',num_threads=1);valid=np.isfinite(g).all(axis=1);r=float(np.corrcoef(g[valid,0],g[valid,1])[0,1]**2)
 z=next(q for q in rows if q['candidate']==x and q['proxy_reference_id']==y);v=float(z['r2_genotype_dosage']);out.append({'candidate':x,'proxy':y,'independent_bed_reader_r2':r,'scan_r2':v,'difference':abs(r-v),'n':int(valid.sum())});assert abs(r-v)<1e-5
json.dump({'method':'Independent bed-reader package genotype decoding and numpy correlation; four selected pairs','pairs':out},open('outputs/stage3_reference_ld_validation.json','w'),indent=2)
print(out)
