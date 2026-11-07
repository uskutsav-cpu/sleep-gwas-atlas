import csv,hashlib,json,math
from pathlib import Path
p=Path(__file__).resolve().parent
read=lambda path:list(csv.DictReader(path.open(),delimiter='\t'))
x=read(Path('/Volumes/Brain6Work/chr7-three-trait-ukb-sensitivity.tsv'));rec=json.loads((p/'matching_receipt.json').read_text());assert len(x)==1055==rec['intersection'];assert hashlib.sha256(Path(rec['matched_input']).read_bytes()).hexdigest()==rec['matched_input_sha256']
def ls(v):
 m=max(v);return m+math.log(math.fsum(math.exp(i-m) for i in v))
def bf(b,s,sd):
 b=float(b);s=float(s);r=sd**2/(sd**2+s*s);return .5*(math.log1p(-r)+r*(b/s)**2)
err=0
for row in read(p/'pairwise.tsv'):
 pre='fi' if row['source']=='meta' else 'ukb_fi';other=row['pair'].lower();a=[bf(z[pre+'_beta'],z[pre+'_se'],.15) for z in x];b=[bf(z[other+'_beta'],z[other+'_se'],.2) for z in x];sa=ls(a);sb=ls(b);sc=ls([i+j for i,j in zip(a,b)]);distinct=sa+sb+math.log1p(-math.exp(sc-sa-sb));w=[0,math.log(1e-4)+sa,math.log(1e-4)+sb,math.log(1e-8)+distinct,math.log(float(row['p12']))+sc];tot=ls(w);err=max(err,max(abs(math.exp(v-tot)-float(row['H'+str(i)])) for i,v in enumerate(w)))
assert err<1e-10
parent=read(p.parent/'clusters.tsv'); current=read(p/'clusters.tsv');cols=['iteration','traits','posterior_prob','regional_prob','candidate_snp','posterior_explained_by_snp','dropped_trait','setting'];assert [{c:r[c] for c in cols} for r in current if r['source']=='meta']==[{c:r[c] for c in cols} for r in parent]
result={'matched_n':len(x),'checked_pairwise_rows':8,'max_posterior_difference':err,'parent_meta_clusters_reproduced':True,'note':'Reproduces arithmetic and source-support consistency, not model calibration or independent biological evidence.'};(p/'verification.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))
