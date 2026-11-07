import csv,hashlib,json,math,pathlib
p=pathlib.Path(__file__).resolve().parent;src=pathlib.Path('/Volumes/Brain6Work/chr7-three-trait-input.tsv')
def read(path):
 with path.open(newline='') as h:return list(csv.DictReader(h,delimiter='\t'))
x=read(src);rec=json.loads((p/'matching_receipt.json').read_text());assert len(x)==rec['matched']==1055 and hashlib.sha256(src.read_bytes()).hexdigest()==rec['external_analysis_input_sha256'] and len({r['rsid'] for r in x})==len(x)
def logsum(xs):
 m=max(xs);return m+math.log(math.fsum(math.exp(v-m) for v in xs))
def abf(b,se,sd):
 b=float(b);se=float(se);r=sd**2/(sd**2+se**2);return .5*(math.log1p(-r)+r*(b/se)**2)
traits={name:[abf(r[name+'_beta'],r[name+'_se'],sd) for r in x] for name,sd in [('fi',.15),('insomnia',.2),('adhd',.2)]};pairs={'FI_insomnia':('fi','insomnia'),'FI_ADHD':('fi','adhd'),'insomnia_ADHD':('insomnia','adhd')};err=0
for row in read(p/'pairwise.tsv'):
 a,b=[traits[n] for n in pairs[row['pair']]];sa=logsum(a);sb=logsum(b);sc=logsum([i+j for i,j in zip(a,b)]);distinct=sa+sb+math.log1p(-math.exp(sc-sa-sb));w=[0,math.log(1e-4)+sa,math.log(1e-4)+sb,math.log(1e-8)+distinct,math.log(float(row['p12']))+sc];total=logsum(w);err=max(err,max(abs(math.exp(v-total)-float(row['H'+str(i)])) for i,v in enumerate(w)))
assert err<1e-10
out={'common_count':len(x),'pairwise_rows_checked':6,'max_absolute_posterior_difference':err,'external_input_sha256':rec['external_analysis_input_sha256'],'note':'Arithmetic and input checks only, not statistical calibration or replication.'};(p/'verification.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out,indent=2))
