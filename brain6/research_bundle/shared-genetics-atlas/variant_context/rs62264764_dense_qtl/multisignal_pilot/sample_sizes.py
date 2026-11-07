import pathlib,csv,statistics,json,collections,hashlib
p=pathlib.Path(__file__).resolve().parent;parent=p.parent;v=list(csv.DictReader((p/'variants.tsv').open(),delimiter='\t'));base=pathlib.Path('/Volumes/Extreme SSD/brain6-work/exploratory-finemap-coloc-v2/inputs');raw={}
for t in ('insomnia','adhd'):
 f=base/('placo_insomnia__adhd_chr3_117142005_118142005__'+t+'.tsv');d=list(csv.DictReader(f.open(),delimiter='\t'));raw[t]={x['SNP']:x for x in d};assert len(raw[t])==len(d)
q=list(csv.DictReader((parent/'Adrenal_Gland_ENSG00000185565.tsv').open(),delimiter='\t'));q={x['variant']:x for x in q}
out=[]
for x in v:
 n1=float(raw['insomnia'][x['rsid']]['N_AUDIT']);n2=float(raw['adhd'][x['rsid']]['N_AUDIT']);nq=float(q[x['variant']]['an'])/2;assert min(n1,n2,nq)>2;out.append({'variant':x['variant'],'qtl_N_from_an':nq,'insomnia_N_AUDIT':n1,'adhd_N_AUDIT':n2})
with (p/'sample_sizes.tsv').open('w',newline='') as f:w=csv.DictWriter(f,list(out[0]),delimiter='\t');w.writeheader();w.writerows(out)
s={}
for k in ('qtl_N_from_an','insomnia_N_AUDIT','adhd_N_AUDIT'):
 a=[x[k] for x in out];s[k]={'min':min(a),'median':statistics.median(a),'max':max(a),'distinct':len(set(a)),'below_max':sum(z<max(a) for z in a)}
s['output_sha256']=hashlib.sha256((p/'sample_sizes.tsv').read_bytes()).hexdigest();(p/'sample_sizes_receipt.json').write_text(json.dumps(s,indent=2)+'\n');print(json.dumps(s,indent=2))
