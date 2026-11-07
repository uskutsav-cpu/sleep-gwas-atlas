import csv,hashlib,json,math,pathlib,collections
p=pathlib.Path(__file__).resolve().parent
def read(q):
 with q.open(newline='') as f:return list(csv.DictReader(f,delimiter='\t'))
a=read(p/'rs62264764.tsv');b=read(p/'rs62264767.tsv');assert len(a)==len(b)==9
keys=lambda r:(r['xQTLtype'],r['context'],r['study'],r['target']);da={keys(r):r for r in a};db={keys(r):r for r in b};assert len(da)==len(db)==9 and set(da)==set(db)
for x,rs,variant in [(a,'rs62264764','chr3:117920728:G:A'),(b,'rs62264767','chr3:117923158:A:C')]:
 for r in x:
  assert r['rsID']==rs and r['variantId']==variant and math.isfinite(float(r['Z'])) and 0<float(r['FDR'])<.05 and float(r['se'])>0
  assert abs(float(r['beta'])/float(r['se'])-float(r['Z']))<1e-5
for k in da:
 for field in ('pvalue','FDR','Z','beta','se'):assert math.isclose(float(da[k][field]),float(db[k][field]),rel_tol=1e-9,abs_tol=1e-12)
gene=json.loads((p/'ensembl_gene.json').read_text());assert gene['id']=='ENSG00000239268' and gene['display_name']=='LINC03051'
result={'counts':dict(collections.Counter(r['xQTLtype'] for r in a)),'same_association_targets_and_statistics':len(da),'variant_effect_alleles':'ALT A for rs62264764; ALT C for rs62264767 per portal convention','ensembl_gene_id':gene['id'],'ensembl_display_name':gene['display_name'],'ensembl_version':gene.get('version'),'sha256':{q.name:hashlib.sha256(q.read_bytes()).hexdigest() for q in (p/'rs62264764.tsv',p/'rs62264767.tsv',p/'ensembl_gene.json')},'notes':'All returned records are selected HMT-significant results, not all tested variants; matching statistics do not prove causal variant or independent replication.'};(p/'verification.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))
