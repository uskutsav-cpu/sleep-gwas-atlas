import collections,csv,hashlib,json,math
from pathlib import Path
p=Path(__file__).resolve().parent
with (p/'rs4916723.tsv').open() as f:rows=list(csv.DictReader(f,delimiter='\t'))
assert len(rows)==16 and len({(r['xQTLtype'],r['context'],r['study'],r['target']) for r in rows})==16
for r in rows:
 assert r['rsID']=='rs4916723' and r['variantId']=='chr5:88558577:A:C' and float(r['se'])>0 and 0<float(r['FDR'])<.05
 assert math.isclose(float(r['beta'])/float(r['se']),float(r['Z']),rel_tol=1e-6,abs_tol=1e-6)
for snp in ('rs2431108','rs77960'):
 text=(p/(snp+'_no_records.html')).read_text();assert 'No records found for' in text and snp in text
for gene,typ,ctx in [('TMEM161B-AS1','snuc-eQTL','opc'),('LINC02060','snuc-eQTL','ast')]:assert any(r['targetGene']==gene and r['xQTLtype']==typ and r['context']==ctx for r in rows)
h=lambda f:hashlib.sha256(f.read_bytes()).hexdigest()
out={'rs4916723_records':len(rows),'counts_by_type':dict(collections.Counter(r['xQTLtype'] for r in rows)),'GRCh38_variant':'chr5:88558577:A:C','effect_allele':'ALT C per portal convention','other_variants':'no records returned in HMT-filtered lookup; not proof of no association or no test','sha256':{f.name:h(f) for f in p.iterdir() if f.suffix in ('.tsv','.html')}}
(p/'verification.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out,indent=2))
