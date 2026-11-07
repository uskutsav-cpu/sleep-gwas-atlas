import csv,gzip,hashlib,json
from pathlib import Path
from collections import Counter
out=Path(__file__).resolve().parent
rank=out.parent/'rs40465_mpra/ranked_screen.tsv'; src=Path('/Volumes/Extreme SSD/brain6-crossdisorder-mpra')
def read(name):
 with gzip.open(src/name,'rt') as f:return list(csv.DictReader(f))
dna=read('dnaMatrix.csv.gz');rna=read('rnaMatrix.csv.gz');exclude=read('excludelist.csv.gz')
a={x['variant']:x for x in dna};b={x['variant']:x for x in rna}; ex={x['variant'] for x in exclude};assert len(a)==len(dna) and len(b)==len(rna)
rows=[]
for r in csv.DictReader(rank.open(),delimiter='\t'):
 if not r.get('grch37'):continue
 key=r['grch37']+'_'+r['grch38'].split(':')[-1].replace('>',':')
 excluded=key.replace(':','_',1) in ex
 rows.append({'rsid':r['rsid'],'region':r['region'],'grch38':r['grch38'],'source_key':key,'dna_present':key in a,'rna_present':key in b,'listed_excluded':excluded})
fields=['rsid','region','grch38','source_key','dna_present','rna_present','listed_excluded']
with (out/'ranked_coverage.tsv').open('w',newline='') as f:
 wr=csv.DictWriter(f,fieldnames=fields,delimiter='\t');wr.writeheader();wr.writerows(rows)
matched=[r for r in rows if r['dna_present'] and r['rna_present']]
assert len(matched)==42
print('matched',len(matched),'excluded',sum(r['listed_excluded'] for r in matched),'rs4916723',next((r for r in rows if r['rsid']=='rs4916723'),None),'region',Counter(r['region'] for r in matched))
h=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
receipt={'ranked_input_sha256':h(rank),'sources':{f.name:h(f) for f in [src/'dnaMatrix.csv.gz',src/'rnaMatrix.csv.gz',src/'excludelist.csv.gz']},'rows_checked':len(rows),'both_count_matrices':len(matched),'among_matches_listed_in_excludelist':sum(r['listed_excluded'] for r in matched),'matches_by_region':dict(Counter(r['region'] for r in matched)),'output_sha256':h(out/'ranked_coverage.tsv'),'note':'presence of a variant in count matrices is not a statement that it passed all publication QC or has a significant allelic effect; do not interpret excludelist alone as QC approval'}
(out/'receipt.json').write_text(json.dumps(receipt,indent=2)+'\n')
