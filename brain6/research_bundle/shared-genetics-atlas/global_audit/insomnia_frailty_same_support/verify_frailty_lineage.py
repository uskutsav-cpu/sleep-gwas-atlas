import csv,gzip,hashlib,json,pathlib,math
p=pathlib.Path(__file__).resolve().parent;base=pathlib.Path('/Volumes/Extreme SSD/Codex Archive/hfrs_fi_ffs_common_support');child=pathlib.Path('/Volumes/Extreme SSD/Codex Archive/hfrs_ukb_sensitivity');root=pathlib.Path('/Volumes/Extreme SSD/sleep-frailty-full-autonomous-v5-20261009/frailty_paper/full_autonomous_research_v5/genetics/new_public_validation');support=json.loads((base/'support.json').read_text());spec=['SNP','A1','A2','N','Z']
def sha(f):
 h=hashlib.sha256()
 with f.open('rb') as s:
  for b in iter(lambda:s.read(1<<20),b''):h.update(b)
 return h.hexdigest()
results={}
for t in ('FI','FFS','HFRS'):
 source=base/f'{t}.sumstats.gz';out=child/f'{t}.sumstats.gz';assert sha(source)==support['output_sha256'][t]
 upstream=root/('munged' if t=='HFRS' else 'common_support')/f'{t}.sumstats.gz';assert sha(upstream)==support['input_sha256'][t]
 with gzip.open(source,'rt') as a,gzip.open(out,'rt') as b:
  ar=csv.DictReader(a,delimiter='\t');br=csv.DictReader(b,delimiter='\t');assert ar.fieldnames==br.fieldnames==spec;want=next(br,None);matched=0;rows=0
  for r in ar:
   rows+=1
   if want and r['SNP']==want['SNP']:
    assert all(r[k]==want[k] for k in ('SNP','A1','A2')) and all(float(r[k])==float(want[k]) for k in ('N','Z')),(t,r,want)
    want=next(br,None);matched+=1
  assert want is None and next(br,None) is None
 assert matched==1134953 and rows==1137927
 results[t]={'parent_path':str(source),'parent_sha256':sha(source),'upstream_path':str(upstream),'upstream_sha256':sha(upstream),'child_path':str(out),'child_sha256':sha(out),'parent_rows':rows,'matched_child_rows':matched,'validation':'ordered exact allele and numeric-value subset; parent and upstream hashes match prior support receipt'}
(p/'frailty_lineage_receipt.json').write_text(json.dumps({'traits':results,'limitations':'Verifies file lineage, not original unharmonized HFRS allele alignment, participant overlap or per-variant sample size.'},indent=2)+'\n');print('verified FI, FFS and HFRS exact subsets and archived upstream file hashes')
