import csv,gzip,hashlib,json,pathlib
p=pathlib.Path(__file__).resolve().parent;parent=pathlib.Path('/Volumes/Extreme SSD/sleep-frailty-full-autonomous-v5-20261009/frailty_paper/full_autonomous_research_v5/genetics/new_public_validation/common_support/S.sumstats.gz');child=pathlib.Path('/Volumes/Extreme SSD/Codex Archive/hfrs_ukb_sensitivity/S.sumstats.gz')
def sha(f):
 h=hashlib.sha256()
 with f.open('rb') as s:
  for b in iter(lambda:s.read(1<<20),b''):h.update(b)
 return h.hexdigest()
with gzip.open(parent,'rt') as a,gzip.open(child,'rt') as b:
 ar=csv.DictReader(a,delimiter='\t');br=csv.DictReader(b,delimiter='\t');expected=['SNP','A1','A2','N','Z'];assert ar.fieldnames==br.fieldnames==expected
 target=next(br,None);n=0;parent_n=0
 for x in ar:
  parent_n+=1
  if target is not None and x['SNP']==target['SNP']:
   assert all(x[k]==target[k] for k in expected),(n,x,target);n+=1;target=next(br,None)
 assert target is None and next(br,None) is None
assert n==1134953 and parent_n==1138262
known=json.loads((p.parents[1]/'placo_diagnostic/diagnostic_parameters.json').read_text())['input_hashes']['S'];assert sha(parent)==known
receipt={'parent_path':str(parent),'parent_sha256':sha(parent),'parent_matches_frozen_R13_source_hash':True,'child_path':str(child),'child_sha256':sha(child),'parent_rows':parent_n,'child_rows':n,'validation':'Child is an ordered exact subset of the pinned R13 insomnia input, including SNP, alleles, N and rounded Z; both gzip files read to EOF.','limitation':'Pin verifies file lineage, not underlying phenotype definition, participant identities or per-variant N.'};(p/'insomnia_lineage_receipt.json').write_text(json.dumps(receipt,indent=2)+'\n');print('verified',n,'child rows as exact ordered subset of',parent_n,'frozen R13 source rows')
