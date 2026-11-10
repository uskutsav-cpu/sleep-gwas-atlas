import pathlib,hashlib,json,tarfile,datetime,glob
root=pathlib.Path('.'); selected=set()
for pattern in ['outputs/stage3*','outputs/stage4*','outputs/shared22_eight_source_evidence_matrix.csv','outputs/shared22_allele_orientation_crosswalk.csv','outputs/integrated_15_candidate_review.csv','outputs/prototype_variants_build_reference.csv','outputs/candidate_trait_direction_alleles.csv','work/stage3*.py','work/stage4*.py','work/stage3*.tsv','work/stage3*.json','work/stage3*aligned.csv']:
 for p in root.glob(pattern):
  if p.is_file() and not p.is_symlink():selected.add(p)
for name in ['work/stage3_remaining_qtl','work/stage3_negr1_qtl','work/stage3_positional_fallback_qtl']:
 d=root/name
 if d.exists():selected.update(p for p in d.rglob('*') if p.is_file() and not p.is_symlink())
def sha(f):
 h=hashlib.sha256()
 with f.open('rb') as x:
  while b:=x.read(1024*1024):h.update(b)
 return h.hexdigest()
files=[{'path':str(p),'bytes':p.stat().st_size,'sha256':sha(p)} for p in sorted(selected)]
manifest={'created_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'scope':'Stage 3 and 4 deliverables, candidate allele crosswalks, accessible molecular inputs and analysis scripts; large external data are referenced by source receipts, not embedded','file_count':len(files),'files':files}
pathlib.Path('outputs/stage5_frozen_input_manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
archive=pathlib.Path('outputs/stage5_frozen_stage3_stage4_inputs.tar.gz')
with tarfile.open(archive,'w:gz') as t:
 for p in sorted(selected):t.add(p,arcname=str(p),recursive=False)
 t.add('outputs/stage5_frozen_input_manifest.json',arcname='stage5_frozen_input_manifest.json')
with tarfile.open(archive,'r:gz') as t:
 members={m.name:m for m in t if m.isfile()};assert len(members)==len(files)+1
 for item in files:
  h=hashlib.sha256();f=t.extractfile(members[item['path']])
  while b:=f.read(1024*1024):h.update(b)
  assert h.hexdigest()==item['sha256'],item['path']
receipt={'archive':str(archive),'archive_bytes':archive.stat().st_size,'archive_sha256':sha(archive),'manifest_sha256':sha(pathlib.Path('outputs/stage5_frozen_input_manifest.json')),'verified_member_count':len(files),'all_member_sha256_verified':True,'large_inputs':'external reference panels and downloaded atlas archives are not embedded; their receipts and source hashes are preserved where available'}
pathlib.Path('outputs/stage5_frozen_input_receipt.json').write_text(json.dumps(receipt,indent=2)+'\n');print(json.dumps(receipt,indent=2))
