import csv,pathlib,hashlib,json,collections,math
p=pathlib.Path(__file__).resolve().parent;root=pathlib.Path('/Volumes/Extreme SSD/sleep-frailty-full-autonomous-v5-20261009/frailty_paper/full_autonomous_research_v5');raw=root/'genetics/new_public_validation/raw/GCST90707291.tsv';filtered=root/'genetics/new_public_validation/filtered/HFRS_HM3.tsv';freeze=root/'research_completion_v2/2026-10-09_genetics_update/manifests/source_freeze/SOURCE_FREEZE_METADATA.json';manifest=root/'research_completion_v2/2026-10-09_genetics_update/manifests/NATIVE_VALIDATION_INPUT_HASHES.json';metadata=json.loads(freeze.read_text())['HFRS'];hashes=json.loads(manifest.read_text())['inputs'];comp=str.maketrans('ACGT','TGCA')
def sha(f):
 h=hashlib.sha256()
 with f.open('rb') as z:
  for b in iter(lambda:z.read(1<<20),b''):h.update(b)
 return h.hexdigest()
assert sha(filtered)==hashes[str(filtered)]
filtered_rows={};duplicate_filtered=set();counts=collections.Counter()
with filtered.open() as f:
 for r in csv.DictReader(f,delimiter='\t'):
  rs=r['SNP'].encode();counts['filtered_rows']+=1
  if rs in filtered_rows:duplicate_filtered.add(rs)
  else:filtered_rows[rs]=(r['A1'],r['A2'],r['BETA'],r['SE'],r['P'],r['N'],r['EAF'],r['CHR'])
assert not duplicate_filtered and len(filtered_rows)==1161680
matches=collections.Counter();good=collections.Counter();bad_examples=[];md5=hashlib.md5();sha_raw=hashlib.sha256();bytes_read=0
with raw.open('rb') as f:
 hdr=f.readline();md5.update(hdr);sha_raw.update(hdr);bytes_read+=len(hdr);assert hdr.decode().strip().split('\t')==['chromosome','base_pair_location','effect_allele','other_allele','beta','standard_error','effect_allele_frequency','p_value','rs_id']
 for line in f:
  md5.update(line);sha_raw.update(line);bytes_read+=len(line);counts['raw_rows']+=1;cols=line.rstrip(b'\r\n').split(b'\t');assert len(cols)==9
  rs=cols[8];target=filtered_rows.get(rs)
  if target is None:continue
  counts['id_rows']+=1;a,b,beta,se,pval,n,eaf,ch=target
  v=[x.decode() for x in cols];chrom,pos,x,y,rawb,rawse,rawf,rawp,_=v
  if chrom!=ch:counts['chrom_mismatch_rows']+=1;continue
  if any(len(q)!=1 or q not in 'ACGT' for q in (a,b,x,y)) or ''.join(sorted(a+b)) in ('AT','CG'):counts['invalid_allele_rows']+=1;continue
  if (x,y)==(a,b):mode='direct';flip=1
  elif (x.translate(comp),y.translate(comp))==(a,b):mode='complement';flip=1
  elif (y,x)==(a,b):mode='swapped';flip=-1
  elif (y.translate(comp),x.translate(comp))==(a,b):mode='complement_swapped';flip=-1
  else:counts['allele_mismatch_rows']+=1;continue
  matches[rs]+=1
  try:
   expected_b=float(rawb)*flip;freq=float(rawf) if flip==1 else 1-float(rawf)
   ok=all(math.isclose(float(u),float(v),rel_tol=1e-10,abs_tol=1e-12) for u,v in [(beta,expected_b),(se,rawse),(pval,rawp),(eaf,freq)]) and float(n)==500737
  except ValueError:ok=False
  if ok:good[rs]+=1;counts[mode]+=1
  else:
   counts['value_mismatch_rows']+=1
   if len(bad_examples)<10:bad_examples.append({'rsid':rs.decode(),'filtered':target,'raw':v,'mode':mode})
assert bytes_read==metadata['bytes'] and md5.hexdigest()==metadata['md5'],(bytes_read,md5.hexdigest(),metadata)
counts['filtered_no_matching_source']=sum(matches[x]==0 for x in filtered_rows);counts['filtered_multiple_matching_source']=sum(matches[x]>1 for x in filtered_rows);counts['filtered_exact_unique_and_values_agree']=sum(matches[x]==1 and good[x]==1 for x in filtered_rows)
missing=[rs for rs in filtered_rows if matches[rs]==0]; selected=set()
with __import__('gzip').open('/Volumes/Extreme SSD/Codex Archive/hfrs_ukb_sensitivity/HFRS.sumstats.gz','rt') as z:
 for r in csv.DictReader(z,delimiter='\t'):
  if r['SNP'].encode() in missing:selected.add(r['SNP'])
result={'missing_filtered_rows':[{'rsid':rs.decode(),'in_final_support':rs.decode() in selected,'filtered':filtered_rows[rs]} for rs in missing],'raw_path':str(raw),'raw_bytes':bytes_read,'raw_md5':md5.hexdigest(),'publisher_md5':metadata['md5'],'raw_sha256':sha_raw.hexdigest(),'filtered_path':str(filtered),'filtered_sha256':sha(filtered),'counts':dict(counts),'value_mismatch_examples':bad_examples,'note':'Source coordinates GRCh38; filtered table omits position. chr/rsID/allele check is not liftover verification; total N used, not per-variant N.'}
(p/'hfrs_raw_alignment_receipt.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result['counts'],indent=2),flush=True)
