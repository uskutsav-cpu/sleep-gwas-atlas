import collections,csv,gzip,hashlib,json,pathlib
from pyliftover import LiftOver
p=pathlib.Path(__file__).resolve().parent;final=pathlib.Path('/Volumes/Extreme SSD/Codex Archive/hfrs_ukb_sensitivity/HFRS.sumstats.gz');ref=pathlib.Path('/Volumes/Extreme SSD/sleep-gwas-atlas-frailty-v1/analysis-workspace/reference/eur_w_ld_chr');chain=pathlib.Path('/Users/swethasunilkumar/Documents/Codex/2026-10-27/loo/work/hg19ToHg38.over.chain.gz');source=pathlib.Path('/Volumes/Extreme SSD/sleep-frailty-full-autonomous-v5-20261009/frailty_paper/full_autonomous_research_v5/genetics/new_public_validation/raw/GCST90707291.tsv');prior=json.loads((p/'hfrs_raw_alignment_receipt.json').read_text());comp=str.maketrans('ACGT','TGCA')
def sha(path):
 h=hashlib.sha256()
 with path.open('rb') as f:
  for b in iter(lambda:f.read(1<<20),b''):h.update(b)
 return h.hexdigest()
selected={}
with gzip.open(final,'rt') as f:
 for x in csv.DictReader(f,delimiter='\t'):
  k=x['SNP'].encode();assert k not in selected;selected[k]=(x['A1'],x['A2'],x['Z'])
assert len(selected)==1134953
reference={};duplicate=set()
for c in range(1,23):
 with gzip.open(ref/f'{c}.l2.ldscore.gz','rt') as f:
  for x in csv.DictReader(f,delimiter='\t'):
   k=x['SNP'].encode()
   if k not in selected:continue
   if k in reference:duplicate.add(k)
   else:reference[k]=(x['CHR'],int(x['BP']))
lift=LiftOver(str(chain));targets={};counts=collections.Counter();examples=collections.defaultdict(list)
def example(kind,rs,detail):
 if len(examples[kind])<10:examples[kind].append({'rsid':rs.decode(),'details':detail})
unresolved=[]
for rs,(a,b,z) in selected.items():
 if rs not in reference:counts['reference_missing']+=1;example('reference_missing',rs,None);unresolved.append({'rsid':rs.decode(),'status':'reference_missing','z':z});continue
 if rs in duplicate:counts['reference_duplicate']+=1;example('reference_duplicate',rs,None);unresolved.append({'rsid':rs.decode(),'status':'reference_duplicate','z':z});continue
 chrom,pos=reference[rs];maps=lift.convert_coordinate('chr'+chrom,pos-1)
 if len(maps)!=1 or maps[0][2] not in ('+','-') or not maps[0][0].startswith('chr'):
  counts['lift_unmapped_or_ambiguous']+=1;example('lift_unmapped_or_ambiguous',rs,{'chr37':chrom,'pos37':pos,'candidates':maps});unresolved.append({'rsid':rs.decode(),'status':'lift_unmapped_or_ambiguous','z':z,'chr37':chrom,'pos37':pos});continue
 chrom38=maps[0][0][3:];pos38=maps[0][1]+1;strand=maps[0][2]
 if strand=='-':a=a.translate(comp);b=b.translate(comp)
 targets[rs]=(chrom38,pos38,frozenset((a,b)),chrom,pos,strand)
counts['unique_lift']=len(targets);hits=collections.Counter();rsrows=collections.Counter();coords=collections.Counter();badalleles=collections.Counter();other_coords=collections.defaultdict(list);md5=hashlib.md5();sha_src=hashlib.sha256();raw_rows=0;size=0
with source.open('rb') as f:
 head=f.readline();md5.update(head);sha_src.update(head);size+=len(head);assert head.decode().strip().split('\t')==['chromosome','base_pair_location','effect_allele','other_allele','beta','standard_error','effect_allele_frequency','p_value','rs_id']
 for line in f:
  md5.update(line);sha_src.update(line);size+=len(line);raw_rows+=1;c=line.rstrip(b'\r\n').split(b'\t');assert len(c)==9
  for rs in set(c[8].split(b',')):
   target=targets.get(rs)
   if target is None:continue
   rsrows[rs]+=1;chrom,pos,alleles,_,_,_=target
   if c[0].decode()!=chrom or int(c[1])!=pos:
    if frozenset((c[2].decode(),c[3].decode()))==alleles:other_coords[rs].append(f'{c[0].decode()}:{c[1].decode()}')
    continue
   coords[rs]+=1
   if frozenset((c[2].decode(),c[3].decode()))!=alleles:badalleles[rs]+=1;continue
   hits[rs]+=1
assert size==prior['raw_bytes'] and md5.hexdigest()==prior['publisher_md5'] and sha_src.hexdigest()==prior['raw_sha256']
for rs,target in targets.items():
 if hits[rs]==1:counts['unique_exact_match']+=1
 elif hits[rs]>1:counts['multiple_exact_match']+=1;example('multiple_exact_match',rs,hits[rs]);unresolved.append({'rsid':rs.decode(),'status':'multiple_exact_match','z':selected[rs][2]})
 elif not rsrows[rs]:counts['rsid_absent']+=1;example('rsid_absent',rs,None);unresolved.append({'rsid':rs.decode(),'status':'rsid_absent','z':selected[rs][2]})
 elif not coords[rs]:
  counts['coordinate_mismatch']+=1;example('coordinate_mismatch',rs,{'lift':target[:2],'source_rsid_rows':rsrows[rs]});counts['coordinate_mismatch_alt_contig' if '_' in target[0] else 'coordinate_mismatch_primary_contig']+=1;unresolved.append({'rsid':rs.decode(),'status':'coordinate_mismatch','z':selected[rs][2],'chr37':target[3],'pos37':target[4],'chr38_lift':target[0],'pos38_lift':target[1],'other_source_coords':';'.join(other_coords[rs])})
 else:counts['allele_mismatch']+=1;example('allele_mismatch',rs,{'lift':target[:2],'coordinate_rows':coords[rs]});unresolved.append({'rsid':rs.decode(),'status':'allele_mismatch','z':selected[rs][2]})
assert sum(counts[x] for x in ('unique_exact_match','multiple_exact_match','rsid_absent','coordinate_mismatch','allele_mismatch'))==len(targets)
assert len(unresolved)==len(selected)-counts['unique_exact_match']
with (p/'hfrs_coordinate_unresolved.tsv').open('w',newline='') as f:
 w=csv.DictWriter(f,fieldnames=['rsid','status','z','chr37','pos37','chr38_lift','pos38_lift','other_source_coords'],delimiter='\t',restval='');w.writeheader();w.writerows(unresolved)
result={'selected':len(selected),'counts':dict(counts),'examples':dict(examples),'raw_source':str(source),'raw_sha256':sha_src.hexdigest(),'raw_md5':md5.hexdigest(),'reference_ldscore_sha256':{str(i):sha(ref/f'{i}.l2.ldscore.gz') for i in range(1,23)},'selected_sha256':sha(final),'chain_path':str(chain),'chain_sha256':sha(chain),'limitation':'Position and allele consistency only, selected rsIDs; does not validate reference sequence, source model, ancestry or sample overlap.'};(p/'hfrs_coordinate_receipt.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result['counts'],indent=2),flush=True)
