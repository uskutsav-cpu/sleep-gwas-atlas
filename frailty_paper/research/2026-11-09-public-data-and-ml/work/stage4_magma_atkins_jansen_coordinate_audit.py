import gzip,collections,json,pathlib
src={'Atkins_FI':('/Volumes/Extreme SSD/sleep-gwas-atlas-frailty-v1/data/raw/frailty.txt.gz',0,1,2),'Jansen_insomnia':('/Volumes/Extreme SSD/sleep-gwas-atlas-frailty-v1/data/raw/insomnia.txt.gz',0,2,3)};want=collections.defaultdict(list)
for name,(path,si,ci,pi) in src.items():
 with gzip.open(path,'rt') as f:
  next(f)
  for l in f:
   a=l.rstrip('\n').split('\t');rs=a[si]
   if rs.startswith('rs') and rs[2:].isdigit() and int(rs[2:])%1000==0:want[rs].append((name,a[ci],a[pi]))
ref=collections.defaultdict(list)
with open('work/stage4_magma_reference/g1000_eur.bim') as f:
 for line in f:
  a=line.split()
  if a[1] in want:ref[a[1]].append((a[0],a[3]))
c=collections.defaultdict(collections.Counter);examples=collections.defaultdict(list)
for rs,uses in want.items():
 for name,ch,pos in uses:
  c[name]['sampled_source_rows']+=1
  if rs not in ref:c[name]['not_directly_in_reference']+=1;continue
  if len(ref[rs])!=1:c[name]['duplicate_reference_id']+=1;continue
  refch,refpos=ref[rs][0];normch='23' if ch=='X' else ch
  if ch in {str(i) for i in range(1,23)}:
   c[name]['autosomal_direct_match_rows']+=1
   if (ch,pos)!=(refch,refpos):c[name]['autosomal_direct_coordinate_mismatch']+=1
  if (normch,pos)==(refch,refpos):c[name]['direct_id_and_coordinate_match_normalized_X']+=1
  else:
   c[name]['direct_id_coordinate_mismatch']+=1
   if len(examples[name])<10:examples[name].append({'rsid':rs,'source':[ch,pos],'reference':ref[rs][0]})
out={'sample':'fixed rsID number divisible by 1000; entire raw source, not restricted to MAGMA filters','results':{k:dict(v) for k,v in c.items()},'examples':examples,'limitations':'direct ID only; not allele, synonym or complete build verification'}
pathlib.Path('outputs/stage4_magma_atkins_jansen_coordinate_sample_audit.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out,indent=2))
