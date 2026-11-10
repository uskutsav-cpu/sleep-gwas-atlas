import csv,gzip,re,collections,openpyxl,json,hashlib
book='work/matoba2026/table3.xlsx';w=openpyxl.load_workbook(book,read_only=True,data_only=True)
s=w['a. Regulatory activity (Veh)'];it=iter(s.iter_rows(values_only=True));h=next(it);loc=collections.defaultdict(set)
for r in it:
 if len(r)>9 and r[1] and r[8] and r[9]:loc[str(r[1])].add((str(r[8]),int(r[9])))
s=w['b. Allelic effect (Veh)'];it=iter(s.iter_rows(values_only=True));h=next(it); effects={}
for r in it:
 if len(r)>7 and r[0] is not None and isinstance(r[7],(float,int)):effects[str(r[0])]=(float(r[7]),str(r[1]),str(r[2]))
seq=collections.defaultdict(lambda:collections.defaultdict(set))
with gzip.open('work/matoba2026/barcode_statistics.gz','rt') as f:
 for r in csv.DictReader(f,delimiter='\t'):
  m=re.search(r'_SNP_(rs\d+)_.*_(A[12])_',r['name'])
  if m and set(r['variant'])<=set('ACGT'):seq[m[1]][m[2]].add(r['variant'])
counts=collections.Counter();good=[]
for rs,(p,ref,alt) in effects.items():
 if len(loc[rs])!=1:counts['location_missing_or_ambiguous']+=1;continue
 chrom=next(iter(loc[rs]))[0];
 if chrom not in ('1','2'):counts['not_heldout_chr']+=1;continue
 if len(seq[rs]['A1'])!=1 or len(seq[rs]['A2'])!=1:counts['sequence_pair_missing_or_ambiguous']+=1;continue
 a=next(iter(seq[rs]['A1']));b=next(iter(seq[rs]['A2']))
 if len(a)!=len(b) or sum(x!=y for x,y in zip(a,b))!=1:counts['not_single_base_pair']+=1;continue
 i=next(i for i in range(len(a)) if a[i]!=b[i])
 if {a[i],b[i]}!={ref,alt}:
  complement=str.maketrans('ACGT','TGCA')
  if {a[i].translate(complement),b[i].translate(complement)}=={ref,alt}: counts['reverse_complement_pair_possible']+=1
  else: counts['unresolved_allele_mismatch']+=1
  continue
 counts['usable_pairs']+=1;counts['detected' if p<.05 else 'not_detected']+=1;good.append(rs)
print(json.dumps({'baseline_reported_effects':len(effects),'exclusions_and_counts':dict(counts),'genomic_build':'not confirmed','warning':'feasibility only; genomic-region and sequence overlap, alleles and licensing remain unverified'},indent=2))
