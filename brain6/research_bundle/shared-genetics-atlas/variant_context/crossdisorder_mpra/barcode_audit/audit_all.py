import csv,gzip,tarfile,re,collections,json,hashlib
from pathlib import Path
root=Path('/Volumes/Extreme SSD/brain6-crossdisorder-mpra');p=Path(__file__).resolve().parent
coverage=[]
with (p.parent/'ranked_coverage.tsv').open() as f:
 for r in csv.DictReader(f,delimiter='\t'):
  if r['dna_present']=='True' and r['rna_present']=='True':coverage.append(r)
want={}
for r in coverage:
 m=re.fullmatch(r'(chr\d+):(\d+)_([ACGT]):([ACGT])',r['source_key']);assert m
 want[r['rsid']]={'chr':m[1],'pos':int(m[2]),'ref':m[3],'alt':m[4]}
assert len(want)==42
by_coord={(v['chr'],v['pos']):(rs,v) for rs,v in want.items()}
assert len(by_coord)==42
mapped=collections.defaultdict(set);barcode_targets=collections.defaultdict(set)
with gzip.open(root/'unique_bcdat.tsv.gz','rt') as f:
 for r in csv.DictReader(f,delimiter='\t'):
  m=re.match(r'^(chr\d+)_(\d+)_.*_([ACGT])$',r['variant'])
  if not m:continue
  match=by_coord.get((m[1],int(m[2])))
  if not match:continue
  rs,v=match
  if m[3] not in (v['ref'],v['alt']) or ':'+v['ref']+':'+v['alt'] not in r['variant']:continue
  tag=(rs,'REF' if m[3]==v['ref'] else 'ALT');mapped[tag].add(r['barcode']);barcode_targets[r['barcode']].add(tag)
observed=collections.defaultdict(set);observed_files=collections.defaultdict(set); reverse=collections.defaultdict(set)
table=str.maketrans('ACGT','TGCA');rc=lambda s:s.translate(table)[::-1]
with tarfile.open(root/'GSE244011_RAW.tar') as t:
 files=t.getmembers();assert len(files)==44
 for m in files:
  with gzip.open(t.extractfile(m),'rt') as f:
   for row in csv.DictReader(f):
    b=row['barcode']
    for tag in barcode_targets.get(b,()):observed[tag].add(b);observed_files[tag].add(m.name)
    for tag in barcode_targets.get(rc(b),()):reverse[tag].add(rc(b))
rows=[]
for rs in sorted(want):
 item={'rsid':rs}
 for a in ('REF','ALT'):
  key=(rs,a);item[a+'_mapped']=len(mapped[key]);item[a+'_observed_exact']=len(observed[key]);item[a+'_files_with_match']=len(observed_files[key]);item[a+'_observed_reverse_complement']=len(reverse[key])
 rows.append(item)
with (p/'coverage_42.tsv').open('w',newline='') as f:
 w=csv.DictWriter(f,fieldnames=list(rows[0]),delimiter='\t');w.writeheader();w.writerows(rows)
summary={'candidates':len(rows),'both_alleles_exact':sum(bool(r['REF_observed_exact'] and r['ALT_observed_exact']) for r in rows),'neither_allele_exact':sum(not r['REF_observed_exact'] and not r['ALT_observed_exact'] for r in rows),'one_allele_only':sum(bool(r['REF_observed_exact'])!=bool(r['ALT_observed_exact']) for r in rows),'any_reverse_complement':sum(bool(r['REF_observed_reverse_complement'] or r['ALT_observed_reverse_complement']) for r in rows),'source_files':len(files),'source_archive_sha256':hashlib.sha256((root/'GSE244011_RAW.tar').read_bytes()).hexdigest(),'barcode_map_sha256':hashlib.sha256((root/'unique_bcdat.tsv.gz').read_bytes()).hexdigest(),'coverage_output_sha256':hashlib.sha256((p/'coverage_42.tsv').read_bytes()).hexdigest()}
(p/'coverage_receipt.json').write_text(json.dumps(summary,indent=2)+'\n');print(json.dumps(summary,indent=2));print('no exact',','.join(r['rsid'] for r in rows if not r['REF_observed_exact'] and not r['ALT_observed_exact']))
