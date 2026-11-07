#!/usr/bin/env python3
import collections,csv,gzip,hashlib,json,re,tarfile
from pathlib import Path
root=Path('/Volumes/Extreme SSD/brain6-crossdisorder-mpra')
out=Path(__file__).resolve().parent
parent=out.parent
variants={'rs185260','rs325485','rs4916723','rs62369151'}
with (parent.parent/'ranked_coverage.tsv').open() as f:
    selected=[r for r in csv.DictReader(f,delimiter='\t') if r['rsid'] in variants]
assert len(selected)==4
want={}
for r in selected:
    m=re.fullmatch(r'(chr\d+):(\d+)_([ACGT]):([ACGT])',r['source_key']);assert m
    want[(m[1],int(m[2]))]=(r['rsid'],m[3],m[4])
assert len(want)==4
mapped=collections.defaultdict(set)
mapping=root/'unique_bcdat.tsv.gz';archive=root/'GSE244011_RAW.tar'
with gzip.open(mapping,'rt') as f:
 for row in csv.DictReader(f,delimiter='\t'):
    m=re.match(r'^(chr\d+)_(\d+)_.*_([ACGT])$',row['variant'])
    if not m or (m[1],int(m[2])) not in want:continue
    rs,ref,alt=want[(m[1],int(m[2]))]
    if m[3] not in (ref,alt) or ':'+ref+':'+alt not in row['variant']:continue
    b=row['barcode'];assert len(b)==20 and set(b)<=set('ACGT')
    mapped[(rs,'REF' if m[3]==ref else 'ALT')].add(b)
assert len(mapped)==8 and all(mapped.values())
neighbors=collections.defaultdict(set)
for tag,barcodes in mapped.items():
 for b in barcodes:
  for i,c in enumerate(b):
   for a in 'ACGT':
    if a!=c:neighbors[b[:i]+a+b[i+1:]].add((tag,b))
observed=collections.defaultdict(set);originals=collections.defaultdict(set)
files=collections.defaultdict(set);reads=collections.Counter();multi_tag=collections.Counter();multi_source=collections.Counter()
with tarfile.open(archive) as tar:
 members=tar.getmembers();assert len(members)==44
 for member in members:
  with gzip.open(tar.extractfile(member),'rt') as f:
   for row in csv.DictReader(f):
    b=row['barcode'];matches=neighbors.get(b)
    if not matches:continue
    tags={tag for tag,_ in matches}
    for tag in tags:
     observed[tag].add(b);files[tag].add(member.name);reads[tag]+=int(row['num_reads'])
     originals[tag].update(source for t,source in matches if t==tag)
     if len(tags)>1:multi_tag[tag]+=1
     if sum(t==tag for t,_ in matches)>1:multi_source[tag]+=1
rows=[]
for tag in sorted(mapped):
 rows.append({'rsid':tag[0],'allele':tag[1],'mapped_labels':len(mapped[tag]),
  'mapped_labels_with_neighbor':len(originals[tag]),'distinct_raw_neighbors':len(observed[tag]),
  'files_with_neighbor':len(files[tag]),'raw_neighbor_reads':reads[tag],
  'multi_tag_file_rows':multi_tag[tag],'multi_original_file_rows':multi_source[tag]})
with (out/'results.tsv').open('w',newline='') as f:
 w=csv.DictWriter(f,fieldnames=list(rows[0]),delimiter='\t',lineterminator='\n');w.writeheader();w.writerows(rows)
def sha(path):
 h=hashlib.sha256()
 with path.open('rb') as f:
  for b in iter(lambda:f.read(1024*1024),b''):h.update(b)
 return h.hexdigest()
prior=json.loads((parent/'coverage_receipt.json').read_text())
assert sha(archive)==prior['source_archive_sha256'] and sha(mapping)==prior['barcode_map_sha256']
(out/'receipt.json').write_text(json.dumps({'source_archive_sha256':sha(archive),'mapping_sha256':sha(mapping),
 'result_sha256':sha(out/'results.tsv'),'files_scanned':len(members),'rows':rows},indent=2)+'\n')
print(json.dumps(rows,indent=2))
