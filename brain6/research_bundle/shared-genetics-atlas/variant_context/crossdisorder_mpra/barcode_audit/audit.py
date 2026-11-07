import csv,gzip,tarfile,collections,re,json,hashlib
from pathlib import Path
root=Path('/Volumes/Extreme SSD/brain6-crossdisorder-mpra');p=Path(__file__).resolve().parent
wanted={'rs4916723':('chr5',87854395,'A','C'),'rs11242522':('chr5',103904914,'T','C')}
barcode_map={};n=collections.Counter();duplicate=collections.Counter()
with gzip.open(root/'unique_bcdat.tsv.gz','rt') as f:
 for r in csv.DictReader(f,delimiter='\t'):
  if r['rsid'] not in wanted:continue
  chrom,pos,ref,alt=wanted[r['rsid']];name=r['variant'];m=re.match(r'^(chr\d+)_(\d+)_.*_([ACGT])$',name)
  if not m or m[1]!=chrom or int(m[2])!=pos or m[3] not in (ref,alt):continue
  tag=(r['rsid'],'REF' if m[3]==ref else 'ALT')
  if r['barcode'] in barcode_map and barcode_map[r['barcode']]!=tag:duplicate[r['barcode']]+=1
  barcode_map[r['barcode']]=tag;n[tag]+=1
assert not duplicate
results=[]
with tarfile.open(root/'GSE244011_RAW.tar') as t:
 members=t.getmembers();assert len(members)==44
 for idx,m in enumerate(members):
  kind='DNA' if idx<22 else 'RNA';rep=idx+1 if idx<22 else idx-21;counts=collections.Counter();detected=collections.Counter()
  with gzip.open(t.extractfile(m),'rt') as f:
   for row in csv.DictReader(f):
    tag=barcode_map.get(row['barcode'])
    if tag:counts[tag]+=int(row['num_reads']);detected[tag]+=1
  for (rs,allele),count in sorted(counts.items()):results.append({'rsid':rs,'allele':allele,'kind':kind,'replicate':rep,'barcode_count':detected[(rs,allele)],'raw_reads':count,'source_file':m.name})
with (p/'raw_barcode_totals.tsv').open('w',newline='') as f:
 w=csv.DictWriter(f,fieldnames=list(results[0]),delimiter='\t');w.writeheader();w.writerows(results)
sha=lambda x:hashlib.sha256((root/x).read_bytes()).hexdigest()
(p/'receipt.json').write_text(json.dumps({'archive_sha256':sha('GSE244011_RAW.tar'),'barcode_map_sha256':sha('unique_bcdat.tsv.gz'),'source_rows':44,'mapped_barcode_labels':{str(k):v for k,v in n.items()},'output_sha256':hashlib.sha256((p/'raw_barcode_totals.tsv').read_bytes()).hexdigest(),'warning':'raw mapped reads and barcodes, not the authors post-QC counts; ordering assumed from file order, no full QC attempted'},indent=2)+'\n');print(n,len(results))
