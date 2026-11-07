import csv,gzip,json,hashlib,math
from pathlib import Path
p=Path(__file__).resolve().parent;m=Path('/Volumes/Extreme SSD/Utsav-Research-Archive/Sleep-GWAS/Sep1-local-dependencies/data');receipt=json.loads((p/'receipt.json').read_text());summary={}
for trait in ('adhd','insomnia'):
 ns={round(x[0],3) for x in receipt[trait]['most_common_n'][:5]}; selected={}
 with gzip.open(m/'munged'/f'{trait}.sumstats.gz','rt') as f:
  for row in csv.DictReader(f,delimiter='\t'):
   try:n=round(float(row['N']),3)
   except (TypeError,ValueError):continue
   if n in ns and n not in selected:selected[n]=row
   if len(selected)==len(ns):break
 assert len(selected)==len(ns)
 ids={x['SNP'] for x in selected.values()};matches={i:[] for i in ids};raw=m/'raw'/f'{trait}.txt.gz';count=0;h=hashlib.sha256()
 with raw.open('rb') as f:
  for block in iter(lambda:f.read(1048576),b''):h.update(block)
 with gzip.open(raw,'rt') as f:
  hdr=next(f).split()
  for line in f:
   count+=1;vals=line.split()
   if len(vals)!=len(hdr):continue
   r=dict(zip(hdr,vals))
   if r['SNP'] in ids:matches[r['SNP']].append(r)
 records=[]
 for n,row in selected.items():
  hits=[x for x in matches[row['SNP']] if {x['A1'],x['A2']}=={row['A1'],row['A2']}];rec={'SNP':row['SNP'],'munged_N':row['N'],'raw_unique_allele_match_count':len(hits)}
  if len(hits)==1:
   hit=hits[0]
   if trait=='adhd':
    ca=int(hit['Nca']);co=int(hit['Nco']);neff=4*ca*co/(ca+co);rec.update(Nca=ca,Nco=co,formula_neff=neff,difference=float(row['N'])-neff)
   else:rec['raw_N']=hit['N']
  records.append(rec)
 summary[trait]={'raw_path':str(raw),'raw_sha256':h.hexdigest(),'raw_rows':count,'records':records}
(p/'raw_sample.json').write_text(json.dumps(summary,indent=2)+'\n');print(json.dumps(summary,indent=2))
