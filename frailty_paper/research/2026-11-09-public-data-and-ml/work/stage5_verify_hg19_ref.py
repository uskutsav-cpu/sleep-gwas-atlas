import csv,collections,subprocess,concurrent.futures,json,pathlib,hashlib,os,datetime
out=pathlib.Path('outputs');q=list(csv.DictReader(open(out/'stage5_queue_snv_reciprocal_liftover.csv')));g=collections.defaultdict(list)
for x in q:
 if x['status']=='unique_positive_reciprocal_mapping':g[(x['chr37'],int(x['pos37'])//100000)].append(x)
dir=pathlib.Path('work/stage5_hg19_reference');dir.mkdir(exist_ok=True)
def fetch(item):
 (chrom,window),v=item;start=min(int(x['pos37'])-1 for x in v);end=max(int(x['pos37']) for x in v);url=f'https://api.genome.ucsc.edu/getData/sequence?genome=hg19;chrom=chr{chrom};start={start};end={end}';target=dir/f'chr{chrom}_{window}.json'
 if target.exists():return target
 r=subprocess.run(['curl','-fsS','--max-time','45',url],capture_output=True,check=True);d=json.loads(r.stdout);assert d['start']==start and d['end']==end and len(d['dna'])==end-start
 tmp=target.with_suffix('.tmp');tmp.write_text(json.dumps({'url':url,'result':d}));os.replace(tmp,target);return target
with concurrent.futures.ThreadPoolExecutor(max_workers=3) as pool:
 for _ in pool.map(fetch,g.items()):pass
rows=[]
for (chrom,window),vs in g.items():
 obj=json.loads((dir/f'chr{chrom}_{window}.json').read_text());d=obj['result'];seq=d['dna'].upper()
 for x in vs:
  base=seq[int(x['pos37'])-1-d['start']];ref=x['source_ref'];alt=x['source_alt'];status='reference_same' if base==ref else 'reference_is_grch38_alt' if base==alt else 'reference_neither_allele'
  rows.append({'variant_GRCh38':x['variant_GRCh38'],'chr37':chrom,'pos37':x['pos37'],'ref38':ref,'alt38':alt,'hg19_reference_base':base,'reference_status':status,'ucsc_source_url':obj['url'],'caveat':'Positive-strand single-base reference comparison; does not independently authenticate original GWAS source genome build or its effect allele convention'})
with open(out/'stage5_reciprocal_snv_hg19_reference_check.csv','w',newline='') as f:w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
report={'date_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'variants':len(rows),'query_intervals':len(g),'status_counts':dict(collections.Counter(x['reference_status'] for x in rows)),'input_sha256':hashlib.sha256((out/'stage5_queue_snv_reciprocal_liftover.csv').read_bytes()).hexdigest(),'source':'https://api.genome.ucsc.edu/getData/sequence','note':'Reference check addresses mapped bases, not model validation or genetic causality.'}
(out/'stage5_reciprocal_snv_hg19_reference_check_audit.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report,indent=2));assert len(rows)==405
