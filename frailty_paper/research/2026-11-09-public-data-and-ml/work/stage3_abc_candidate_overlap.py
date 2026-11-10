import csv,gzip,pathlib,hashlib,json,collections,time
p=pathlib.Path('work/stage3_region_cache/ABC_Nasser2021_predictions_GRCh37.txt.gz');m=list(csv.DictReader(open('outputs/shared22_eight_source_evidence_matrix.csv')));x=next(x for x in csv.DictReader(open('outputs/integrated_15_candidate_review.csv')) if x['variant']=='rs323509');markers=[{'variant':r['variant'],'chr':'chr'+r['chr_GRCh37_assumed'],'pos':int(r['position_GRCh37_assumed'])} for r in m]+[{'variant':'rs323509','chr':'chr'+x['GRCh37_chr'],'pos':int(x['GRCh37_pos'])}];assert len(markers)==23
by=collections.defaultdict(list)
for x in markers:by[x['chr']].append(x)
h=hashlib.sha256()
with p.open('rb') as f:
 while b:=f.read(8*1024*1024):h.update(b)
rows=0;out=[];start=time.time()
with gzip.open(p,'rt') as f:
 header=f.readline().rstrip('\n').split('\t'); assert header[0:3]==['chr','start','end'] and 'ABC.Score' in header and 'CellType' in header
 for line in f:
  rows+=1;a=line.rstrip('\n').split('\t');match=by.get(a[0])
  if not match:continue
  lo=int(a[1]);hi=int(a[2]);
  for z in match:
   if lo<=z['pos']-1<hi:out.append({'candidate':z['variant'],'candidate_GRCh37_1based':z['chr']+':'+str(z['pos']),**dict(zip(header,a))})
 if rows%1000000==0:print('rows read',rows,flush=True)
outfile=pathlib.Path('outputs/stage3_abc_candidate_regulatory_links.csv');fields=['candidate','candidate_GRCh37_1based']+header
with outfile.open('w',newline='') as f:w=csv.DictWriter(f,fieldnames=fields);w.writeheader();w.writerows(out)
counts=collections.Counter(x['candidate'] for x in out);cls=collections.Counter(x['class'] for x in out)
pathlib.Path('outputs/stage3_abc_candidate_regulatory_receipt.json').write_text(json.dumps({'source':'https://mitra.stanford.edu/engreitz/oak/public/Nasser2021/AllPredictions.AvgHiC.ABC0.015.minus150.ForABCPaperV3.txt.gz','source_paper':'https://doi.org/10.1038/s41586-021-03446-x','genome_build':'GRCh37 per AnnotSV data-source documentation; candidate input GRCh37 assumed','local_source_sha256':h.hexdigest(),'source_bytes':p.stat().st_size,'provider_checksum_verified':False,'gzip_read_to_end':True,'source_rows':rows,'matched_prediction_rows':len(out),'candidates_with_matches':len(counts),'rows_by_candidate':dict(counts),'rows_by_class':dict(cls),'definition':'1-based candidate coordinate converted to 0-based and intersected source 0-based half-open intervals; exact position only; no LD proxy expansion, allele-specific effect or causal evidence','seconds':round(time.time()-start,1)},indent=2)+'\n')
print('source rows',rows,'matched',len(out),'candidates',len(counts),'classes',dict(cls),flush=True)
