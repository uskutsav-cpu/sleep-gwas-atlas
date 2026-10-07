from pathlib import Path
import json,hashlib,datetime,subprocess,gzip,csv,shutil
root=Path.cwd(); out=root/'brain6/translational_psychiatry_research_v1/research_v1/discovery';out.mkdir(parents=True,exist_ok=True)
raw=root/'work/research-discovery-raw';raw.mkdir(exist_ok=True)
sources=[dict(trait='insomnia',url='https://vu.data.surf.nl/index.php/s/ACjPo9hYZPy2zxl/download?path=%2F&files=Insomnia_sumstats_Jansenetal.txt.gz',bytes=314778585,sha256='32848cd92a6324c9048cad4a804cf1546299af1d0281981c6ca43c2046a9065b',access='PUBLIC_CURRENT_HTTP200_SOURCE_LICENSE_REVIEW_REQUIRED_BEFORE_REDISTRIBUTION'),dict(trait='adhd',url='https://ndownloader.figshare.com/files/40036684',bytes=200588408,sha256='c58a96031ec44b1edba81c91d603037a2efd03fcd946531424866e3f5f1d40c5',access='FIGSHARE_CC_BY_4.0_PUBLIC')]
windows={'A':(11,112459489,114257728),'B':(6,100630147,102636772),'C':(5,103447968,104447968)}
manifest={'frozen_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'purpose':'Recover exact discovery beta/SE/alleles/nativeN to determine admissibility of genotype-backed local follow-up; not independent replication','source_identity_from':'Historical frozen source ledger, current HTTPcontent length and Figshare officialmetadata','planned_bytes':sum(x['bytes'] for x in sources),'disk_free_before':shutil.disk_usage(root).free,'max_total_raw_bytes':600000000,'windows_grch37':windows,'sources':sources,'expected_peak_RAM_bytes':100000000,'stop_conditions':['SHA256orlength mismatch','free disk <700MBbeforeacquisition','source access condition change','regionalinputinterpretation hold untilsourcecontractreview'],'script_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
p=out/'acquisition_manifest.json'
if p.exists():raise SystemExit('Do not overwrite acquisition freeze')
p.write_text(json.dumps(manifest,indent=2)+'\n')
receipts=[]
for src in sources:
 path=raw/(src['trait']+'.txt.gz'); print('Acquiring',src['trait'],flush=True)
 subprocess.run(['curl','--fail','--silent','--show-error','--location','--retry','2','--max-time','900','--max-filesize',str(src['bytes']),src['url'],'-o',str(path)],check=True)
 h=hashlib.file_digest(path.open('rb'),'sha256').hexdigest()
 if path.stat().st_size!=src['bytes'] or h!=src['sha256']:raise ValueError('Exact source mismatch')
 handles={};writers={};counts={k:0 for k in windows}
 with gzip.open(path,'rt') as f:
  header=next(f).split();idxchr=header.index('CHR');idxbp=header.index('BP'); print(src['trait'],'header',header,flush=True)
  for key in windows:
   handles[key]=gzip.open(out/(src['trait']+'_'+key+'.raw.tsv.gz'),'wt');writers[key]=csv.writer(handles[key],delimiter='\t');writers[key].writerow(header)
  total=0
  for line in f:
   total+=1;row=line.split()
   try:chr=int(row[idxchr]);bp=int(row[idxbp])
   except (ValueError,IndexError):continue
   for key,(chrom,left,right) in windows.items():
    if chr==chrom and left<=bp<=right:writers[key].writerow(row);counts[key]+=1
  for handle in handles.values():handle.close()
 rec={'trait':src['trait'],'source_sha256':h,'bytes':path.stat().st_size,'rows_scanned':total,'regional_counts':counts,'regional_sha256':{key:hashlib.file_digest((out/(src['trait']+'_'+key+'.raw.tsv.gz')).open('rb'),'sha256').hexdigest() for key in windows},'finished_utc':datetime.datetime.now(datetime.timezone.utc).isoformat()};receipts.append(rec)
 (out/'acquisition_receipts.json').write_text(json.dumps(receipts,indent=2)+'\n');print(rec,flush=True)
