"""Exact-allele permissible AVI annotation; reuse verified project score rows."""
from pathlib import Path
import sys,csv,json,hashlib,collections,concurrent.futures,threading,time,math,argparse
ROOT=Path(__file__).resolve().parents[2];WORK=Path(__file__).resolve().parent;OUT=ROOT/'outputs/sleep_frailty_ai_v3';PREV=ROOT/'outputs/sleep_frailty_ai_v2'
sys.path.insert(0,str(ROOT/'work/stage5_env/lib/python3.13/site-packages'));import pysam
URL='https://grr.iossifovlab.com/hg38/scores/AVI/alphagenome_variant_impact_score_snvs.tsv.gz';INDEX=Path('/Volumes/Extreme SSD/codex-sleep-frailty-work/stage4_atlas_cache/avi_snvs.tsv.gz.tbi')
assert hashlib.md5(INDEX.read_bytes()).hexdigest()=='911ee04b6e8fba1882465221b5a14f2c'
manifest=WORK/'avi_query_pairs.csv';rows=list(csv.DictReader(manifest.open()));wanted={r['variant'] for r in rows};known={};prior_files=[]
parser=argparse.ArgumentParser();parser.add_argument('--workers',type=int,default=4);args=parser.parse_args();assert 1<=args.workers<=8
for stem in ['avi','external_avi']:
    path=PREV/(stem+'_exact_allele_scores.csv');audit=json.loads((PREV/(stem+'_retrieval_audit.json')).read_text());assert audit['output_sha256']==hashlib.sha256(path.read_bytes()).hexdigest() and not audit['failed_groups']
    prior_files.append({'file':str(path),'sha256':audit['output_sha256']})
    for r in csv.DictReader(path.open()):
        if r['variant_GRCh38'] in wanted and r['status']=='exact_allele_match':known[r['variant_GRCh38']]={'avi_raw':float(r['avi_raw']),'avi_phred':float(r['avi_phred'])}
groups=collections.defaultdict(set)
for r in rows:
    if r['variant'] in known:continue
    groups[(r['chromosome'],(int(r['position_hg38'])-1)//16384)].add((int(r['position_hg38']),r['genome_reference_allele'],r['tested_alternate_allele']))
cache=WORK/'avi_chunks';cache.mkdir(exist_ok=True);local=threading.local()
def fetch(item):
    (chrom,window),keys=item;keys=sorted(keys);sig=hashlib.sha256(json.dumps([URL,chrom,keys]).encode()).hexdigest();path=cache/f'{chrom}_{window}.json'
    if path.exists():
        saved=json.loads(path.read_text())
        if saved.get('signature')==sig and saved.get('complete'):return saved
    for attempt in range(3):
        try:
            if not hasattr(local,'tabix'):local.tabix=pysam.TabixFile(URL,index=str(INDEX))
            sought=set(keys);found={}
            for line in local.tabix.fetch(chrom,min(k[0] for k in keys)-1,max(k[0] for k in keys)):
                f=line.split('\t');key=(int(f[1]),f[2],f[3])
                if key in sought:
                    score={'avi_raw':float(f[4]),'avi_phred':float(f[5])};assert all(math.isfinite(v) for v in score.values())
                    vk=f'{chrom}:{key[0]}:{key[1]}:{key[2]}';assert vk not in found or found[vk]==score;found[vk]=score
            result={'signature':sig,'complete':True,'found':found,'requested':[f'{chrom}:{p}:{a}:{b}' for p,a,b in keys]};tmp=path.with_suffix('.tmp');tmp.write_text(json.dumps(result,separators=(',',':'))+'\n');tmp.replace(path);return result
        except Exception:
            if hasattr(local,'tabix'):local.tabix.close();del local.tabix
            if attempt==2:raise
            time.sleep(attempt+1)
items=sorted(groups.items(),key=lambda x:(int(x[0][0].removeprefix('chr')),x[0][1]));errors=[];all_found=dict(known);completed=set(known)
print('Unique requested',len(wanted),'reused exact scores',len(known),'new query groups',len(items),'workers',args.workers,flush=True)
with concurrent.futures.ThreadPoolExecutor(max_workers=args.workers) as pool:
    futures={pool.submit(fetch,item):item[0] for item in items}
    for n,f in enumerate(concurrent.futures.as_completed(futures),1):
        try:
            r=f.result();all_found.update(r['found']);completed.update(r['requested'])
        except Exception as e:errors.append({'group':futures[f],'error':repr(e)})
        if n%100==0 or n==len(items):print('AVI groups',n,'/',len(items),'errors',len(errors),flush=True)
target=OUT/'avi_exact_allele_scores.csv'
with target.open('w',newline='') as f:
    w=csv.DictWriter(f,fieldnames=['variant_GRCh38','avi_raw','avi_phred','status']);w.writeheader()
    for key in sorted(wanted):
        status='exact_allele_match' if key in all_found else 'not_found' if key in completed else 'query_failed'
        w.writerow({'variant_GRCh38':key,**all_found.get(key,{}),'status':status})
audit={'source':URL,'official_artifact_designation':'AVI SNV scores, Permissive Use Downloadable Artifacts, observed 2026-11-06','terms':'https://deepmind.google.com/science/alphagenome/terms','downloads':'https://deepmind.google.com/science/alphagenome/downloads','genome_build':'GRCh38; source Siraj alleles uniquely reciprocal-mapped from GRCh37','index_md5':hashlib.md5(INDEX.read_bytes()).hexdigest(),'whole_data_table_publisher_md5_not_independently_verified':'d7753d6a3aba87a0b90ecb9df3b789e3','manifest_sha256':hashlib.sha256(manifest.read_bytes()).hexdigest(),'prior_score_files':prior_files,'requested':len(wanted),'reused_exact_scores':len(known),'exact_matches':len(all_found),'unavailable':len(wanted-completed) + len(completed-set(all_found)),'failed_groups':errors,'query_groups':len(items),'maximum_query_workers':args.workers,'output_sha256':hashlib.sha256(target.read_bytes()).hexdigest(),'training_usage':'Only designated downloaded AVI scores. No detailed API predictions, AVI feature breakdowns, or other restricted output derivatives.','missingness':'Unmatched variants remain unavailable, never zero biological effect; audit possible reference/build differences','clinical_use':False}
(OUT/'avi_retrieval_audit.json').write_text(json.dumps(audit,indent=2)+'\n');print(json.dumps(audit,indent=2),flush=True)
if errors:sys.exit(1)
