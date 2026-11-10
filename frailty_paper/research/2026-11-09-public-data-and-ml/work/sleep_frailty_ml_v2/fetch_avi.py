"""Resumable exact-allele AVI annotation, without downloading the 88 GB table."""
from pathlib import Path
import sys,csv,json,hashlib,collections,concurrent.futures,threading,time,argparse

ROOT=Path(__file__).resolve().parents[2];WORK=Path(__file__).resolve().parent;OUT=ROOT/'outputs/sleep_frailty_ai_v2'
# Existing official PyPI pysam installation, same CPython ABI; no credentials used.
sys.path.insert(0,str(ROOT/'work/stage5_env/lib/python3.13/site-packages'))
import pysam

URL='https://grr.iossifovlab.com/hg38/scores/AVI/alphagenome_variant_impact_score_snvs.tsv.gz'
INDEX=Path('/Volumes/Extreme SSD/codex-sleep-frailty-work/stage4_atlas_cache/avi_snvs.tsv.gz.tbi')
assert hashlib.md5(INDEX.read_bytes()).hexdigest()=='911ee04b6e8fba1882465221b5a14f2c'
p=argparse.ArgumentParser();p.add_argument('--limit',type=int,default=0);p.add_argument('--workers',type=int,default=4);p.add_argument('--pairs-csv');p.add_argument('--prefix',default='avi');args=p.parse_args()
assert args.prefix.replace('_','').isalnum()
cache=WORK/(args.prefix+'_chunks');cache.mkdir(exist_ok=True);OUT.mkdir(exist_ok=True)
groups=collections.defaultdict(set)
if not args.pairs_csv:
    with (ROOT/'outputs/sleep_frailty_ai_v1/split_manifest.csv').open() as f:
        for r in csv.DictReader(f):
            chrom,pos,ref,alt=r['variant'].split('-');groups[('chr'+chrom,(int(pos)-1)//16384)].add((int(pos),ref,alt))
with (Path(args.pairs_csv) if args.pairs_csv else ROOT/'work/sleep_frailty_ml/project_sequence_pairs.csv').open() as f:
    for r in csv.DictReader(f):
        chrom=r['chromosome'];chrom=chrom if chrom.startswith('chr') else 'chr'+chrom
        pos=int(r['position_hg38']);groups[(chrom,(pos-1)//16384)].add((pos,r['genome_reference_allele'],r['tested_alternate_allele']))
local=threading.local()
def fetch(item):
    (chrom,window),wanted=item;wanted=sorted(wanted)
    signature=hashlib.sha256(json.dumps([URL,chrom,wanted]).encode()).hexdigest()
    path=cache/f'{chrom}_{window}.json'
    if path.exists():
        old=json.loads(path.read_text())
        if old.get('signature')==signature and old.get('complete'):return old
    for attempt in range(3):
        try:
            if not hasattr(local,'tabix'):local.tabix=pysam.TabixFile(URL,index=str(INDEX))
            target=set(wanted);found={};unexpected_reference=set()
            for line in local.tabix.fetch(chrom,min(k[0] for k in wanted)-1,max(k[0] for k in wanted)):
                f=line.split('\t');key=(int(f[1]),f[2],f[3])
                if key in target:
                    name=f'{chrom}:{key[0]}:{key[1]}:{key[2]}'
                    value={'avi_raw':float(f[4]),'avi_phred':float(f[5])}
                    assert name not in found or found[name]==value
                    found[name]=value
            data={'chromosome':chrom,'signature':signature,'source':URL,'complete':True,'requested':[f'{chrom}:{p}:{r}:{a}' for p,r,a in wanted],'found':found}
            tmp=path.with_suffix('.tmp');tmp.write_text(json.dumps(data,separators=(',',':'))+'\n');tmp.replace(path)
            return data
        except Exception:
            if hasattr(local,'tabix'):local.tabix.close();del local.tabix
            if attempt==2:raise
            time.sleep(attempt+1)
items=sorted(groups.items(),key=lambda x:(int(x[0][0].removeprefix('chr')),x[0][1]))
if args.limit:items=items[:args.limit]
print('Query groups',len(items),'workers',args.workers,flush=True)
results=[];errors=[]
with concurrent.futures.ThreadPoolExecutor(max_workers=args.workers) as pool:
    futures={pool.submit(fetch,x):x[0] for x in items}
    for i,future in enumerate(concurrent.futures.as_completed(futures),1):
        try:results.append(future.result())
        except Exception as e:errors.append({'group':futures[future],'error':repr(e)})
        if i%50==0 or i==len(items):print('Finished groups',i,'/',len(items),'errors',len(errors),flush=True)
all_found={k:v for r in results for k,v in r['found'].items()};requested=sorted({k for r in results for k in r['requested']})
target=OUT/(args.prefix+('_partial.csv' if args.limit else '_exact_allele_scores.csv'))
with target.open('w',newline='') as f:
    writer=csv.DictWriter(f,fieldnames=['variant_GRCh38','avi_raw','avi_phred','status']);writer.writeheader()
    for k in requested:writer.writerow({'variant_GRCh38':k,**all_found.get(k,{}),'status':'exact_allele_match' if k in all_found else 'not_found'})
audit={'source':URL,'official_artifact_designation':'AVI SNV scores, Permissive Use Downloadable artifacts section, observed 2026-11-06','official_downloads':'https://deepmind.google.com/science/alphagenome/downloads','index_md5':hashlib.md5(INDEX.read_bytes()).hexdigest(),'whole_table_md5_publisher_reported_not_verified':'d7753d6a3aba87a0b90ecb9df3b789e3','query_groups':len(items),'failed_groups':errors,'requested_unique_alleles_in_completed_groups':len(requested),'exact_matches':len(all_found),'missing':len(requested)-len(all_found),'output_sha256':hashlib.sha256(target.read_bytes()).hexdigest(),'partial_run':bool(args.limit),'missing_interpretation':'Unavailable score, not zero effect','generation_scope':'Only downloaded AVI scores, not detailed API outputs or feature breakdowns'}
(OUT/(args.prefix+('_partial_audit.json' if args.limit else '_retrieval_audit.json'))).write_text(json.dumps(audit,indent=2)+'\n')
print(json.dumps(audit,indent=2),flush=True)
if errors:sys.exit(1)
