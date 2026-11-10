"""Bounded public-score/native-sequence retrieval, with resumable exact signatures."""
import sys,argparse,subprocess,urllib.request,time,concurrent.futures
from common import *
def one(r):
    sys.path.insert(0,str(ROOT/'work/stage5_env/lib/python3.13/site-packages'));import pysam
    assert hashlib.md5(INDEX.read_bytes()).hexdigest()=='911ee04b6e8fba1882465221b5a14f2c';c=r['chromosome'];p=r['position_hg38'];start=p-150;end=start+300;url=f'https://api.genome.ucsc.edu/getData/sequence?genome=hg38&chrom={c}&start={start}&end={end}';signature=sha(DATA/'teacher_regions.json')+'_'+r['anchor_id'];cache=DATA/'teacher_cache';cache.mkdir(exist_ok=True);target=cache/(r['anchor_id']+'.json')
    if target.exists():
        old=json.loads(target.read_text())
        if old['signature']==signature and old['complete']:return old
    request=urllib.request.Request(url,headers={'User-Agent':'SleepFrailtyAcademicResearch/1.0'});raw=urllib.request.urlopen(request,timeout=35).read();d=json.loads(raw);seq=d['dna'].upper();assert d['chrom']==c and d['start']==start and d['end']==end and len(seq)==300
    with pysam.TabixFile(AVI_URL,index=str(INDEX)) as tab:
        scores={}
        for line in tab.fetch(c,p-16,p+15):
            f=line.split('\t');key=(int(f[1]),f[2],f[3]);value=[float(f[4]),float(f[5])];assert np.isfinite(value).all();assert key not in scores or scores[key]==value;scores[key]=value
    pairs=[];missing=[]
    for offset in range(-15,16):
        j=offset+15;ref=seq[j:j+270];position=p+offset
        if set(ref)-set('ACGT'):missing.append({'position':position,'reason':'non_ACGT_native_context'});continue
        base=ref[134]
        for alt in 'ACGT':
            if alt==base:continue
            key=(position,base,alt)
            if key not in scores:missing.append({'variant':f'{c}:{position}:{base}:{alt}','reason':'not_exactly_matched_in_AVI'});continue
            pairs.append({'variant_GRCh38':f'{c}:{position}:{base}:{alt}','offset':offset,'reference_base':base,'alternate_base':alt,'avi_raw':scores[key][0],'avi_phred':scores[key][1]})
    result={'anchor':r,'signature':signature,'complete':True,'native_sequence_300bp':seq,'UCSC_source_url':url,'UCSC_payload_sha256':hashlib.sha256(raw).hexdigest(),'AVI_source':AVI_URL,'AVI_index_md5':hashlib.md5(INDEX.read_bytes()).hexdigest(),'pairs':pairs,'unavailable':missing,'interpretation':'Synthetic single-base substitutions with teacher-predicted scores, not measured experimental effects'};write_json(target,result);return result
def main():
    parser=argparse.ArgumentParser();parser.add_argument('--one');parser.add_argument('--stage',choices=['pilot','expanded','reserved'],default='pilot');parser.add_argument('--workers',type=int,default=4);args=parser.parse_args();rows=json.loads((DATA/'teacher_regions.json').read_text());assert sha(DATA/'teacher_regions.json')==json.loads((OUT/'teacher_region_plan.json').read_text())['anchor_manifest_sha256']
    if args.one:
        one(next(r for r in rows if r['anchor_id']==args.one));return
    selected=[r for r in rows if r['pilot']] if args.stage=='pilot' else [r for r in rows if (r['split']=='reserved_teacher_test')==(args.stage=='reserved')];log=[]
    def fetch(r):
        last=None
        for attempt in range(3):
            try:
                subprocess.run([sys.executable,str(Path(__file__).resolve()),'--one',r['anchor_id']],timeout=100,check=True,capture_output=True,text=True);path=DATA/'teacher_cache'/(r['anchor_id']+'.json');d=json.loads(path.read_text());assert d['complete'] and d['signature']==sha(DATA/'teacher_regions.json')+'_'+r['anchor_id'];return {'anchor_id':r['anchor_id'],'cache_sha256':sha(path),'pairs':len(d['pairs']),'unavailable':len(d['unavailable'])}
            except Exception as e:last=repr(e);time.sleep(attempt+1)
        raise RuntimeError(r['anchor_id']+': '+str(last))
    errors=[]
    with concurrent.futures.ThreadPoolExecutor(max_workers=args.workers) as pool:
        futures={pool.submit(fetch,r):r for r in selected}
        for n,f in enumerate(concurrent.futures.as_completed(futures),1):
            try:log.append(f.result())
            except Exception as e:errors.append({'anchor_id':futures[f]['anchor_id'],'error':repr(e)})
            if n%20==0 or n==len(selected):print('Teacher collection',args.stage,n,'/',len(selected),'errors',len(errors),flush=True)
    write_json(OUT/('teacher_'+args.stage+'_retrieval.json'),{'selected_anchors':len(selected),'completed_anchors':len(log),'new_exact_teacher_pairs':sum(r['pairs'] for r in log),'unavailable':sum(r['unavailable'] for r in log),'chunks':log,'failed':errors,'training_targets_are_AlphaGenome_predictions':True,'reserved_teacher_test_collected':args.stage=='reserved','manifest_sha256':sha(DATA/'teacher_regions.json')});assert not errors
if __name__=='__main__':main()
