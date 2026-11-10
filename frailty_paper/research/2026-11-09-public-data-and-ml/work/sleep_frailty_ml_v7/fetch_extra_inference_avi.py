"""Exact public AVI values for expanded project/candidate-source inputs only.

The fitted-model training/test AVI arrays are asserted byte-identical.
"""
import sys,sqlite3,threading,concurrent.futures,collections
from common import *
from fetch_avi import URL,INDEX_MD5
sys.path.insert(0,str(ROOT/'work/stage5_env/lib/python3.13/site-packages'))
import pysam
def main():
    wanted={r['variant_GRCh38'] for r in json.loads((DATA/'all_project_variant_metadata.json').read_text()) if not r['unscorable_reason']}
    for r in csv.DictReader((DATA/'neuroimmune_physical_design_verified.csv').open()):wanted.add(r['variant_GRCh38'])
    original=json.loads((OUT/'AVI_collection_receipt.json').read_text())
    for name,h in original['outputs'].items():assert sha(DATA/name)==h
    index=SOURCES/'AVI_SNVS.tsv.gz.tbi';assert hashlib.md5(index.read_bytes()).hexdigest()==INDEX_MD5;db=sqlite3.connect(DATA/'avi_exact_scores.sqlite');db.execute('PRAGMA journal_mode=WAL');present={r[0] for r in db.execute('SELECT variant FROM scores')};groups=collections.defaultdict(set)
    for key in wanted-present:
        c,p,a,b=key.split(':');p=int(p);assert len(a)==len(b)==1 and a!=b;groups[(c,(p-1)//16384)].add((p,a,b))
    local=threading.local()
    def query(item):
        (c,w),keys=item;error=None
        for attempt in range(3):
            try:
                if not hasattr(local,'tab'):local.tab=pysam.TabixFile(URL,index=str(index))
                values={}
                if c in local.tab.contigs:
                    for line in local.tab.fetch(c,min(x[0] for x in keys)-1,max(x[0] for x in keys)):
                        f=line.split('\t');key=(int(f[1]),f[2],f[3])
                        if key in keys:
                            v=(float(f[4]),float(f[5]));assert np.isfinite(v).all() and v[1]>=0;assert key not in values or values[key]==v;values[key]=v
                return [(f'{c}:{p}:{a}:{b}',*values.get((p,a,b),(None,None)),int((p,a,b) in values),'additional_exact_public_AVI_inference_query') for p,a,b in keys]
            except Exception as e:
                error=repr(e)
                if hasattr(local,'tab'):
                    try:local.tab.close()
                    except Exception:pass
                    del local.tab
                time.sleep(attempt+1)
        raise RuntimeError(str((c,w))+':'+str(error))
    errors=[];done=0;started=time.time()
    with concurrent.futures.ThreadPoolExecutor(4) as pool:
        futures={pool.submit(query,item):item[0] for item in groups.items()}
        for future in concurrent.futures.as_completed(futures):
            try:db.executemany('INSERT OR REPLACE INTO scores VALUES (?,?,?,?,?)',future.result())
            except Exception as e:errors.append({'window':list(futures[future]),'error':repr(e)})
            done+=1
            if done%25==0:db.commit()
            if done%50==0 or done==len(groups):write_json(OUT/'extra_AVI_progress.json',{'done':done,'total':len(groups),'errors':errors,'seconds':time.time()-started});print('EXTRA EXACT AVI',done,'/',len(groups),flush=True)
    db.commit();assert not errors,errors;rows=[]
    for key in sorted(wanted):
        r=db.execute('SELECT raw,phred,available FROM scores WHERE variant=?',(key,)).fetchone();assert r is not None;rows.append({'variant_GRCh38':key,'raw':r[0] if r[2] else '','phred':r[1] if r[2] else '','available':r[2]})
    path=DATA/'extra_exact_AVI_scores.csv'
    with path.open('w') as f:w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
    db.close()
    for name,h in original['outputs'].items():assert sha(DATA/name)==h
    write_json(OUT/'extra_AVI_inference_receipt.json',{'requested_unique_variants':len(wanted),'new_query_windows':len(groups),'available':sum(r['available'] for r in rows),'missing':sum(not r['available'] for r in rows),'source':URL,'index_MD5':INDEX_MD5,'maximum_parallel_queries':4,'output':str(path),'sha256':sha(path),'original_training_and_reserved_AVI_arrays_unchanged':True,'purpose':'Additional genomic-inference and candidate-validation inputs only. No model fitting or outcome-based selection. Missing scores are unavailable, never biological zero. AVI is a permissible composite, not independent experimental truth.','whole_remote_data_MD5_independently_verified':False});print('EXTRA PUBLIC INFERENCE AVI COMPLETE',len(wanted),flush=True)
if __name__=='__main__':main()
