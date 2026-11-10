"""Bounded exact-allele AVI collection; one SQLite cache avoids tiny-file disk waste."""
import sys,gzip,sqlite3,concurrent.futures,threading,datetime
from common import *
sys.path.insert(0,str(ROOT/'work/stage5_env/lib/python3.13/site-packages'))
import pysam
URL='https://grr.iossifovlab.com/hg38/scores/AVI/alphagenome_variant_impact_score_snvs.tsv.gz'
INDEX_MD5='911ee04b6e8fba1882465221b5a14f2c'
def main():
    assert not (OUT/'reserved_evaluation_started.json').exists()
    bank_audit=json.loads((OUT/'development_data_audit.json').read_text());assert sha(DATA/'sequence_bank.csv.gz')==bank_audit['outputs']['sequence_bank.csv.gz']
    with gzip.open(DATA/'sequence_bank.csv.gz','rt') as f:bank_keys=[r['variant_GRCh38'] for r in csv.DictReader(f)]
    test=list(csv.DictReader((DATA/'reserved_physical_design.csv').open()));test_keys=[r['variant_GRCh38'] for r in test];wanted=set(bank_keys)|set(test_keys)
    index=SOURCES/'AVI_SNVS.tsv.gz.tbi'
    if not index.exists():fetch(URL+'.tbi',index,10_000_000)
    assert hashlib.md5(index.read_bytes()).hexdigest()==INDEX_MD5,'AVI source index changed; investigate before reusing old scores'
    db=sqlite3.connect(DATA/'avi_exact_scores.sqlite');db.execute('PRAGMA journal_mode=WAL');db.execute('CREATE TABLE IF NOT EXISTS scores (variant TEXT PRIMARY KEY, raw REAL, phred REAL, available INTEGER NOT NULL, source TEXT NOT NULL)');db.execute('CREATE TABLE IF NOT EXISTS completed_queries (chrom TEXT, window INTEGER, signature TEXT, PRIMARY KEY(chrom,window))')
    prior=ROOT/'outputs/sleep_frailty_ai_v3/avi_exact_allele_scores.csv';receipt=json.loads((prior.parent/'avi_retrieval_audit.json').read_text());assert sha(prior)==receipt['output_sha256'] and receipt['source']==URL and receipt['index_md5']==INDEX_MD5
    imported=0
    with prior.open() as f:
        for r in csv.DictReader(f):
            if r['variant_GRCh38'] in wanted and r['status']=='exact_allele_match':db.execute('INSERT OR IGNORE INTO scores VALUES (?,?,?,?,?)',(r['variant_GRCh38'],float(r['avi_raw']),float(r['avi_phred']),1,'byte_verified_V3_source'));imported+=1
    db.commit();existing={r[0] for r in db.execute('SELECT variant FROM scores')};groups=collections_defaultdict_set()
    for key in sorted(wanted-existing):
        c,p,a,b=key.split(':');p=int(p);assert len(a)==len(b)==1 and a!=b;groups[(c,(p-1)//16384)].add((p,a,b))
    local=threading.local()
    def query(item):
        (c,w),keys=item;keys=set(keys);last=None
        for attempt in range(3):
            try:
                if not hasattr(local,'tab'):local.tab=pysam.TabixFile(URL,index=str(index))
                found={}
                for line in local.tab.fetch(c,min(k[0] for k in keys)-1,max(k[0] for k in keys)):
                    f=line.split('\t');key=(int(f[1]),f[2],f[3])
                    if key not in keys:continue
                    value=(float(f[4]),float(f[5]));assert np.isfinite(value).all() and value[1]>=0;assert key not in found or found[key]==value;found[key]=value
                rows=[]
                for p,a,b in sorted(keys):
                    v=found.get((p,a,b));rows.append((f'{c}:{p}:{a}:{b}',v[0] if v else None,v[1] if v else None,int(v is not None),'current_exact_public_AVI_query'))
                return c,w,rows
            except Exception as e:
                last=repr(e)
                if hasattr(local,'tab'):
                    try:local.tab.close()
                    except Exception:pass
                    del local.tab
                time.sleep(attempt+1)
        raise RuntimeError(str((c,w))+': '+str(last))
    # Validate the old score cache against fresh same-source queries, using
    # already observed development variants and no experimental outcomes.
    probes=list(db.execute('SELECT variant,raw,phred FROM scores WHERE available=1 LIMIT 8'));probe_errors=[]
    for key,raw,phred in probes:
        c,p,a,b=key.split(':');p=int(p);_,_,rows=query(((c,(p-1)//16384),{(p,a,b)}));r=rows[0];assert r[3]==1;err=max(abs(raw-r[1]),abs(phred-r[2]));assert err<1e-9;probe_errors.append(err)
    write_json(OUT/'AVI_collection_plan.json',{'source':URL,'index_md5':INDEX_MD5,'source_index_sha256':sha(index),'requested_unique_variants':len(wanted),'bank_sequence_rows':len(bank_keys),'reserved_design_rows':len(test_keys),'reused_prior_rows':imported,'remaining_query_windows':len(groups),'maximum_parallel_queries':4,'source_cache_current_query_max_errors':probe_errors,'permissive_designation_receipt_sha256':sha(OUT/'AVI_permissive_artifact_receipt.json'),'prior_source_receipt_sha256':sha(prior.parent/'avi_retrieval_audit.json'),'scope':'Only permissible downloaded composite AVI SNV input values. No RNA effects are training targets here and no detailed API tracks/feature importance derivatives are used. Missing scores are unavailable, never zero effects. Native-context AVI does not prove effect in an engineered reporter.'})
    failures=[];done=0;started=time.time();signature=sha(DATA/'sequence_bank.csv.gz')+'_'+sha(DATA/'reserved_physical_design.csv')
    with concurrent.futures.ThreadPoolExecutor(4) as pool:
        futures={pool.submit(query,item):item[0] for item in groups.items()}
        for future in concurrent.futures.as_completed(futures):
            try:
                c,w,rows=future.result();db.executemany('INSERT OR REPLACE INTO scores VALUES (?,?,?,?,?)',rows);db.execute('INSERT OR REPLACE INTO completed_queries VALUES (?,?,?)',(c,w,signature))
            except Exception as e:failures.append({'query':list(futures[future]),'error':repr(e)})
            done+=1
            if done%25==0 or done==len(groups):db.commit()
            if done%200==0 or done==len(groups):
                progress={'done_windows':done,'total_windows':len(groups),'failures':failures,'seconds':time.time()-started,'signature':signature};write_json(OUT/'AVI_collection_progress.json',progress);print('V7 exact AVI windows',done,'/',len(groups),'failures',len(failures),flush=True)
    db.commit();assert not failures,failures
    values={key:(raw,phred,available) for key,raw,phred,available in db.execute('SELECT variant,raw,phred,available FROM scores')};assert wanted<=set(values)
    for name,keys in [('avi_bank.npy',bank_keys),('avi_reserved.npy',test_keys)]:
        x=np.array([[values[k][0] if values[k][2] else np.nan,values[k][1] if values[k][2] else np.nan,values[k][2]] for k in keys],dtype='float32');np.save(DATA/name,x)
    db.close();write_json(OUT/'AVI_collection_receipt.json',{'requested_unique_variants':len(wanted),'available':sum(values[k][2] for k in wanted),'missing':sum(not values[k][2] for k in wanted),'failures':failures,'source':URL,'index_md5':INDEX_MD5,'whole_remote_data_table_MD5_not_independently_verified':True,'outputs':{n:sha(DATA/n) for n in ['avi_bank.npy','avi_reserved.npy']},'source_plan_sha256':sha(OUT/'AVI_collection_plan.json'),'matching':'Exact chromosome,1basedGRCh38position,genomeREF andALT, not position-only or reference-allele inference.','purpose':'Matched added-input arm; not experimental truth, predicted gene-expression direction or clinical risk.'});print('AVI MATCHED INPUT COLLECTION COMPLETE',len(wanted),flush=True)
def collections_defaultdict_set():
    import collections
    return collections.defaultdict(set)
if __name__=='__main__':main()
