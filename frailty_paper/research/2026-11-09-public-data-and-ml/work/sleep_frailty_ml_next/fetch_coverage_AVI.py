"""Exact-allele, resumable public AVI inputs; never outcome labels or teacher targets."""
import sys,sqlite3,threading,concurrent.futures,time
from pathlib import Path
ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT/'work/sleep_frailty_ml_v9'))
from v9_common import *
sys.path.append(str(ROOT/'work/stage5_env/lib/python3.13/site-packages'))
import pysam
DEST=ROOT/'outputs/sleep_frailty_ml_next'
URL='https://grr.iossifovlab.com/hg38/scores/AVI/alphagenome_variant_impact_score_snvs.tsv.gz'

def main():
    assert_previous_sources();index=V7/'sources/AVI_SNVS.tsv.gz.tbi';assert hashlib.md5(index.read_bytes()).hexdigest()=='911ee04b6e8fba1882465221b5a14f2c'
    wanted=set()
    for name in ['OL27_development_sequence_bank.csv.gz','reserved_physical_sequence_bank.csv.gz']:wanted.update(r['variant_GRCh38'] for r in csv_rows(DEST/name))
    db=sqlite3.connect(DEST/'coverage_AVI_exact_scores.sqlite');db.execute('CREATE TABLE IF NOT EXISTS scores (variant TEXT PRIMARY KEY, raw REAL, phred REAL, available INTEGER, source TEXT)');db.commit()
    present={r[0] for r in db.execute('SELECT variant FROM scores')};copied=0
    for path in [V7/'avi_exact_scores.sqlite',V8/'fresh_AVI_exact_scores.sqlite',DATA/'extra_AVI_exact_scores.sqlite']:
        if not path.exists():continue
        old=sqlite3.connect(path.as_uri()+'?mode=ro',uri=True)
        for key in sorted(wanted-present):
            row=old.execute('SELECT variant,raw,phred,available,source FROM scores WHERE variant=?',(key,)).fetchone()
            if row is not None:db.execute('INSERT INTO scores VALUES (?,?,?,?,?)',row);present.add(key);copied+=1
        old.close()
    db.commit();groups=collections.defaultdict(set)
    for key in wanted-present:
        c,p,a,b=key.split(':');p=int(p);assert a in 'ACGT' and b in 'ACGT' and a!=b;groups[(c,(p-1)//16384)].add((p,a,b))
    local=threading.local()
    def query(item):
        (c,window),keys=item
        for attempt in range(3):
            try:
                if not hasattr(local,'tab'):local.tab=pysam.TabixFile(URL,index=str(index))
                found={}
                for line in local.tab.fetch(c,min(k[0] for k in keys)-1,max(k[0] for k in keys)):
                    f=line.split('\t');key=(int(f[1]),f[2],f[3])
                    if key in keys:
                        value=(float(f[4]),float(f[5]));assert np.isfinite(value).all() and value[1]>=0;assert key not in found or found[key]==value;found[key]=value
                return [(f'{c}:{p}:{a}:{b}',*found.get((p,a,b),(None,None)),int((p,a,b) in found),'exact_public_permissive_AVI_coverage_query') for p,a,b in keys]
            except Exception:
                if hasattr(local,'tab'):local.tab.close();del local.tab
                if attempt==2:raise
    started=time.monotonic();done=0;failures=[]
    print('EXACT AVI INPUT QUERY START',len(wanted),'alleles',len(groups),'remote windows',copied,'copied',flush=True)
    with concurrent.futures.ThreadPoolExecutor(max_workers=2) as pool:
        pending={pool.submit(query,item):item[0] for item in sorted(groups.items())}
        for future in concurrent.futures.as_completed(pending):
            try:db.executemany('INSERT OR REPLACE INTO scores VALUES (?,?,?,?,?)',future.result())
            except Exception as error:failures.append({'window':list(pending[future]),'error':repr(error)})
            done+=1
            if done%25==0 or done==len(groups):
                db.commit();progress={'windows_done':done,'windows_total':len(groups),'failed_windows':len(failures),'cached_alleles':db.execute('SELECT count(*) FROM scores').fetchone()[0],'elapsed_seconds':time.monotonic()-started,'PID':os.getpid()};write_json(DEST/'coverage_AVI_query_progress.json',progress);print('EXACT AVI INPUT QUERY',json.dumps(progress),flush=True)
    db.commit();write_json(DEST/'coverage_AVI_query_attempt_receipt.json',{'windows':len(groups),'failures':failures,'elapsed_seconds':time.monotonic()-started})
    assert not failures,failures[:5]
    rows=[]
    for key in sorted(wanted):
        value=db.execute('SELECT raw,phred,available FROM scores WHERE variant=?',(key,)).fetchone();assert value is not None
        rows.append({'variant_GRCh38':key,'raw':value[0] if value[2] else '','phred':value[1] if value[2] else '','available':value[2]})
    from audit_and_lock_library_measurements import write_gzip_csv
    write_gzip_csv(DEST/'coverage_AVI_comparison_inputs.csv.gz',rows);db.close()
    # Old banks are reused only in their original exact bank order.
    from feature_banks import avi_bank
    old=avi_bank();expansion=json.loads((DEST/'coverage_expansion_receipt.json').read_text());assert old.shape[0]==expansion['appended_bank_start']
    sequence=list(csv_rows(DEST/'OL27_development_sequence_bank.csv.gz'));mapping={r['variant_GRCh38']:r for r in rows};extra=np.full((len(sequence),3),np.nan,dtype='float32')
    for i,r in enumerate(sequence):
        score=mapping[r['variant_GRCh38']];extra[i]=[float(score['raw']) if score['raw'] else np.nan,float(score['phred']) if score['phred'] else np.nan,int(score['available'])]
    np.save(DEST/'coverage_appended_AVI.npy',extra)
    write_json(DEST/'coverage_AVI_input_receipt.json',{'passed':True,'source':URL,'index_MD5':hashlib.md5(index.read_bytes()).hexdigest(),'requested_alleles':len(wanted),'available':sum(r['available'] for r in rows),
              'comparison_inputs_sha256':sha(DEST/'coverage_AVI_comparison_inputs.csv.gz'),'appended_bank_sha256':sha(DEST/'coverage_appended_AVI.npy'),
              'official_permissive_artifact_receipt_sha256':sha(V7/'outputs/AVI_permissive_artifact_receipt.json'),
              'use':'Permissive composite AVI scores as three added inputs only. Not restricted detailed AlphaGenome outputs; not experimental labels. Gains cannot be attributed uniquely to AlphaGenome.',
              'no_outcomes_or_model_performance_used':True,'missing_policy':'Keep missing flag; training-only median fill at inference. Missing is not zero biological effect.',
              'code_sha256':sha(Path(__file__)),'elapsed_seconds':time.monotonic()-started});assert_previous_sources();print('EXACT AVI INPUTS COMPLETE',flush=True)

if __name__=='__main__':main()
