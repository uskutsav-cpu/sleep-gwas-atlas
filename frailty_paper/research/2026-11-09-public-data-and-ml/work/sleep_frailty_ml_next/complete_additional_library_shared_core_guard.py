"""Connect exact shared DNA before exclusions using compact, collision-checked hashes."""
import sys,time,hashlib,collections
from pathlib import Path
ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT/'work/sleep_frailty_ml_v9'))
from v9_common import *

SCREEN=ROOT/'outputs/sleep_frailty_ai_v9_validation/next_source_screen'
DEST=SCREEN/'Siraj_unused_libraries'


def main():
    started=time.monotonic();source=DEST/'full_reference_verified_with_spatial_roles.csv';rows=list(csv_rows(source));parents=list(range(max(int(r['component_root']) for r in rows)+1));case_root=np.array([int(r['component_root']) for r in rows],dtype='int32')
    def find(i):
        while parents[i]!=i:parents[i]=parents[parents[i]];i=parents[i]
        return i
    def union(a,b):
        a,b=find(a),find(b)
        if a!=b:parents[max(a,b)]=min(a,b)
    def core_at(record):
        r=rows[int(record['case'])];field='alternate' if record['allele'] else 'reference';k=int(record['start']);s=r[field][k:k+148]
        return rc(s) if record['reverse'] else s
    count=sum(4*(min(int(r['changed_position']),len(r['reference'])-148)-max(0,int(r['changed_position'])-147)+1) for r in rows)
    dtype=[('h','<u8'),('case','<u4'),('start','u1'),('allele','u1'),('reverse','u1')];records=np.empty(count,dtype=dtype);prefixes=set();cursor=0
    for i,r in enumerate(rows):
        j=int(r['changed_position']);lo=max(0,j-147);hi=min(j,len(r['reference'])-148)
        for allele,field in enumerate(['reference','alternate']):
            for k in range(lo,hi+1):
                value=r[field][k:k+148]
                for reverse,s in enumerate([value,rc(value)]):
                    h=int.from_bytes(hashlib.blake2b(s.encode(),digest_size=8).digest(),'big');records[cursor]=(h,i,k,allele,reverse);prefixes.add(s[:12]);cursor+=1
    assert cursor==count;order=np.argsort(records['h'],kind='stable');records=records[order];del order;keys=np.ascontiguousarray(records['h'])
    boundary=(keys[1:]==keys[:-1])&(case_root[records['case'][1:]]!=case_root[records['case'][:-1]]);candidates=np.unique(keys[1:][boundary]);joins=0
    for h in candidates:
        lo,hi=np.searchsorted(keys,[h,h],side='left')[0],np.searchsorted(keys,h,side='right');seen={}
        for record in records[lo:hi]:
            s=core_at(record);root=int(case_root[int(record['case'])])
            if s in seen:
                if find(root)!=find(seen[s]):joins+=1;union(root,seen[s])
            else:seen[s]=root
    print('NEW LIBRARY SHARED-CORE GRAPH',count,'cores',joins,'cross-region joins',flush=True)
    histories=[V6/'sequence_bank.csv.gz',V7/'sequence_bank.csv.gz',DATA/'extra_sequence_bank.csv.gz',V8/'matoba_full_verified_physical_design.csv',SCREEN/'Tardaguila_full_reference_verified_design.csv'];shared_bad=set();history_summary=[];cached={}
    def flush(batch):
        if not batch:return 0
        hashes=np.array([int.from_bytes(hashlib.blake2b(s.encode(),digest_size=8).digest(),'big') for s in batch],dtype='uint64');positions=np.searchsorted(keys,hashes);valid=positions<len(keys);ids=np.flatnonzero(valid);ids=ids[keys[positions[ids]]==hashes[ids]];matches=0
        for i in ids:
            h=int(hashes[i]);s=batch[i];cachekey=(h,s)
            if cachekey not in cached:
                lo=np.searchsorted(keys,np.uint64(h),side='left');hi=np.searchsorted(keys,np.uint64(h),side='right');found=set()
                for record in records[lo:hi]:
                    if core_at(record)==s:found.add(find(int(case_root[int(record['case'])])))
                cached[cachekey]=found
            found=cached[cachekey]
            if found:shared_bad.update(found);matches+=1
        return matches
    for path in histories:
        batch=[];windows=0;matches=0
        for r in csv_rows(path):
            j=int(r['changed_position']);length=len(r['reference'])
            if length<148:continue
            lo=max(0,j-147);hi=min(j,length-148)
            for field in ['reference','alternate']:
                s=r[field]
                for k in range(lo,hi+1):
                    windows+=1
                    if s[k:k+12] in prefixes:batch.append(s[k:k+148])
                    if len(batch)>=8192:matches+=flush(batch);batch=[]
        matches+=flush(batch);history_summary.append({'path':str(path),'sha256':sha(path),'variant_containing_148bp_windows_examined':windows,'literal_sequence_verified_match_windows':matches});print('NEW LIBRARY HISTORICAL CORE SCAN',path.parent.name,windows,matches,flush=True)
    spatial_bad={find(int(r['component_root'])) for r in rows if r['after_spatial_exclusion']=='False'};reserved={find(int(r['component_root'])) for r in rows if r['library']=='OL33'};output=[]
    for r in rows:
        root=find(int(r['component_root']));why=[]
        if root in spatial_bad:why.append('spatial_history_overlap')
        if root in shared_bad:why.append('literal_shared_148bp_history')
        role='excluded_supervised_history' if why else 'reserved_candidate' if root in reserved else 'development_candidate'
        output.append({**r,'joint_component_after_exact_shared_DNA':root,'phase_after_supervised_history_and_shared_core':role,'supervised_history_exclusion_reason':'|'.join(why),'pretraining_and_assay_audits_still_required':True})
    write_csv(DEST/'full_reference_and_shared_core_roles.csv',output)
    summary={'passed':True,'input_designs':len(rows),'candidate_actual_core_records':count,'hash_algorithm':'BLAKE2b 64-bit for indexing only; every hash match is verified by full 148bp literal DNA before joining or exclusion. Both orientations included.',
             'cross_spatial_component_DNA_joins':joins,'histories':history_summary,'by_library':{lib:dict(collections.Counter(r['phase_after_supervised_history_and_shared_core'] for r in output if r['library']==lib)) for lib in ['OL27','OL33']},
             'shared_contexts_connected_before_exclusions':True,'all_OL33_components_held_out_of_development':True,'new_source_outcomes_or_model_performance_used':False,'final_new_test_eligibility_complete':False,
             'remaining':['Pretraining exposure and source assay/target availability audits','Freeze source roles and next-model/uncertainty protocol before performance evaluation'],'elapsed_seconds':time.monotonic()-started}
    write_json(DEST/'shared_core_and_historical_overlap_audit.json',summary);assert_previous_sources();print(json.dumps(summary),flush=True)


if __name__=='__main__':main()
