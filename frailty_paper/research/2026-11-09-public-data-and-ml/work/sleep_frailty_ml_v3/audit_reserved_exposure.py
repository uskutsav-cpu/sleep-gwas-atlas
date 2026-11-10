"""Conservative prior-project region exclusion before any reserved inference.

Use only coordinates/alleles from previously analyzed molecular studies, never
reserved-model scores or outcome-dependent eligibility. Development fits stay
unchanged. Original full reserved results remain a labelled comparison.
"""
import bisect,collections,datetime
from common import *
def main():
    assert_no_reserved_evaluation();sal=WORK/'salomon_manifest_local.csv';sir=WORK/'siraj_manifest_local.csv';rr=list(csv.DictReader(sir.open()));nsal=sum(1 for _ in csv.DictReader(sal.open()));sources={'previous_Salomon_study':sal,'previous_Retallick_study':ROOT/'work/sleep_frailty_ml_v2/external_verified_all_local.csv'};excluded=set();details={}
    for name,path in sources.items():
        positions=collections.defaultdict(list)
        for r in csv.DictReader(path.open()):positions[int(r['chromosome'].removeprefix('chr'))].append(int(r['position_hg38']))
        for c in positions:positions[c].sort()
        close=[];exact=[]
        for j,r in enumerate(rr):
            if r['split']!='test':continue
            c=int(r['chromosome']);p=int(r['position_hg38']);a=positions[c];k=bisect.bisect_left(a,p-1000);l=bisect.bisect_left(a,p)
            if k<len(a) and a[k]<=p+1000:close.append(j)
            if l<len(a) and a[l]==p:exact.append(j)
        excluded.update(close);details[name]={'source_file':str(path),'sha256':sha(path),'within_1000bp_count':len(close),'same_coordinate_count':len(exact),'within_1000bp_Siraj_local_indices':close}
    original=[i for i,r in enumerate(rr) if r['split']=='test'];retained=sorted(set(original)-excluded);output={'created_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'source_manifest_sha256':{'salomon':sha(sal),'siraj':sha(sir)},'sources':details,'rule':'Exclude Siraj reserved rows within1000bp inclusive of any coordinate in either previously analyzed Salomon or Retallick molecular study. Original project queue neighborhood exclusions already applied by preparation.','original_Siraj_reserved_n':len(original),'additional_excluded_n':len(excluded),'strict_Siraj_reserved_n':len(retained),'excluded_Siraj_local_indices':sorted(excluded),'excluded_global_row_indices':[nsal+i for i in sorted(excluded)],'retained_Siraj_global_row_indices':[nsal+i for i in retained],'previously_inspected_project_queue_exclusion':'Already applied in original data preparation, not duplicated here','adaptation_disclosure':'Pre-evaluation refinement after development model fitting had begun. Eligibility uses previous-study coordinates only, before any reserved model predictions. Preserve the original full Siraj set as a labelled comparison. No development rows, model targets, folds or fitted weights changed.','genome_build':'GRCh38 coordinates; prior Retallick all-variant build/reference verification already completed and source manifest hash retained','not_established':'No guarantee of absence from original foundation-model sequence training, study publications or external resources','reserved_labels_used_for_exclusion':False,'reserved_predictions_evaluated':False}
    path=OUT/'reserved_exposure_audit.json';assert not path.exists();path.write_text(json.dumps(output,indent=2)+'\n');print('Strict reserved Siraj set',len(retained),'excluded',len(excluded),'before prediction',flush=True)
if __name__=='__main__':main()
