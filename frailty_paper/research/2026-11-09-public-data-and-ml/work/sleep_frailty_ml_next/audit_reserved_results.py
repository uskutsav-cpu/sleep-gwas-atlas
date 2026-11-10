"""Separate source parser, rank calculation and cluster replay of the fresh result."""
import sys,time
from pathlib import Path
ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT/'work/sleep_frailty_ml_v9'))
from v9_common import *
DEST=ROOT/'outputs/sleep_frailty_ml_next';FINAL=DEST/'validation';SCREEN=ROOT/'outputs/sleep_frailty_ai_v9_validation/next_source_screen'

def table(path,delimiter):
    opening=gzip.open if str(path).endswith('.gz') else open
    with opening(path,'rt') as handle:
        reader=csv.reader(handle,delimiter=delimiter);columns=next(reader);index={name:i for i,name in enumerate(columns)}
        for row in reader:yield index,row

def ranks(v):
    order=np.argsort(v,kind='stable');z=v[order];starts=np.r_[0,np.flatnonzero(z[1:]!=z[:-1])+1];counts=np.diff(np.r_[starts,len(v)]);out=np.empty(len(v));out[order]=np.repeat(starts+(counts-1)*.5,counts);return out

def rank_correlation(a,b):
    if len(a)<3 or np.ptp(a)==0 or np.ptp(b)==0:return None
    x,y=ranks(a),ranks(b);x-=x.mean();y-=y.mean();return float(np.dot(x,y)/np.sqrt(np.dot(x,x)*np.dot(y,y)))

def measure(y,p,lbase,cbase,line):
    errors=[];lb=[];cb=[];rho=[]
    for label in sorted(set(line)):
        mask=line==label;errors.append(sum(abs(y[mask]-p[mask]))/mask.sum());lb.append(sum(abs(y[mask]-lbase[mask]))/mask.sum());cb.append(sum(abs(y[mask]-cbase[mask]))/mask.sum());r=rank_correlation(y[mask],p[mask])
        if r is not None:rho.append(r)
    error=float(np.mean(errors));return error,1-error/float(np.mean(lb)),1-error/float(np.mean(cb)),float(np.mean(rho)) if rho else None

def main():
    start=time.monotonic();assert_previous_sources();result=json.loads((FINAL/'reserved_outcome_evaluation.json').read_text());before=FINAL/'before_reserved_outcomes_lock.json';assert sha(before)==result['before_outcome_lock_sha256'];lock=json.loads(before.read_text())
    for path,h in lock['files'].items():assert sha(path)==h,path
    truths={}
    for columns,row in table(SCREEN/'Tardaguila2026_TableS3.tsv','\t'):
        get=lambda name:row[columns[name]];truths[('Tardaguila2026',get('SNP'),get('comparison'),get('window'),get('Cell_Type'))]=abs(float(get('log2Skew')))
    retrieval=json.loads((SCREEN/'Siraj_unused_libraries/targeted_retrieval_receipt.json').read_text())
    for item in retrieval['files']:
        for columns,row in table(ROOT/item['path'],'\t'):
            get=lambda name:row[columns[name]]
            if get('comb')==get('SNP')+'_center_fwd_ref':
                try:value=abs(float(get('Log2Skew')))
                except ValueError:continue
                truths[('Siraj2026',item['library'],'chr'+get('SNP').removeprefix('chr'),item['cell'])]=value
    rows=list(csv_rows(DEST/'reserved_observation_metadata.csv.gz'));physical={int(r['bank_id']):r for r in csv_rows(DEST/'reserved_physical_sequence_bank.csv.gz')};joined=list(csv_rows(FINAL/'reserved_joined_measurements_and_predictions.csv.gz'));assert len(rows)==len(joined)
    y=[];clusters=[]
    for i,r in enumerate(rows):
        key=('Tardaguila2026',r['source_SNP'],r['comparison'],r['window'],r['source_cell']) if r['source']=='Tardaguila2026' else ('Siraj2026',r['library'],r['original_variant_GRCh37'],r['source_cell']);value=truths[key];assert abs(float(joined[i]['absolute_original_Log2Skew'])-value)<1e-12;y.append(value);s=physical[int(r['bank_index'])];cluster=s['chromosome']+':'+str((int(s['position_hg38'])-1)//1000000);assert joined[i]['genomic_1Mb_window']==cluster;clusters.append(cluster)
    y=np.array(y);clusters=np.array(clusters);lines=np.array([r['cell_line'] for r in rows]);source=np.array([r['source'] for r in rows]);library=np.array([r['library'] for r in rows]);masks={'blood_K562':(source=='Tardaguila2026')&(lines=='K562'),'OL33_three_known_lines':(source=='Siraj2026')&(library=='OL33')};predictions={};audited={};protocol=json.loads((FINAL/'final_settings_lock.json').read_text())['new_test_protocol']
    for arm in ['sequence','AVI']:
        with np.load(FINAL/(arm+'_reserved_predictions.npz')) as z:predictions[arm]={k:z[k] for k in z.files}
        p=predictions[arm];unknown=np.array([r['exact_line_known']!='True' for r in rows]);assert not p['pretest_policy_accepted'][unknown].any();audited[arm]={}
        for name,mask in masks.items():
            ix=np.flatnonzero(mask);yy=y[ix];pp=p['magnitude'][ix];lb=p['baseline'][ix];cb=p['cell_baseline'][ix];ss=lines[ix];gg=clusters[ix];actual=measure(yy,pp,lb,cb,ss);reported=result['primary_accuracy'][arm][name]
            for value,key in zip(actual,['macro_MAE','relative_MAE_gain_vs_line_median','relative_MAE_gain_vs_cell_median','macro_within_cell_line_Spearman']):
                if value is None:assert reported[key] is None
                else:assert abs(value-reported[key])<1e-7,(arm,name,key,value,reported[key])
            coverage=float(np.mean(abs(yy-pp)<=p['interval_halfwidth'][ix]));assert abs(coverage-result['confidence'][arm][name]['observed_coverage'])<1e-12
            # Independently replay the fixed-window paired bootstrap, including ranks.
            windows=sorted(set(gg));parts=[np.flatnonzero(gg==w) for w in windows];rng=np.random.default_rng(protocol['bootstrap_seed']);samples=[[],[],[],[]]
            for iteration in range(protocol['bootstrap_replicates']):
                indices=np.concatenate([parts[j] for j in rng.integers(len(parts),size=len(parts))]);v=measure(yy[indices],pp[indices],lb[indices],cb[indices],ss[indices]);samples[0].append(v[1]);samples[1].append(v[2]);samples[2].append(v[3])
                if arm=='AVI':
                    reference=predictions['sequence']['magnitude'][ix][indices];base_error=measure(yy[indices],reference,lb[indices],cb[indices],ss[indices])[0];samples[3].append(1-v[0]/base_error)
            errors=[]
            for values,key in zip(samples,['gain_vs_line_95CI','gain_vs_cell_95CI','macro_within_line_Spearman_95CI','gain_vs_sequence_95CI']):
                if not values:continue
                values=[v for v in values if v is not None]
                if len(values)<100:assert reported[key] is None;continue
                ci=np.quantile(values,[.025,.975]);error=float(np.max(abs(ci-np.asarray(reported[key]))));assert error<1e-7,(arm,name,key,error);errors.append(error)
            audited[arm][name]={'source_targets_reconstructed':len(ix),'macro_errors_and_ranks_recomputed':True,'coverage':coverage,'cluster_CI_maximum_error':max(errors)};print('INDEPENDENT RESERVED RESULT VERIFIED',arm,name,flush=True)
    paired_audits={}
    with np.load(FINAL/'reserved_measured_pair_identities.npz') as z:pair_banks={k:z[k] for k in z.files}
    for arm,p in predictions.items():
        paired_audits[arm]={}
        for name,mask in masks.items():
            paired_audits[arm][name]={}
            for label,pairs in pair_banks.items():
                pairs=pairs[np.all(mask[pairs],axis=1)];reported=result['paired_difference_results'][arm][name][label];assert len(pairs)==reported['pairs']
                if not len(pairs):continue
                a,b=pairs[:,0],pairs[:,1];dy=y[a]-y[b];dp=p['magnitude'][a]-p['magnitude'][b];db=p['baseline'][a]-p['baseline'][b];keys=np.array([lines[i]+'|'+lines[j] for i,j in pairs]);base=float(abs(dy-db).mean());errors=[]
                values={'difference_MAE':float(abs(dy-dp).mean()),'difference_Spearman':rank_correlation(dy,dp),'relative_gain_vs_training_median_difference':1-float(abs(dy-dp).mean())/base if base else None,'ordering_accuracy_nonzero_measured_differences':float(np.mean(np.sign(dy[dy!=0])==np.sign(dp[dy!=0]))) if (dy!=0).any() else None}
                within=[]
                for key in sorted(set(keys)):
                    ix=keys==key
                    if ix.sum()>=30:
                        r=rank_correlation(dy[ix],dp[ix])
                        if r is not None:within.append(r)
                values['macro_within_fixed_cell_pair_Spearman']=float(np.mean(within)) if within else None
                for key,value in values.items():
                    if value is None:assert reported[key] is None
                    else:assert abs(value-reported[key])<1e-7,(arm,name,label,key)
                paired_audits[arm][name][label]={'pairs':len(pairs),'point_errors_ordering_and_ranks_independently_recomputed':True}
    write_json(FINAL/'reserved_results_independent_audit.json',{'passed':True,'all_source_targets_reconstructed':len(y),'primary_results':audited,'paired_difference_point_results':paired_audits,'pair_CIs_independently_replayed':False,'all_frozen_pretest_hashes_unchanged':True,'unknown_lines_pretest_accepted':0,
              'original_result_sha256':sha(FINAL/'reserved_outcome_evaluation.json'),'code_sha256':sha(Path(__file__)),'scope':'Reconstruction and computational audit; does not improve adverse scientific results.','elapsed_seconds':time.monotonic()-start});assert_previous_sources();print('RESERVED RESULTS INDEPENDENT AUDIT COMPLETE',flush=True)

if __name__=='__main__':main()
