"""Versioned continuous test after the canonical public SNP output proved missing."""
import subprocess,datetime
from common import *
from scipy.stats import spearmanr
from additional_reserved_checks import WeightedRanks,correlation,summarize
def main():
    original=json.loads((OUT/'final_model_lock.json').read_text())
    for path,h in original['files'].items():assert sha(path)==h,path
    for path,h in original['model_checkpoints'].items():assert sha(path)==h,path
    audit=json.loads((OUT/'kosicki_continuous_recovery_audit.json').read_text());assert audit['source_pvalues_recovered'] is False
    path=OUT/'kosicki_reconstructed_continuous_effects.csv';assert sha(path)==audit['reconstructed_effects_sha256']
    critical=[path,WORK/'ADDITIONAL_VALIDATION_PROTOCOL.txt',WORK/'recover_neuronal_measurements.py',WORK/'export_neuronal_replicates.R',WORK/'evaluate_continuous_recovery.py',WORK/'check_continuous_effects.R',OUT/'kosicki_blinded_predictions.npz',OUT/'kosicki_continuous_recovery_audit.json']
    lock=OUT/'continuous_recovery_evaluation_lock.json'
    if not lock.exists():write_json(lock,{'utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'original_model_lock_sha256':sha(OUT/'final_model_lock.json'),'files':{str(p):sha(p) for p in critical},'model_effect_comparisons_already_performed':False,'qualification':'Recovery endpoint fixed after discovering missing canonical SNP results. Original predictions and model weights were fixed before any measured outcomes. This is explicitly a versioned continuous recovery, not recovered FDR validation.'})
    for p,h in json.loads(lock.read_text())['files'].items():assert sha(p)==h,p
    rows=list(csv.DictReader(path.open()));ix=np.array([int(r['design_index']) for r in rows]);effect=np.array([float(r['absolute_mean_zscore_difference']) for r in rows]);logeffect=np.abs(np.array([float(r['mean_log2_ratio_difference']) for r in rows]));tiles=np.array([r['reference_name'] for r in rows]);blocks=np.array([r['chromosome']+':'+str(int(r['position_hg38'])//1000000) for r in rows])
    with np.load(OUT/'kosicki_blinded_predictions.npz') as z:
        assert z['name'][ix].tolist()==[r['name'] for r in rows];models={k.removeprefix('model::'):z[k][ix] for k in z.files if k.startswith('model::')}
    unique,group=np.unique(blocks,return_inverse=True);nblocks=len(unique);er=WeightedRanks(effect);pr={n:WeightedRanks(p) for n,p in models.items()};local={};results={}
    for name,p in models.items():
        values=[];groups=[]
        for tile in np.unique(tiles):
            ti=np.flatnonzero(tiles==tile)
            if len(ti)>=10 and len(np.unique(effect[ti]))>1 and len(np.unique(p[ti]))>1:values.append(float(spearmanr(effect[ti],p[ti]).statistic));groups.append(int(np.bincount(group[ti],minlength=nblocks).argmax()))
        local[name]=(np.array(values),np.array(groups,dtype='int32'));results[name]={'effect_magnitude_Spearman':float(spearmanr(effect,p).statistic),'within_reference_tile_mean_Spearman':float(np.mean(values)),'within_reference_tile_groups':len(values),'absolute_log2_ratio_Spearman_sensitivity':float(spearmanr(logeffect,p).statistic)}
    draw={n:[] for n in models};ldraw={n:[] for n in models};rng=np.random.default_rng(20261107)
    for b in range(1000):
        counts=np.bincount(rng.integers(0,nblocks,size=nblocks),minlength=nblocks);w=counts[group].astype('float64');ranks=er.ranks(w)
        for n in models:
            value=correlation(ranks,pr[n].ranks(w),w);draw[n].append(value);lv,lg=local[n];lw=counts[lg];ldraw[n].append(float(np.dot(lw,lv)/lw.sum()) if lw.sum() else np.nan)
            if b<3 and n=='primary':
                repeated=np.repeat(np.arange(len(effect)),w.astype('int32'));assert abs(value-float(spearmanr(effect[repeated],models[n][repeated]).statistic))<1e-12
        if b%250==0:print('Continuous frozen test block draws',b,flush=True)
    src=OUT/'kosicki_continuous_R_input.csv';dest=OUT/'kosicki_continuous_R_metrics.csv'
    with src.open('w',newline='') as f:w=csv.writer(f);w.writerow(['row_id','effect','tile',*models]);w.writerows(zip(ix,effect,tiles,*models.values()))
    subprocess.run(['Rscript',str(WORK/'check_continuous_effects.R'),str(src),str(dest)],check=True,env=r_environment())
    checks=[]
    for r in csv.DictReader(dest.open()):
        result=results[r['model']];error=max(abs(float(r['Spearman'])-result['effect_magnitude_Spearman']),abs(float(r['within_tile_Spearman'])-result['within_reference_tile_mean_Spearman']));assert error<1e-12;assert int(r['tiles'])==result['within_reference_tile_groups'];checks.append({'model':r['model'],'two_metrics_max_error':error})
    paired={}
    for a,b in [('primary',n) for n in models if n!='primary']+[('contrast_hard','contrast_standard'),('contrast_standard','contrast_neural_only')]:paired[a+' minus '+b]={'global_Spearman':summarize(np.array(draw[a])-np.array(draw[b])),'within_tile_mean_Spearman':summarize(np.array(ldraw[a])-np.array(ldraw[b]))}
    report={'study':'Kosicki2025','eligible_pairs':len(rows),'unique_genomic_contrasts':len({r['variant_GRCh38'] for r in rows}),'reference_tiles':len(set(tiles)),'genomic_1Mb_blocks':nblocks,'models':results,'block_intervals':{n:{'global_Spearman':summarize(draw[n]),'within_tile_mean_Spearman':summarize(ldraw[n])} for n in models},'paired_differences':paired,'independent_R_checks':checks,'recovery_lock_sha256':sha(lock),'original_model_lock_sha256':sha(OUT/'final_model_lock.json'),'source_audit':'kosicki_continuous_recovery_audit.json','canonical_binary_labels_available':False,'uncertainty_label_coverage_validated':False,'scope':'Recovered actual continuous experimental effects. No fabricated p-values or binary labels. One published supplementary tile name differs; exact genome variant, physical sequence and rounded measured effect agree. Both neuronal resources use WTC11-derived NGN2 neurons; independent donors are not established. No clinical/genetic-score/causal validation.'}
    write_json(OUT/'kosicki_continuous_results.json',report);print('CONTINUOUS RECOVERY TEST AND INDEPENDENT R CHECKS COMPLETE',json.dumps(results['primary']),flush=True)
if __name__=='__main__':main()
