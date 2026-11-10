"""Recover actual continuous effects, never invent missing original p-values."""
import subprocess,collections,openpyxl
from common import *
def main():
    lock=json.loads((OUT/'final_model_lock.json').read_text())
    for path,h in lock['files'].items():assert sha(path)==h,path
    assert (OUT/'final_test_evaluation_started.json').exists()
    p=DATA/'reserved_outcome_source/tested_counts.rds';entry=next(r for r in json.loads((SOURCES/'kosicki_sealed_outcome_manifest.json').read_text()) if r['member']=='mpravista/data/input/tested_counts.rds');assert sha(p)==entry['sha256'];export=DATA/'reserved_outcome_source/published_replicate_readouts.csv'
    if not export.exists():subprocess.run(['Rscript',str(WORK/'export_neuronal_replicates.R'),str(p),str(export),str(OUT/'final_model_lock.json')],check=True,env=r_environment())
    byname=collections.defaultdict(dict)
    with export.open() as f:
        for r in csv.DictReader(f):
            rep=int(r['rep']);assert rep in [1,2,3] and rep not in byname[r['name']];values=[float(r[k]) for k in ['zscore_neg','log2.ratio','BC_count']];assert np.isfinite(values).all();byname[r['name']][rep]=values
    reference=list(csv.DictReader((DATA/'reserved_outcome_source/final_tests_atomic.csv').open()));ref_errors=[]
    for r in reference:
        if r['label_class']!='reference' or r['keep']!='TRUE':continue
        assert r['name'] in byname and set(byname[r['name']])=={1,2,3};mean=np.mean([v[0] for v in byname[r['name']].values()]);ref_errors.append(abs(mean-float(r['zscore_neg'])))
    assert ref_errors and max(ref_errors)<1e-10,max(ref_errors)
    design=list(csv.DictReader((DATA/'kosicki_blind_design.csv').open()));kept=[];excluded=collections.Counter()
    for i,r in enumerate(design):
        ref=byname.get(r['reference_name']);alt=byname.get(r['name'])
        if ref is None or alt is None:excluded['missing_reference_or_alternate_measurement']+=1;continue
        if set(ref)!={1,2,3} or set(alt)!={1,2,3}:excluded['not_all3replicates']+=1;continue
        if min(v[2] for v in list(ref.values())+list(alt.values()))<15:excluded['reference_or_alternate_fewer_than15barcodes']+=1;continue
        delta=np.array([alt[j][0]-ref[j][0] for j in [1,2,3]]);logdelta=np.array([alt[j][1]-ref[j][1] for j in [1,2,3]]);mean=float(delta.mean());kept.append({'design_index':i,**r,'mean_zscore_difference':mean,'absolute_mean_zscore_difference':abs(mean),'mean_log2_ratio_difference':float(logdelta.mean()),'rep1_zscore_difference':float(delta[0]),'rep2_zscore_difference':float(delta[1]),'rep3_zscore_difference':float(delta[2]),'minimum_barcode_count':min(v[2] for v in list(ref.values())+list(alt.values()))})
    book=openpyxl.load_workbook(SOURCES/'recovery_supplement.xlsx',read_only=True,data_only=True);pub=[row for row in book['5_GWAS'].values if row[0] and str(row[0]).startswith('chr')];checks=[]
    for row in pub:
        key=str(row[0]).replace('-',':');tile=str(row[1]);values=[r for r in kept if r['variant_GRCh38']==key];exact=[r for r in values if r['reference_name'].split('|')[0]==tile];values=exact or values;assert len(values)==1,(key,tile,len(values));value=values[0]['mean_zscore_difference'];published=float(row[5]);error=abs(value-published);checks.append({'variant':key,'published_tile':tile,'sequence_verified_tile':values[0]['reference_name'].split('|')[0],'tile_name_matches':tile==values[0]['reference_name'].split('|')[0],'reconstructed':value,'published_rounded':published,'absolute_error':error,'matches_to_published_precision':bool(error<=.005000001)})
    print('Published hit comparisons',checks,flush=True)
    assert len(checks)==7 and all(r['matches_to_published_precision'] for r in checks),checks
    target=OUT/'kosicki_reconstructed_continuous_effects.csv'
    with target.open('w',newline='') as f:w=csv.DictWriter(f,fieldnames=list(kept[0]));w.writeheader();w.writerows(kept)
    audit={'source':'Kosicki2025publishedtested_counts.rds','rows_in_replicate_table':sum(len(v) for v in byname.values()),'QC_passed_reference_means_verified':len(ref_errors),'reference_mean_maximum_error':max(ref_errors),'eligible_paired_sequences':len(kept),'unique_genomic_contrasts':len({r['variant_GRCh38'] for r in kept}),'exclusions':dict(excluded),'published_GWAS_hit_effect_checks':checks,'source_RDS_sha256':sha(p),'replicate_export_sha256':sha(export),'reconstructed_effects_sha256':sha(target),'original_model_lock_sha256':sha(OUT/'final_model_lock.json'),'recovery_protocol_sha256':sha(WORK/'ADDITIONAL_VALIDATION_PROTOCOL.txt'),'source_pvalues_recovered':False,'source_significance_calls_invented':False,'scope':'Actual continuous paired RNA/DNAreadout recovery from3publishedreplicates. Canonical source per-variant p/qmeasurements remain absent. No new model choice, fitting, recalibration or model-performance calculation performed here.'}
    write_json(OUT/'kosicki_continuous_recovery_audit.json',audit);print('CONTINUOUS MEASUREMENTS RECOVERED',len(kept),'pairs;',len(ref_errors),'reference means and7publishedhit effects verified',flush=True)
if __name__=='__main__':main()
