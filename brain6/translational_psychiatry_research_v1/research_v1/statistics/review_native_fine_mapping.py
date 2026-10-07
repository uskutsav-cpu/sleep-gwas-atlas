"""Post-result independent native-model posterior review. No source GWAS rows exported."""
import json, math
import numpy as np
import pandas as pd
from scipy.special import logsumexp
from calibration import HERE, ROOT, sha
from review_parent_analyses import exclusive_pair_logsum

LD=HERE.parent/'ld';LOCAL=ROOT/'work/ld_genotypes_research_v1/C_finemap_inputs'

def main():
    protocol=LD/'FINEMAP_PRE_FIT_PROTOCOL.md'
    expected=(LD/'FINEMAP_PRE_FIT_PROTOCOL.sha256').read_text().split()[0]
    assert sha(protocol)==expected
    manifest=json.loads((LD/'finemap_input_manifest.json').read_text());assert manifest['protocol_sha256']==expected
    for row in manifest['local_only_input_files'].values():assert sha(row['path'])==row['sha256']
    assert sha(ROOT/'work/ld_genotypes_research_v1/chr5_signed_LD.npz')==manifest['LD_full_sha256']
    sample=pd.read_csv(LD/'chr5_sample_manifest.tsv',sep='\t');assert len(sample)==503
    variants=pd.read_csv(LD/'chr5_variants.tsv',sep='\t');assert (variants.n_reference==503).all() and not variants.variant_key.duplicated().any()
    raw=np.load(ROOT/'work/ld_genotypes_research_v1/chr5_signed_LD.npz')
    shared=pd.read_csv(LOCAL/'insomnia_C_shared.tsv',sep='\t')['variant_key'].tolist()
    index={k:i for i,k in enumerate(raw['variant_keys'])};ix=[index[k] for k in shared]
    full=raw['R'];small=np.fromfile(LOCAL/'shared_signed_LD.float64.bin',dtype='<f8').reshape((2185,2185),order='F')
    LDerror=float(np.max(np.abs(full[np.ix_(ix,ix)]-small)));assert LDerror<1e-12
    native=pd.read_csv(LD/'native_coloc_all_signal_pairs.tsv',sep='\t');rows=[];assert len(native)==12
    for model,g in native.groupby('model',sort=False):
        bfpath=LOCAL/f'{model}_native_BFs_LOCAL_ONLY.tsv';bf=pd.read_csv(bfpath,sep='\t');assert bf.variant_key.tolist()==shared
        l1=bf.logBF1.to_numpy();l2=bf.logBF2.to_numpy()
        base=np.array([0,math.log(1e-4)+logsumexp(l1),math.log(1e-4)+logsumexp(l2),math.log(1e-4)+math.log(1e-4)+exclusive_pair_logsum(l1,l2),logsumexp(l1+l2)])
        primary=g[g.p12==1e-5].iloc[0]
        ppprimary=np.array([primary[f'PP.H{i}.abf'] for i in range(5)])
        for r in g.itertuples(index=False,name=None):
            row=dict(zip(g.columns,r));terms=base.copy();terms[4]+=math.log(row['p12']);pp=np.exp(terms-logsumexp(terms))
            prior=ppprimary.copy();prior[4]*=row['p12']/1e-5;prior/=prior.sum()
            reported=np.array([row[f'PP.H{i}.abf'] for i in range(5)])
            err=float(np.max(np.abs(pp-reported)));rerr=float(np.max(np.abs(prior-reported)));total=float(reported.sum())
            ratio=reported[4]/(reported[3]+reported[4])
            assert np.all((reported>=0)&(reported<=1)) and abs(total-1)<1e-12 and err<1e-10 and rerr<1e-10
            rows.append({'model':model,'p12':row['p12'],'n_snps':len(bf),'native_H4':reported[4],'explicit_configuration_H4':pp[4],'all_hypotheses_sum':total,'max_configuration_error':err,'max_prior_reweighting_error':rerr,'conditional_H4_H3H4_ratio_error':abs(ratio-row['H4_over_H3_plus_H4']),'local_BF_sha256':sha(bfpath),'classification':'POST_RESULT_SAME_SOURCE_NUMERICAL_REVIEW','independent_replication':False})
    pd.DataFrame(rows).to_csv(HERE/'native_coloc_configuration_and_prior_checks.tsv',sep='\t',index=False,float_format='%.17g')
    checks=pd.read_csv(HERE/'native_fine_mapping_RDS_checks.tsv',sep='\t');assert len(checks)==8
    summary={'classification':'POST_RESULT_SEPARATE_AUTOMATED_NATIVE_FINE_MAPPING_REVIEW','frozen_protocol_sha256':expected,'native_fits_reviewed_from_RDS':8,'credible_sets_reviewed':8,'signal_pair_posterior_conditions_reviewed':12,'full_to_shared_signed_LD_max_error':LDerror,'max_configuration_posterior_error':max(r['max_configuration_error'] for r in rows),'max_prior_reweighting_error':max(r['max_prior_reweighting_error'] for r in rows),'minimum_conditional_CS_coverage':float(checks.independent_conditional_coverage.min()),'minimum_CS_purity':float(checks.independent_purity_min.min()),'max_PIP_error':float(checks.max_PIP_error.max()),'max_CS_alpha_error':float(checks.max_CS_alpha_error.max()),'max_purity_error':float(checks.max_purity_error.max()),'max_ELBO_roundtrip_error':float(checks.max_ELBO_roundtrip_error.max()),'minimum_ELBO_step':float(checks.min_ELBO_change.min()),'reference_samples':503,'shared_variants':2185,'primary_H4':float(native[(native.model=='PRIMARY_L10_MEDIAN_NEFF')&(native.p12==1e-5)]['PP.H4.abf'].iloc[0]),'primary_low_prior_H4':float(native[(native.model=='PRIMARY_L10_MEDIAN_NEFF')&(native.p12==1e-6)]['PP.H4.abf'].iloc[0]),'independent_two_trait_replications':0,'new_GWAS_fits_run_by_reviewer':0}
    (HERE/'native_fine_mapping_review_summary.json').write_text(json.dumps(summary,indent=2)+'\n');print(json.dumps(summary,indent=2))

if __name__=='__main__':main()
