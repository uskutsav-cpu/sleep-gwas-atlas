"""Sequential, resumable clump -> intact locus -> signed LD -> SuSiE -> coloc.

A missing/failing locus is accounted for explicitly. Successful siblings can be
preserved, but a failed or underpowered unit is never promoted as a null result.
This orchestration has only software tests here; the native engines require the
user's verified files, LD reference, R packages and PLINK installation.
"""
from __future__ import annotations
from pathlib import Path
from .io import (ContractError,require,read_json,read_tsv,write_json,write_tsv,
                 file_record,sha256,json_hash,safe_id,check_hash)
from .artifacts import verify_artifact,verify_current_artifact,transaction
from .pairs import verify_lock
from .pipeline import create_job
from .executor import run_job
from .ld import map_loci,prepare_locus


def _settings_file(path,value):
    if path.exists():require(read_json(path)==value,'Changed native settings; start a new reviewed run')
    else:write_json(path,value)
    return path


def run_followup(settings_path, *, allow_synthetic=False):
    c=read_json(settings_path)
    require(c.get('reviewed') is True,'Deep follow-up requires reviewed native settings')
    lock=verify_lock(c['pair_lock']);synthetic=lock['synthetic']
    require(not synthetic or allow_synthetic,'Synthetic run requires explicit opt-in')
    require(c['pair_id'] in {p['pair_id'] for p in lock['pairs']},'Pair is not selected')
    require(c.get('method_compatibility_reviewed') is True,'Verify the actual discovery method, including legacy PLACO vs PLACO+')
    index=Path(c['pair_index']).resolve();collated=Path(c['collated']).resolve()
    ir=verify_artifact(index);cr=verify_artifact(collated)
    require(ir['synthetic']==cr['synthetic']==synthetic,'Synthetic data contamination')
    require(cr.get('scientific_status')=='PASS','Cannot follow up rejected pleiotropy outputs')
    require(c.get('discovery_method') in {'PLACO','PLACO_PLUS'},'Record actual discovery method')
    require(c.get('ancestry')=='EUR' and c.get('genome_build') in {'GRCh37','GRCh38'},'Explicit reference metadata required')
    require(0<=c['maximum_locus_failure_rate']<=1,'Invalid failure-rate cap')
    root=Path(c['output_root']).resolve();root.mkdir(parents=True,exist_ok=True)
    specs=root/'job_specs';specs.mkdir(exist_ok=True)
    source_files=[settings_path,c['pair_lock'],index/'receipt.json',collated/'receipt.json',c['ld_blocks']]
    source_files += [c['reference_prefix']+s for s in ['.bed','.bim','.fam']]
    signature={'settings':file_record(settings_path),'inputs':[file_record(p) for p in source_files]}
    binding=root/'followup.lock.json'
    _settings_file(binding,signature)
    if (root/'summary').exists():return verify_current_artifact(root/'summary')
    def native(name,method,settings,outputs):
        sp=_settings_file(specs/(name+'.settings.json'),settings)
        jp=create_job(sp,method=method,job_id=name,inputs={'reviewed_followup':settings_path},
             outputs={**outputs,'semantic_status':'status.json'},root=specs,reviewed=True,
             synthetic=synthetic,minimum_free_bytes=c.get('minimum_free_bytes',2**30),threads=1)
        return run_job(jp,root,allow_synthetic=synthetic)
    # New-family threshold; the original legacy output and threshold are untouched.
    clump={**c['clumping'],'mode':'clump','reference_prefix':c['reference_prefix'],
           'associations':str(collated/'results.tsv.gz'),'input_receipt':str(collated),
           'p_column':'P_PLACO','p1':lock['genomewide_family_threshold'],'threads':1}
    native('clump','clump',clump,{'leads':'leads.tsv'})
    loci_path=root/'loci.tsv'
    loci=map_loci(root/'clump/leads.tsv',c['ld_blocks'],None)
    if loci_path.exists():
        require(list(read_tsv(loci_path))==[{k:str(v) for k,v in r.items()} for r in loci],
                'Derived locus table was modified')
    else:
        write_tsv(loci_path,['locus_id','CHR','START','STOP','lead_snp','lead_p','clump_leads','n_clump_leads'],loci)
    summaries=[];failures=0
    for locus in loci:
        lid=safe_id(locus['locus_id']);state={'pair_id':c['pair_id'],'locus_id':lid,
           'status':'UNRESOLVED','trait_h4':'NA','pip_max':'NA','n_signal_comparisons':0,
           'method':c['discovery_method'],'reason':''}
        try:
            ld={'mode':'ld','plink':c['clumping']['plink'],'reference_prefix':c['reference_prefix'],
                'reference_ancestry':c['ancestry'],'genome_build':c['genome_build'],
                'chromosome':int(locus['CHR']),'start':int(locus['START']),'stop':int(locus['STOP']),
                'max_working_bytes':c['max_working_bytes'],'threads':1}
            native(lid+'_ld','ld',ld,{'matrix':'signed.ld','bim':'locus.bim'})
            lm={'kind':'signed_r','counted_allele':'BIM_A1','synthetic':synthetic,
                'ancestry':c['ancestry'],'genome_build':c['genome_build'],
                'matrix_path':str(root/(lid+'_ld')/'signed.ld'),
                'bim_path':str(root/(lid+'_ld')/'locus.bim')}
            lm['matrix_sha256']=sha256(lm['matrix_path']);lm['bim_sha256']=sha256(lm['bim_path'])
            lmp=_settings_file(specs/(lid+'.ld.json'),lm)
            prepared=root/(lid+'_aligned')
            if prepared.exists():verify_current_artifact(prepared)
            else:prepare_locus(index,lmp,locus,root,lid+'_aligned',synthetic=synthetic,
                 max_working_bytes=c['max_working_bytes'],min_coverage=c['minimum_ld_coverage'])
            fits=[]
            for number in [1,2]:
                sn=f'{lid}_susie{number}'
                sc={**c['susie_common'],**c['trait'+str(number)],'trait_number':number,
                    'locus_file':str(prepared/'locus.tsv'),'ld_file':str(prepared/'ld.tsv')}
                native(sn,'susie',sc,{'fit':'fit.rds','pip':'pip.tsv','cs':'credible_sets.tsv','diagnostics':'diagnostics.json'})
                fits.append(root/sn/'fit.rds')
            co={**c['coloc'],'fit1':str(fits[0]),'fit2':str(fits[1])}
            receipt=native(lid+'_coloc','coloc',co,{'coloc':'coloc.tsv','fit':'coloc.rds'})
            state['status']=receipt['scientific_status']
            ps=[float(r['PIP']) for n in [1,2] for r in read_tsv(root/f'{lid}_susie{n}/pip.tsv')]
            state['pip_max']=max(ps) if ps else 'NA'
            if state['status']=='PASS':
                rs=list(read_tsv(root/(lid+'_coloc')/'coloc.tsv'))
                state['trait_h4']=max(float(r['PP.H4.abf']) for r in rs)
                state['n_signal_comparisons']=len(rs)
                state['reason']='Maximum H4 is descriptive only; inspect ALL signal pairs and prior sensitivity'
            else:state['reason']='Insufficient signal-level evidence; not a negative result'
        except (ContractError,FileNotFoundError,KeyError,ValueError) as exc:
            failures+=1;state['status']='FAILED_OR_BLOCKED';state['reason']=str(exc)
        summaries.append(state)
    status='NO_SIGNAL' if not loci else ('FAILED_QC_NOT_CONSUMED' if failures/len(loci)>c['maximum_locus_failure_rate']
            else 'INSUFFICIENT_EVIDENCE' if any(s['status']!='PASS' for s in summaries) else 'PASS')
    # Summary completion describes evidence accounting, not a mechanism finding.
    with transaction(root,'summary',stage='deep_followup_summary',inputs=source_files+[loci_path],
                     parameters={'pair_id':c['pair_id']},synthetic=synthetic) as (work,meta):
        columns=['pair_id','locus_id','status','trait_h4','pip_max','n_signal_comparisons','method','reason']
        write_tsv(work/'loci_evidence.tsv',columns,summaries)
        meta['scientific_status']=status
        write_json(work/'summary.json',{'status':status,'loci':len(loci),'failures':failures,
            'insufficient':sum(s['status']=='INSUFFICIENT_EVIDENCE' for s in summaries),
            'interpretation':'No causal-mechanism claim. Maximum H4 is a descriptive signal-prioritization summary, not multiplicity-adjusted evidence.'})
    return verify_current_artifact(root/'summary')
