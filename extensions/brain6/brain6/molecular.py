"""Reviewed molecular-QTL follow-up, preserving every attempted feature/locus.

Runs the native SuSiE/coloc adapters on complete cis-QTL summary statistics.
A high posterior is not itself a validated mechanism or a multiplicity-adjusted
claim. All signal-pair tables and prior sensitivity files remain in native jobs.
"""
from __future__ import annotations
from pathlib import Path
from .io import (read_json, read_tsv, require, write_json, write_tsv,
                 file_record, check_hash, safe_id, ContractError)
from .artifacts import transaction, verify_artifact, verify_current_artifact
from .pairs import verify_lock
from .followup import _settings_file
from .pipeline import create_job
from .executor import run_job
from .qtl import locus as qtl_locus
from .ld import index_pair, prepare_locus


def run(settings_path, *, allow_synthetic=False):
    c=read_json(settings_path)
    require(c.get('reviewed') is True and c.get('feature_family_frozen') is True,
            'Freeze feature/tissue/locus queries before inspecting QTL colocalization')
    lock=verify_lock(c['pair_lock']);synthetic=lock['synthetic']
    require(not synthetic or allow_synthetic,'Synthetic molecular work requires opt-in')
    selected={p['pair_id']:p for p in lock['pairs']}
    queries=c.get('queries',[])
    require(queries and len({q['query_id'] for q in queries})==len(queries),'Empty/duplicate molecular query list')
    # Preflight identity and source binding for every query before execution.
    inputs=[settings_path,c['pair_lock']];seen=set()
    for q in queries:
        safe_id(q['query_id']);require(q['pair_id'] in selected,'Unselected molecular pair')
        require(q['gwas_trait'] in {selected[q['pair_id']]['sleep_trait'],selected[q['pair_id']]['disease_trait']},
                'QTL query GWAS trait must belong to the selected pair')
        require(q.get('feature_id') and q.get('selection_reason'),'Explain feature selection')
        identity=(q['pair_id'],q['gwas_trait'],q['qtl_index'],q['feature_id'],str(q['region']))
        require(identity not in seen,'Duplicated feature/locus query');seen.add(identity)
        gr=verify_artifact(q['gwas']);qr=verify_artifact(q['qtl_index'])
        require(gr['synthetic']==qr['synthetic']==synthetic,'Synthetic molecular contamination')
        require(gr.get('source',{}).get('trait_id')==q['gwas_trait'],'GWAS trait identity mismatch')
        lm=read_json(q['ld_manifest'])
        check_hash(lm['matrix_path'],lm['matrix_sha256']);check_hash(lm['bim_path'],lm['bim_sha256'])
        inputs += [Path(q['gwas'])/'receipt.json',Path(q['gwas'])/'variants.sqlite',
                   Path(q['qtl_index'])/'receipt.json',Path(q['qtl_index'])/'qtl.sqlite',
                   q['ld_manifest'],lm['matrix_path'],lm['bim_path']]
    root=Path(c['output_root']).resolve();root.mkdir(parents=True,exist_ok=True)
    specs=root/'job_specs';specs.mkdir(exist_ok=True)
    _settings_file(root/'molecular.lock.json',{'inputs':[file_record(p) for p in inputs]})
    if (root/'summary').exists():return verify_current_artifact(root/'summary')
    def prepared(name,fn):
        target=root/name
        if target.exists():verify_current_artifact(target)
        else:fn()
        return target
    def native(name,method,settings,outputs):
        sp=_settings_file(specs/(name+'.settings.json'),settings)
        jp=create_job(sp,method=method,job_id=name,inputs={'frozen_molecular_settings':settings_path},
             outputs={**outputs,'semantic_status':'status.json'},root=specs,reviewed=True,
             synthetic=synthetic,minimum_free_bytes=c.get('minimum_free_bytes',2**30),threads=1)
        return run_job(jp,root,allow_synthetic=synthetic)
    rows=[];failures=0;successful=[]
    for q in queries:
        ident=q['query_id'];region={**q['region'],'locus_id':ident}
        row={'query_id':ident,'pair_id':q['pair_id'],'gwas_trait':q['gwas_trait'],
             'feature_id':q['feature_id'],'status':'NOT_RUN','max_h4':'NA','signal_comparisons':0,
             'coloc_artifact':'','reason':''}
        try:
            pp=prepared(ident+'_pair',lambda:qtl_locus(q['gwas'],q['qtl_index'],q['feature_id'],region,
                 root,ident+'_pair',minimum_variants=c.get('minimum_variants',20),
                 minimum_overlap=c.get('minimum_qtl_overlap',.8)))
            ix=prepared(ident+'_index',lambda:index_pair(pp,root,ident+'_index'))
            aligned=prepared(ident+'_aligned',lambda:prepare_locus(ix,q['ld_manifest'],region,
                 root,ident+'_aligned',synthetic=synthetic,max_working_bytes=c['max_working_bytes'],
                 min_coverage=c['minimum_ld_coverage']))
            fit_names=[]
            for number,key in ((1,'gwas_susie'),(2,'qtl_susie')):
                name=ident+f'_susie{number}'
                sc={**c['susie_common'],**q[key],'trait_number':number,
                    'locus_file':str(aligned/'locus.tsv'),'ld_file':str(aligned/'ld.tsv')}
                native(name,'susie',sc,{'fit':'fit.rds','pip':'pip.tsv','cs':'credible_sets.tsv','diagnostics':'diagnostics.json'})
                fit_names.append(root/name/'fit.rds')
            name=ident+'_coloc'
            cr=native(name,'coloc',{**c['coloc'],'fit1':str(fit_names[0]),'fit2':str(fit_names[1])},
                      {'coloc':'coloc.tsv','fit':'coloc.rds'})
            row.update(status=cr['scientific_status'],coloc_artifact=str(root/name))
            successful.append(root/name/'receipt.json')
            if row['status']=='PASS':
                signals=list(read_tsv(root/name/'coloc.tsv',['PP.H4.abf']))
                require(signals,'PASS molecular coloc must contain signal comparisons')
                values=[float(r['PP.H4.abf']) for r in signals]
                require(all(0<=v<=1 for v in values),'Invalid coloc posterior')
                row.update(max_h4=max(values),signal_comparisons=len(signals))
            row['reason']='Descriptive maximum only; inspect all signal pairs, H3, and prior sensitivity'
        except (ContractError,FileNotFoundError,KeyError,ValueError) as exc:
            failures+=1;row.update(status='FAILED_OR_BLOCKED',reason=str(exc))
        rows.append(row)
    with transaction(root,'summary',stage='molecular_followup',inputs=inputs+successful,
                     parameters={'queries':len(rows),'feature_family_frozen':True},synthetic=synthetic) as (work,meta):
        write_tsv(work/'molecular_evidence.tsv',list(rows[0]),rows)
        # Completion of all queries is separate from evidence of a mechanism.
        state='PASS' if not failures else 'INSUFFICIENT_EVIDENCE'
        meta['scientific_status']=state
        write_json(work/'status.json',{'status':state,'attempted':len(rows),'failed_or_blocked':failures,
          'insufficient':sum(r['status']=='INSUFFICIENT_EVIDENCE' for r in rows),
          'note':'QTL colocalization is not a mediation test. No automatically confirmed gene, cell, or causal mechanism.'})
    return verify_current_artifact(root/'summary')
