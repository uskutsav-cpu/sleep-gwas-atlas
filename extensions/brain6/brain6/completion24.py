"""Twenty-four-step computational coverage and fail-closed release packaging.

A checklist is not a substitute for native evidence. Every accepted unit binds
an exact artifact fingerprint, actual method, explicit scope and a reviewable
source contract. Post-discovery locus/feature scopes stay unresolved until
expanded from a separately hashed manifest; one locus cannot certify a family.
"""
from __future__ import annotations
import collections
import shutil
from pathlib import Path
from .pairs import verify_lock
from .artifacts import transaction, verify_artifact, verify_current_artifact
from .io import read_json, write_json, read_tsv, write_tsv, require, json_hash, file_record, check_hash, safe_id, safe_write_path

STEPS = {
 1:('environment','GLOBAL'),2:('pair_selection','GLOBAL'),3:('dense_gwas_qc','TRAIT'),
 4:('harmonization','PAIR'),5:('global_baseline','PAIR'),6:('replication','PAIR'),
 7:('local_architecture','PAIR'),8:('genomewide_pleiotropy','PAIR'),9:('independent_loci','PAIR'),
 10:('signed_ld','LOCUS'),11:('fine_mapping','LOCUS_TRAIT'),12:('trait_colocalization','LOCUS'),
 13:('molecular_qtl','FEATURE'),14:('regulatory_annotation','LOCUS'),15:('cell_enrichment','PAIR'),
 16:('pathway_enrichment','PAIR'),17:('bidirectional_mr','DIRECTION'),18:('genomic_sem','GLOBAL'),
 19:('conjunction_fdr','PAIR'),20:('cross_disorder','GLOBAL'),21:('robustness','PAIR'),
 22:('provenance','GLOBAL'),23:('regression_validation','GLOBAL'),24:('release','ACTION')}

# Structural allowlist: a frozen plan may narrow these identities, not label an
# unrelated QC table as a native analysis. External methods require a dedicated
# reviewed importer/validator rather than editing a receipt's status.
ACCEPTED_STAGES = {
 1: {'native_environment_validation'}, 2: set(), 3: {'normalize'},
 4: {'join_pair'}, 5: {'native_job'}, 6: {'replication_family'},
 7: {'local_rg_family'}, 8: {'collate_placo','joint_placo_family'},
 9: {'native_job'}, 10: {'native_job','prepare_locus'},
 11: {'native_job'}, 12: {'native_job'}, 13: {'native_job'},
 14: {'cell_annotation_scores'}, 15: {'cell_enrichment'},
 16: {'pathway_exact_stratified'}, 17: {'native_job'},
 18: {'native_job'}, 19: {'native_job'},
 20: {'cross_disorder_block_comparison'}, 21: {'robustness_family'},
 22: {'artifact_provenance_audit'}, 23: {'regression_validation'}
}
ACCEPTED_METHODS = {5:{'rg'},9:{'clump'},10:{'ld'},11:{'susie'},
 12:{'coloc'},13:{'coloc'},17:{'mr'},18:{'genomicsem'},19:{'pleiofdr'}}

SCOPED={10,11,12,13,14}
COMPLETE={'PASS','NO_SIGNAL'}
NONNUMERIC={'NO_ELIGIBLE_PAIR','NO_VALID_REPLICATION','UNDERPOWERED','DISCORDANT',
            'NO_SIGNAL','NOT_APPLICABLE','INSUFFICIENT_EVIDENCE'}


def draft(pair_lock, out):
    lock=verify_lock(pair_lock)
    pairs=[r['pair_id'] for r in lock['pairs']]
    traits=sorted({r[k] for r in lock['pairs'] for k in ['sleep_trait','disease_trait']})
    requirements=[]
    for step,(stage,scope) in STEPS.items():
        if step==24:continue
        units=['GLOBAL'] if scope=='GLOBAL' else traits if scope=='TRAIT' else [p+'__'+d for p in pairs for d in ['forward','reverse']] if scope=='DIRECTION' else pairs
        for unit in units:
            requirements.append({'id':f's{step:02d}.{unit}','step':step,'stage':stage,'unit':unit,
              'expanded_scope_required':step in SCOPED,'scope_manifest':None,'bindings':[],
              'accepted_artifact_stages':sorted(ACCEPTED_STAGES[step]), 'accepted_methods':sorted(ACCEPTED_METHODS.get(step,set())), 'required':True,
              'absence_reason':'Not bound to production evidence',
              'expected_inputs':[]})
    value={'schema_version':1,'reviewed':False,'scientific_review_by':None,
      'pair_lock':file_record(pair_lock),'synthetic':lock['synthetic'],
      'requirements':requirements,'no_eligible_disorders':[r['disease_trait'] for r in lock.get('disorder_decisions',[]) if r['role']=='NO_ELIGIBLE_PAIR'],
      'release_exports':[], 'note':'DRAFT ONLY. Define accepted method identities and locus/feature scopes before results. Do not remove failed requirements after seeing results.'}
    write_json(out,value);return value


def _bound_file(record):
    require(isinstance(record,dict) and record.get('path') and record.get('sha256'),'Expected pinned file record')
    check_hash(record['path'],record['sha256'])
    if 'bytes' in record:require(Path(record['path']).stat().st_size==record['bytes'],'Bound input size changed')
    return Path(record['path'])


def validate_plan(c):
    require(c.get('schema_version')==1,'Unknown 24-step plan schema')
    _bound_file(c['pair_lock']);lock=verify_lock(c['pair_lock']['path'])
    require(lock['synthetic']==c.get('synthetic',False),'Synthetic pair-lock mismatch')
    for key in ['atlas','config','decisions']:_bound_file(lock[key])
    req=c['requirements'];ids=[r['id'] for r in req]
    require(req and len(ids)==len(set(ids)),'Empty or duplicate requirement IDs')
    for r in req:
        safe_id(r['id']);require(r['step'] in range(1,24),'Invalid computational step')
        require(r['stage']==STEPS[r['step']][0],'Mismatched step/stage identity')
        require(isinstance(r['required'],bool),'Required flag must be explicit')
        require(set(r.get('accepted_artifact_stages',[]))<=ACCEPTED_STAGES[r['step']],
                'Wrong artifact class for this computational step')
        require(set(r.get('accepted_methods',[]))<=ACCEPTED_METHODS.get(r['step'],set()),
                'Wrong native method for this computational step')
        require(bool(r.get('expanded_scope_required'))==(r['step'] in SCOPED),
                'Cannot disable required locus/feature scope expansion')
        if not r['required']:
            require(r.get('prespecified_omission_reason') and r.get('omission_reviewed_before_results') is True,
                    'Omitting a failed method after seeing results is forbidden')
    # All required scopes from the base plan must be represented; additional
    # sensitivity requirements are permitted but no disease/pair disappears.
    pairs=[r['pair_id'] for r in lock['pairs']]
    traits=sorted({r[k] for r in lock['pairs'] for k in ['sleep_trait','disease_trait']})
    seen={(r['step'],r['unit']) for r in req}
    for step,(stage,scope) in STEPS.items():
        if step==24:continue
        units=['GLOBAL'] if scope=='GLOBAL' else traits if scope=='TRAIT' else [p+'__'+d for p in pairs for d in ['forward','reverse']] if scope=='DIRECTION' else pairs
        require({(step,u) for u in units}<=seen,f'Missing planned units in step {step}')
    return lock


def audit(manifest, root, name, *, allow_synthetic=False):
    c=read_json(manifest);lock=validate_plan(c);synthetic=lock['synthetic']
    require(not synthetic or allow_synthetic,'Synthetic audit requires explicit flag')
    inputs=[manifest,c['pair_lock']['path']];rows=[]
    reviewed=c.get('reviewed') is True and bool(c.get('scientific_review_by'))
    for req in c['requirements']:
        row={'id':req['id'],'step':req['step'],'stage':req['stage'],'unit':req['unit'],
             'status':'MISSING','scientific_status':'UNRESOLVED','expected_units':1,'bound_units':0,
             'reason':req.get('absence_reason',''),'artifacts':''}
        if not req['required']:
            row.update(status='PRESPECIFIED_OUT_OF_SCOPE',scientific_status='NOT_APPLICABLE',reason=req['prespecified_omission_reason']);rows.append(row);continue
        if req['step']==2:
            row.update(status='VERIFIED',scientific_status='PASS',bound_units=1,reason='Content-verified reviewed pair lock including no-eligible decisions')
            rows.append(row);continue
        expected=['WHOLE_SCOPE']
        if req.get('expanded_scope_required'):
            scope=req.get('scope_manifest')
            if not scope:
                row.update(status='UNEXPANDED_SCOPE',reason='Actual locus/feature universe is not bound');rows.append(row);continue
            sp=_bound_file(scope);inputs.append(sp)
            entries=list(read_tsv(sp,['unit_id']))
            expected=[x['unit_id'] for x in entries]
            require(len(expected)==len(set(expected)),'Duplicate expanded locus/feature unit')
            # A genuinely empty selected locus set requires the discovery
            # artifact and NO_SIGNAL state, not an arbitrary empty TSV.
            if not expected:
                no=req.get('no_signal_evidence')
                if no:
                    rr=verify_artifact(no['path'],no['fingerprint']);inputs.append(Path(no['path'])/'receipt.json')
                    require(rr['synthetic']==synthetic and rr.get('scientific_status')=='NO_SIGNAL'
                            and rr['stage'] in {'collate_placo','joint_placo_family','native_job'},
                            'Invalid no-signal applicability evidence')
                    if rr['stage']=='native_job':require(rr['parameters']['method']=='clump','Empty scope must come from locus discovery')
                    row.update(status='VERIFIED',scientific_status='NOT_APPLICABLE_NO_SIGNAL',expected_units=0,bound_units=0,reason='Empty planned scope bound to verified no-signal discovery')
                else:row.update(status='EMPTY_SCOPE_UNJUSTIFIED',reason='No upstream NO_SIGNAL evidence')
                rows.append(row);continue
        row['expected_units']=len(expected)
        bindings=req.get('bindings',[]);units=[b['unit_id'] for b in bindings]
        require(len(units)==len(set(units)) and set(units)<=set(expected),'Duplicate/extra evidence unit')
        require(len({(b['path'],b['fingerprint']) for b in bindings})==len(bindings),'The same artifact cannot stand in for distinct expanded units')
        if not bindings:rows.append(row);continue
        accepted=req.get('accepted_artifact_stages',[])
        require(accepted,'Freeze accepted artifact types; arbitrary JSON is not production evidence')
        states=[];paths=[]
        for binding in bindings:
            rr=verify_artifact(binding['path'],binding['fingerprint'])
            require(rr.get('synthetic')==synthetic,'Synthetic artifact cannot enter empirical completion')
            if not synthetic:
                require(req.get('expected_inputs'),
                        'Production evidence must pin the unit-specific input bytes, not just a passing receipt')
            require(rr['stage'] in accepted,'Artifact method/stage differs from frozen requirement')
            if rr['stage']=='native_job':
                require(rr['parameters']['method'] in req.get('accepted_methods',[]),'Unexpected native method')
                if req['step']==18:
                    setting_record=rr['parameters'].get('inputs',{}).get('settings')
                    require(setting_record,'Genomic SEM model evidence must bind its settings')
                    settings=_bound_file(setting_record);inputs.append(settings)
                    require(read_json(settings).get('mode')=='model',
                            'Covariance preparation is not a validated Genomic SEM model')
            pinned={r['sha256'] for r in rr['inputs']}
            for expected_input in req.get('expected_inputs',[]):
                ip=_bound_file(expected_input);inputs.append(ip)
                require(expected_input['sha256'] in pinned,'Artifact was not computed with required input bytes')
            inputs.append(Path(binding['path'])/'receipt.json');paths.append(binding['path'])
            ss=rr.get('scientific_status','UNRESOLVED');states.append(ss)
        row['bound_units']=len(bindings);row['artifacts']=';'.join(paths)
        if len(bindings)!=len(expected):row.update(status='PARTIAL_SCOPE',reason='One passing unit cannot represent all planned loci/features')
        elif all(s in COMPLETE for s in states):row.update(status='VERIFIED',scientific_status='PASS',reason='Every declared unit has pinned, QC-accepted evidence')
        elif all(s in COMPLETE|NONNUMERIC for s in states):
            row.update(status='TERMINAL_LIMITATION',scientific_status='INSUFFICIENT_EVIDENCE',reason='Executed evidence remains underpowered, unavailable, discordant or unresolved; not a validated scientific finish')
        else:row.update(status='FAILED_OR_UNRESOLVED',reason='At least one artifact failed scientific acceptance')
        rows.append(row)
    resolved=sum(r['status'] in {'VERIFIED','PRESPECIFIED_OUT_OF_SCOPE'} for r in rows)
    ready=reviewed and resolved==len(rows)
    with transaction(root,name,stage='computational_24_audit',inputs=inputs,parameters=c,synthetic=synthetic) as (work,meta):
        write_tsv(work/'coverage.tsv',list(rows[0]),rows)
        summary=[]
        for step in range(1,25):
            group=[r for r in rows if r['step']==step]
            state=('READY_TO_SEAL' if ready else 'BLOCKED') if step==24 else 'VERIFIED' if group and all(r['status'] in {'VERIFIED','PRESPECIFIED_OUT_OF_SCOPE'} for r in group) else 'INCOMPLETE'
            summary.append({'step':step,'stage':STEPS[step][0],'status':state,'requirements':len(group),
                            'resolved':sum(r['status'] in {'VERIFIED','PRESPECIFIED_OUT_OF_SCOPE'} for r in group)})
        write_tsv(work/'steps_1_24.tsv',list(summary[0]),summary)
        meta['scientific_status']='PASS' if ready else 'INSUFFICIENT_EVIDENCE'
        write_json(work/'status.json',{'status':'READY_TO_SEAL' if ready else 'INCOMPLETE','reviewed':reviewed,
             'requirements':len(rows),'resolved':resolved,'unresolved':len(rows)-resolved,'synthetic':synthetic,
             'publication_ready':False,'note':'Completeness is scoped to this frozen plan; no study causality or wet-lab validation is certified.'})
    return Path(root)/name


def seal(manifest, audit_dir, root, name, *, allow_synthetic=False):
    c=read_json(manifest);lock=validate_plan(c)
    require(c.get('reviewed') is True and c.get('scientific_review_by'),'Final plan review is required')
    require(not lock['synthetic'] or allow_synthetic,'Synthetic results cannot become an empirical release')
    rr=verify_current_artifact(audit_dir)
    require(rr['stage']=='computational_24_audit','Release needs the actual 24-step audit')
    from .provenance24 import verify_graph
    verify_graph([{'path':str(audit_dir),'fingerprint':rr['fingerprint']}],synthetic=lock['synthetic'])
    require(rr['synthetic']==lock['synthetic'] and rr.get('scientific_status')=='PASS','Computational plan is incomplete; release refused')
    require(any(r['sha256']==file_record(manifest)['sha256'] for r in rr['inputs']),'Audit belongs to a different plan')
    exports=c.get('release_exports',[]);require(exports,'No licensed summary artifacts selected for release')
    destinations={'receipt.json','RELEASE.json','COMPUTATIONAL_COVERAGE.tsv'};bound_receipts={Path(b['path']).resolve() for r in c['requirements'] for b in r.get('bindings',[])}
    inputs=[manifest,Path(audit_dir)/'receipt.json'];verified=[]
    for e in exports:
        require(e.get('redistribution_reviewed') is True and e.get('content_class')=='aggregate_summary','Raw/individual GWAS data must not be bundled automatically')
        source=Path(e['artifact']).resolve();require(source in bound_receipts,'Release file is outside audited artifact inventory')
        er=verify_artifact(source,e['fingerprint']);p=safe_write_path(source,e['file'])
        require(any(o['path']==e['file'] for o in er['outputs']),'Export is not a sealed output')
        require(p.suffix.lower() in {'.tsv','.csv','.json','.txt','.md'},'Only portable aggregate text summaries may be bundled')
        require(p.stat().st_size<=20*1024**2,'Large export needs separate data repository/licensing review')
        destination=e['destination'];safe_write_path(root,destination)
        require(destination not in destinations,'Release output collision');destinations.add(destination)
        inputs.append(p);verified.append((p,destination))
    with transaction(root,name,stage='computational_release',inputs=inputs,parameters=c,synthetic=lock['synthetic']) as (work,meta):
        for p,dest in verified:
            q=safe_write_path(work,dest);q.parent.mkdir(parents=True,exist_ok=True);shutil.copyfile(p,q)
        shutil.copyfile(Path(audit_dir)/'coverage.tsv',work/'COMPUTATIONAL_COVERAGE.tsv')
        write_json(work/'RELEASE.json',{'computational_scope_complete':True,'synthetic':lock['synthetic'],
            'pair_lock_sha256':lock['lock_sha256'],'plan_sha256':file_record(manifest)['sha256'],
            'publication_ready':False,'science_claim':'No automatic causal-mechanism certification'})
        meta['scientific_status']='PASS'
    return Path(root)/name
