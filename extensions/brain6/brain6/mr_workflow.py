"""Exposure-only MR preparation and family-wide result accounting.

This module prepares *inputs* to native TwoSampleMR. It does not turn PLACO
hits into instruments, infer cohort independence, or assert causality.
"""
from __future__ import annotations
import math
import sqlite3
from pathlib import Path
from .artifacts import transaction, verify_artifact
from .gwas import FIELDS, align
from .io import check_hash, read_json, read_tsv, require, write_json, write_tsv
from .pairs import verify_lock
from .stats import bh


def candidates(exposure, root, name, *, instrument_p=5e-8, minimum_f=10.0):
    exposure = Path(exposure).resolve()
    rec = verify_artifact(exposure)
    require(rec.get('stage') == 'normalize', 'Use normalized full exposure GWAS, not a PLACO pair')
    require(0 < instrument_p < 1 and math.isfinite(minimum_f) and minimum_f > 0,
            'Invalid exposure instrument thresholds')
    params = {'instrument_p': instrument_p, 'minimum_f': minimum_f,
              'exposure_only': True, 'exposure': str(exposure)}
    with transaction(root, name, stage='mr_candidates',
                     inputs=[exposure/'receipt.json', exposure/'variants.sqlite'],
                     parameters=params, synthetic=rec['synthetic']) as (work, meta):
        db = sqlite3.connect(f'{(exposure/"variants.sqlite").as_uri()}?mode=ro', uri=True)
        counts = {'scanned': 0, 'p_pass': 0, 'retained': 0}
        def rows():
            for r in db.execute('SELECT snp,chr,bp,a1,a2,beta,se,p,n,eaf,info FROM variants ORDER BY chr,bp,snp'):
                counts['scanned'] += 1
                if r[7] >= instrument_p:
                    continue
                counts['p_pass'] += 1
                if (r[5]/r[6])**2 < minimum_f:
                    continue
                counts['retained'] += 1
                yield dict(zip(FIELDS, ['NA' if x is None else x for x in r]))
        try:
            write_tsv(work/'candidates.tsv', FIELDS, rows())
        finally:
            db.close()
        meta.update(scientific_status='PASS' if counts['retained'] else 'NO_SIGNAL',
                    source=rec['source'])
        write_json(work/'status.json', {'status': meta['scientific_status'], **counts,
                   'instrument_selection': 'EXPOSURE_P_AND_F_ONLY',
                   'outcome_used_for_selection': False})
    return Path(root)/name


def harmonize(candidates_dir, clump_dir, outcome_dir, root, name, *, max_eaf_difference=.15, max_instrument_r2=.001, minimum_clump_kb=10000):
    """Extract outcome statistics only AFTER exposure-only native LD clumping.

    Binds the clumping job's hashed settings to this candidate table. Every lost
    instrument is retained in attrition.tsv. Outcome P is never a filter.
    """
    cand, clump, outcome = [Path(x).resolve() for x in (candidates_dir, clump_dir, outcome_dir)]
    cr, lr, yr = [verify_artifact(x) for x in (cand, clump, outcome)]
    require(cr['stage'] == 'mr_candidates' and yr['stage'] == 'normalize', 'Unexpected MR input artifacts')
    require(cr['synthetic'] == lr['synthetic'] == yr['synthetic'], 'Synthetic/empirical MR mix')
    require(lr['stage'] == 'native_job' and lr.get('scientific_status') in {'PASS','NO_SIGNAL'},
            'Clumping must be an accepted native job')
    require(lr['parameters'].get('method') == 'clump', 'Receipt is not an LD-clumping job')
    settings_rec = lr['parameters']['inputs']['settings']
    check_hash(settings_rec['path'], settings_rec['sha256'])
    sc = read_json(settings_rec['path'])
    require(Path(sc['associations']).resolve() == cand/'candidates.tsv' and sc['p_column'] == 'P',
            'Clumping must use the exposure-only candidate table, not PLACO or outcome P')
    require(sc['p1'] == cr['parameters']['instrument_p'], 'Exposure/clump P threshold mismatch')
    require(0<max_instrument_r2<1 and minimum_clump_kb>0,'Invalid MR LD policy')
    require(0<float(sc.get('r2',1))<=max_instrument_r2 and sc.get('kb',0)>=minimum_clump_kb,
            'MR LD clumping does not meet the instrument independence policy')
    for key in ('genome_build','ancestry'):
        require(cr['source'][key] == yr['source'][key], 'MR metadata mismatch: '+key)
    require(0 <= max_eaf_difference <= 1, 'Invalid EAF tolerance')
    exposure_rows = {r['SNP']: r for r in read_tsv(cand/'candidates.tsv', FIELDS)}
    leads = list(read_tsv(clump/'leads.tsv', ['SNP']))
    require(len({r['SNP'] for r in leads}) == len(leads), 'Duplicate clumped MR instrument')
    with transaction(root, name, stage='mr_harmonization',
                     inputs=[cand/'receipt.json', cand/'candidates.tsv', clump/'receipt.json',
                             clump/'leads.tsv', outcome/'receipt.json', outcome/'variants.sqlite'],
                     parameters={'max_eaf_difference': max_eaf_difference,'max_instrument_r2':max_instrument_r2,
                                 'minimum_clump_kb':minimum_clump_kb},
                     synthetic=cr['synthetic']) as (work, meta):
        db = sqlite3.connect(f'{(outcome/"variants.sqlite").as_uri()}?mode=ro', uri=True)
        keep_x, keep_y, attrition = [], [], []
        try:
            for lead in leads:
                snp = lead['SNP']
                require(snp in exposure_rows, 'Clumped SNP was not an eligible exposure instrument')
                x = exposure_rows[snp]
                values = db.execute('SELECT snp,chr,bp,a1,a2,beta,se,p,n,eaf,info FROM variants WHERE snp=?', (snp,)).fetchone()
                reason = 'RETAINED'
                y = None
                if values is None:
                    reason = 'MISSING_OUTCOME'
                else:
                    y = dict(zip(FIELDS, values))
                    if (int(x['CHR']), int(x['BP'])) != (y['CHR'], y['BP']):
                        reason = 'COORDINATE_MISMATCH'
                    else:
                        orient = align(x['A1'],x['A2'],y['A1'],y['A2'])
                        if orient is None:
                            reason = 'ALLELE_MISMATCH'
                        else:
                            sign, _ = orient
                            y.update(BETA=sign*y['BETA'], A1=x['A1'], A2=x['A2'])
                            if y['EAF'] is not None and sign == -1:
                                y['EAF'] = 1-y['EAF']
                            if x['EAF'] != 'NA' and y['EAF'] is not None and abs(float(x['EAF'])-y['EAF']) > max_eaf_difference:
                                reason = 'FREQUENCY_MISMATCH'
                attrition.append({'SNP':snp,'status':reason})
                if reason == 'RETAINED':
                    keep_x.append(x)
                    keep_y.append({k:'NA' if v is None else v for k,v in y.items()})
        finally:
            db.close()
        write_tsv(work/'exposure.tsv', FIELDS, keep_x)
        write_tsv(work/'outcome.tsv', FIELDS, keep_y)
        write_tsv(work/'attrition.tsv', ['SNP','status'], attrition)
        state = 'PASS' if len(keep_x) >= 3 else 'INSUFFICIENT_EVIDENCE'
        meta.update(scientific_status=state, exposure_source=cr['source'], outcome_source=yr['source'])
        write_json(work/'status.json', {'status':state,'clumped_instruments':len(leads),
                'aligned_instruments':len(keep_x),'outcome_p_filter':False,
                'note':'Inputs to native MR, not a causal result. Run reverse direction using independently exposure-selected instruments.'})
    return Path(root)/name


def collate(manifest_path, root, name):
    """Keep both directions of every frozen pair in the IVW testing family.

    The manifest maps pair_id + FORWARD/REVERSE to native MR artifact paths.
    Missing/failed jobs occupy family slots but retain NA results.
    """
    c = read_json(manifest_path)
    lock = verify_lock(c['pair_lock'])
    require(c.get('reviewed') is True, 'Review MR result bindings')
    expected = {(r['pair_id'],d) for r in lock['pairs'] for d in ('FORWARD','REVERSE')}
    pair_map={r['pair_id']:r for r in lock['pairs']}
    supplied = {}
    for row in c.get('results',[]):
        key = row['pair_id'], row['direction']
        require(key in expected and key not in supplied, 'Unexpected/duplicate MR unit')
        supplied[key] = row
    rows = []
    inputs = [manifest_path,c['pair_lock']]
    for pair,direction in sorted(expected):
        row = dict(pair_id=pair,direction=direction,status='NOT_RUN',beta=None,se=None,p=None,n_instruments=0,reason='No bound result')
        bound = supplied.get((pair,direction))
        if bound is not None:
            if bound.get('path'):
                path = Path(bound['path'])
                r = verify_artifact(path, bound.get('fingerprint'))
                require(r['synthetic'] == lock['synthetic'], 'Synthetic MR contamination')
                require(r['stage'] == 'native_job' and r['parameters'].get('method') == 'mr',
                        'Expected a native MR result, not diagnostic IVW arithmetic')
                require(r.get('scientific_status') in {'PASS','INSUFFICIENT_EVIDENCE'}, 'Rejected MR artifact')
                if not lock['synthetic']:
                    require(bound.get('harmonization_dir'),'Bind the MR harmonization receipt to verify pair and direction')
                    hd=Path(bound['harmonization_dir']).resolve();hr=verify_artifact(hd)
                    require(hr['stage']=='mr_harmonization' and not hr['synthetic'],'Invalid empirical MR preparation')
                    expected_pair=pair_map[pair]
                    et,ot=expected_pair['sleep_trait'],expected_pair['disease_trait']
                    if direction=='REVERSE':et,ot=ot,et
                    require(hr['exposure_source']['trait_id']==et and hr['outcome_source']['trait_id']==ot,
                            'Native MR input traits do not match the frozen pair/direction')
                    sr=r['parameters']['inputs']['settings'];check_hash(sr['path'],sr['sha256'])
                    nc=read_json(sr['path'])
                    require(Path(nc['exposure']).resolve()==hd/'exposure.tsv' and Path(nc['outcome']).resolve()==hd/'outcome.tsv',
                            'Native MR job did not use the bound harmonized exposure/outcome')
                    inputs.append(hd/'receipt.json')
                inputs += [path/'receipt.json',path/'mr.tsv']
                row.update(status=r['scientific_status'],reason='')
                if row['status'] == 'PASS':
                    vals=[x for x in read_tsv(path/'mr.tsv',['method','b','se','pval','nsnp'])
                          if x['method'] == 'Inverse variance weighted']
                    require(len(vals)==1, 'Exactly one primary IVW estimate required per direction')
                    v=vals[0]; b,se,p=map(float,(v['b'],v['se'],v['pval']))
                    require(all(map(math.isfinite,(b,se,p))) and se>0 and 0<=p<=1,'Invalid IVW result')
                    require(float(v['nsnp']).is_integer() and int(v['nsnp'])>=3,'Primary IVW needs at least three retained instruments')
                    row.update(beta=b,se=se,p=p,n_instruments=int(v['nsnp']))
            else:
                require(bound.get('status') in {'BLOCKED_BY_DATA','BLOCKED_BY_SOFTWARE','FAILED_QC','UNDERPOWERED'},
                        'Missing MR artifact needs an explicit blocker')
                require(bound.get('reason'),'Missing MR blocker reason')
                row.update(status=bound['status'],reason=bound['reason'])
        rows.append(row)
    qs = bh([r['p'] for r in rows]) if rows else []
    for r,q in zip(rows,qs):
        r['family_fdr']=q
        r['interpretation']='NO_DIRECTIONAL_RESULT' if q is None else ('DIRECTIONAL_ASSOCIATION_REQUIRES_IV_REVIEW' if q<.05 else 'NO_FDR_SIGNIFICANT_DIRECTIONAL_EVIDENCE')
    with transaction(root,name,stage='mr_family_summary',inputs=inputs,
                     parameters={'expected_tests':len(expected)},synthetic=lock['synthetic']) as (work,meta):
        columns=['pair_id','direction','status','beta','se','p','family_fdr','n_instruments','reason','interpretation']
        write_tsv(work/'mr_family.tsv',columns,({k:'NA' if v is None else v for k,v in r.items()} for r in rows))
        done=sum(r['p'] is not None for r in rows)
        meta['scientific_status']='PASS' if done == len(rows) and rows else 'INSUFFICIENT_EVIDENCE'
        write_json(work/'status.json',{'status':meta['scientific_status'],'expected_tests':len(expected),
                   'tested':done,'unresolved':len(rows)-done,'scope':'PRIMARY_IVW_FAMILY_ONLY',
                   'note':'Method agreement and nonsignificant Egger tests do not establish valid IV assumptions.'})
    return Path(root)/name
