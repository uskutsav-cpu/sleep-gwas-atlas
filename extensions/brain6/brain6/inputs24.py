"""Inventory original source contracts and actual host bytes, without assuming downloaded=usable."""
from __future__ import annotations
import re
from pathlib import Path
from .io import read_tsv, write_tsv, write_json, sha256, require
from .artifacts import transaction
from .pairs import BRAIN, SLEEP


def inventory(repo, root, name, *, traits=None, hash_dense=False):
    repo=Path(repo).resolve()
    panel=repo/'config/analysis_panel.tsv';registry=repo/'config/public_gwas_sources.tsv';schema=repo/'config/gwas_schemas.tsv'
    p={r['trait_id']:r for r in read_tsv(panel,['trait_id','source_id','raw_file','build'])}
    sources={r['source_id']:r for r in read_tsv(registry,['source_id','download_url','archive_name','archive_sha256','archive_bytes'])}
    schemas={(r['trait_id'],r['source_id']):r for r in read_tsv(schema,['trait_id','source_id','effect_convention','sample_size'])}
    needed=sorted(traits or (set(BRAIN)|set(SLEEP)))
    require(set(needed)<=set(p),'Unknown source traits')
    rows=[];inputs=[panel,registry,schema]
    for t in needed:
        row=p[t];source=sources[row['source_id']];sc=schemas[(t,row['source_id'])]
        # Input symlinks may point at a read-only shared data store. Unlike
        # output links they are legal, but targets must actually be present.
        raw=repo/'data/raw'/row['raw_file']
        dense=repo/'data/harmonized'/(t+'.harmonized.tsv.gz')
        archive=repo/'data/raw/.archives'/source['archive_name']
        state='MISSING_DENSE_BYTES';dh=''; ah=''
        if dense.is_file():
            state='PRESENT_NOT_ROW_VALIDATED'
            if hash_dense:dh=sha256(dense);inputs.append(dense);state='HASHED_NOT_ROW_VALIDATED'
        if archive.is_file() and hash_dense:
            ah=sha256(archive);inputs.append(archive)
            require(ah==source['archive_sha256'] and archive.stat().st_size==int(source['archive_bytes']),f'Raw archive mismatch: {t}')
        issues=[]
        if sc['chromosome']=='ABSENT' or sc['position']=='ABSENT':issues.append('AUDITED_VARIANT_COORDINATE_MAP_REQUIRED')
        effect=sc['effect_convention']
        if 'lmm' in (effect+' '+sc.get('notes','')).lower():issues.append('BINARY_LMM_SCALE_NOT_AUTOMATIC_LOG_ODDS')
        if 'NEFF' in sc['sample_size'].upper():issues.append('SOURCE_SPECIFIC_EFFECTIVE_N_CONVERSION')
        if row['build'] not in {'hg19','GRCh37'}:issues.append('BUILD_RECONCILIATION_REQUIRED')
        if sc.get('info')=='ABSENT':issues.append('NO_ROW_WISE_INFO_IN_SOURCE')
        rows.append({'trait_id':t,'source_id':row['source_id'],'source_build':row['build'],
          'raw_path':str(raw),'raw_present':raw.is_file(),'dense_path':str(dense),'dense_state':state,
          'dense_bytes':dense.stat().st_size if dense.is_file() else 0,'dense_sha256':dh,
          'archive_path':str(archive),'archive_present':archive.is_file(),'expected_archive_bytes':source['archive_bytes'],
          'expected_archive_sha256':source['archive_sha256'],'observed_archive_sha256':ah,
          'effect_convention':effect,'sample_size_convention':sc['sample_size'],
          'download_url':source['download_url'],'source_notes':source.get('notes',''),
          'required_reviews':';'.join(issues)})
    with transaction(root,name,stage='source_inventory',inputs=inputs,parameters={'traits':needed,'hash_dense':hash_dense},synthetic=False) as (work,meta):
        write_tsv(work/'sources.tsv',list(rows[0]),rows)
        meta['scientific_status']='INSUFFICIENT_EVIDENCE'
        write_json(work/'summary.json',{'traits':len(rows),'dense_present':sum(r['dense_bytes']>0 for r in rows),
            'raw_present':sum(r['raw_present'] for r in rows),'raw_archives_hash_verified':sum(bool(r['observed_archive_sha256']) for r in rows),
            'production_ready':False,'rule':'Existence or an archive hash is not a substitute for schema, effect-scale, variant mapping, ancestry, and row-level QC. Existing dense manifests must be rebound to these exact bytes.'})
    return Path(root)/name
