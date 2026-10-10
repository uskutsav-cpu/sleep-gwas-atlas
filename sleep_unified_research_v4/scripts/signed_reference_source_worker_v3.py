"""Authenticated finite-reference G2--G5 only. No GWAS, simulation or fit."""
import argparse
from collections import Counter
import csv
import fcntl
import gzip
import hashlib
import json
import os
from pathlib import Path
import sqlite3
import subprocess
import sys
import time

from signed_reference_source_common_v3 import (load_frozen, save, binary_new,
    streamed_new, Guard, regular_sha, assert_no_termination,record_output,INTENDED_OUTPUT_SHA256,verify_outputs,identity_gate)
from canonical_calibration_common_v4_8 import validate_inherited_lock
from signed_reference_compact_v3 import (StaticWriter, GenotypeWriter, static_rows, genotype_rows,
    CAPACITY, SCHEMA, SCHEMA_SHA256)
from signed_reference_decode_v2 import (hashes, regular, extract_strict66,
    decode_a1, qc_column, standardized, allele_sign, frozen_control_indices,
    operator, close_array, BED_HEADER, BYTES_PER_SNP)


def line(value):
    return (json.dumps(value,sort_keys=True,separators=(',',':'),allow_nan=False)+'\n').encode()


def consume_jsonl(path,expected,guard):
    if regular_sha(path)!=expected:raise RuntimeError('SOURCE_ORDERED_INPUT_CHANGED_BEFORE_PARSE')
    with regular(path).open() as f:
        for text in f:
            guard();yield json.loads(text)
    if regular_sha(path)!=expected:raise RuntimeError('SOURCE_ORDERED_INPUT_CHANGED_AFTER_PARSE')


def curl_file(plan,url,out,head=False):
    # The sole source worker owns this synchronous curl child in its own group.
    # No resume, retry, redirect, disabled TLS, source mirror or automatic restart.
    cap=1048576 if head or url.endswith('/8292725') else plan['candidate']['bytes']
    command=[plan['curl'],'-q','--fail','--silent','--show-error','--proto','=https',
        '--proto-redir','=https','--connect-timeout','30','--max-time','6900',
        '--max-filesize',str(cap)]
    if head:command+=['--head']
    command+=['--output',str(out),'--write-out','%{http_code}\n%{url_effective}\n',url]
    result=subprocess.run(command,capture_output=True,check=False)
    expected=hashlib.sha256(result.stdout).hexdigest()
    binary_new(str(out)+'.curl_status',result.stdout)
    if result.returncode or result.stdout.decode().splitlines()!=['200',url]:
        raise RuntimeError('EXACT_TLS_SOURCE_CURL_FAILED_NO_RETRY: '+str(result.returncode))
    if regular_sha(str(out)+'.curl_status')!=expected:raise RuntimeError('CURL_STATUS_CHANGED')
    with regular(out).open('rb') as file:os.fsync(file.fileno())
    fd=os.open(Path(out).parent,os.O_RDONLY)
    try:os.fsync(fd)
    finally:os.close(fd)
    identity=hashes(out);record_output(out,identity['sha256']);return identity


def acquire(plan,root,guard):
    candidate=plan['candidate'];meta=root/'qc/fresh_zenodo_metadata.json';head=root/'qc/fresh_archive_HEAD.txt'
    guard(force=True);mh=curl_file(plan,'https://zenodo.org/api/records/8292725',meta)
    m=load_frozen(meta,mh['sha256']);matches=[f for f in m.get('files',[]) if f.get('key')==candidate['file_name']]
    if len(matches)!=1 or m.get('id')!=8292725 or m.get('doi')!='10.5281/zenodo.8292725':raise RuntimeError('OFFICIAL_RECORD_EXACT_VERSION_FAILED')
    f=matches[0]
    if (f.get('id'),f.get('size'),f.get('checksum'),f.get('links',{}).get('self'))!=(candidate['file_id'],candidate['bytes'],candidate['upstream_md5'],candidate['content_url']):raise RuntimeError('OFFICIAL_ARCHIVE_FILE_ID_SIZE_MD5_URL_CHANGED')
    if m.get('metadata',{}).get('access_right')!='open' or m['metadata'].get('license',{}).get('id')!='cc-by-4.0':raise RuntimeError('SOURCE_ACCESS_OR_LICENSE_CHANGED')
    guard(force=True);hh=curl_file(plan,candidate['content_url'],head,True)
    raw=regular(head).read_text();sizes=[int(s.split(':',1)[1].strip()) for s in raw.splitlines() if s.lower().startswith('content-length:')]
    if sizes!=[candidate['bytes']] or regular_sha(head)!=hh['sha256']:raise RuntimeError('EXACT_HEAD_CONTENT_LENGTH_REQUIRED')
    archive=root/'archive/1000G_Phase3_plinkfiles.tgz'
    guard(force=True);identity=curl_file(plan,candidate['content_url'],archive)
    guard(force=True)
    if identity['bytes']!=candidate['bytes'] or 'md5:'+identity['md5']!=candidate['upstream_md5']:raise RuntimeError('UPSTREAM_ARCHIVE_FULL_SIZE_MD5_FAILED')
    save(root/'qc/G2_authentication.json',dict(status='AUTHENTICATED_NEW_CANDIDATE_NOT_OLD_COMPACT_IDENTITY',
        metadata_sha256=mh['sha256'],HEAD_sha256=hh['sha256'],archive_identity=identity,
        G1_old_local_continuity='NOT_REQUESTED_NOT_RUN_NO_AUTHENTICATION_CREDIT',candidate=candidate))
    return archive,identity


def configure_db(path):
    db=sqlite3.connect(path)
    db.execute('PRAGMA journal_mode=DELETE');db.execute('PRAGMA temp_store=MEMORY')
    db.execute('PRAGMA cache_size=-32768')
    return db


def source_format(plan,root,extracted,guard):
    files=extracted['files'];panel=Path(plan['authoritative_panel_path'])
    if regular_sha(panel)!=plan['reference_body_sha256'][str(panel)]:raise RuntimeError('AUTHORITATIVE_PANEL_CHANGED')
    with panel.open() as f:
        reader=csv.DictReader(f,delimiter='\t')
        if reader.fieldnames!=['sample','pop','super_pop','gender']:raise RuntimeError('AUTHORITATIVE_PANEL_SCHEMA')
        eur=[r['sample'] for r in reader if r['super_pop']=='EUR']
    if len(eur)!=503 or len(set(eur))!=503:raise RuntimeError('EXACT503_EUR_PANEL_REQUIRED')
    path=root/'extracted/source_index.sqlite';db=configure_db(path)
    db.execute('CREATE TABLE source(SNP TEXT PRIMARY KEY,CHR INT,BP INT,A1 TEXT,A2 TEXT,row0 INT,offset INT)')
    first_fam=None;order=None;counts={};allbed={}
    try:
        for chrom in range(1,23):
            guard(force=True);fam=regular(files[f'1000G.EUR.QC.{chrom}.fam']['path'])
            if regular_sha(fam)!=files[fam.name]['sha256']:raise RuntimeError('FAM_CHANGED_BEFORE_PARSE')
            payload=fam.read_bytes()
            if first_fam is None:
                first_fam=payload;rows=[x.split() for x in payload.decode().splitlines()]
                if len(rows)!=503 or any(len(r)!=6 for r in rows):raise RuntimeError('EXACT503_SIXFIELD_FAM_REQUIRED')
                order=[[r[0],r[1]] for r in rows]
                if len({r[1] for r in order})!=503 or {r[1] for r in order}!=set(eur):raise RuntimeError('FAM_AUTHORITATIVE_EUR_SAMPLE_SET_FAILED')
                binary_new(root/'qc/private_sample_order.json',line(order))
            elif payload!=first_fam:raise RuntimeError('ALL22_FAM_BYTE_IDENTICAL_REQUIRED')
            if regular_sha(fam)!=files[fam.name]['sha256']:raise RuntimeError('FAM_CHANGED_AFTER_PARSE')
            bim=regular(files[f'1000G.EUR.QC.{chrom}.bim']['path']);bed=regular(files[f'1000G.EUR.QC.{chrom}.bed']['path'])
            for p in [bim,bed]:
                if hashes(p,guard)['sha256']!=files[p.name]['sha256']:raise RuntimeError('SOURCE_MEMBER_CHANGED_BEFORE_FORMAT_PARSE')
            n=0;previous=0
            with bim.open() as f:
                for text in f:
                    guard();r=text.split()
                    if len(r)!=6 or r[0]!=str(chrom) or not r[3].isdigit() or int(r[3])<=0 or int(r[3])<previous or not r[1]:raise RuntimeError('BIM_SIXFIELD_CHR_BP_OR_ORDER_FAILED')
                    previous=int(r[3]);db.execute('INSERT INTO source VALUES(?,?,?,?,?,?,?)',(r[1],chrom,previous,r[4],r[5],n,3+n*126));n+=1
            db.commit()
            with bed.open('rb') as f:
                if f.read(3)!=BED_HEADER or bed.stat().st_size!=3+n*126:raise RuntimeError('BED_EXACT_VARIANT_MAJOR_SIZE_REQUIRED')
                while True:
                    guard();b=f.read(126*2048)
                    if not b:break
                    if len(b)%126 or any(x & 0xC0 for x in b[125::126]):raise RuntimeError('BED_COMPLETE_ROWS_OR_PADDING_FAILED')
            for p in [bim,bed]:
                if hashes(p,guard)['sha256']!=files[p.name]['sha256']:raise RuntimeError('SOURCE_MEMBER_CHANGED_DURING_FORMAT_CHECK')
            counts[str(chrom)]=n;allbed[str(chrom)]=str(bed)
        if regular_sha(panel)!=plan['reference_body_sha256'][str(panel)]:raise RuntimeError('PANEL_CHANGED_AFTER_PARSE')
    finally:db.close()
    index_sha=regular_sha(path)
    record_output(path,index_sha)
    save(root/'qc/G3_format.json',dict(status='EXACT_SOURCE_FORMAT_PASS',source_index_sha256=index_sha,
        chromosome_variant_counts=counts,sample_order_sha256=hashlib.sha256(line(order)).hexdigest(),
        samples=503,all22_fam_byte_identical=True,padding_all_rows_checked=True,
        full_archive_build='UNRESOLVED_NOT_ASSERTED',selected_GRCh37_compatibility_required_before_QC=True))
    return path,index_sha,allbed


def prepare_static(plan,root,source_index,index_sha,guard):
    # Entire static universe is retained and sealed before genotype QC.
    source=sqlite3.connect('file:'+str(source_index)+'?mode=ro',uri=True)
    source.execute('PRAGMA cache_size=-32768');p=root/'qc/reference_index.sqlite';db=configure_db(p)
    db.execute('PRAGMA journal_mode=OFF');db.execute('PRAGMA max_page_count='+str(CAPACITY['coordinate_SQLite_max_pages']))
    if db.execute('PRAGMA page_size').fetchone()[0]!=CAPACITY['coordinate_SQLite_page_bytes']:raise RuntimeError('COMPACT_SQLITE_EXACT_PAGE_SIZE')
    db.execute('CREATE TABLE coord(SNP TEXT PRIMARY KEY,CHR INT,BP INT,A1 TEXT,A2 TEXT)')
    coordinate=plan['coordinate_map_path']
    if regular_sha(coordinate)!=plan['reference_body_sha256'][coordinate]:raise RuntimeError('COORDINATE_MAP_CHANGED_BEFORE_PARSE')
    with gzip.open(regular(coordinate),'rt') as f:
        reader=csv.DictReader(f,delimiter='\t')
        if reader.fieldnames!=['SNP','CHR','BP','A1','A2']:raise RuntimeError('PINNED_COORDINATE_SCHEMA_REQUIRED')
        for r in reader:
            guard();db.execute('INSERT INTO coord VALUES(?,?,?,?,?)',(r['SNP'],int(r['CHR']),int(r['BP']),r['A1'],r['A2']))
    db.commit()
    if regular_sha(coordinate)!=plan['reference_body_sha256'][coordinate]:raise RuntimeError('COORDINATE_MAP_CHANGED_AFTER_PARSE')
    db.execute('CREATE TABLE master(SNP TEXT PRIMARY KEY)');exclusions=Counter();rank=0
    schema_path=root/'qc/compact_schema.json';schema_sha=save(schema_path,SCHEMA);load_frozen(schema_path,schema_sha)
    writer=StaticWriter(root/'qc',plan['expected_master_reference_rows'])
    try:
        for chrom in range(1,23):
            path=plan['reference_files_numeric_order'][chrom-1]
            if regular_sha(path)!=plan['reference_body_sha256'][path]:raise RuntimeError('MASTER_REFERENCE_CHANGED_BEFORE_PARSE')
            with gzip.open(regular(path),'rt') as f:
                header=f.readline().split()
                if not all(k in header for k in ['CHR','SNP','BP']):raise RuntimeError('MASTER_REFERENCE_SCHEMA')
                previous=0
                for text in f:
                    guard();tokens=text.split();r=dict(zip(header,tokens))
                    if len(tokens)!=len(header) or r['CHR']!=str(chrom) or not r['BP'].isdigit() or int(r['BP'])<previous or int(r['BP'])<=0:raise RuntimeError('MASTER_REFERENCE_ORDER_OR_COORDINATE')
                    previous=int(r['BP']);db.execute('INSERT INTO master VALUES(?)',(r['SNP'],))
                    item=dict(SNP=r['SNP'],CHR=chrom,BP=previous,master_rank0=rank,master_rank1=rank+1);rank+=1
                    c=db.execute('SELECT CHR,BP,A1,A2 FROM coord WHERE SNP=?',(r['SNP'],)).fetchone()
                    s=source.execute('SELECT CHR,BP,A1,A2,row0,offset FROM source WHERE SNP=?',(r['SNP'],)).fetchone()
                    reason=None
                    if chrom==6 and 25000000<=previous<=34000000:reason='MHC_INCLUSIVE'
                    elif c is None:reason='MISSING_COORDINATE_MAP_ID'
                    elif (c[0],c[1])!=(chrom,previous):reason='REFERENCE_COORDINATE_MAP_CONFLICT'
                    elif s is None:reason='MISSING_SOURCE_ID'
                    elif (s[0],s[1])!=(chrom,previous):reason='SOURCE_GRCh37_COORDINATE_CONFLICT'
                    else:
                        sign,reason=allele_sign(c[2:4],s[2:4])
                        if reason is None:item.update(target_A1=c[2],target_A2=c[3],source_A1=s[2],source_A2=s[3],sign=sign,source_row0=s[4],BED_byte_offset=s[5],source_member=f'1000G.EUR.QC.{chrom}.bed')
                    if reason:item['exclusion']=reason;exclusions[reason]+=1
                    writer.append(item)
            if regular_sha(path)!=plan['reference_body_sha256'][path]:raise RuntimeError('MASTER_REFERENCE_CHANGED_AFTER_PARSE')
            db.commit()
        binding=writer.close()
    except BaseException:
        writer.abort();raise
    finally:source.close();db.close()
    with p.open('rb') as file:os.fsync(file.fileno())
    if regular_sha(source_index)!=index_sha:raise RuntimeError('STATIC_INPUT_CHANGED')
    record_output(p,regular_sha(p))
    if rank!=plan['expected_master_reference_rows'] or binding['static_rows']==0:raise RuntimeError('MASTER_UNIVERSE_OR_EMPTY_STATIC_SELECTION')
    # Exact canonical original stream identities are reconstructed on each read.
    save(root/'qc/G4_static_seal.json',dict(status='FROZEN_BEFORE_GENOTYPE_QC',master_reference_rows=rank,
      static_rows=binding['static_rows'],static_ordered_sha256=binding['logical_sha256']['static_selection'],
      exclusion_counts=dict(exclusions),static_exclusions_sha256=binding['logical_sha256']['static_exclusions'],
      compact_static=binding,compact_schema_file_sha256=schema_sha,
      compatibility='EXACT_SELECTED_ID_CHR_BP_ORDERED_OR_SWAPPED_ALLELES_TO_PINNED_GRCh37;FULL_SOURCE_BUILD_UNRESOLVED'))
    return binding




def read_dose(rows,bedpaths,np):
    if not 0<len(rows)<=2048 or len({r['CHR'] for r in rows})!=1:raise RuntimeError('SINGLE_CHROMOSOME_BOUNDED_READ_REQUIRED')
    path=regular(bedpaths[str(rows[0]['CHR'])]);payload=bytearray()
    with path.open('rb') as f:
        for r in rows:
            if r['BED_byte_offset']!=3+r['source_row0']*126:raise RuntimeError('EXACT_SOURCE_OFFSET_FAILED')
            f.seek(r['BED_byte_offset']);b=f.read(126)
            if len(b)!=126:raise RuntimeError('SOURCE_BED_SHORT_SELECTED_ROW')
            payload.extend(b)
    return decode_a1(bytes(payload),np)


def genotype_qc(root,static_binding,bedpaths,np,guard):
    excludes=Counter();chromcounts=Counter();index=0
    writer=GenotypeWriter(root/'qc',static_binding)
    def one_group(group):
        nonlocal index
        doses=read_dose(group,bedpaths,np)
        for i,row in enumerate(group):
            q=qc_column(doses[:,i],np);item={**row,**q}
            if q['exclusion']:excludes[q['exclusion']]+=1
            else:
                standardized(doses[:,i],row['sign'],np);item['J_index0']=index;index+=1;chromcounts[str(row['CHR'])]+=1
            writer.append(item)
    try:
        group=[]
        for row in static_rows(root/'qc',static_binding,guard):
            if group and (len(group)==2048 or row['CHR']!=group[0]['CHR']):one_group(group);group=[]
            group.append(row)
        if group:one_group(group)
        binding=writer.close()
    except BaseException:
        writer.abort();raise
    if binding['J_rows']!=index:raise RuntimeError('INTENDED_J_COUNT_CHANGED')
    if any(chromcounts[str(c)]<192 for c in range(1,23)):raise RuntimeError('G5_INCOMPLETE_SOURCE_CHROMOSOME_FEWER_THAN192')
    save(root/'qc/G4_J_seal.json',dict(status='REFERENCE_ONLY_J_FROZEN',J_count=index,J_ordered_sha256=binding['logical_sha256']['J_ordered'],
      chromosome_J_counts=dict(chromcounts),genotype_exclusion_counts=dict(excludes),genotype_exclusions_sha256=binding['logical_sha256']['genotype_exclusions'],
      compact_genotype=binding,no_MAF_threshold=True,no_GWAS_outcomes=True,no_native_mask_selection=True))
    return binding,chromcounts




def actual_controls(plan,root,static_binding,qc_binding,bedpaths,np,guard,chromcounts):
    jsha=qc_binding['logical_sha256']['J_ordered']
    controls=root/'qc/controls';controls.mkdir();records=[]
    # Only metadata rows, never full chromosome dense genotypes, are retained.
    frozen={str(c):frozen_control_indices(chromcounts[str(c)]) for c in range(1,23)}
    selection=controls/'fixed_G5_selections.json';selection_sha=save(selection,dict(J_sha256=jsha,indices=frozen,sign_flip='ZERO_BASED_EVERY_THIRD_SELECTED_COLUMN'))
    load_frozen(selection,selection_sha)
    def read_g(rows):
        dose=read_dose(rows,bedpaths,np)
        return np.column_stack([standardized(dose[:,i],r['sign'],np) for i,r in enumerate(rows)])
    for chrom in range(1,23):
        guard(force=True);chrrows=[r for r in genotype_rows(root/'qc',static_binding,qc_binding,guard) if r['CHR']==chrom];choice=frozen[str(chrom)];indices=sorted(set(sum(choice.values(),[])));selected=[chrrows[i] for i in indices]
        window=controls/f'chr{chrom}.bed';payload=bytearray(BED_HEADER)
        with regular(bedpaths[str(chrom)]).open('rb') as f:
            for r in selected:f.seek(r['BED_byte_offset']);payload.extend(f.read(126))
        wsha=binary_new(window,bytes(payload));m=controls/f'chr{chrom}.metadata.json';output=controls/f'chr{chrom}.B.npz';receipt=controls/f'chr{chrom}.B.receipt.json'
        spec=dict(samples=503,rows=selected,bed_path=str(window),bed_sha256=wsha,maintained_reader_path=plan['maintained_decoderB_path'],output_path=str(output),receipt_path=str(receipt))
        msha=save(m,spec);load_frozen(m,msha)
        command=[plan['decoderB_python'],'-B',plan['decoderB_wrapper'],'--metadata',str(m),'--metadata-sha256',msha]
        with (controls/f'chr{chrom}.B.log').open('xb') as log:
            rc=subprocess.run(command,stdout=log,stderr=subprocess.STDOUT,check=False).returncode
        guard(force=True)
        if rc:raise RuntimeError('INDEPENDENT_MAINTAINED_SOURCE_DECODER_FAILED: '+str(chrom))
        rsha=regular_sha(receipt);br=load_frozen(receipt,rsha)
        if br['metadata_sha256']!=msha or br['source_window_sha256']!=wsha or br['maintained_reader_sha256']!=plan['maintained_decoderB_sha256'] or br['kept_source_rows']!=list(range(len(selected))) or br['rows']!=selected:raise RuntimeError('INDEPENDENT_DECODER_B_RECEIPT_BINDING')
        osha=br['output_sha256']
        if regular_sha(output)!=osha:raise RuntimeError('DECODER_B_ARRAY_CHANGED_BEFORE_PARSE')
        with np.load(regular(output),allow_pickle=False) as b:
            raw=b['A2_dose'];bG=b['standardized_A2']
        if regular_sha(output)!=osha or regular_sha(receipt)!=rsha:raise RuntimeError('DECODER_B_OUTPUT_CHANGED_AFTER_PARSE')
        record_output(receipt,rsha);record_output(output,osha);record_output(controls/f'chr{chrom}.B.log',regular_sha(controls/f'chr{chrom}.B.log'))
        A=read_dose(selected,bedpaths,np);converted=np.where(raw==9,-1,2-raw)
        if not np.array_equal(A,converted):raise RuntimeError('INDEPENDENT_DISCRETE503_SOURCE_DECODE_DIFFERS')
        G=read_g(selected);sign=np.array([r['sign'] for r in selected]);B=-bG*sign
        if not close_array(G,B,np):raise RuntimeError('INDEPENDENT_STANDARDIZED_SIGNED_SOURCE_DECODE_DIFFERS')
        checked=[]
        for name,ix in {**choice,'combined_distant_and_windows':indices}.items():
            rows=[chrrows[i] for i in ix];g=read_g(rows)
            independent_g=B[:,[indices.index(i) for i in ix]];R=independent_g.T@independent_g/503
            v=np.array([((r['master_rank0']*17+3)%23-11)/11 for r in rows],dtype=np.float64);want=R@v
            if not np.any(v):raise RuntimeError('G5_INCOMPLETE_DETERMINISTIC_WINDOW_PROBE_ZERO_NO_REPLACEMENT')
            for chunk in [1,7,64,2048]:
                for rev in [False,True]:
                    result,_=operator(rows,read_g,v,np,chunk,rev,guard)
                    if not close_array(result,want,np):raise RuntimeError('DIRECT_VS_TWO_PASS_OPERATOR_FAILED')
            d=np.where(np.arange(len(rows))%3==0,-1.,1.)
            flipped=read_g([{**r,'sign':r['sign']*int(d[i])} for i,r in enumerate(rows)]);Rnew=flipped.T@flipped/503
            if not close_array(flipped,g*d,np) or not close_array(Rnew,d[:,None]*R*d[None,:],np):raise RuntimeError('PREDETERMINED_SIGNED_R_NOT_SQUARED_R_FAILED')
            if np.linalg.matrix_rank(g)>502 or not np.all(np.isfinite(R)) or np.any(np.abs(np.diag(R)-1)>1e-10):raise RuntimeError('FINITE_REFERENCE_CENTERED_GRAM_RANK_OR_DIAGONAL')
            checked.append(name)
        v=np.array([((r['master_rank0']*17+3)%23-11)/11 for r in chrrows],dtype=np.float64)
        w=np.array([((r['master_rank0']*19+7)%29-14)/14 for r in chrrows],dtype=np.float64)
        if not np.any(v) or not np.any(w):raise RuntimeError('G5_INCOMPLETE_GLOBAL_PROBE_ZERO_NO_REPLACEMENT')
        rv,gv=operator(chrrows,read_g,v,np,2048,False,guard);rw,gw=operator(chrrows,read_g,w,np,2048,True,guard)
        if not close_array(np.array([v@rv]),np.array([gv@gv/503]),np) or not close_array(np.array([w@rv]),np.array([v@rw]),np):raise RuntimeError('FULL_CHROMOSOME_ENERGY_OR_SYMMETRY_FAILED')
        records.append(dict(CHR=chrom,J_count=len(chrrows),independent_rows=len(selected),source_window_sha256=wsha,B_array_sha256=osha,B_receipt_sha256=rsha,
            controls=checked,full_chromosome_energy_symmetry=True,between_storage_window_correlations_retained=True))
    # Every reconstruction checks immutable physical and original logical hashes.
    save(root/'qc/G5_operator_controls.json',dict(status='ACTUAL_AUTHENTICATED_SOURCE_CONTROLS_PASS',fixed_selection_sha256=selection_sha,chromosomes=records,
      maintained_B_padding_independent_rejection=False,zero_variance_fallback_not_credited=True,
      finite_reference_503_conditional_only=True,population_LD_or_actual_GWAS_mask_adequacy_established=False,S_simulations_run=0,C_calibration_credit=0,estimator_calls=0,biological41Cov_release_admitted=False))


def run(plan_path,plan_sha,started,admission_path,admission_sha,heavy_fd,raw_fd):
    plan=load_frozen(plan_path,plan_sha)
    validate_inherited_lock(heavy_fd)
    raw=regular(plan['raw_transfer_family_mutex']);stat=os.fstat(raw_fd);current=raw.stat()
    if (stat.st_dev,stat.st_ino)!=(current.st_dev,current.st_ino):raise RuntimeError('SOURCE_INHERITED_RAW_TRANSFER_LOCK_IDENTITY_DIFFERS')
    fcntl.flock(raw_fd,fcntl.LOCK_EX|fcntl.LOCK_NB)
    identity_gate(plan,plan_path,admission_path,plan_sha,admission_sha)
    import numpy as np
    if sys.version_info[:3]!=(3,11,11) or np.__version__!='1.26.4':raise RuntimeError('SOURCE_DECLARED_A_RUNTIME_VERSION_DIFFERS')
    root=Path(plan['output_root']);guard=Guard(plan,started)
    archive,identity=acquire(plan,root,guard)
    guard(force=True);extracted=extract_strict66(archive,root/'extracted/members',plan['resource_contract']['extracted_all_partial_and_failed_cap_bytes'],guard)
    save(root/'qc/archive_members.json',extracted)
    source,index_sha,bedpaths=source_format(plan,root,extracted,guard)
    static_binding=prepare_static(plan,root,source,index_sha,guard)
    qc_binding,counts=genotype_qc(root,static_binding,bedpaths,np,guard)
    actual_controls(plan,root,static_binding,qc_binding,bedpaths,np,guard,counts)
    jsha=qc_binding['logical_sha256']['J_ordered']
    # Every input source member/archive rehashed after the entire source stage.
    for item in extracted['files'].values():
        if hashes(item['path'],guard)['sha256']!=item['sha256']:raise RuntimeError('SOURCE_MEMBER_CHANGED_AFTER_CONTROLS')
        record_output(item['path'],item['sha256'])
    if hashes(archive,guard)!=identity:raise RuntimeError('AUTHENTICATED_ARCHIVE_CHANGED_AFTER_USE')
    guard(force=True)
    verify_outputs(INTENDED_OUTPUT_SHA256)
    save(root/'qc/source_worker_result.json',dict(status='SOURCE_G2_G5_PASS_PRIVATE_PROVISIONAL',plan_sha256=plan_sha,
        archive_identity=identity,J_ordered_sha256=jsha,chromosome_J_counts=dict(counts),
        source_index_sha256=index_sha,full_source_build_unresolved=True,old_compact_identity_not_claimed=True,
        ordinary_output_sha256=dict(INTENDED_OUTPUT_SHA256),
        S_simulations_run=0,C_calibration_credit=0,estimator_calls=0,biological41Cov_release_admitted=False))
    guard(force=True)


if __name__=='__main__':
    p=argparse.ArgumentParser(allow_abbrev=False);p.add_argument('--plan',required=True);p.add_argument('--plan-sha256',required=True);p.add_argument('--started-monotonic',type=float,required=True)
    p.add_argument('--admission',required=True);p.add_argument('--admission-sha256',required=True);p.add_argument('--heavy-lock-fd',type=int,required=True);p.add_argument('--raw-lock-fd',type=int,required=True);a=p.parse_args()
    run(a.plan,a.plan_sha256,a.started_monotonic,a.admission,a.admission_sha256,a.heavy_lock_fd,a.raw_lock_fd)
