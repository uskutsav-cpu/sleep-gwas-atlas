"""Stdlib provenance gates for unchanged historical extension preprocessing."""
import csv
import gzip
import hashlib
import json
import math
import os
import plistlib
import subprocess
from pathlib import Path


def sha(path):
    h=hashlib.sha256()
    with Path(path).open('rb') as f:
        for b in iter(lambda:f.read(1<<20),b''):h.update(b)
    return h.hexdigest()


def body_hashes(path):
    h=hashlib.sha256();m=hashlib.md5();n=0
    with Path(path).open('rb') as f:
        for b in iter(lambda:f.read(1<<20),b''):h.update(b);m.update(b);n+=len(b)
    return dict(sha256=h.hexdigest(),md5=m.hexdigest(),bytes=n)


def write_new(path,data):
    payload=json.dumps(data,indent=2,allow_nan=False)+'\n'
    with Path(path).open('x') as f:f.write(payload);f.flush();os.fsync(f.fileno())


def check_bindings(plan,include_archived=False):
    for path,digest in plan['dependencies_sha256'].items():
        if path in plan['archived_input_sha256'] and not include_archived:continue
        if sha(path)!=digest:raise RuntimeError('FROZEN_DEPENDENCY_CHANGED: '+path)


def acquisition_receipt_gate(plan,member,expected_receipt_sha):
    path=Path(member['acquisition_receipt'])
    if sha(path)!=expected_receipt_sha:raise RuntimeError('ACQUISITION_RECEIPT_CHANGED')
    r=json.loads(path.read_text());m=member['acquisition_member'];a=plan['acquisition_operational_identity']
    if r['status']!='EXACT_IMMUTABLE_SOURCE_ACQUIRED' or r['member']!=m or r['plan_sha256']!=plan['acquisition_plan_sha256']:
        raise RuntimeError('RAW_SOURCE_NOT_EXACT_SUCCESSFUL_CHECKPOINT')
    if r['returncode']!=0 or r['stop_reason'] is not None or not r['teardown']['teardown_verified'] or r['teardown'].get('remaining_group_members',[]) or r['teardown'].get('cleanup_error'):
        raise RuntimeError('ACQUISITION_EXIT_OR_TEARDOWN_UNCLEARED')
    if not r['post_cleanup_hash_resource_identity_gates_pass']:raise RuntimeError('ACQUISITION_FINAL_SEAL_UNCLEARED')
    expected=(m['expected_bytes'],m['expected_md5'],m['expected_sha256'])
    if (r['actual_size'],r['actual_md5'],r['actual_sha256'])!=expected or (r['final_seal_md5'],r['final_seal_sha256'])!=expected[1:]:raise RuntimeError('ACQUISITION_SOURCE_IDENTITY_DIFFERS')
    offset=a['resume_prefix_bytes'] if m['index']==1 else 0
    if r['resume_offset']!=offset:raise RuntimeError('ACQUISITION_RESUME_OFFSET_DIFFERS')
    logbase=Path(a['source_folder'])/'logs'/m['extension_trait_id']
    command=[a['curl'],'-q','--fail','--location','--max-redirs','5','--proto','=https','--proto-redir','=https','--tlsv1.2',
             '--max-time',str(a['per_body_seconds_limit']),'--max-filesize',str(m['expected_bytes']),
             '--dump-header',str(logbase)+'.headers.txt','--output',m['body_path']+'.partial']
    if offset:command+=['--continue-at',str(offset)]
    command.append(m['url'])
    if r['command']!=command:raise RuntimeError('ACQUISITION_EXACT_COMMAND_DIFFERS')
    h=r['observed_headers']
    if h['status']!=(206 if offset else 200) or h['x-amz-version-id']!=m['expected_s3_version_id'] or h['etag']!=m['expected_etag'] or int(h['content-length'])!=m['expected_bytes']-offset:
        raise RuntimeError('ACQUISITION_HTTP_IDENTITY_DIFFERS')
    if offset and h['content-range']!='bytes '+str(offset)+'-'+str(m['expected_bytes']-1)+'/'+str(m['expected_bytes']):raise RuntimeError('ACQUISITION_EXACT_RANGE_DIFFERS')
    for state in [r['resource_before'],r['resource_after']]:
        if state['internal_free_bytes']<plan['guard']['internal_floor_bytes'] or state['ssd_free_bytes']<plan['guard']['SSD_floor_bytes']:raise RuntimeError('ACQUISITION_RESOURCE_FLOORS_NOT_PRESERVED')
    body=Path(m['body_path'])
    if body.is_symlink() or not body.is_file() or body.stat().st_size!=m['expected_bytes']:raise RuntimeError('RAW_BODY_MISSING_OR_NOT_REGULAR_EXACT_SIZE')
    return r


def compare_munged(new,archived,expected_rows):
    """Compare entire decompressed streams and every schema field; always reach EOF."""
    digests=[hashlib.sha256(),hashlib.sha256()];rows=[0,0];counts=[{'finite_N_Z':0,'missing_or_nonfinite_N_Z':0},{'finite_N_Z':0,'missing_or_nonfinite_N_Z':0}]
    unequal=0;first=[];header=None
    with gzip.open(new,'rb') as a,gzip.open(archived,'rb') as b:
        heads=[a.readline(),b.readline()]
        if heads[0]!=heads[1] or heads[0].rstrip(b'\r\n')!=b'SNP\tA1\tA2\tZ\tN':raise RuntimeError('MUNGED_CANONICAL_HEADER_DIFFERS')
        header=heads[0].decode().rstrip('\r\n')
        for h,x in zip(digests,heads):h.update(x)
        while True:
            lines=[a.readline(),b.readline()]
            if not any(lines):break
            parsed=[]
            for i,line in enumerate(lines):
                if not line:parsed.append(None);continue
                rows[i]+=1;digests[i].update(line)
                fields=line.rstrip(b'\r\n').split(b'\t')
                if len(fields)!=5:raise RuntimeError('MUNGED_ROW_WIDTH_DIFFERS')
                parsed.append(fields)
                try:finite=all(math.isfinite(float(x)) for x in fields[3:5])
                except ValueError:finite=False
                counts[i]['finite_N_Z' if finite else 'missing_or_nonfinite_N_Z']+=1
            if lines[0]!=lines[1]:
                unequal+=1
                if len(first)<10:first.append(dict(row=max(rows)+1,changed_columns=[n for i,n in enumerate(['SNP','A1','A2','Z','N']) if parsed[0] is None or parsed[1] is None or parsed[0][i]!=parsed[1][i]],literal_line_or_order_changed=True))
    return dict(status='EXACT_FULL_DECOMPRESSED_TEMPLATE_MATCH' if unequal==0 and rows==[expected_rows,expected_rows] else 'DECOMPRESSED_TEMPLATE_MISMATCH_PRESERVED',
                header=header,rows_new=rows[0],rows_archived=rows[1],expected_template_rows=expected_rows,
                unequal_rows=unequal,first_mismatch_diagnostics=first,counts_new=counts[0],counts_archived=counts[1],
                decompressed_sha256_new=digests[0].hexdigest(),decompressed_sha256_archived=digests[1].hexdigest(),
                entire_ordered_SNP_allele_N_Z_and_literal_missingness_compared=True,both_gzip_CRC_and_EOF_verified=True,
                compressed_sha256_new=sha(new),compressed_sha256_archived=sha(archived),
                compressed_byte_identity_required=False,numerical_tolerance_for_munged_comparison=None)


def harmonized_receipt_gate(plan,member):
    r=json.loads(Path(member['harmonization_receipt']).read_text());v=r['source_verification'];m=member['acquisition_member']
    if r['extension_trait_id']!=m['extension_trait_id'] or r['source_url']!='LOCAL_FILE' or r['panel_path']!=plan['panel'] or r['source_file_name']!=m['filename']:
        raise RuntimeError('NEW_HARMONIZATION_RECEIPT_SOURCE_OR_PANEL_DIFFERS')
    if v['verification_status']!='PASS' or v['source']!=m['body_path'] or (v['observed_size_bytes'],v['observed_md5'],v['observed_sha256'])!=(m['expected_bytes'],m['expected_md5'],m['expected_sha256']):
        raise RuntimeError('FULL_RAW_EOF_HASHES_FROM_HARMONIZER_DIFFER')
    if r['harmonized_output']!=member['harmonized'] or r['harmonization_qc']!=member['harmonization_qc'] or r['harmonized_output_sha256']!=sha(member['harmonized']):raise RuntimeError('HARMONIZED_CURRENT_CONTENT_DIFFERS')
    with Path(member['harmonization_qc']).open() as f:qc={x['metric']:x['count'] for x in csv.DictReader(f,delimiter='\t')}
    for key,expected in [('reference_sha256',plan['reference_sha256']),('source_sha256',m['expected_sha256']),('source_md5',m['expected_md5']),('output_sha256',r['harmonized_output_sha256'])]:
        if qc.get(key)!=expected:raise RuntimeError('NEW_HARMONIZATION_QC_HASH_DIFFERS: '+key)
    if int(qc['source_size_bytes'])!=m['expected_bytes'] or int(qc['output_rows'])<=0:raise RuntimeError('HARMONIZATION_QC_INVALID_SIZE_OR_EMPTY')
    return dict(receipt=r,new_filter_counts={k:int(v) for k,v in qc.items() if v.isdigit() and not k.endswith('_sha256')},
                original_filter_count_concordance='NOT_ESTABLISHED_ORIGINAL_QC_TSV_UNAVAILABLE',
                original_harmonized_compressed_SHA_not_a_reproduction_target=True)


def physical_mount():
    info=plistlib.loads(subprocess.run(['/usr/sbin/diskutil','info','-plist','/Volumes/Extreme SSD'],capture_output=True,check=True).stdout)
    if not os.path.ismount('/Volumes/Extreme SSD') or info.get('VolumeUUID','').upper()!='77FD98CC-B09E-3DAA-ACC0-82FF2B776B28' or info.get('ReadOnlyVolume',False):
        raise RuntimeError('EXACT_WRITABLE_PHYSICAL_EXTREME_SSD_NOT_PRESENT')
    return dict(mount='/Volumes/Extreme SSD',VolumeUUID=info['VolumeUUID'].upper(),read_only=False)
