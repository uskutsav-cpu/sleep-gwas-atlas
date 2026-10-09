#!/usr/bin/env python3
"""Independent bounded compressed/content audit; no native helper invocation."""
import argparse
from datetime import datetime, timezone
import gzip
import hashlib
import json
from pathlib import Path
import resource
import shutil
import struct
import subprocess
import sys
import time
import zlib

HERE=Path(__file__).resolve().parent
PACKAGE=HERE.parent
ROOT=PACKAGE.parent
SSD=Path('/Volumes/Extreme SSD/sleep-unified-research-v1/recovery-2026-10-09/independent_cancer_content_review_v3')
PLAN=HERE/'cancer_materialization_independent_resource_plan_v3.json'
TRAITS=['breast_cancer','ovarian_cancer']
BUFFER=65536
BEGAN=time.monotonic()
LAST_GUARD=0.
PEAK=0


def save(path,value):
    with path.open('x') as f:json.dump(value,f,indent=2);f.write('\n')


def guard(force=False):
    global LAST_GUARD,PEAK
    now=time.monotonic()
    if not force and now-LAST_GUARD<1:return
    peak=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    PEAK=max(PEAK,peak if sys.platform=='darwin' else peak*1024)
    if PEAK>128*1024**2 or shutil.disk_usage(ROOT).free<128*1024**2 or shutil.disk_usage(SSD.parent).free<1024**3 or now-BEGAN>1800:
        raise RuntimeError('INDEPENDENT_CONTENT_REVIEW_RESOURCE_GUARD')
    LAST_GUARD=now


def file_digest(path):
    before=path.stat();digest=hashlib.sha256()
    with path.open('rb') as f:
        while True:
            block=f.read(BUFFER)
            if not block:break
            digest.update(block);guard()
    after=path.stat()
    return {'path':str(path),'bytes':after.st_size,'sha256':digest.hexdigest(),
            'size_mtime_inode_stable':(before.st_size,before.st_mtime_ns,before.st_ino)==(after.st_size,after.st_mtime_ns,after.st_ino)}


def freeze():
    if shutil.disk_usage(ROOT).free<256*1024**2 or shutil.disk_usage(SSD.parent).free<1024**3:
        raise RuntimeError('REVIEW_LAUNCH_RESOURCE_GATE')
    if PLAN.exists() or SSD.exists():raise RuntimeError('PRIOR_REVIEW_NAMESPACE_PRESERVED')
    SSD.mkdir();(SSD/'tmp').mkdir()
    native_plans={t:PACKAGE/'manifests'/(t+'_materialization_plan_v3.json') for t in TRAITS}
    plan={'frozen_utc':datetime.now(timezone.utc).isoformat(),'review_code_sha256':file_digest(Path(__file__))['sha256'],
          'native_plan_hashes':{t:file_digest(p)['sha256'] for t,p in native_plans.items()},
          'native_plan_paths':{t:str(p) for t,p in native_plans.items()},
          'traits':TRAITS,'workers':1,'network_bytes':0,'hash_and_inflate_read_buffer_bytes':BUFFER,
          'internal_launch_minimum_bytes':256*1024**2,'internal_free_at_freeze_bytes':shutil.disk_usage(ROOT).free,
          'internal_emergency_floor_bytes':128*1024**2,'observed_RSS_limit_bytes':128*1024**2,
          'SSD_runtime_reserve_bytes':1024**3,'maximum_seconds':1800,'maximum_new_retained_review_bytes':8*1024**2,
          'TMPDIR':str(SSD/'tmp'),'resource_basis':'Hashes and gzip content use 64 KiB buffers; no row dictionaries, reference loads, native helper or estimator.',
          'source_and_outputs_read_only':True,'scope':'Original helper/code and receipt binding; exact old/new gzip content and framing comparison only',
          'complete_native_QC_or_biological_findings_admitted':False,'h2_or_rg_estimated':False}
    save(PLAN,plan)
    print(json.dumps({'status':'REVIEW_RESOURCE_PLAN_FROZEN','sha256':file_digest(PLAN)['sha256']}),flush=True)


def gzip_content(path):
    before=path.stat();digest=hashlib.sha256();size=0;newlines=0;crc=0;last=b''
    with gzip.open(path,'rb') as f:
        while True:
            block=f.read(BUFFER)
            if not block:break
            digest.update(block);size+=len(block);newlines+=block.count(b'\n');crc=zlib.crc32(block,crc);last=block[-1:];guard()
    after=path.stat()
    return {'decompressed_sha256':digest.hexdigest(),'decompressed_bytes':size,
            'newline_count_including_header':newlines,'full_gzip_crc_pass':True,
            'actual_line_count_including_header':newlines+int(size>0 and last!=b'\n'),
            'ends_with_newline':last==b'\n','independent_CRC32':crc & 0xffffffff,
            'size_mtime_inode_stable':(before.st_size,before.st_mtime_ns,before.st_ino)==(after.st_size,after.st_mtime_ns,after.st_ino)}


def gzip_framing(path):
    # Hash the exact compressed bytes after the first complete gzip header.
    # Equality here proves every non-initial-header byte is identical, even
    # without assuming that gzip_content's decoded equality implies this.
    with path.open('rb') as f:
        fixed=f.read(10)
        if len(fixed)!=10 or fixed[:3]!=b'\x1f\x8b\x08' or fixed[3]&0xe0:raise RuntimeError('INVALID_GZIP_HEADER')
        flags=fixed[3];metadata={'flags':flags,'mtime':int.from_bytes(fixed[4:8],'little'),'xfl':fixed[8],'os':fixed[9]}
        extra=b'';filename=b'';comment=b'';header_crc=b''
        if flags&4:
            width=f.read(2)
            if len(width)!=2:raise RuntimeError('TRUNCATED_GZIP_EXTRA')
            extra=f.read(int.from_bytes(width,'little'))
            if len(extra)!=int.from_bytes(width,'little'):raise RuntimeError('TRUNCATED_GZIP_EXTRA')
        def terminated():
            value=bytearray()
            while True:
                b=f.read(1)
                if not b:raise RuntimeError('TRUNCATED_GZIP_HEADER_STRING')
                if b==b'\0':return bytes(value)
                value.extend(b)
                if len(value)>1024**2:raise RuntimeError('GZIP_HEADER_RESOURCE_CAP')
        if flags&8:filename=terminated()
        if flags&16:comment=terminated()
        if flags&2:
            header_crc=f.read(2)
            if len(header_crc)!=2:raise RuntimeError('TRUNCATED_HEADER_CRC')
        offset=f.tell();total=path.stat().st_size
        if total-offset<8:raise RuntimeError('MISSING_GZIP_TRAILER')
        payload=hashlib.sha256();remaining=total-offset-8
        while remaining:
            b=f.read(min(BUFFER,remaining))
            if not b:raise RuntimeError('TRUNCATED_GZIP_PAYLOAD')
            payload.update(b);remaining-=len(b);guard()
        trailer=f.read(8)
        metadata.update({'extra_hex':extra.hex(),'filename_latin1':filename.decode('latin1'),'comment_latin1':comment.decode('latin1'),
                         'header_crc_hex':header_crc.hex(),'header_bytes':offset})
        return {'header':metadata,'nonheader_deflate_bytes':total-offset-8,
                'nonheader_deflate_sha256':payload.hexdigest(),'trailer_hex':trailer.hex(),
                'trailer_CRC32':struct.unpack('<I',trailer[:4])[0],'trailer_ISIZE':struct.unpack('<I',trailer[4:])[0]}


def verify():
    if not PLAN.exists():raise RuntimeError('FROZEN_REVIEW_PLAN_REQUIRED')
    plan=json.loads(PLAN.read_text());plan_before=file_digest(PLAN)
    if plan['review_code_sha256']!=file_digest(Path(__file__))['sha256']:raise RuntimeError('REVIEW_CODE_CHANGED_AFTER_FREEZE')
    if shutil.disk_usage(ROOT).free<plan['internal_launch_minimum_bytes']:raise RuntimeError('REVIEW_LAUNCH_RESOURCE_GATE')
    guard(True)
    checks=[];results={}
    def check(name,actual,expected):checks.append({'name':name,'pass':actual==expected,'actual':actual,'expected':expected})
    for trait in TRAITS:
        p=Path(plan['native_plan_paths'][trait]);native=json.loads(p.read_text())
        r=PACKAGE/'logs'/(trait+'_materialization_receipt_v3.json')
        if not r.exists():raise RuntimeError('WAIT_FOR_PARENT_COMPLETION_SIGNAL_AND_RECEIPT')
        receipt=json.loads(r.read_text());receipt_hash=file_digest(r)
        check(trait+'_native_success',receipt['status'],'NATIVE_RAW_MATERIALIZATION_EXACT_CONTENT')
        check(trait+'_native_plan_frozen_at_review_freeze',file_digest(p)['sha256'],plan['native_plan_hashes'][trait])
        check(trait+'_native_plan_receipt_binding',file_digest(p)['sha256'],receipt['plan_sha256'])
        check(trait+'_no_outcomes_or_full_QC_admitted',not receipt['h2_or_rg_estimated'] and not receipt['full_QC_chain_reproduced'],True)
        paths={k:Path(v) for k,v in native['bound_paths'].items()}
        before={k:file_digest(v) for k,v in paths.items()}
        check(trait+'_all_actual_dependencies_match_frozen_plan',{k:v['sha256'] for k,v in before.items()},native['inputs_sha256_before'])
        check(trait+'_actual_dependencies_match_execution_before_after',{k:v['sha256'] for k,v in before.items()},receipt['inputs_sha256_before'])
        check(trait+'_native_before_after_hashes_match',receipt['inputs_sha256_before'],receipt['inputs_sha256_after'])
        git={}
        for rel,registered in native['unchanged_original_main_helper_Git_blobs'].items():
            original=subprocess.check_output(['git','rev-parse','659d01cf:'+rel],cwd=ROOT,text=True).strip()
            head=subprocess.check_output(['git','rev-parse','HEAD:'+rel],cwd=ROOT,text=True).strip()
            disk=subprocess.check_output(['git','hash-object',str(ROOT/rel)],cwd=ROOT,text=True).strip()
            git[rel]={'original_main_blob':original,'HEAD_blob':head,'actual_disk_blob':disk}
            check(trait+'_unchanged_original_git_'+rel,[original,head,disk],[registered]*3)
        acquisition=json.loads(paths['acquisition_receipt'].read_text())
        check(trait+'_source_registry_acquisition_full_SHA',before['original_container']['sha256'],acquisition['expected_sha256'])
        check(trait+'_source_registry_acquisition_exact_bytes',before['original_container']['bytes'],acquisition['expected_bytes'])
        old=paths['historical_raw'];new=Path(receipt['output_path'])
        old_compressed=before['historical_raw'];new_compressed=file_digest(new)
        check(trait+'_new_compressed_SHA_receipt',new_compressed['sha256'],receipt['output_sha256'])
        equal_compressed=old_compressed['sha256']==new_compressed['sha256']
        check(trait+'_gzip_byte_equality_claim',equal_compressed,receipt['gzip_bytes_identical'])
        old_decoded=gzip_content(old);new_decoded=gzip_content(new)
        for label,decoded in [('historical',old_decoded),('native',new_decoded)]:
            reduced={k:decoded[k] for k in ['decompressed_sha256','decompressed_bytes','newline_count_including_header','full_gzip_crc_pass']}
            check(trait+'_'+label+'_decoded_receipt_binding',reduced,receipt[label+'_decompressed'])
            check(trait+'_'+label+'_newline_terminated',decoded['ends_with_newline'],True)
        check(trait+'_independent_decoded_byte_SHA_count_identity',old_decoded,new_decoded)
        old_frame=gzip_framing(old);new_frame=gzip_framing(new)
        for label,decoded,frame in [('historical',old_decoded,old_frame),('native',new_decoded,new_frame)]:
            check(trait+'_'+label+'_trailer_CRC_matches_independent_content_CRC',frame['trailer_CRC32'],decoded['independent_CRC32'])
            check(trait+'_'+label+'_trailer_ISIZE_matches_bytes_mod_2to32',frame['trailer_ISIZE'],decoded['decompressed_bytes']%(2**32))
        nonheader_equal=all(old_frame[k]==new_frame[k] for k in ['nonheader_deflate_bytes','nonheader_deflate_sha256','trailer_hex'])
        header_differences=[k for k in old_frame['header'] if old_frame['header'][k]!=new_frame['header'][k]]
        metadata_only=nonheader_equal and all(k in ['mtime','filename_latin1','header_bytes'] for k in header_differences)
        after={k:file_digest(v) for k,v in paths.items()}
        check(trait+'_all_bound_dependencies_unchanged_during_independent_review',{k:v['sha256'] for k,v in before.items()},{k:v['sha256'] for k,v in after.items()})
        check(trait+'_new_output_unchanged_after_review',file_digest(new)['sha256'],new_compressed['sha256'])
        check(trait+'_native_receipt_unchanged_after_review',file_digest(r)['sha256'],receipt_hash['sha256'])
        check(trait+'_all_file_stat_stability',all(v['size_mtime_inode_stable'] for v in list(before.values())+list(after.values())+[new_compressed]),True)
        check(trait+'_native_log_binding',file_digest(PACKAGE/'logs'/(trait+'_native_materialization_v3.log'))['sha256'],receipt['log_sha256'])
        results[trait]={'native_plan_sha256':file_digest(p)['sha256'],'native_receipt':receipt_hash,
                       'input_dependency_hashes_before':before,'input_dependency_hashes_after':after,'original_git_blobs':git,
                       'historical_compressed':old_compressed,'new_compressed':new_compressed,
                       'historical_content':old_decoded,'new_content':new_decoded,
                       'historical_gzip_framing':old_frame,'new_gzip_framing':new_frame,
                       'full_compressed_bytes_identical':equal_compressed,'decompressed_bytes_identical':old_decoded==new_decoded,
                       'all_bytes_after_initial_gzip_header_identical':nonheader_equal,
                       'initial_header_changed_fields':header_differences,
                       'only_initial_filename_time_metadata_differs':metadata_only,
                       'h2_or_rg_estimated':False,'full_QC_chain_reproduced':False}
        print(json.dumps({'trait':trait,'decoded_SHA':new_decoded['decompressed_sha256'],'decoded_bytes':new_decoded['decompressed_bytes'],
                          'lines':new_decoded['actual_line_count_including_header'],'compressed_equal':equal_compressed,'metadata_only_difference':metadata_only}),flush=True)
    check('independent_review_resource_plan_unchanged',file_digest(PLAN)['sha256'],plan_before['sha256'])
    guard(True)
    retained=sum(x.stat().st_size for x in SSD.rglob('*') if x.is_file())
    check('review_SSD_retained_within_plan',retained<=plan['maximum_new_retained_review_bytes'],True)
    result={'completed_utc':datetime.now(timezone.utc).isoformat(),'status':'INDEPENDENT_CANCER_CONTENT_STREAM_PASS' if all(x['pass'] for x in checks) else 'INDEPENDENT_CANCER_CONTENT_STREAM_FAILURE',
            'review_code':file_digest(Path(__file__)),'review_resource_plan_sha256':plan_before['sha256'],
            'check_count':len(checks),'checks':checks,'traits':results,'peak_observed_RSS_bytes':PEAK,
            'elapsed_seconds':time.monotonic()-BEGAN,'review_retained_SSD_bytes':retained,
            'original_inputs_modified_by_checker':False,'native_materializer_launched_by_checker':False,
            'h2_or_rg_estimated':False,'pair_tests_admitted':0,'complete_native_reproduction':False,
            'limits':['Exact compressed identity is distinguished from exact decompressed content; metadata-only attribution requires independently identical compressed payload and trailer.',
                      'Content equality covers source-specific raw materialization; no independent biological allele/effect/phenotype interpretation, harmonization, munging, h2/rg or 396/1200-family reproduction is certified.',
                      'Original helper Git identity and input before/after hashes support observed code/input binding. They cannot prove absence of transient edits or validate every source-specific scientific assumption.',
                      'Observed resource thresholds are sampled; no large reference dictionary or native analysis is loaded.']}
    save(HERE/'cancer_materialization_stream_review_v3.json',result)
    print(json.dumps({'status':result['status'],'checks':len(checks),'RSS':PEAK,'elapsed':result['elapsed_seconds']}),flush=True)
    if not all(x['pass'] for x in checks):raise SystemExit(1)


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--freeze',action='store_true');args=parser.parse_args()
    if args.freeze:freeze()
    else:verify()
