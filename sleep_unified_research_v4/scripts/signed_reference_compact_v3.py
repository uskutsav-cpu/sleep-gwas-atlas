"""Lossless fixed-schema private retention; no SNP removal or float rounding."""
import hashlib
import json
import os
from pathlib import Path
import struct
from signed_reference_source_common_v3 import regular_sha, record_output
from signed_reference_decode_v2 import regular

MASTER = struct.Struct('<BQQIQb4sB')
QC = struct.Struct('<IHHHIB')
HEADER = struct.Struct('<8sQ')
MASTER_MAGIC = b'SRMSTR02'
QC_MAGIC = b'SRGENQ02'
SENTINEL = (1 << 64) - 1
STATIC_REASONS = [None, 'MHC_INCLUSIVE', 'MISSING_COORDINATE_MAP_ID',
                  'REFERENCE_COORDINATE_MAP_CONFLICT', 'MISSING_SOURCE_ID',
                  'SOURCE_GRCh37_COORDINATE_CONFLICT', 'INVALID_TARGET_ALLELES',
                  'PALINDROMIC', 'ALLELE_CONFLICT_NO_COMPLEMENT_RESCUE']
QC_REASONS = [None, 'ALL_MISSING', 'FEWER_THAN_TWO_OBSERVED',
              'EXACT_ZERO_VARIANCE', 'MISSING_GT25_OF503']
CAPACITY = dict(master_rows=1290028, master_record_bytes=MASTER.size,
    genotype_record_bytes=QC.size, header_bytes=HEADER.size,
    SNP_UTF8_pool_cap_bytes=64*(1<<20), coordinate_SQLite_page_bytes=4096,
    coordinate_SQLite_max_pages=32768, coordinate_SQLite_cap_bytes=128*(1<<20),
    control_component_cap_bytes=64*(1<<20), companion_component_cap_bytes=32*(1<<20),
    other_operational_metadata_cap_bytes=32*(1<<20), failure_metadata_reserve_bytes=16*(1<<20))
CAPACITY['maximum_master_and_genotype_fixed_bytes'] = 1290028*(MASTER.size+QC.size)+2*HEADER.size
CAPACITY['maximum_planned_QC_bytes'] = CAPACITY['maximum_master_and_genotype_fixed_bytes'] + sum(
    CAPACITY[k] for k in ['SNP_UTF8_pool_cap_bytes', 'coordinate_SQLite_cap_bytes',
    'control_component_cap_bytes', 'companion_component_cap_bytes',
    'other_operational_metadata_cap_bytes', 'failure_metadata_reserve_bytes'])
CAPACITY['unchanged_QC_cap_bytes'] = 512*(1<<20)
CAPACITY['structural_margin_bytes'] = CAPACITY['unchanged_QC_cap_bytes']-CAPACITY['maximum_planned_QC_bytes']
SCHEMA = dict(version=2, byte_order='little', master_format=MASTER.format,
    master_fields=['CHR_u8','BP_u64','SNP_pool_offset_u64','SNP_UTF8_length_u32',
                   'source_row0_u64_or_sentinel','sign_i8','target_A1_A2_source_A1_A2_ASCII4','static_reason_u8'],
    master_rank0='exact sequential original master row index', master_rank1='master_rank0+1',
    source_member='1000G.EUR.QC.{CHR}.bed', BED_byte_offset='3+126*source_row0',
    genotype_format=QC.format,
    genotype_fields=['original_master_rank0_u32','observed_count_u16','sum_dose_u16',
                     'sum_dose2_u16','integer_variance_numerator_u32','QC_reason_u8'],
    missing_count='503-observed_count', J_index0='sequential eligible QC row count; never master rank',
    binary64_frequency='sum_dose/(2*observed_count) if observed_count else None; original expression',
    binary64_MAF='min(f,1-f) if f is not None else None; original expression',
    static_reasons=STATIC_REASONS, genotype_reasons=QC_REASONS,
    schema_reconstruction='exact original static/J/exclusion dictionaries and canonical JSON bytes',
    capacity=CAPACITY)
SCHEMA_SHA256 = hashlib.sha256(json.dumps(SCHEMA,sort_keys=True,separators=(',',':')).encode()).hexdigest()


def canonical_line(row):
    return (json.dumps(row,sort_keys=True,separators=(',',':'),allow_nan=False)+'\n').encode()


class IntendedWriter:
    def __init__(self,path,cap):
        self.path=Path(path); self.cap=cap; self.total=0; self.sha=hashlib.sha256()
        if any(p.is_symlink() for p in self.path.parents):raise RuntimeError('COMPACT_PARENT_SYMLINK')
        self.file=self.path.open('xb')
    def write(self,payload):
        if self.total+len(payload)>self.cap:raise RuntimeError('COMPACT_COMPLETE_RETENTION_CAP_NO_TRUNCATION')
        self.sha.update(payload);self.file.write(payload);self.total+=len(payload)
    def close(self):
        self.file.flush();os.fsync(self.file.fileno());self.file.close()
        fd=os.open(self.path.parent,os.O_RDONLY)
        try:os.fsync(fd)
        finally:os.close(fd)
        expected=self.sha.hexdigest()
        if regular_sha(self.path)!=expected:raise RuntimeError('COMPACT_INTENDED_STREAM_CHANGED')
        record_output(self.path,expected);return expected
    def abort(self):
        if not self.file.closed:self.file.flush();self.file.close()


class StaticWriter:
    def __init__(self,root,expected_rows):
        root=Path(root);self.expected=expected_rows;self.count=0;self.selected=0
        self.master=IntendedWriter(root/'master_static.bin',HEADER.size+expected_rows*MASTER.size)
        self.pool=IntendedWriter(root/'master_SNP_utf8.bin',CAPACITY['SNP_UTF8_pool_cap_bytes'])
        self.master.write(HEADER.pack(MASTER_MAGIC,expected_rows))
        self.logical={k:hashlib.sha256() for k in ['static_selection','static_exclusions']}
    def append(self,row):
        rank=self.count
        if row['master_rank0']!=rank or row['master_rank1']!=rank+1 or rank>=self.expected:raise RuntimeError('COMPACT_ORIGINAL_MASTER_RANK')
        text=row['SNP'].encode('utf8');offset=self.pool.total;reason=row.get('exclusion')
        code=STATIC_REASONS.index(reason)
        if code:
            source=SENTINEL;sign=0;alleles=b'\0'*4;role='static_exclusions'
        else:
            source=row['source_row0'];sign=row['sign'];alleles=''.join(row[k] for k in ['target_A1','target_A2','source_A1','source_A2']).encode('ascii');role='static_selection'
            if row['BED_byte_offset']!=3+126*source or row['source_member']!=f"1000G.EUR.QC.{row['CHR']}.bed" or sign not in [-1,1] or len(alleles)!=4:raise RuntimeError('COMPACT_SELECTED_SOURCE_FIELDS')
        if not 1<=row['CHR']<=22 or not 0<row['BP']<1<<64 or not text or len(text)>=1<<32 or not 0<=source<=SENTINEL:raise RuntimeError('COMPACT_LOSSLESS_FIELD_RANGE_ENTIRE_STAGE_FAIL')
        payload=MASTER.pack(row['CHR'],row['BP'],offset,len(text),source,sign,alleles,code)
        self.pool.write(text);self.master.write(payload);self.logical[role].update(canonical_line(row))
        self.count+=1;self.selected+=code==0
    def close(self):
        if self.count!=self.expected:raise RuntimeError('COMPACT_ALL_MASTER_ROWS_REQUIRED')
        return dict(schema_sha256=SCHEMA_SHA256,master_sha256=self.master.close(),SNP_pool_sha256=self.pool.close(),
                    master_rows=self.count,static_rows=self.selected,logical_sha256={k:v.hexdigest() for k,v in self.logical.items()})
    def abort(self):self.master.abort();self.pool.abort()


def static_rows(root,binding,guard=lambda:None,selected=True):
    root=Path(root);master=root/'master_static.bin';pool=root/'master_SNP_utf8.bin'
    if binding.get('schema_sha256')!=SCHEMA_SHA256:raise RuntimeError('COMPACT_EXACT_SCHEMA_REQUIRED')
    expected={str(master):binding['master_sha256'],str(pool):binding['SNP_pool_sha256']}
    for p,h in expected.items():
        if regular_sha(p)!=h:raise RuntimeError('COMPACT_STATIC_CHANGED_BEFORE_PARSE')
    logical={k:hashlib.sha256() for k in ['static_selection','static_exclusions']};count=0;sel=0;offset_want=0
    with regular(master).open('rb') as f,regular(pool).open('rb') as ids:
        if f.read(HEADER.size)!=HEADER.pack(MASTER_MAGIC,binding['master_rows']):raise RuntimeError('COMPACT_MASTER_HEADER')
        for rank in range(binding['master_rows']):
            guard();b=f.read(MASTER.size)
            if len(b)!=MASTER.size:raise RuntimeError('COMPACT_MASTER_SHORT_RECORD')
            chrom,bp,offset,length,source,sign,alleles,code=MASTER.unpack(b)
            if not 1<=chrom<=22 or not bp or offset!=offset_want or not length or code>=len(STATIC_REASONS):raise RuntimeError('COMPACT_MASTER_SCHEMA_VALUES')
            text=ids.read(length)
            if len(text)!=length:raise RuntimeError('COMPACT_ID_POOL_SHORT')
            offset_want+=length;row=dict(SNP=text.decode('utf8'),CHR=chrom,BP=bp,master_rank0=rank,master_rank1=rank+1)
            if code:
                if source!=SENTINEL or sign or alleles!=b'\0'*4:raise RuntimeError('COMPACT_EXCLUDED_FIELD_SENTINELS')
                row['exclusion']=STATIC_REASONS[code];role='static_exclusions'
            else:
                if source==SENTINEL or sign not in [-1,1] or any(x not in b'ACGT' for x in alleles):raise RuntimeError('COMPACT_SELECTED_FIELD_SCHEMA')
                row.update(target_A1=chr(alleles[0]),target_A2=chr(alleles[1]),source_A1=chr(alleles[2]),source_A2=chr(alleles[3]),
                    sign=sign,source_row0=source,BED_byte_offset=3+126*source,source_member=f'1000G.EUR.QC.{chrom}.bed');role='static_selection';sel+=1
            logical[role].update(canonical_line(row));count+=1
            if (selected and not code) or (not selected and code):yield row
        if f.read(1) or ids.read(1):raise RuntimeError('COMPACT_MASTER_OR_POOL_TRAILING')
    if count!=binding['master_rows'] or sel!=binding['static_rows'] or {k:v.hexdigest() for k,v in logical.items()}!=binding['logical_sha256']:raise RuntimeError('COMPACT_STATIC_LOGICAL_STREAM_DIFFERS')
    for p,h in expected.items():
        if regular_sha(p)!=h:raise RuntimeError('COMPACT_STATIC_CHANGED_AFTER_PARSE')


class GenotypeWriter:
    def __init__(self,root,static_binding):
        self.expected=static_binding['static_rows'];self.count=0;self.J=0
        self.writer=IntendedWriter(Path(root)/'genotype_QC.bin',HEADER.size+self.expected*QC.size)
        self.writer.write(HEADER.pack(QC_MAGIC,self.expected));self.previous=-1
        self.logical={k:hashlib.sha256() for k in ['J_ordered','genotype_exclusions']}
    def append(self,item):
        if item['master_rank0']<=self.previous or self.count>=self.expected:raise RuntimeError('COMPACT_QC_MASTER_ORDER')
        q=qc_values(item['observed_count'],item['sum_dose'],item['sum_dose2'],item['integer_variance_numerator'])
        if any(item[k]!=v for k,v in q.items()):raise RuntimeError('COMPACT_EXACT_SUFFICIENT_STATISTICS_OR_FLOAT_RECONSTRUCTION')
        code=QC_REASONS.index(item['exclusion']);role='genotype_exclusions' if code else 'J_ordered'
        if not code:
            if item['J_index0']!=self.J:raise RuntimeError('COMPACT_EXACT_J_INDEX')
            self.J+=1
        self.writer.write(QC.pack(item['master_rank0'],item['observed_count'],item['sum_dose'],item['sum_dose2'],item['integer_variance_numerator'],code))
        self.logical[role].update(canonical_line(item));self.previous=item['master_rank0'];self.count+=1
    def close(self):
        if self.count!=self.expected:raise RuntimeError('COMPACT_ALL_STATIC_QC_ROWS_REQUIRED')
        return dict(schema_sha256=SCHEMA_SHA256,genotype_QC_sha256=self.writer.close(),QC_rows=self.count,J_rows=self.J,
                    logical_sha256={k:v.hexdigest() for k,v in self.logical.items()})
    def abort(self):self.writer.abort()


def qc_values(n,total,total2,numerator):
    if not 0<=n<=503 or not 0<=total<=2*n or not 0<=total2<=4*n or numerator!=n*total2-total*total or numerator<0:raise RuntimeError('COMPACT_EXACT_INTEGER_QC_SCHEMA')
    miss=503-n
    reason='ALL_MISSING' if n==0 else 'FEWER_THAN_TWO_OBSERVED' if n<2 else 'EXACT_ZERO_VARIANCE' if numerator<=0 else 'MISSING_GT25_OF503' if miss>25 else None
    f=total/(2*n) if n else None
    return dict(observed_count=n,missing_count=miss,sum_dose=total,sum_dose2=total2,integer_variance_numerator=numerator,
                A1_frequency=f,MAF=min(f,1-f) if f is not None else None,exclusion=reason)


def genotype_rows(root,static_binding,binding,guard=lambda:None,eligible=True):
    path=Path(root)/'genotype_QC.bin';expected=binding['genotype_QC_sha256']
    if binding.get('schema_sha256')!=SCHEMA_SHA256 or regular_sha(path)!=expected:raise RuntimeError('COMPACT_QC_CHANGED_BEFORE_PARSE_OR_SCHEMA')
    logical={k:hashlib.sha256() for k in ['J_ordered','genotype_exclusions']};count=0;index=0
    with regular(path).open('rb') as f:
        if f.read(HEADER.size)!=HEADER.pack(QC_MAGIC,static_binding['static_rows']):raise RuntimeError('COMPACT_QC_HEADER')
        for row in static_rows(root,static_binding,guard):
            guard();b=f.read(QC.size)
            if len(b)!=QC.size:raise RuntimeError('COMPACT_QC_SHORT_RECORD')
            rank,n,total,total2,num,code=QC.unpack(b)
            if rank!=row['master_rank0'] or code>=len(QC_REASONS):raise RuntimeError('COMPACT_QC_EXACT_MASTER_JOIN')
            q=qc_values(n,total,total2,num)
            if q['exclusion']!=QC_REASONS[code]:raise RuntimeError('COMPACT_QC_REASON_MISMATCH')
            item={**row,**q};role='genotype_exclusions' if code else 'J_ordered'
            if not code:item['J_index0']=index;index+=1
            logical[role].update(canonical_line(item));count+=1
            if (eligible and not code) or (not eligible and code):yield item
        if f.read(1):raise RuntimeError('COMPACT_QC_TRAILING')
    if count!=binding['QC_rows'] or index!=binding['J_rows'] or {k:v.hexdigest() for k,v in logical.items()}!=binding['logical_sha256']:raise RuntimeError('COMPACT_QC_LOGICAL_STREAM_DIFFERS')
    if regular_sha(path)!=expected:raise RuntimeError('COMPACT_QC_CHANGED_AFTER_PARSE')
