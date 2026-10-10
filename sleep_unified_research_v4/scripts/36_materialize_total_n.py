#!/usr/bin/env python3
"""Result-free, streaming N-only derivatives for the frozen MS/melanoma diagnostic."""
import argparse
import datetime
import gzip
import hashlib
import json
import math
import os
import resource
import shutil
import sys
from pathlib import Path

MANIFEST_SHA='3ce7b50bb4e1289ef5a69d8e46dcef759963703ef7b317e156bdee0d4c40e5e1'
SSD=Path('/Volumes/Extreme SSD/sleep-unified-research-v1/research-completion-2026-10-09/sleep_unified_research_v4/sensitivities')
ROOT=Path(__file__).resolve().parents[2]
CONFIG={'ms':376169,'melanoma':290130}
MISSING={b'',b'NA',b'NaN',b'nan',b'NAN',b'#NA',b'.'}


def sha(p):
    h=hashlib.sha256()
    with Path(p).open('rb') as f:
        for b in iter(lambda:f.read(1<<20),b''):h.update(b)
    return h.hexdigest()


def guard():
    internal=shutil.disk_usage(ROOT).free
    output=shutil.disk_usage(SSD).free
    rss=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    if sys.platform!='darwin':rss*=1024
    if internal<3*(1<<30) or output<5*(1<<30) or rss>128*(1<<20):
        raise RuntimeError('DERIVATIVE_RESOURCE_GUARD: '+str((internal,output,rss)))
    return dict(internal_free_bytes=internal,SSD_free_bytes=output,self_peak_RSS_bytes=rss)


def transform(line,total):
    ending=b'\r\n' if line.endswith(b'\r\n') else b'\n' if line.endswith(b'\n') else b''
    body=line[:-len(ending)] if ending else line
    fields=body.split(b'\t')
    if len(fields)!=5:raise ValueError('expected exactly five fields')
    n=fields[-1]
    try:finite=math.isfinite(float(n))
    except ValueError:
        if n not in MISSING:raise ValueError('unexpected invalid N token')
        finite=False
    if finite and float(n)<=0:raise ValueError('nonpositive original N')
    prefix=b'\t'.join(fields[:4])+b'\t'
    out=prefix+(str(total).encode() if finite else n)+ending
    canonical=prefix+(b'<FINITE_N>' if finite else n)+ending
    return out,canonical,finite


class DigestSink:
    def __init__(self):self.h=hashlib.sha256();self.count=0
    def write(self,b):self.h.update(b);self.count+=len(b);return len(b)
    def flush(self):pass


def materialize(source,dest,total):
    before=sha(source)
    counts=dict(template_rows=0,finite_N_changed_rows=0,nonfinite_N_preserved_rows=0)
    invariant=hashlib.sha256()
    initial=guard()
    partial=dest.with_name(dest.name+'.partial')
    with source.open('rb') as raw_in,gzip.GzipFile(fileobj=raw_in,mode='rb') as src,partial.open('xb') as raw_out:
        with gzip.GzipFile(filename='',fileobj=raw_out,mode='wb',compresslevel=6,mtime=0) as out:
            header=src.readline()
            if header.rstrip(b'\r\n')!=b'SNP\tA1\tA2\tZ\tN':raise ValueError('source header mismatch')
            out.write(header);invariant.update(header)
            for line in src:
                transformed,canonical,finite=transform(line,total)
                out.write(transformed);invariant.update(canonical)
                counts['template_rows']+=1
                counts['finite_N_changed_rows' if finite else 'nonfinite_N_preserved_rows']+=1
                if counts['template_rows']%50000==0:guard()
    if counts['template_rows']!=1217311:raise ValueError('template cardinality changed')
    # Independently compare complete decompressed streams, including every row and literal missing token.
    verified=0
    out_invariant=hashlib.sha256()
    with gzip.open(source,'rb') as src,gzip.open(partial,'rb') as out:
        header=src.readline()
        if out.readline()!=header:raise ValueError('header mutated')
        out_invariant.update(header)
        for original in src:
            actual=out.readline()
            expected,canonical,finite=transform(original,total)
            if actual!=expected:raise ValueError('full-stream row differs at '+str(verified+2))
            _,actual_canonical,actual_finite=transform(actual,total)
            if actual_finite!=finite or actual_canonical!=canonical:raise ValueError('missingness/invariant differs')
            out_invariant.update(actual_canonical);verified+=1
            if verified%50000==0:guard()
        if out.read(1):raise ValueError('extra derivative rows')
    if invariant.hexdigest()!=out_invariant.hexdigest():raise ValueError('invariant hash mismatch')
    # A second deterministic compression goes directly to a hash sink, avoiding a second stored derivative.
    sink=DigestSink()
    with gzip.open(source,'rb') as src,gzip.GzipFile(filename='',fileobj=sink,mode='wb',compresslevel=6,mtime=0) as out:
        out.write(src.readline())
        for i,line in enumerate(src,1):
            out.write(transform(line,total)[0])
            if i%50000==0:guard()
    digest=sha(partial)
    if sink.h.hexdigest()!=digest or sink.count!=partial.stat().st_size:raise ValueError('gzip determinism failed')
    after=sha(source)
    if before!=after:raise ValueError('original source mutated')
    os.rename(partial,dest)
    return dict(source=str(source),source_sha256_before=before,source_sha256_after=after,
                derivative=str(dest),derivative_bytes=dest.stat().st_size,derivative_sha256=digest,
                new_finite_N=total,counts=counts,full_stream_verified_rows=verified,
                ordered_SNP_A1_A2_Z_missing_N_invariant_sha256=invariant.hexdigest(),
                original_and_derivative_invariant_match=True,literal_missingness_preserved=True,
                gzip=dict(filename='',mtime=0,compresslevel=6,independent_second_compression_matches=True),
                resource_initial=initial,resource_final=guard())


def main():
    parser=argparse.ArgumentParser();parser.add_argument('--manifest',type=Path,default=ROOT/'sleep_unified_research_v4/manifests/frozen_estimator_sensitivity_members_v4.json');args=parser.parse_args()
    if sha(args.manifest)!=MANIFEST_SHA:raise ValueError('frozen member manifest changed')
    m=json.loads(args.manifest.read_text())
    for path,expected in m['files'].items():
        if sha(path)!=expected:raise ValueError('frozen dependency changed: '+path)
    (SSD/'total_N_inputs').mkdir(parents=True,exist_ok=True);(SSD/'proofs').mkdir(exist_ok=True)
    receipt=SSD/'proofs/total_N_derivatives_v1.json'
    if receipt.exists():raise FileExistsError('derivative receipt already sealed')
    results=[]
    for trait,total in CONFIG.items():
        source=next(Path(p) for p in m['source_hashes'] if Path(p).name==trait+'.sumstats.gz')
        dest=SSD/'total_N_inputs'/(trait+'.sumstats.gz')
        if dest.exists() or dest.with_name(dest.name+'.partial').exists():raise FileExistsError('no overwrite of derivative or partial')
        if sha(source)!=m['source_hashes'][str(source)]:raise ValueError('original source SHA differs')
        results.append(materialize(source,dest,total))
    record=dict(schema='frozen_sensitivity_total_N_source_proof_v1',completed_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),
                status='RESULT_FREE_ONLY_FINITE_N_CHANGED',manifest_sha256=MANIFEST_SHA,worker_sha256=sha(Path(__file__)),
                worker_python=sys.version,results=results,estimator_calls=0,scientific_filter_changes=0,
                memory_limit_bytes=128*(1<<20),guard_internal_bytes=3*(1<<30),guard_output_bytes=5*(1<<30))
    with receipt.open('x') as f:json.dump(record,f,indent=2,allow_nan=False);f.write('\n')
    print(json.dumps(dict(receipt=str(receipt),receipt_sha256=sha(receipt),status=record['status'],results=results),indent=2))


if __name__=='__main__':main()
