#!/usr/bin/env python3
"""Hash recovered core sources and compare LD files to the pinned official tar.

No source file is changed. Expensive checks are sequential and resumable only
when path, size and nanosecond mtime still match an earlier check.
"""
import argparse
import csv
import hashlib
import json
from pathlib import Path
import tarfile
import urllib.request
from datetime import datetime, timezone

from importlib.util import spec_from_file_location, module_from_spec
spec = spec_from_file_location('recovery', Path(__file__).with_name('00_recover_inventory.py'))
recovery = module_from_spec(spec); spec.loader.exec_module(recovery)

REFERENCE_URL = 'https://zenodo.org/records/8182036/files/eur_w_ld_chr.tar.gz?download=1'
REFERENCE_SHA256 = '9537f00eb0d163a935aaa2cf04b358b7cf21852279b9c7925802526f6060b069'
REF = recovery.FAILED / 'work/track_b_completion/local_dependency_copies/ref'

def qc_fields(path):
    output = {}
    for line in path.read_text().splitlines():
        parts = line.split('\t')
        if len(parts) == 2: output[parts[0]] = parts[1]
    return output

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--repo', type=Path, default=Path(__file__).resolve().parents[2])
    parser.add_argument('--output', type=Path, default=Path(__file__).resolve().parents[1])
    parser.add_argument('--reference-tar', type=Path)
    args = parser.parse_args(); repo, out = args.repo, args.output
    ledger = out / 'tables/native_input_hash_checks.tsv'
    fields = ['kind','trait_id','path','bytes','mtime_ns','expected_sha256','actual_sha256','status','expected_hash_evidence']
    prior = recovery.read_tsv(ledger) if ledger.exists() else []
    cache = {r['path']:r for r in prior}
    rows = []
    def check(path, expected, kind, trait='', evidence=''):
        if not path.is_file():
            r = dict(kind=kind,trait_id=trait,path=str(path),bytes='',mtime_ns='',expected_sha256=expected,
                     actual_sha256='',status='MISSING',expected_hash_evidence=evidence)
        else:
            st = path.stat(); old = cache.get(str(path),{})
            if old.get('bytes') == str(st.st_size) and old.get('mtime_ns') == str(st.st_mtime_ns):
                actual = old['actual_sha256']
            else: actual = recovery.sha256(path)
            status = ('MATCH_EXPECTED_SHA256' if actual == expected else 'HASH_MISMATCH') if expected else 'HASH_RECORDED_NO_HISTORICAL_EXPECTED_HASH'
            r = dict(kind=kind,trait_id=trait,path=str(path),bytes=st.st_size,mtime_ns=st.st_mtime_ns,
                     expected_sha256=expected,actual_sha256=actual,status=status,expected_hash_evidence=evidence)
        rows.append(r)
        recovery.write_tsv(ledger, rows, fields)
        print(f"{kind} {trait or path.name}: {r['status']}", flush=True)
        return r

    # Verify the reference by decompressed tar member bytes, not filename alone.
    archive = args.reference_tar or repo.parent / 'reference_check/eur_w_ld_chr.tar.gz'
    archive.parent.mkdir(parents=True, exist_ok=True)
    if not archive.exists():
        with urllib.request.urlopen(REFERENCE_URL, timeout=90) as src, archive.open('xb') as dst:
            while chunk := src.read(1024*1024): dst.write(chunk)
    reference = check(archive, REFERENCE_SHA256, 'official_reference_archive', evidence=REFERENCE_URL)
    if reference['status'] != 'MATCH_EXPECTED_SHA256': raise SystemExit('REFERENCE_TAR_HASH_FAILED')
    with tarfile.open(archive, 'r:gz') as tar:
        for member in tar.getmembers():
            if not member.isfile(): continue
            rel = Path(member.name)
            if rel.parts[0] != 'eur_w_ld_chr': continue
            # Extra archived *_old files are not consumed or promoted.
            with tar.extractfile(member) as f:
                expected = hashlib.sha256(f.read()).hexdigest()
            check(REF / rel, expected, 'official_ld_reference_member', evidence=REFERENCE_SHA256+':'+member.name)
    check(REF/'hm3_grch37_variant_map.tsv.gz','6775a7a0d3ca90dc74e472180b1d77103bc238129c4f969358f46307e5c306b4','variant_map',evidence='docs/reference_panel_provenance.md')
    check(REF/'hg38ToHg19.over.chain.gz','14a712e8e147d9fc8e9d87d51977b46f6f8ddb93efbe5d0843d86b6205f587b1','liftover_chain',evidence='config/liftover_plans.tsv')
    panel = recovery.read_tsv(repo / 'config/analysis_panel.tsv')
    for trait in panel:
        tid = trait['trait_id']
        qc = recovery.DEPS/'data/harmonized'/f'{tid}.qc.txt'
        metadata = qc_fields(qc) if qc.exists() else {}
        expected = metadata.get('prefilter_source_sha256', metadata.get('infile_sha256',''))
        raw = recovery.DEPS/'data/raw'/trait['raw_file']
        if not raw.exists() and (recovery.RAW/'raw'/trait['raw_file']).exists():
            raw = recovery.RAW/'raw'/trait['raw_file']
        check(raw,expected,'core_raw',tid,str(qc))
        if metadata.get('prefilter_source_sha256'):
            check(recovery.DEPS/'data/harmonized/.prefilter'/Path(metadata['infile']).name,
                  metadata['infile_sha256'],'core_prefiltered',tid,str(qc))
        check(recovery.DEPS/'data/munged'/f'{tid}.sumstats.gz','','core_munged',tid)
        check(recovery.DEPS/'data/harmonized'/f'{tid}.harmonized.tsv.gz','','core_harmonized',tid)
    # Exact archives cover source release identity separately from derived raw bytes.
    seen = set()
    for source in recovery.read_tsv(repo/'config/public_gwas_sources.tsv'):
        path = recovery.DEPS/'data/raw/.archives'/source['archive_name']
        if str(path) in seen: continue
        seen.add(str(path))
        check(path,source['archive_sha256'],'core_source_archive',source['trait_ids'],'config/public_gwas_sources.tsv')

    manifest = recovery.read_tsv(out/'SSD_INPUT_MANIFEST.tsv')
    checked = {r['path']:r for r in rows}
    for r in manifest:
        c = checked.get(r['path'])
        if c:
            r.update(expected_sha256=c['expected_sha256'],actual_sha256=c['actual_sha256'],verification_status=c['status'])
            if c['kind'].startswith('core_'):
                trait = next((t for t in panel if t['trait_id']==c['trait_id']),{})
                r.update(phenotype_identity=c['trait_id'], source_release=trait.get('dataset_version','UNKNOWN'),
                         genome_build=trait.get('build','UNKNOWN'),ancestry=trait.get('ancestry','UNKNOWN'))
    recovery.write_tsv(out/'SSD_INPUT_MANIFEST.tsv',manifest,list(manifest[0]))
    summary = {'completed_utc':datetime.now(timezone.utc).isoformat(),'checks':len(rows),
               'matched':sum(r['status']=='MATCH_EXPECTED_SHA256' for r in rows),
               'mismatches':[r for r in rows if r['status']=='HASH_MISMATCH'],
               'missing':[r for r in rows if r['status']=='MISSING'],
               'core_munged_present':sum(r['kind']=='core_munged' and r['status']!='MISSING' for r in rows),
               'source_input_files_modified':False}
    (out/'logs/native_input_verification_receipt_v1.json').write_text(json.dumps(summary,indent=2)+'\n')
    print(json.dumps(summary,indent=2))

if __name__ == '__main__': main()
