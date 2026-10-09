#!/usr/bin/env python3
"""Combine immutable discovery inventory and separate verification ledgers."""
from collections import Counter
import csv
from datetime import datetime,timezone
import hashlib
import json
from pathlib import Path
import shutil

PACKAGE=Path(__file__).resolve().parents[1]
ROOT=PACKAGE.parent
FIELDS=['path','bytes','expected_sha256','actual_sha256','source_release','genome_build','ancestry','phenotype_identity','effect_coding','permitted_use','required_downstream_tasks','verification_status','asset_kind','expected_hash_evidence']

def read(p):
    with p.open(newline='') as f:return list(csv.DictReader(f,delimiter='\t'))
def write(p,rows):
    with p.open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=FIELDS,delimiter='\t',extrasaction='ignore');w.writeheader();w.writerows(rows)
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()

def main():
    initial=PACKAGE/'manifests/ssd_inventory_initial_v1.tsv'
    if not initial.exists():shutil.copyfile(PACKAGE/'SSD_INPUT_MANIFEST.tsv',initial)
    inventory={r['path']:r for r in read(initial)}
    core={r['trait_id']:r for r in read(ROOT/'config/analysis_panel.tsv')}
    schemas={r['trait_id']:r for r in read(ROOT/'config/gwas_schemas.tsv')}
    ext={r['extension_trait_id']:r for r in read(ROOT/'discovery_extension/config/candidate_traits.tsv')}
    rep={r['replication_source_id']:r for r in read(ROOT/'discovery_extension/config/replication_manifest.tsv') if r['source_curation_status']=='COMPLETE_BEFORE_RESULTS'}
    def add(path,values,trait=None,kind=''):
        base=inventory.get(path,{'path':path})
        base.update(values);base['asset_kind']=kind or base.get('asset_kind','source_asset')
        base.setdefault('permitted_use','RESEARCH_ONLY; SOURCE_TERMS_AND_ATTRIBUTION_REQUIRED; NO_RAW_REDISTRIBUTION')
        base.setdefault('required_downstream_tasks','SOURCE_QC_AND_NATIVE_REPRODUCTION_PENDING')
        if trait in core:
            c=core[trait];base.update(source_release=c['dataset_version'],genome_build=c['build'],ancestry=c['ancestry'],phenotype_identity=trait,effect_coding=schemas[trait]['effect_convention'])
            if kind in ['core_munged','core_harmonized']:base['effect_coding']='Z_OR_BETA_FOR_REPORTED_A1; PER_SOURCE_ORIENTATION_CHAIN_PENDING_REPLAY'
        for key in ['source_release','genome_build','ancestry','phenotype_identity','effect_coding']:
            base.setdefault(key,'UNKNOWN_SOURCE_METADATA_REVIEW_REQUIRED')
        inventory[path]=base
    checks=PACKAGE/'tables/native_input_hash_checks.tsv'
    for r in read(checks):
        metadata={k:r.get(k,'') for k in ['bytes','expected_sha256','actual_sha256','expected_hash_evidence']}
        metadata['verification_status']=r['status']
        add(r['path'],metadata,r['trait_id'],r['kind'])
    recovered=PACKAGE/'tables/archived_extension_recovery.tsv'
    for r in read(recovered):
        relative=r['archive_member'].split('sleep-gwas-atlas/',1)[1]
        is_munged=relative.endswith('.sumstats.gz')
        trait=Path(relative).name.replace('.sumstats.gz','').replace('.munge.log','')
        v={k:r.get(k,'') for k in ['bytes','expected_sha256','actual_sha256','expected_hash_evidence']}
        v['verification_status']=r['status']
        if trait in ext:
            c=ext[trait];v.update(source_release=c['source']+';'+c['study_accession'],genome_build=c['build'],ancestry=c['ancestry'],phenotype_identity=c['phenotype_name'],effect_coding='SIGNED_Z_FOR_A1; N_EFFECTIVE_FOR_BINARY' if is_munged else 'NOT_APPLICABLE_LOG_OR_QC')
        elif trait in rep:
            c=rep[trait];v.update(source_release=c['replication_study_accession'],genome_build=c['build'],ancestry=c['ancestry'],phenotype_identity=c['replication_phenotype_definition'],effect_coding='SIGNED_Z_FOR_A1; N_EFFECTIVE' if is_munged else 'NOT_APPLICABLE_LOG_OR_QC')
        elif '/logs/' in relative or '/qc/' in relative:
            v.update(source_release='ORIGINAL_ARCHIVED_RUN',genome_build='NOT_APPLICABLE_LOG_OR_QC',ancestry='EUR',phenotype_identity=Path(relative).stem,effect_coding='NOT_APPLICABLE_LOG_OR_QC')
        else:
            v.update(source_release='PINNED_PANUKBB_HAPMAP3_REFERENCE',genome_build='GRCh37',ancestry='EUR',phenotype_identity='HAPMAP3_VARIANT_ALLELE_MAP',effect_coding='REFERENCE_ALLELES; NOT_A_GWAS_EFFECT')
        add(r['path'],v,kind='recovered_munged' if is_munged else 'recovered_reference_log_or_QC')
    acquisition=PACKAGE/'logs/mvp_insomnia_acquisition_receipt_v1.json'
    c=json.loads(acquisition.read_text())
    add(c['path'],dict(bytes=c['actual_bytes'],expected_sha256='',actual_sha256=c['actual_sha256'],
        expected_hash_evidence='UPSTREAM_MD5='+c['expected_md5']+'; logs/mvp_insomnia_acquisition_receipt_v1.json',
        verification_status='UPSTREAM_MD5_AND_BYTES_MATCH; CURRENT_SHA256_PINNED',
        source_release='MVP_GCST90475826_2025_03_19',genome_build='GRCh38',ancestry='MVP_EUR',
        phenotype_identity='CLINICAL_INSOMNIA_PHECODE_327.4',effect_coding='LOG_OR_FOR_EFFECT_ALLELE; CI_DERIVED_SE_REQUIRES_AUDIT',
        required_downstream_tasks='SCHEMA_N_ALLELE_BUILD_QC_H2_POWER_AND_PHENOTYPE_TRANSPORT_GATES'),kind='new_conditionally_admitted_raw')
    rows=[inventory[p] for p in sorted(inventory)]
    write(PACKAGE/'SSD_INPUT_MANIFEST.tsv',rows)
    receipt={'completed_utc':datetime.now(timezone.utc).isoformat(),'initial_inventory_rows':len(read(initial)),
       'unified_rows':len(rows),'hash_status_counts':dict(Counter(r['verification_status'] for r in rows)),
       'source_ledgers_sha256':{str(p.relative_to(PACKAGE)):sha(p) for p in [initial,checks,recovered,acquisition]},
       'manifest_sha256':sha(PACKAGE/'SSD_INPUT_MANIFEST.tsv'),
       'unhashed_entries_not_assumed_to_match':True,'original_SSD_files_changed':False}
    (PACKAGE/'logs/unified_inventory_receipt_v1.json').write_text(json.dumps(receipt,indent=2)+'\n')
    print(json.dumps(receipt,indent=2))
if __name__=='__main__':main()
