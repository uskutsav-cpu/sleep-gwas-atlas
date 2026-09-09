"""Read-only full traversal of archived LAVA shards; never repurpose HM3 data as dense GWAS."""
from __future__ import annotations
import collections
import gzip
import hashlib
import math
from pathlib import Path
from .artifacts import transaction
from .io import read_json, write_json, write_tsv, require, sha256, safe_write_path


def audit_shards(repo, root, name, *, trait=None):
    repo=Path(repo).resolve()
    manifest=repo/'results/track_b/lava_chromosome_inputs.provenance.json'
    c=read_json(manifest)
    require(c['schema_version']=='track-b-lava-chromosome-inputs.1','Unknown shard provenance schema')
    require(c['chromosomes']==list(range(1,23)),'Expected 22 autosomes')
    require(trait is None or trait in c['trait_order'],'Trait not in shard manifest')
    shards=[r for r in c['shards'] if trait is None or r['phenotype']==trait]
    traits={r['phenotype'] for r in shards}
    require(len(shards)==len(traits)*22,'Incomplete chromosome/trait family')
    require(len({(r['phenotype'],r['chromosome']) for r in shards})==len(shards),'Duplicate shard entry')
    out=[];inputs=[manifest]
    for record in shards:
        path=safe_write_path(repo,record['path'])
        require(path.is_file(),f'Missing actual shard bytes: {path}')
        inputs.append(path)
        require(path.stat().st_size==record['bytes'] and sha256(path)==record['sha256'],f'Shard checksum mismatch: {path}')
        n=missing=invalid=palindromic=0; min_n=math.inf; max_n=0.;z2sum=0.;valid=0
        seen=set();dups=0;decomp=hashlib.sha256()
        with gzip.open(path,'rb') as f:
            header=f.readline();decomp.update(header)
            require(header.rstrip(b'\r\n').decode().split('\t')==['SNP','A1','A2','Z','N'],'Unexpected LDSC shard columns')
            for raw in f:
                decomp.update(raw);n+=1
                values=raw.rstrip(b'\r\n').decode().split('\t')
                require(len(values)==5,'Malformed shard field count')
                snp,a1,a2,z,ns=values
                require(bool(snp),'Empty SNP ID')
                key=snp.lower();dups+=key in seen;seen.add(key)
                if any(v in {'','NA','nan','NaN'} for v in [a1,a2,z,ns]):missing+=1;continue
                try:z=float(z);ns=float(ns)
                except ValueError:invalid+=1;continue
                if not(math.isfinite(z) and math.isfinite(ns) and ns>0 and a1 in 'ACGT' and len(a1)==1 and a2 in 'ACGT' and len(a2)==1 and a1!=a2):
                    invalid+=1;continue
                palindromic+=a1+a2 in {'AT','TA','CG','GC'}
                valid+=1;z2sum+=z*z;min_n=min(min_n,ns);max_n=max(max_n,ns)
        require(n==record['data_row_count'],'Recorded and traversed row counts differ')
        require(dups==record['duplicate_snp_row_count']==0,'Duplicate SNP in shard')
        require(decomp.hexdigest()==record['decompressed_sha256'],'Decompressed byte identity mismatch')
        out.append({'trait':record['phenotype'],'chromosome':record['chromosome'],'rows':n,'valid_rows':valid,
             'missing_statistics_rows':missing,'invalid_statistics_rows':invalid,'palindromic_rows':palindromic,
             'mean_z_squared':z2sum/valid if valid else 'NA','minimum_n':min_n if valid else 'NA',
             'maximum_n':max_n if valid else 'NA','source_path':record['source_path'],
             'scope':'MUNGED_REFERENCE_INTERSECTION_NOT_DENSE_GWAS','sha256':record['sha256']})
    with transaction(root,name,stage='archived_shard_audit',inputs=inputs,parameters={'trait':trait},synthetic=False) as (work,meta):
        write_tsv(work/'shards.tsv',list(out[0]),out)
        status='PASS' if not any(r['invalid_statistics_rows'] for r in out) else 'FAILED_QC_NOT_CONSUMED'
        meta['scientific_status']=status
        write_json(work/'summary.json',{'status':status,'shards':len(out),'traits':len(traits),
             'rows':sum(r['rows'] for r in out),'valid_rows':sum(r['valid_rows'] for r in out),
             'missing_statistics_rows':sum(r['missing_statistics_rows'] for r in out),
             'invalid_statistics_rows':sum(r['invalid_statistics_rows'] for r in out),
             'raw_dense_gwas_verified':False,'new_local_rg_or_placo_inference':False,
             'limitation':'Verifies only the retained reference-intersected LDSC-format shards. Upstream whole-genome files and reference matrices were not reconstructed. Do not use these for dense PLACO or fine-mapping.'})
    return Path(root)/name
