"""Explicit read-only import of an already completed legacy PLACO result."""
from pathlib import Path
from .artifacts import transaction
from .io import read_json,read_tsv,require,check_hash,write_tsv,write_json


def import_legacy(manifest,root,relative):
    c=read_json(manifest)
    require(c.get("reviewed") is True,"Legacy source/column mapping must be reviewed")
    require(c.get("pair_id")=="insomnia__adhd","This import preserves the existing ADHD Pair B")
    require(c.get("status")=="COMPLETE_WITH_HITS" or c.get("status")=="COMPLETE_NO_HITS","Legacy terminal state unverified")
    require(c.get("method") in {"PLACO","PLACO_PLUS"},"Record actual legacy method, do not assume")
    require(c.get("input_fingerprint") and c.get("run_fingerprint"),"Missing legacy provenance")
    check_hash(c["path"],c["sha256"])
    columns=c["column_map"]
    fields=["SNP","CHR","BP","Z1","Z2","P_PLACO","status"]
    require(set(columns)==set(fields),"Map all required legacy fields, including failures")
    with transaction(root,relative,stage="import_legacy",inputs=[manifest,c["path"]],
                     parameters=c,synthetic=False) as (work,meta):
        n=write_tsv(work/"chunk.tsv",fields,
                    ({name:r[col] for name,col in columns.items()} for r in read_tsv(c["path"],columns.values())))
        require(n==c["expected_rows"],"Legacy output row count mismatch")
        meta["scientific_status"]="PASS"
        meta["legacy_method"]=c["method"]
        meta["reuse_is_not_replication"]=True
        write_json(work/"status.json",{"status":"PASS","imported_rows":n,"source_unmodified":True,
                   "note":"No new PLACO computation. New-family reporting does not create an independent replication."})
    return Path(root)/relative
