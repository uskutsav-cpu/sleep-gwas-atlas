"""A small synchronous DAG runner; Snakemake may schedule the same task graph.

All dependencies are file artifacts. A finished process is never silently
promoted to a causal result. The pair lock and runtime are immutable per run.
"""
from __future__ import annotations
import json
import os
import sys
from pathlib import Path
from .artifacts import verify_artifact, verify_current_artifact
from .io import (ContractError, check_hash, file_record, json_hash, read_json,
                 require, safe_id, safe_write_path, write_json)
from .pairs import verify_lock


def sorted_tasks(tasks):
    ids=[t["id"] for t in tasks]
    require(len(ids)==len(set(ids)),"Duplicate task IDs")
    for t in tasks:
        safe_id(t["id"])
        require(all(d in ids for d in t.get("depends_on",[])),"Unknown DAG dependency")
    pending={t["id"]:t for t in tasks}; result=[]; done=set()
    while pending:
        ready=sorted(k for k,t in pending.items() if set(t.get("depends_on",[]))<=done)
        require(ready,"Task graph contains a cycle")
        for k in ready:result.append(pending.pop(k));done.add(k)
    return result


def expand(value, artifacts):
    """Only complete @task/path strings resolve; no arbitrary evaluation."""
    if isinstance(value,dict):return {k:expand(v,artifacts) for k,v in value.items()}
    if isinstance(value,list):return [expand(v,artifacts) for v in value]
    if isinstance(value,str) and value.startswith("@"):
        parts=value[1:].split("/",1); name=parts[0]
        require(name in artifacts,f"Unresolved artifact reference: {value}")
        return str(safe_write_path(artifacts[name],parts[1])) if len(parts)==2 else str(artifacts[name])
    return value


def create_job(settings_file, *, method, job_id, inputs, outputs, root,
               reviewed=False, synthetic=False, minimum_free_bytes=2**30, threads=1):
    """Register every referenced input, including package wrapper source files."""
    base=Path(__file__).resolve().parents[1]
    script_dir=base/"scripts"
    if not (script_dir/"common.R").exists():script_dir=Path(__file__).resolve().parent/"resources"
    native_python={"munge","h2","rg","clump","ld","pleiofdr"}
    if method in native_python:
        script=script_dir/"native_cli.py"
        argv=[sys.executable,"{input:adapter}","{input:settings}"]
    else:
        require(method in {"placo","susie","coloc","lava","genomicsem","mr"},"Unsupported R adapter")
        script=script_dir/f"{method}.R"
        argv=["Rscript","{input:adapter}","{input:settings}"]
    from .bindings import discover_inputs
    records=discover_inputs(method,read_json(settings_file))
    records.update({name:file_record(path) for name,path in inputs.items()})
    records.update(adapter=file_record(script),settings=file_record(settings_file))
    if method not in native_python:records["common"]=file_record(script_dir/"common.R")
    job={"schema_version":1,"job_id":safe_id(job_id),"method":method,"reviewed":reviewed,
         "synthetic":synthetic,"argv":argv,"inputs":records,"outputs":outputs,
         "minimum_free_bytes":minimum_free_bytes,"threads":threads}
    output=safe_write_path(root,job_id+".job.json")
    if output.exists():require(read_json(output)==job,"Job configuration changed; create a new run")
    else:write_json(output,job)
    return output


def run_pipeline(runtime_path, *, only=None, allow_synthetic=False):
    runtime=read_json(runtime_path)
    require(runtime.get("schema_version")==1,"Unknown runtime schema")
    require(runtime.get("reviewed") is True,"Runtime bindings and task settings require review")
    synthetic=runtime.get("synthetic",False)
    require(not synthetic or allow_synthetic,"Synthetic pipeline requires explicit flag")
    lock=verify_lock(runtime["pair_lock"])
    require(lock["synthetic"]==synthetic,"Pair-lock synthetic mismatch")
    root=Path(runtime["output_root"]).resolve();root.mkdir(parents=True,exist_ok=True)
    tasks=sorted_tasks(runtime["tasks"])
    require(bool(tasks),"Empty workflow is not a completed project")
    signature={"runtime_sha256":file_record(runtime_path)["sha256"],"pair_lock":lock["lock_sha256"]}
    binding=root/"runtime.lock.json"
    if binding.exists():require(read_json(binding)==signature,"Runtime was changed within an existing run")
    else:write_json(binding,signature)
    artifacts={t["id"]:safe_write_path(root,t["id"]) for t in tasks}
    allowed_ids={r["pair_id"] for r in lock["pairs"]}
    if only is not None:require(only in artifacts,"Unknown task ID")
    for task in tasks:
        if only is not None and task["id"]!=only:continue
        for dep in task.get("depends_on",[]):
            receipt=verify_artifact(artifacts[dep])
            require(task["kind"]=="campaign_audit" or receipt.get("scientific_status","PASS") in {"PASS","NO_SIGNAL"},
                    f"Rejected/insufficient dependency: {dep}")
        # Validate dependency output bytes on every resume. Inputs too, where still mounted.
        if artifacts[task["id"]].exists():
            receipt=verify_current_artifact(artifacts[task["id"]])
            for rec in receipt.get("inputs",[]):check_hash(rec["path"],rec["sha256"])
            continue
        if task.get("pair_id"):
            require(task["pair_id"] in allowed_ids,"Task pair is outside frozen family")
            if task["kind"]=="native" and task.get("method")=="placo" and task["pair_id"]=="insomnia__adhd":
                raise ContractError("Legacy Pair B is protected: import its verified canonical outputs instead of rerunning PLACO")
        p=expand(task.get("parameters",{}),artifacts)
        kind=task["kind"]
        if kind=="normalize":
            from .gwas import normalize
            normalize(p["source"],root,task["id"],synthetic=synthetic)
        elif kind=="join":
            from .gwas import join_pair
            join_pair(p.pop("left"),p.pop("right"),root,task["id"],synthetic=synthetic,**p)
        elif kind=="index":
            from .ld import index_pair
            index_pair(p["pair_dir"],root,task["id"])
        elif kind=="split":
            from .pleiotropy import split_pair
            split_pair(p["pair_dir"],root,task["id"],chunk_rows=p.get("chunk_rows",20000))
        elif kind=="locus":
            from .ld import prepare_locus
            prepare_locus(p.pop("pair_index"),p.pop("ld_manifest"),p.pop("locus"),root,task["id"],synthetic=synthetic,**p)
        elif kind=="collate_placo":
            from .pleiotropy import collate_placo
            collate_placo(p.pop("manifest"),root,task["id"],synthetic=synthetic,**p)
        elif kind=="mr_candidates":
            from .mr_workflow import candidates
            candidates(p.pop("exposure"),root,task["id"],**p)
        elif kind=="mr_harmonize":
            from .mr_workflow import harmonize
            harmonize(p.pop("candidates_dir"),p.pop("clump_dir"),p.pop("outcome_dir"),root,task["id"],**p)
        elif kind=="mr_family":
            from .mr_workflow import collate
            collate(p["manifest"],root,task["id"])
        elif kind=="replication_family":
            from .replication import evaluate
            evaluate(p["manifest"],root,task["id"])
        elif kind=="qtl_index":
            from .qtl import index
            index(p["source"],root,task["id"],synthetic=synthetic)
        elif kind=="qtl_locus":
            from .qtl import locus
            locus(p.pop("gwas_dir"),p.pop("qtl_dir"),p.pop("feature"),p.pop("region"),root,task["id"],**p)
        elif kind=="annotation_scores":
            from .annotations import score
            score(p["variants"],p["annotations"],root,task["id"],synthetic=synthetic)
        elif kind=="campaign_audit":
            from .campaign import audit
            audit(p["manifest"],root,task["id"])
        elif kind in {"joint_placo_family","local_rg_family","pathway_family","block_comparison","robustness_family","cell_family"}:
            if kind=="joint_placo_family":
                from .family24 import joint_placo as fn
            elif kind=="local_rg_family":
                from .family24 import local_family as fn
            elif kind=="pathway_family":
                from .pathway24 import run as fn
            elif kind=="block_comparison":
                from .compare24 import compare as fn
            elif kind=="cell_family":
                from .cell24 import run as fn
            else:
                from .robustness24 import summarize as fn
            require(read_json(p["manifest"]).get("synthetic",False)==synthetic,"Synthetic manifest mismatch")
            fn(p["manifest"],root,task["id"])
        elif kind=="native":
            from .executor import run_job
            jobdir=root/"job_specs";jobdir.mkdir(exist_ok=True)
            settings=jobdir/(task["id"]+".settings.json")
            if settings.exists():require(read_json(settings)==p,"Native settings changed")
            else:write_json(settings,p)
            jobpath=create_job(settings,method=task["method"],job_id=task["id"],
                    inputs=expand(task.get("input_files",{}),artifacts),outputs=task["outputs"],root=jobdir,
                    reviewed=True,synthetic=synthetic,minimum_free_bytes=task.get("minimum_free_bytes",2**30),
                    threads=task.get("threads",1))
            run_job(jobpath,root,allow_synthetic=synthetic)
        else:raise ContractError(f"Unknown task kind: {kind}")
    return {k:str(v) for k,v in artifacts.items()}
