"""Inspect a checkout without changing files or claiming to inspect a missing repo."""
from __future__ import annotations
import importlib.util
import platform
import shutil
import subprocess
from pathlib import Path
from .io import file_record,write_json,require


def audit_repo(repo,out):
    root=Path(repo).resolve()
    require(root.is_dir(),f"Checkout directory is missing: {root}")
    known=["config/analysis_panel.tsv","config/analysis_panel.lock.json", "results/tables/rg_matrix.tsv",
           "results/tables/h2_summary.tsv","results/track_b/pleiotropy/results/placo/A.full.tsv.gz",
           "results/track_b/pleiotropy/results/placo/B.full.tsv.gz",
           "results/track_b/pleiotropy/results/placo/CONTROL.full.tsv.gz"]
    result={"repo":str(root),"read_only":True,"files":{},"git":{},"environment":{}}
    for name in known:
        p=root/name
        result["files"][name]=file_record(p) if p.is_file() else {"status":"NOT_FOUND"}
    if (root/".git").exists():
        for name,args in [("commit",["rev-parse","HEAD"]),("branch",["branch","--show-current"]),
                          ("working_tree",["status","--short"])]:
            r=subprocess.run(["git","-C",str(root),*args],capture_output=True,text=True,check=False)
            result["git"][name]=r.stdout.strip() if r.returncode==0 else "UNAVAILABLE"
    env=result["environment"]
    env.update(python=platform.python_version(),platform=platform.platform(),free_bytes=shutil.disk_usage(root).free)
    for command in ["Rscript","plink","plink2","matlab","snakemake"]:env[command]=shutil.which(command)
    for package in ["numpy","scipy","matplotlib"]:env[package]=importlib.util.find_spec(package) is not None
    result["scientific_status"]="NOT_REEVALUATED_BY_ENVIRONMENT_AUDIT"
    write_json(out,result)
    return result
