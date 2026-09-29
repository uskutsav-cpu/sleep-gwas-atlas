"""Command line. Errors are explicit and nonzero; no implicit network runs."""
from __future__ import annotations
import argparse
import json
import sys
from pathlib import Path
from .io import ContractError,read_json


def parser():
    p=argparse.ArgumentParser(prog="brain6",description="Auditable six-brain-disorder sleep-GWAS extension")
    sub=p.add_subparsers(dest="command",required=True)
    a=sub.add_parser("demo",help="Run artificial software integration example, not GWAS inference")
    a.add_argument("--out",required=True);a.add_argument("--no-figure",action="store_true")
    a=sub.add_parser("audit",help="Read-only audit of actual local checkout")
    a.add_argument("--repo",required=True);a.add_argument("--out",required=True)
    a=sub.add_parser("rank",help="Rank the 72 measured brain-subset rows using existing family FDR")
    a.add_argument("--matrix",required=True);a.add_argument("--config",required=True);a.add_argument("--out",required=True)
    a=sub.add_parser("freeze",help="Freeze explicitly reviewed decisions; never invent pairs")
    for k in ["matrix","config","decisions","out","reviewer"]:a.add_argument("--"+k,required=True)
    a.add_argument("--synthetic",action="store_true")
    a=sub.add_parser("verify-lock");a.add_argument("path")
    a=sub.add_parser("verify-artifact");a.add_argument("path")
    a=sub.add_parser("download")
    for k in ["url","out","sha256"]:a.add_argument("--"+k,required=True)
    a.add_argument("--bytes",required=True,type=int)
    a=sub.add_parser("normalize")
    a.add_argument("--source",required=True);a.add_argument("--root",required=True);a.add_argument("--name",required=True)
    a.add_argument("--synthetic",action="store_true")
    a=sub.add_parser("join")
    for k in ["left","right","root","name"]:a.add_argument("--"+k,required=True)
    a.add_argument("--min-variants",type=int,default=1000000)
    a.add_argument("--min-overlap",type=float,default=.8)
    a.add_argument("--synthetic",action="store_true")
    a=sub.add_parser("index-pair")
    for k in ["pair","root","name"]:a.add_argument("--"+k,required=True)
    a=sub.add_parser("split-pair")
    for k in ["pair","root","name"]:a.add_argument("--"+k,required=True)
    a.add_argument("--chunk-rows",type=int,default=20000)
    a=sub.add_parser("collate-placo")
    for k in ["manifest","root","name"]:a.add_argument("--"+k,required=True)
    a.add_argument("--expected-rows",type=int,required=True);a.add_argument("--n-pairs",type=int,required=True)
    a.add_argument("--max-failure-rate",type=float,default=.001);a.add_argument("--synthetic",action="store_true")
    a=sub.add_parser("map-loci")
    for k in ["leads","blocks","out"]:a.add_argument("--"+k,required=True)
    a=sub.add_parser("prepare-locus")
    for k in ["index","ld-manifest","locus","root","name"]:a.add_argument("--"+k,required=True)
    a.add_argument("--min-coverage",type=float,default=.8)
    a.add_argument("--max-memory-gib",type=float,default=4)
    a.add_argument("--synthetic",action="store_true")
    a=sub.add_parser("run-job")
    a.add_argument("--job",required=True);a.add_argument("--root",required=True);a.add_argument("--synthetic",action="store_true")
    a=sub.add_parser("make-job")
    for k in ["settings","method","name","bindings","out-root"]:a.add_argument("--"+k,required=True)
    a.add_argument("--reviewed",action="store_true");a.add_argument("--synthetic",action="store_true")
    a=sub.add_parser("pipeline")
    a.add_argument("--runtime",required=True);a.add_argument("--only");a.add_argument("--synthetic",action="store_true")
    a=sub.add_parser("import-legacy")
    for k in ["manifest","root","name"]:a.add_argument("--"+k,required=True)
    a=sub.add_parser("integrate")
    for k in ["loci","molecular","out"]:a.add_argument("--"+k,required=True)
    a=sub.add_parser("cell-enrichment")
    a.add_argument("--table",required=True);a.add_argument("--out",required=True)
    a.add_argument("--permutations",type=int,default=10000);a.add_argument("--seed",type=int,default=20260908)
    a=sub.add_parser("cross-disorder")
    a.add_argument("--evidence",required=True);a.add_argument("--out",required=True)
    a=sub.add_parser("report");a.add_argument("--root",required=True);a.add_argument("--out",required=True)
    a=sub.add_parser("forest");a.add_argument("--ranking",required=True);a.add_argument("--out",required=True)
    a=sub.add_parser("heatmap");a.add_argument("--ranking",required=True);a.add_argument("--out",required=True)
    a=sub.add_parser("manhattan");a.add_argument("--results",required=True);a.add_argument("--out",required=True)
    a.add_argument("--threshold",type=float);a.add_argument("--bin-bp",type=int,default=100000)
    a=sub.add_parser("qq");a.add_argument("--sqlite",required=True);a.add_argument("--out",required=True)
    a=sub.add_parser("plan-placo")
    for k in ["pair","chunks","pair-lock","pair-id","source","source-sha256","out"]:a.add_argument("--"+k,required=True)
    a.add_argument("--reviewed",action="store_true")
    a=sub.add_parser("execute-placo")
    a.add_argument("--plan",required=True);a.add_argument("--synthetic",action="store_true")
    a=sub.add_parser("assemble-preparation")
    for k in ["pair-lock","sources","output-root","runtime-out"]:a.add_argument("--"+k,required=True)
    a.add_argument("--reviewed",action="store_true")
    a=sub.add_parser("deep-followup");a.add_argument("--settings",required=True);a.add_argument("--synthetic",action="store_true")
    a=sub.add_parser("mr-candidates",help="Select exposure-only MR candidates before native LD clumping")
    for k in ["exposure","root","name"]:a.add_argument("--"+k,required=True)
    a.add_argument("--instrument-p",type=float,default=5e-8);a.add_argument("--minimum-f",type=float,default=10)
    a=sub.add_parser("mr-harmonize",help="Extract outcomes after verified exposure-only LD clumping")
    for k in ["candidates","clump","outcome","root","name"]:a.add_argument("--"+k,required=True)
    a.add_argument("--max-instrument-r2",type=float,default=.001)
    a.add_argument("--minimum-clump-kb",type=int,default=10000)
    for cmd in ["mr-family","replication-family","campaign-audit"]:
        a=sub.add_parser(cmd)
        for k in ["manifest","root","name"]:a.add_argument("--"+k,required=True)
    a=sub.add_parser("campaign-template")
    a.add_argument("--pair-lock",required=True);a.add_argument("--out",required=True)
    a=sub.add_parser("qtl-index")
    for k in ["source","root","name"]:a.add_argument("--"+k,required=True)
    a.add_argument("--synthetic",action="store_true")
    a=sub.add_parser("qtl-locus")
    for k in ["gwas","qtl","feature","region","root","name"]:a.add_argument("--"+k,required=True)
    a.add_argument("--minimum-variants",type=int,default=20)
    a=sub.add_parser("annotation-scores")
    for k in ["variants","annotations","root","name"]:a.add_argument("--"+k,required=True)
    a.add_argument("--synthetic",action="store_true")
    a=sub.add_parser("molecular-followup")
    a.add_argument("--settings",required=True);a.add_argument("--synthetic",action="store_true")
    for cmd in ["joint-placo-family","local-family","pathways24","compare24","robustness24","provenance24","cell-family24","validate24"]:
        a=sub.add_parser(cmd)
        for key in ["manifest","root","name"]:a.add_argument("--"+key,required=True)
    a=sub.add_parser("doctor24");a.add_argument("--out",required=True);a.add_argument("--expected")
    for cmd in ["atlas-review","covariance-review","source-inventory","shard-audit"]:
        a=sub.add_parser(cmd)
        for key in ["repo","root","name"]:a.add_argument("--"+key,required=True)
        if cmd=="atlas-review":a.add_argument("--archive-date")
        if cmd=="source-inventory":a.add_argument("--hash-dense",action="store_true")
        if cmd=="shard-audit":a.add_argument("--trait")
    a=sub.add_parser("plan24");a.add_argument("--pair-lock",required=True);a.add_argument("--out",required=True)
    a=sub.add_parser("audit24")
    for key in ["manifest","root","name"]:a.add_argument("--"+key,required=True)
    a.add_argument("--synthetic",action="store_true")
    a=sub.add_parser("seal24")
    for key in ["manifest","audit","root","name"]:a.add_argument("--"+key,required=True)
    a.add_argument("--synthetic",action="store_true")
    a=sub.add_parser("host-audit24")
    a.add_argument("--repo",required=True);a.add_argument("--out",required=True)
    a.add_argument("--archive-date");a.add_argument("--hash-dense",action="store_true");a.add_argument("--scan-shards",action="store_true")
    return p


def main(argv=None):
    a=parser().parse_args(argv)
    try:
        if a.command=="mr-candidates":
            from .mr_workflow import candidates
            result=str(candidates(a.exposure,a.root,a.name,instrument_p=a.instrument_p,minimum_f=a.minimum_f))
        elif a.command=="mr-harmonize":
            from .mr_workflow import harmonize
            result=str(harmonize(a.candidates,a.clump,a.outcome,a.root,a.name,max_instrument_r2=a.max_instrument_r2,minimum_clump_kb=a.minimum_clump_kb))
        elif a.command=="mr-family":
            from .mr_workflow import collate
            result=str(collate(a.manifest,a.root,a.name))
        elif a.command=="replication-family":
            from .replication import evaluate
            result=str(evaluate(a.manifest,a.root,a.name))
        elif a.command=="campaign-template":
            from .campaign import template
            result=template(a.pair_lock,a.out)
        elif a.command=="campaign-audit":
            from .campaign import audit
            result=str(audit(a.manifest,a.root,a.name))
        elif a.command=="qtl-index":
            from .qtl import index
            result=str(index(a.source,a.root,a.name,synthetic=a.synthetic))
        elif a.command=="qtl-locus":
            from .qtl import locus
            result=str(locus(a.gwas,a.qtl,a.feature,read_json(a.region),a.root,a.name,minimum_variants=a.minimum_variants))
        elif a.command=="annotation-scores":
            from .annotations import score
            result=str(score(a.variants,a.annotations,a.root,a.name,synthetic=a.synthetic))
        elif a.command=="demo":
            from .demo import run_demo
            result=run_demo(a.out,make_figure=not a.no_figure)
        elif a.command=="audit":
            from .audit import audit_repo
            result=audit_repo(a.repo,a.out)
        elif a.command=="rank":
            from .pairs import write_ranking
            result={"ranked_rows":len(write_ranking(a.matrix,a.config,a.out)),"out":a.out}
        elif a.command=="freeze":
            from .pairs import freeze_pairs
            result=freeze_pairs(a.matrix,a.config,a.decisions,a.out,a.reviewer,synthetic=a.synthetic)
        elif a.command=="verify-lock":
            from .pairs import verify_lock
            result=verify_lock(a.path)
        elif a.command=="verify-artifact":
            from .artifacts import verify_artifact
            result=verify_artifact(a.path)
        elif a.command=="download":
            from .download import download
            result=str(download(a.url,a.out,a.sha256,a.bytes))
        elif a.command=="normalize":
            from .gwas import normalize
            result=str(normalize(a.source,a.root,a.name,synthetic=a.synthetic))
        elif a.command=="join":
            from .gwas import join_pair
            result=str(join_pair(a.left,a.right,a.root,a.name,min_variants=a.min_variants,min_overlap=a.min_overlap,synthetic=a.synthetic))
        elif a.command=="index-pair":
            from .ld import index_pair
            result=str(index_pair(a.pair,a.root,a.name))
        elif a.command=="split-pair":
            from .pleiotropy import split_pair
            result=str(split_pair(a.pair,a.root,a.name,chunk_rows=a.chunk_rows))
        elif a.command=="collate-placo":
            from .pleiotropy import collate_placo
            result=str(collate_placo(a.manifest,a.root,a.name,expected_rows=a.expected_rows,n_pairs=a.n_pairs,max_failure_rate=a.max_failure_rate,synthetic=a.synthetic))
        elif a.command=="map-loci":
            from .ld import map_loci
            result=map_loci(a.leads,a.blocks,a.out)
        elif a.command=="prepare-locus":
            from .ld import prepare_locus
            result=str(prepare_locus(a.index,a.ld_manifest,read_json(a.locus),a.root,a.name,
                       max_working_bytes=int(a.max_memory_gib*2**30),min_coverage=a.min_coverage,synthetic=a.synthetic))
        elif a.command=="run-job":
            from .executor import run_job
            result=run_job(a.job,a.root,allow_synthetic=a.synthetic)
        elif a.command=="make-job":
            from .pipeline import create_job
            bindings=read_json(a.bindings)
            result=str(create_job(a.settings,method=a.method,job_id=a.name,inputs=bindings["inputs"],
                       outputs=bindings["outputs"],root=a.out_root,reviewed=a.reviewed,synthetic=a.synthetic))
        elif a.command=="pipeline":
            from .pipeline import run_pipeline
            result=run_pipeline(a.runtime,only=a.only,allow_synthetic=a.synthetic)
        elif a.command=="import-legacy":
            from .legacy import import_legacy
            result=str(import_legacy(a.manifest,a.root,a.name))
        elif a.command=="integrate":
            from .evidence import integrate
            result={"rows":len(integrate(a.loci,a.molecular,a.out))}
        elif a.command=="cell-enrichment":
            from .evidence import cell_enrichment
            result=cell_enrichment(a.table,a.out,permutations=a.permutations,seed=a.seed)
        elif a.command=="cross-disorder":
            from .evidence import cross_disorder
            result=cross_disorder(a.evidence,a.out)
        elif a.command=="report":
            from .report import report
            result=report(a.root,a.out)
        elif a.command=="forest":
            from .report import forest
            forest(a.ranking,a.out);result={"out":a.out}
        elif a.command=="heatmap":
            from .plots import heatmap
            heatmap(a.ranking,a.out);result={"out":a.out}
        elif a.command=="manhattan":
            from .plots import manhattan
            result=manhattan(a.results,a.out,threshold=a.threshold,bin_bp=a.bin_bp)
        elif a.command=="qq":
            from .plots import qq
            qq(a.sqlite,a.out);result={"out":a.out}
        elif a.command=="plan-placo":
            from .planning import placo_jobs
            result=str(placo_jobs(a.pair,a.chunks,a.pair_lock,a.pair_id,a.source,a.source_sha256,a.out,reviewed=a.reviewed))
        elif a.command=="execute-placo":
            from .planning import execute_placo_plan
            result=execute_placo_plan(a.plan,allow_synthetic=a.synthetic)
        elif a.command=="assemble-preparation":
            from .assembly import assemble_preparation
            result=assemble_preparation(a.pair_lock,a.sources,a.output_root,a.runtime_out,reviewed=a.reviewed)
        elif a.command=="molecular-followup":
            from .molecular import run
            result=run(a.settings,allow_synthetic=a.synthetic)
        elif a.command=="deep-followup":
            from .followup import run_followup
            result=run_followup(a.settings,allow_synthetic=a.synthetic)
        elif a.command in {"joint-placo-family","local-family","pathways24","compare24","robustness24","provenance24","cell-family24","validate24"}:
            if a.command=="joint-placo-family":
                from .family24 import joint_placo as fn
            elif a.command=="local-family":
                from .family24 import local_family as fn
            elif a.command=="pathways24":
                from .pathway24 import run as fn
            elif a.command=="compare24":
                from .compare24 import compare as fn
            elif a.command=="robustness24":
                from .robustness24 import summarize as fn
            elif a.command=="cell-family24":
                from .cell24 import run as fn
            elif a.command=="validate24":
                from .validation24 import run as fn
            else:
                from .provenance24 import audit as fn
            result=fn(a.manifest,a.root,a.name)
        elif a.command=="doctor24":
            from .doctor import inspect_environment
            result=inspect_environment(a.out,expected=a.expected)
        elif a.command=="atlas-review":
            from .atlas_review import review
            result=review(a.repo,a.root,a.name,archive_date=a.archive_date)
        elif a.command=="covariance-review":
            from .atlas_review import inspect_covariance
            result=inspect_covariance(a.repo,a.root,a.name)
        elif a.command=="source-inventory":
            from .inputs24 import inventory
            result=inventory(a.repo,a.root,a.name,hash_dense=a.hash_dense)
        elif a.command=="shard-audit":
            from .shards24 import audit_shards
            result=audit_shards(a.repo,a.root,a.name,trait=a.trait)
        elif a.command=="plan24":
            from .completion24 import draft
            result=draft(a.pair_lock,a.out)
        elif a.command=="audit24":
            from .completion24 import audit
            result=audit(a.manifest,a.root,a.name,allow_synthetic=a.synthetic)
        elif a.command=="seal24":
            from .completion24 import seal
            result=seal(a.manifest,a.audit,a.root,a.name,allow_synthetic=a.synthetic)
        elif a.command=="host-audit24":
            from .host24 import inspect
            result=inspect(a.repo,a.out,archive_date=a.archive_date,hash_dense=a.hash_dense,scan_shards=a.scan_shards)
        else:raise ContractError("Unknown command")
        print(json.dumps(result,indent=2,default=str,allow_nan=False))
        return 0
    except (ContractError,OSError,KeyError,ValueError) as e:
        print(f"brain6: {type(e).__name__}: {e}",file=sys.stderr)
        return 2

if __name__=="__main__":raise SystemExit(main())
