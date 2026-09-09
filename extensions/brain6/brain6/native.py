"""Command-line adapters for LDSC and PLINK, with explicit input contracts."""
from __future__ import annotations
import csv
import math
import os
import re
import shutil
import subprocess
from pathlib import Path
from .artifacts import verify_artifact
from .io import (ContractError, read_json, read_tsv, require, write_json, write_tsv)


NUMBER = r"[-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][-+]?\d+)?"


def parse_h2(text: str) -> dict:
    estimate = re.search(r"Total (Observed|Liability) scale h2:\s*("+NUMBER+r")\s*\(("+NUMBER+r")\)",text)
    intercept = re.search(r"Intercept:\s*("+NUMBER+r")\s*\(("+NUMBER+r")\)",text)
    require(estimate is not None and intercept is not None, "Incomplete LDSC h2 log")
    h2, se, ic, ise = float(estimate[2]), float(estimate[3]), float(intercept[1]), float(intercept[2])
    require(all(math.isfinite(x) for x in [h2,se,ic,ise]) and se>0 and ise>=0, "Invalid LDSC h2 result")
    return dict(scale=estimate[1],h2=h2,se=se,z=h2/se,intercept=ic,intercept_se=ise)


def parse_rg(text: str) -> list[dict]:
    lines = text.splitlines()
    fields = None; rows = []
    for line in lines:
        tokens = line.split()
        if tokens[:2] == ["p1","p2"] and {"rg","se","p"} <= set(tokens):
            fields = tokens; continue
        if fields is not None:
            if len(tokens) != len(fields):
                if rows: break
                continue
            raw = dict(zip(fields,tokens))
            try:
                row = {"p1":raw["p1"],"p2":raw["p2"], **{k:float(raw[k]) for k in ["rg","se","p"]}}
            except ValueError:
                continue
            require(all(math.isfinite(row[k]) for k in ["rg","se","p"]), "Nonfinite LDSC rg")
            require(abs(row["rg"])<=1 and row["se"]>0 and 0<=row["p"]<=1, "Nonphysical LDSC rg result")
            if "gcov_int" in raw:
                row["gcov_int"] = float(raw["gcov_int"])
            rows.append(row)
    require(rows, "No complete LDSC genetic-correlation summary")
    return rows


def _run(argv):
    require(shutil.which(str(argv[0])) is not None, f"Missing executable: {argv[0]}")
    result = subprocess.run([str(x) for x in argv], check=False)
    require(result.returncode==0, f"Native command failed ({result.returncode}): {argv[0]}")


def _ldsc_base(c, out):
    return [c["python"],c["ldsc_script"],"--ref-ld-chr",c["ref_ld_chr"],
            "--w-ld-chr",c["weights_chr"],"--out",str(out/"ldsc")]


def execute_adapter(settings):
    c=read_json(settings)
    root=os.environ.get("BRAIN6_JOB_DIR")
    require(root, "Use brain6 run-job for synchronous execution and output protection")
    out=Path(root)
    mode=c["mode"]
    if mode=="pleiofdr":
        from .matlab import run_pleiofdr
        return run_pleiofdr(c,out)
    if mode=="munge":
        argv=[c["python"],c["munge_script"],"--sumstats",c["sumstats"],
              "--merge-alleles",c["hapmap3"],"--snp","SNP","--a1","A1","--a2","A2",
              "--p","P","--N-col","N","--signed-sumstats","BETA,0","--out",str(out/"munged")]
        _run(argv)
        n=sum(1 for _ in read_tsv(out/"munged.sumstats.gz",["SNP","A1","A2","Z","N"]))
        require(n>=c["min_variants"], "Munged SNP count below policy")
        status={"status":"PASS","n_variants":n}
    elif mode=="h2":
        argv=_ldsc_base(c,out)+["--h2",c["sumstats"]]
        if c.get("population_prevalence") is not None:
            require(c.get("prevalence_citation") and 0<c["population_prevalence"]<1
                    and 0<c.get("sample_prevalence",0)<1,"Invalid/unsourced prevalence")
            argv += ["--pop-prev",str(c["population_prevalence"]),"--samp-prev",str(c["sample_prevalence"])]
        _run(argv)
        h=parse_h2((out/"ldsc.log").read_text())
        passed=(0<h["h2"]<=1 and h["z"]>=c["min_h2_z"] and
                c["intercept_min"]<=h["intercept"]<=c["intercept_max"])
        status={"status":"PASS" if passed else "FAILED_QC",**h}
        write_tsv(out/"h2.tsv",list(h),[h])
    elif mode=="rg":
        require(c.get("h2_qc_verified") is True,"Both h2/QC results must be verified")
        argv=_ldsc_base(c,out)+["--rg",c["sumstats1"]+","+c["sumstats2"]]
        # No --no-intercept or hardcoded overlap intercept shortcuts.
        _run(argv)
        rows=parse_rg((out/"ldsc.log").read_text())
        require(len(rows)==1,"Expected exactly one pair per job")
        write_tsv(out/"rg.tsv",list(rows[0]),rows)
        status={"status":"PASS","pairs":1,"note":"Raw P; original family FDR is not replaced"}
    elif mode=="clump":
        if c.get("input_receipt"):
            receipt=verify_artifact(c["input_receipt"])
            require(receipt.get("scientific_status")=="PASS", "Cannot clump rejected PLACO output")
        require(0<c["r2"]<1 and c["kb"]>0 and 0<c["p1"]<=1,"Invalid clump thresholds")
        n=0
        def assoc():
            nonlocal n
            for r in read_tsv(c["associations"],["SNP",c["p_column"]]):
                if r[c["p_column"]] in {"NA","","nan"}: continue
                p=float(r[c["p_column"]])
                require(math.isfinite(p) and 0<=p<=1,"Invalid clumping P")
                n += p <= c["p1"]
                yield {"SNP":r["SNP"],"P":p}
        write_tsv(out/"association.tsv",["SNP","P"],assoc())
        if n:
            _run([c["plink"],"--bfile",c["reference_prefix"],"--clump",str(out/"association.tsv"),
                  "--clump-p1",str(c["p1"]),"--clump-p2",str(c.get("p2",1)),
                  "--clump-r2",str(c["r2"]),"--clump-kb",str(c["kb"]),
                  "--out",str(out/"clump"),"--threads",str(c.get("threads",1))])
            require((out/"clump.clumped").exists(),"No clump output despite significant inputs; check reference coverage")
            with (out/"clump.clumped").open() as f:
                lines=[x.split() for x in f if x.strip()]
            require(lines and {"SNP","CHR","BP","P"}<=set(lines[0]),"Unexpected PLINK clump output schema")
            leads=[]
            for values in lines[1:]:
                row=dict(zip(lines[0],values))
                leads.append({k:row[k] for k in ["SNP","CHR","BP","P"]})
            require(leads,"No usable reference-supported clumps")
        else: leads=[]
        write_tsv(out/"leads.tsv",["SNP","CHR","BP","P"],leads)
        status={"status":"PASS" if leads else "NO_SIGNAL","clump_leads":len(leads),
                "note":"Map clump leads into intact predeclared LD blocks before counting loci"}
    elif mode=="ld":
        require(c.get("reference_ancestry")=="EUR" and c.get("genome_build") in {"GRCh37","GRCh38"},"Reference metadata required")
        require(c["start"]>=1 and c["stop"]>=c["start"],"Invalid locus bounds")
        # First materialize a locus-specific BIM that defines the actual SNP order.
        _run([c["plink"],"--bfile",c["reference_prefix"],"--chr",str(c["chromosome"]),
              "--from-bp",str(c["start"]),"--to-bp",str(c["stop"]),"--keep-allele-order",
              "--make-bed","--out",str(out/"locus")])
        n=sum(1 for _ in (out/"locus.bim").open())
        require(8*n*n*8<=c["max_working_bytes"],"BLOCKED_BY_COMPUTE; no locus thinning")
        _run([c["plink"],"--bfile",str(out/"locus"),"--keep-allele-order","--r","square",
              "--out",str(out/"signed"),"--threads",str(c.get("threads",1))])
        status={"status":"PASS","n_snps":n,"matrix_kind":"signed_r","counted_allele":"BIM_A1"}
    else: raise ContractError(f"Unknown native CLI adapter: {mode}")
    write_json(out/"status.json",status)
