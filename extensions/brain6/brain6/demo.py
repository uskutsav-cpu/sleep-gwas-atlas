"""Small synthetic integration exercise. NEVER empirical project evidence."""
from __future__ import annotations
import csv
import json
import math
from pathlib import Path
import numpy as np
from .artifacts import transaction,verify_artifact
from .gwas import normalize,join_pair,FIELDS
from .ld import index_pair
from .pairs import BRAIN,SLEEP,freeze_pairs,write_ranking
from .pleiotropy import split_pair
from .report import report,forest
from .io import require,write_json,write_tsv,sha256
from .stats import two_sided_normal,ivw


def run_demo(out, *, make_figure=True):
    out=Path(out).resolve()
    require(not out.exists(),"Choose a new empty demo output directory")
    out.mkdir(parents=True)
    (out/"SYNTHETIC_ONLY.txt").write_text("All values in this directory are artificial software-test inputs.\nNo LDSC, PLACO, SuSiE, coloc, LAVA, GenomicSEM, or empirical MR was run.\n")
    rng=np.random.default_rng(20260908)
    inputs=out/"inputs";inputs.mkdir()
    n=1100
    z=rng.normal(size=n)
    for number in [1,2]:
        rows=[]
        for i in range(n):
            a1,a2="A","C"
            beta=(z[i]+.2*(number==2))*0.03
            eaf=.3
            if number==2 and i%3==0:
                a1,a2=a2,a1;beta=-beta;eaf=1-eaf
            rows.append(dict(zip(FIELDS,[f"rs{i+1}",i%22+1,(i//22+1)*10000,a1,a2,
                         beta,.03,two_sided_normal(beta/.03),10000,eaf,.99])))
        # Two identical IDs are intentionally quarantined, not silently kept.
        rows += [dict(rows[0])]
        path=inputs/f"trait{number}.tsv"
        write_tsv(path,FIELDS,rows)
        source={"trait_id":f"synthetic{number}","path":str(path),"sha256":sha256(path),
                "study_id":"SYNTHETIC","phenotype_definition":"Artificial test trait",
                "genome_build":"GRCh37","ancestry":"EUR","effect_scale":"beta",
                "n_semantics":"total","synthetic":True,"column_map":{k:k for k in FIELDS},
                "qc":{"min_maf":.01,"require_eaf":True,"min_n":100}}
        spec=inputs/f"source{number}.json";write_json(spec,source)
        normalize(spec,out,f"normalized{number}",synthetic=True)
    join_pair(out/"normalized1",out/"normalized2",out,"joined",min_variants=100,
              min_overlap=.99,synthetic=True)
    index_pair(out/"joined",out,"indexed")
    split_pair(out/"joined",out,"chunks",chunk_rows=200)
    diseases=BRAIN+[f"synthetic_external{i}" for i in range(27)]
    rows=[]
    for disease in diseases:
        for sleep in SLEEP:
            rg=.12+.02*rng.random();se=.025
            rows.append(dict(sleep_trait=sleep,disease_trait=disease,rg=rg,se=se,
                       p=two_sided_normal(rg/se),fdr=.001+.01*rng.random(),analysis_status="PRIMARY"))
    matrix=inputs/"synthetic_matrix.tsv";write_tsv(matrix,list(rows[0]),rows)
    cfg={"synthetic":True,"brain_traits":BRAIN,"sleep_traits":SLEEP,"expected_matrix_rows":396,
         "matrix_qc_column":"analysis_status","primary_qc_values":["PRIMARY"],
         "native_policy":{"pleiotropy":"PLACO_PLUS"},"experiment":"SYNTHETIC_DEMO"}
    cf=inputs/"config.json";write_json(cf,cfg)
    ranking=out/"SYNTHETIC_ranking.tsv"
    ranked=write_ranking(matrix,cf,ranking)
    decisions=[]
    for b in BRAIN:
        choice=next(r for r in ranked if r["disease_trait"]==b and r["selection_rank"]==1)
        decisions.append(dict(disease_trait=b,sleep_trait=choice["sleep_trait"],role="PRIMARY",reason="Synthetic test selection"))
    dp=inputs/"decisions.tsv";write_tsv(dp,list(decisions[0]),decisions)
    freeze_pairs(matrix,cf,dp,out/"SYNTHETIC_pair_lock.json","synthetic-test-review",synthetic=True)
    bx=np.array([.1,.2,.15,.12]);by=2*bx
    diagnostic=ivw(bx,np.ones(4)*.01,by,np.ones(4)*.03)
    diagnostic["synthetic"]=True;write_json(out/"SYNTHETIC_ivw_diagnostic.json",diagnostic)
    if make_figure:forest(ranking,out/"SYNTHETIC_forest.png")
    report(out,out/"SYNTHETIC_report.md")
    result={"synthetic":True,"normalized_variants_per_trait":1099,"joined_variants":1099,
            "ranked_brain_pairs":72,"atlas_rows":396,"native_genetics_methods_run":[],
            "message":"Software integration demonstration only; no empirical discoveries."}
    write_json(out/"demo_summary.json",result)
    return result
