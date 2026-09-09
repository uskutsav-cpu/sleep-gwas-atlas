#!/usr/bin/env python3
"""Reproducible artificial integration exercise; no native inference is simulated."""
import argparse
from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from brain6.demo import run_demo
from brain6.io import read_json,read_tsv,write_json,write_tsv,sha256
from brain6.gwas import FIELDS
from brain6.qtl import REQUIRED,index,locus
from brain6.ld import index_pair
from brain6.mr_workflow import candidates,collate
from brain6.replication import evaluate
from brain6.annotations import score
from brain6.evidence import cell_enrichment
from brain6.campaign import template,audit
from brain6.plots import heatmap


def run(out):
    out=Path(out).resolve()
    base=run_demo(out)
    inp=out/'inputs';lock=out/'SYNTHETIC_pair_lock.json'
    # QTL effects are artificial and never described as genetic findings.
    q=[]
    for r in read_tsv(inp/'trait1.tsv'):
        if int(r['CHR'])==1:
            q.append(dict(zip(REQUIRED,['SYNTHETIC_GENE',r['CHR'],r['BP'],r['A2'],r['A1'],
               -float(r['BETA']),.01,r['P'],500,1-float(r['EAF'])])))
    qp=inp/'SYNTHETIC_qtl.tsv';write_tsv(qp,REQUIRED,q)
    qs=dict(path=str(qp),sha256=sha256(qp),synthetic=True,study_id='SYNTHETIC',genome_build='GRCh37',
          ancestry='EUR',molecular_trait='expression',tissue='SYNTHETIC_CORTEX',source_uri='synthetic://qtl',
          citation='Artificial software fixture',n_semantics='total',complete_cis_summary=True,
          column_map={k:k for k in REQUIRED})
    spec=inp/'SYNTHETIC_qtl.json';write_json(spec,qs)
    qi=index(spec,out,'qtl_index',synthetic=True)
    loc=locus(out/'normalized1',qi,'SYNTHETIC_GENE',dict(CHR=1,START=1,STOP=1000000),out,'qtl_locus')
    index_pair(loc,out,'qtl_pair_index')
    candidates(out/'normalized1',out,'mr_candidates')
    mf=inp/'SYNTHETIC_mr_family.json';write_json(mf,dict(reviewed=True,pair_lock=str(lock),results=[]))
    collate(mf,out,'mr_family_NOT_RUN')
    rf=inp/'SYNTHETIC_replication.json';write_json(rf,dict(reviewed=True,pair_lock=str(lock),results=[]))
    evaluate(rf,out,'replication_NOT_RUN')
    # 30 independent artificial blocks, one artificial cell annotation, reviewed
    # stratum labels for testing only. This is not a scientific enrichment study.
    bed=inp/'SYNTHETIC_cell.bed';bed.write_text(''.join(f'chr1\t{i*10000-1}\t{i*10000+100}\n' for i in range(1,16)))
    an=inp/'SYNTHETIC_annotations.json';write_json(an,dict(reviewed=True,synthetic=True,genome_build='GRCh37',
        independent_loci_reviewed=True,matched_controls_reviewed=True,annotations=[dict(cell_type='SYNTHETIC_CELL',
        path=str(bed),sha256=sha256(bed),genome_build='GRCh37',coordinate_system='BED0',
        citation='Artificial software fixture',source_uri='synthetic://cell')]))
    vs=[dict(pair_id='SYNTHETIC_PAIR',locus_id=f'SYNTHETIC_L{i}',SNP=f'rs{i}',CHR=1,BP=i*10000,
       PIP=1,selected=int(i%3==0),stratum='artificial_stratum') for i in range(1,31)]
    vp=inp/'SYNTHETIC_pip.tsv';write_tsv(vp,list(vs[0]),vs)
    scores=score(vp,an,out,'cell_scores',synthetic=True)
    cell_enrichment(scores/'scores.tsv',out/'SYNTHETIC_cell_enrichment.tsv',permutations=1000)
    ct=inp/'SYNTHETIC_campaign.json';template(lock,ct);cov=audit(ct,out,'campaign_INCOMPLETE')
    heatmap(out/'SYNTHETIC_ranking.tsv',out/'SYNTHETIC_heatmap.png')
    summary={**base,'qtl_local_variants':len(list(read_tsv(loc/'pair.tsv.gz'))),
      'mr_family_expected_tests':read_json(out/'mr_family_NOT_RUN/status.json')['expected_tests'],
      'mr_family_tested':0,'campaign':read_json(cov/'status.json'),
      'native_genetics_methods_run':[],
      'interpretation':'Only Python data preparation, accounting and artificial annotation testing ran. No empirical inference.'}
    write_json(out/'V02_DEMO_SUMMARY.json',summary)
    return summary

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--out',required=True);a=p.parse_args()
    import json
    print(json.dumps(run(a.out),indent=2))
