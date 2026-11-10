#!/usr/bin/env python3
"""Apply the published PsychAD model quality criteria to project molecular candidate genes."""
from pathlib import Path
import sqlite3,pandas as pd
root=Path(__file__).resolve().parents[1];o=root/'outputs';p=pd.read_csv(o/'stage5_variant_gene_tissue_evidence_atlas.csv',usecols=['gene_id','gene_symbol']).dropna().drop_duplicates('gene_id');rows=[]
for cell in ['Endo','Mural']:
 db=sqlite3.connect(f'file:{root}/work/PsychAD_snTWAS_EUR_Class_{cell}.db?mode=ro',uri=True)
 x=pd.read_sql_query('select gene,gene_name,n_snps_in_model,"pred.perf.R2" as r2,"pred.perf.pval" as pvalue,"pred.perf.qval" as qvalue from extra',db);db.close();assert x.gene.is_unique
 x=p.merge(x,left_on='gene_id',right_on='gene',how='left',validate='one_to_one');x.insert(0,'cell_class',cell);x['passes_paper_filter']=(x.n_snps_in_model>0)&(x.r2>=.01)&(x.qvalue<=.05);rows.append(x)
y=pd.concat(rows);y.to_csv(o/'stage6_psychad_priority_gene_model_qc.csv',index=False);print(y.groupby('cell_class').agg(matched=('gene','count'),pass_count=('passes_paper_filter','sum')).to_string())
