import math
from pathlib import Path
import pytest
from brain6.io import write_json,write_tsv,sha256
from brain6.gwas import FIELDS
from brain6.pairs import BRAIN,SLEEP

@pytest.fixture
def source_factory(tmp_path):
    counter=[0]
    def make(rows=None,**changes):
        counter[0]+=1;i=counter[0]
        if rows is None:
            rows=[dict(zip(FIELDS,[f"rs{j}",1,j*100,"A","C",.1,.02,.00001,10000,.3,.99])) for j in range(1,11)]
        path=tmp_path/f"source{i}.tsv";write_tsv(path,FIELDS,rows)
        s={"trait_id":f"trait{i}","path":str(path),"sha256":sha256(path),"study_id":"SYNTHETIC",
           "phenotype_definition":"Artificial unit-test phenotype","ancestry":"EUR","genome_build":"GRCh37",
           "effect_scale":"beta","n_semantics":"total","synthetic":True,
           "column_map":{k:k for k in FIELDS},"qc":{"require_eaf":True,"min_n":100}}
        s.update(changes)
        config=tmp_path/f"spec{i}.json";write_json(config,s)
        return config,s,rows
    return make

@pytest.fixture
def atlas_factory(tmp_path):
    def make():
        rows=[]
        for b in BRAIN+[f"other{i}" for i in range(27)]:
            for i,s in enumerate(SLEEP):
                rows.append(dict(sleep_trait=s,disease_trait=b,rg=.1+i*.01,se=.02,p=.0001,fdr=.001+i*.001,analysis_status="PRIMARY"))
        matrix=tmp_path/"matrix.tsv";write_tsv(matrix,list(rows[0]),rows)
        cfg={"synthetic":True,"expected_matrix_rows":396,"matrix_qc_column":"analysis_status",
             "primary_qc_values":["PRIMARY"],"native_policy":{"pleiotropy":"PLACO_PLUS"}}
        config=tmp_path/"cfg.json";write_json(config,cfg)
        decisions=[dict(disease_trait=b,sleep_trait="insomnia",role="PRIMARY",reason="Synthetic fixture") for b in BRAIN]
        d=tmp_path/"decisions.tsv";write_tsv(d,list(decisions[0]),decisions)
        return matrix,config,d,rows,cfg,decisions
    return make
