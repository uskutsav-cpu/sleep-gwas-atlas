#!/usr/bin/env python3
"""Create a target-gene MetaXcan-compatible model and covariance slice."""
from pathlib import Path
import gzip,sqlite3
root=Path(__file__).resolve().parents[1];w=root/'work';p=w/'stage6_psychad_rbm6_model.db';p.unlink(missing_ok=True)
c=sqlite3.connect(p);c.execute('ATTACH DATABASE ? AS source',(str(w/'PsychAD_snTWAS_EUR_Class_Mural.db'),))
for t in ['construction','sample_info']:c.execute(f'create table {t} as select * from source.{t}')
for t in ['extra','extra_full','weights']:c.execute(f"create table {t} as select * from source.{t} where gene='ENSG00000004534'")
c.execute('alter table extra add column genename text');c.execute('alter table extra add column `n.snps.in.model` integer');c.execute('update extra set genename=gene_name, `n.snps.in.model`=n_snps_in_model');c.commit();c.close();n=0
with gzip.open(w/'PsychAD_snTWAS_EUR_Class_Mural.db.txt.gz','rt') as f, gzip.open(w/'stage6_psychad_rbm6_cov.txt.gz','wt') as o:
 o.write(next(f))
 for l in f:
  if l.startswith('ENSG00000004534\t'):o.write(l);n+=1
assert n==496;print('covariance rows',n)
