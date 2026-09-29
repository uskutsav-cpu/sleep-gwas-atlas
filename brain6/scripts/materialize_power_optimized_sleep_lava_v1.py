#!/usr/bin/env python3
"""Materialize one Brain6 power-screen trait against frozen LAVA loci/reference."""
from __future__ import annotations
import csv,gzip,hashlib,io,json,math,sqlite3
from pathlib import Path
ROOT=Path(__file__).resolve().parents[2]
BASE=ROOT/'brain6/results/power_optimized_sensitivity_v1'
NORMAL=BASE/'normalized/sleep_duration_continuous_dashti_2019'
REFERENCE=ROOT/'ref/lava/ukb_v1.1'
LOCUS=ROOT/'ref/lava/blocks_s2500_m25_f1_w200.GRCh37_hg19.locfile'
FAMILY=ROOT/'brain6/config/lava_family_canonical_v3.json'
OUT=BASE/'lava_inputs_v1'
INFO_FIELDS=('SNP','CHR','POS','A1','A2','NOBS','MISS','FREQ','NCORRS')
def sha(p):
 h=hashlib.sha256()
 with p.open('rb') as f:
  for b in iter(lambda:f.read(1024*1024),b''):h.update(b)
 return h.hexdigest()
def reference_ids(path,chrom):
 with path.open(newline='') as f:
  r=csv.DictReader(f,delimiter='\t');
  if tuple(r.fieldnames or ())!=INFO_FIELDS:raise ValueError('Reference index schema mismatch')
  out=set()
  for row in r:
   if int(row['CHR'])!=chrom:raise ValueError('Cross chromosome reference record')
   snp=row['SNP'].lower()
   if not snp or snp in out:raise ValueError('Empty/duplicate reference ID')
   out.add(snp)
 return out
def main():
 if OUT.exists():raise FileExistsError(f'Refusing to overwrite {OUT}')
 lock=json.loads(FAMILY.read_text())
 if sha(LOCUS)!=lock['locus_definition']['sha256']:raise ValueError('Frozen locus hash mismatch')
 loci=[]
 with LOCUS.open() as f:
  header=f.readline().split()
  if header!=['LOC','CHR','START','STOP']:raise ValueError('Locus schema mismatch')
  for l in f:
   if l.strip():
    a=l.split(); loci.append({'LOC':a[0],'CHR':int(a[1]),'START':int(a[2]),'STOP':int(a[3])})
 if len(loci)!=2495:raise ValueError('Expected 2,495 frozen loci')
 ref_prov=json.loads((REFERENCE/'reference.provenance.json').read_text())
 if ref_prov.get('verification')!='SHA-256 verified after official HTTPS acquisition':raise ValueError('Reference seal missing')
 conn=sqlite3.connect(f'file:{NORMAL/"variants.sqlite"}?mode=ro',uri=True)
 OUT.mkdir(parents=True)
 record=[]; coverage=[]; genome_overlap=0; total_source=0
 try:
  for chrom in range(1,23):
   refpath=REFERENCE/f'lava-ukb-v1.1_chr{chrom}.info'; ids=reference_ids(refpath,chrom)
   srcids={str(r[0]).lower() for r in conn.execute('select snp from variants where chr=?',(chrom,))}; genome_overlap+=len(ids & srcids); total_source+=len(srcids)
   for loc in (x for x in loci if x['CHR']==chrom):
    d=OUT/f"locus_{loc['LOC']}"; d.mkdir()
    sp=d/'sleep_duration_continuous_dashti_2019.sumstats.tsv.gz'
    raw=sp.open('xb')
    with gzip.GzipFile(filename='',mode='wb',compresslevel=6,fileobj=raw,mtime=0) as gz, io.TextIOWrapper(gz,encoding='utf-8',newline='') as t:
     w=csv.writer(t,delimiter='\t',lineterminator='\n'); w.writerow(('SNP','A1','A2','Z','N'))
     n=0
     for snp,a1,a2,beta,se,N in conn.execute('select snp,a1,a2,beta,se,n from variants where chr=? and bp>=? and bp<=? order by bp,snp',(chrom,loc['START'],loc['STOP'])):
      snp=str(snp).lower()
      if snp not in ids:continue
      beta,se,N=map(float,(beta,se,N)); z=beta/se
      if not all(map(math.isfinite,(z,N))) or N<=0:raise ValueError('Invalid z/N')
      w.writerow((snp,a1,a2,format(z,'.17g'),format(N,'.17g'))); n+=1
    info=d/'input_info.tsv'
    with info.open('x',newline='') as f:
     w=csv.writer(f,delimiter='\t',lineterminator='\n');w.writerow(('phenotype','cases','controls','filename'));w.writerow(('sleep_duration_continuous_dashti_2019',446118,0,sp.name))
    coverage.append({'locus_id':loc['LOC'],'chromosome':chrom,'reference_shared_variants':n})
    record.append({'locus_id':loc['LOC'],'input_info_sha256':sha(info),'sumstats_sha256':sha(sp),'sumstats_rows':n})
 finally:conn.close()
 with (OUT/'reference_coverage.tsv').open('x',newline='') as f:
  w=csv.DictWriter(f,fieldnames=('locus_id','chromosome','reference_shared_variants'),delimiter='\t',lineterminator='\n');w.writeheader();w.writerows(coverage)
 p={'schema_version':1,'analysis_id':'brain6_power_optimized_sensitivity_v1_longsleep_trait_only_screen','scope':'Brain6-only diagnostic inputs; no FI outputs; no canonical v3 modification','source_receipt_sha256':sha(NORMAL/'receipt.json'),'source_variants_after_qc':6344850,'frozen_reference_unique_overlap':genome_overlap,'source_unique_variants':total_source,'reference_overlap_fraction_of_source':genome_overlap/total_source,'loci':len(loci),'locus_definition_sha256':sha(LOCUS),'reference_provenance_sha256':sha(REFERENCE/'reference.provenance.json'),'reference_coverage_rows':len(coverage),'records':record}
 p['coverage_tsv_sha256']=sha(OUT/'reference_coverage.tsv')
 with (OUT/'materialization.provenance.json').open('x') as f:json.dump(p,f,indent=2,sort_keys=True);f.write('\n')
 print(json.dumps({k:p[k] for k in ('analysis_id','source_unique_variants','frozen_reference_unique_overlap','reference_overlap_fraction_of_source','loci','coverage_tsv_sha256')},sort_keys=True))
if __name__=='__main__':main()
