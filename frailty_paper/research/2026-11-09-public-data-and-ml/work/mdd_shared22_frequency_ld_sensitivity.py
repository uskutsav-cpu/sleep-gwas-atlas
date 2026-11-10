#!/usr/bin/env python3
"""Explore 22-count marginal-Z sensitivity; no participant outcomes."""
import csv, hashlib, json, math
from pathlib import Path
import numpy as np
import pysam

root=Path(__file__).resolve().parent.parent;out=root/'outputs'
score_file=out/'unvalidated_22_shared_direction_score.csv';source=root/'work/mdd_symptoms/candidate_hits_private.json'
with score_file.open() as f:score=list(csv.DictReader(f))
assert len(score)==22 and len({r['variant'] for r in score})==22
hits={(x['chr'],x['pos'],frozenset(x['pair'])):x for x in json.loads(source.read_text())}
with (root/'work/1kg_phase3_panel.tsv').open() as f:panel=list(csv.DictReader(f,delimiter='\t'))
samples=sorted(r['sample'] for r in panel if r['super_pop']=='EUR');assert len(samples)==503
handles={}; cols=[]
for r in score:
 c=r['chr_GRCh37_assumed'];pos=int(r['position_GRCh37_assumed'])
 if c not in handles:
  handles[c]=pysam.VariantFile(f'https://hgdownload.soe.ucsc.edu/gbdb/hg19/1000Genomes/phase3/ALL.chr{c}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz')
 found=[v for v in handles[c].fetch(c,pos-1,pos) if v.pos==pos and v.id==r['variant']]
 if len(found)!=1 or len(found[0].alleles)!=2 or set(found[0].alleles)!={r['counted_allele'],r['other_allele']}:
  raise ValueError(f'ambiguous reference {r["variant"]}')
 rec=found[0]; a=rec.alleles.index(r['counted_allele']);values=[]
 for person in samples:
  gt=rec.samples[person].get('GT')
  if gt is None or len(gt)!=2 or any(x is None for x in gt):raise ValueError('missing reference genotype')
  values.append(sum(x==a for x in gt))
 cols.append(values)
genotypes=np.asarray(cols,dtype=float).T
R=np.corrcoef(genotypes,rowvar=False);chrom=np.array([r['chr_GRCh37_assumed'] for r in score]);R=np.where(chrom[:,None]==chrom[None,:],R,0.)
old=json.loads((out/'stage6_mdd_22_weighted_sensitivity.json').read_text())['results'];result=[]
for trait in ('insomnia','fatigue'):
 z=[];sd=[];ses=[]
 for r in score:
  key=(r['chr_GRCh37_assumed'],int(r['position_GRCh37_assumed']),frozenset((r['counted_allele'],r['other_allele'])))
  h=hits[key]['hits'][trait][0];n1=float(h['ncase']);n0=float(h['ncontrol']);p=(n1*float(h['fcas'])+n0*float(h['fcon']))/(n1+n0)
  assert 0<p<1 and float(h['se'])>0
  ses.append(float(h['se']));sd.append(math.sqrt(2*p*(1-p)));z.append(float(h['beta'])/float(h['se'])*(1 if h['a1']==r['counted_allele'] else -1))
 a=np.asarray(sd);z=np.asarray(z);v=float(a@z/np.sqrt(a@R@a));previous=next(x['weighted_z'] for x in old if x['trait']==trait)
 precision=np.asarray(ses)*a;inverse_se=1/np.asarray(ses);v2=float(inverse_se@z/np.sqrt(inverse_se@R@inverse_se))
 result.append({'trait':trait,'target_frequency_approx_z':v,'inverse_se_approx_z':v2,
                'se_times_target_hwe_sd_cv':float(precision.std(ddof=1)/precision.mean()),
                'prior_reference_approx_z':previous,'nominal_p':math.erfc(abs(v)/math.sqrt(2))})
receipt={'score_sha256':hashlib.sha256(score_file.read_bytes()).hexdigest(),'source_aggregate_sha256':hashlib.sha256(source.read_bytes()).hexdigest(),
         'n':503,'markers':22,'results':result,'method':'case-count-weighted frequencies and HWE scaling or inverse-SE weights; external EUR within-chromosome LD, zero between chromosomes',
         'limitations':'Postselected approximation only; meta-analysis frequencies, per-variant samples, allele coding, target LD and sample independence uncertain; no score-level prediction or calibrated test. Remote VCF full checksum not verified.'}
(out/'stage6_mdd_shared22_target_freq_ld_sensitivity.json').write_text(json.dumps(receipt,indent=2)+'\n');print(json.dumps(result,indent=2))
