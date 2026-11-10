#!/usr/bin/env python3
"""Collapse repeated gene/trait rows but preserve source, cell and assay values."""
import argparse,csv,collections,hashlib,pathlib,gzip
p=argparse.ArgumentParser();p.add_argument('--gtex',required=True);p.add_argument('--traits',required=True);p.add_argument('--store',required=True);a=p.parse_args();root=pathlib.Path('outputs');store=pathlib.Path(a.store)
with (root/'stage5_reciprocal_snv_hg19_reference_check.csv').open() as f:
 mapping={f"chr{x['chr37']}:{x['pos37']}:{x['ref38']}:{x['alt38']}":x['variant_GRCh38'] for x in csv.DictReader(f) if x['reference_status']=='reference_same'}
with (root/'stage5_indel_hg19_mapping.csv').open() as f:
 for x in csv.DictReader(f):
  if x['status']=='mapped_reference_checked_normalized':mapping[f"chr{x['hg19_chrom']}:{x['hg19_position']}:{x['hg19_ref']}:{x['hg19_alt']}"]=x['variant_GRCh38']
with (root/'stage5_variant_multitrait_gwas_matrix.csv').open() as f:gw={x['variant_GRCh38_ref_alt']:x for x in csv.DictReader(f)}
keys=['variant','cell_type','library','best_library','type','log2FC','logPadj_BF','Log2Skew','Skew_logPadj','active_screen','emvar_screen'];groups={}
for file,source in [(a.gtex,'GTEx'),(a.traits,'complex_trait')]:
 with open(file) as f:
  for r in csv.DictReader(f):
   variant=mapping[r['variant']];group=(variant,source,*[r[k] for k in keys]);g=groups.setdefault(group,{'variant_GRCh38':variant,'source':source,**{k:r[k] for k in keys},'source_record_count':0,'source_contexts':set()});g['source_record_count']+=1;g['source_contexts'].add('|'.join(r.get(k,'') for k in ('gene','tissue','trait')))
fields=['variant_GRCh38','source',*keys,'source_record_count','source_contexts','model_DNASE_track','model_DNASE_raw','model_DNASE_quantile','candidate_gene_context','molecular_pip_max','Jansen_insomnia_p','FinnGen_R13_insomnia_p','Mak_FI_p','Mak_FP_p','Ye_FFS_p','caveat'];counts=collections.Counter()
with (root/'stage5_siraj2026_allcell_reporter_model_comparison.csv').open('w',newline='') as out:
 w=csv.DictWriter(out,fieldnames=fields);w.writeheader()
 for g in groups.values():
  variant=g['variant_GRCh38'];biosample={'SKNSH':'SK-N-SH','HEPG2':'HepG2','K562':'K562','HCT116':'HCT116','A549':'A549','GM12878':'GM12878'}.get(g['cell_type'])
  with gzip.open(store/(hashlib.sha256(variant.encode()).hexdigest()[:20]+'.csv.gz'),'rt') as f:
   tracks=[r for r in csv.DictReader(f) if r['output_type']=='DNASE' and r['biosample_name']==biosample]
  if not tracks:tracks=[None]
  g['source_contexts']=';'.join(sorted(g['source_contexts']));gg=gw[variant]
  for track in tracks:
   r={**g,'model_DNASE_track':track['track_name'] if track else '', 'model_DNASE_raw':track['raw_score'] if track else '', 'model_DNASE_quantile':track['quantile_score'] if track else '', 'candidate_gene_context':gg['candidate_gene_context'],'molecular_pip_max':gg['molecular_pip_max_across_sets'],**{t+'_p':gg[t+'_p'] for t in ['Jansen_insomnia','FinnGen_R13_insomnia','Mak_FI','Mak_FP','Ye_FFS']},'caveat':'Source records share cohorts/loci and may be duplicated; model and reporter assays measure different features; allele direction, exact sample identity and training independence need verification; no causal claim'};w.writerow(r);counts['rows']+=1;counts['no_dnase']+=track is None
print('groups',len(groups),'variants',len({x['variant_GRCh38'] for x in groups.values()}),'emvar variants',len({x['variant_GRCh38'] for x in groups.values() if x['emvar_screen']=='TRUE'}),'rows',counts['rows'],'without matched DNase',counts['no_dnase'])
