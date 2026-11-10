#!/usr/bin/env python3
"""Collect same-allele T-cell reporter and model scores without asserting direction."""
import argparse,csv,gzip,pathlib,hashlib
from collections import defaultdict,Counter
p=argparse.ArgumentParser();p.add_argument('--store',required=True);a=p.parse_args();root=pathlib.Path('outputs');store=pathlib.Path(a.store)
with (root/'stage5_ho2025_mpra_exact_alleles.csv').open(newline='') as f:mpra=list(csv.DictReader(f))
with (root/'stage5_variant_multitrait_gwas_matrix.csv').open(newline='') as f:gw={r['variant_GRCh38_ref_alt']:r for r in csv.DictReader(f)}
index=defaultdict(list)
for x in mpra:index[x['variant_key']].append(x)
fields=['variant_key','rsid','source_sheet','mpra_sig','mpra_log2skew','mpra_skew_logFDR','mpra_source_row','model_output','model_scorer','model_biosample','model_assay','model_track','model_track_strand','model_data_source','model_raw_alt_minus_ref','model_quantile','candidate_gene_context','molecular_pip_max','Atkins_FI_p','Jansen_insomnia_p','Mak_FI_p','Mak_FP_p','Ye_FFS_p','FinnGen_R13_insomnia_p','interpretation']
n=Counter()
with (root/'stage5_ho2025_mpra_model_comparison.csv').open('w',newline='') as out:
 w=csv.DictWriter(out,fieldnames=fields);w.writeheader()
 for key,items in index.items():
  path=store/(hashlib.sha256(key.encode()).hexdigest()[:20]+'.csv.gz')
  with gzip.open(path,'rt',newline='') as f:
   for s in csv.DictReader(f):
    if s['output_type'] not in ('ATAC','DNASE','CAGE','PROCAP'):continue
    for x in items:
     if x['source_sheet'].startswith('1.'):
      if s['biosample_name'] not in ('CD4-positive, alpha-beta T cell','naive thymus-derived CD4-positive, alpha-beta T cell'):continue
     elif 'jurkat' not in s['biosample_name'].lower():continue
     g=gw[key];row={'variant_key':key,'rsid':x['rsid'],'source_sheet':x['source_sheet'],'mpra_sig':x['mpra_sig'],'mpra_log2skew':x['Log2Skew'],'mpra_skew_logFDR':x['Skew_logFDR'],'mpra_source_row':x['source_row'],'model_output':s['output_type'],'model_scorer':s['variant_scorer'],'model_biosample':s['biosample_name'],'model_assay':s['Assay title'],'model_track':s['track_name'],'model_track_strand':s['track_strand'],'model_data_source':s['data_source'],'model_raw_alt_minus_ref':s['raw_score'],'model_quantile':s['quantile_score'],'candidate_gene_context':g['candidate_gene_context'],'molecular_pip_max':g['molecular_pip_max_across_sets'],'interpretation':('Model CAGE ontology says naive CD4, but published source accession is labeled memory conventional CD4; verify model version and source; ' if s['output_type']=='CAGE' and s['biosample_name']=='naive thymus-derived CD4-positive, alpha-beta T cell' else '')+'Reporter activity and chromatin/TSS scores are distinct assays; cells and stimulation differ; confirm reporter allele orientation and training independence; no causal or concordance claim'}
     for t in ['Atkins_FI','Jansen_insomnia','Mak_FI','Mak_FP','Ye_FFS','FinnGen_R13_insomnia']:row[t+'_p']=g[t+'_p']
     w.writerow(row);n[x['source_sheet']]+=1
print(dict(n))
