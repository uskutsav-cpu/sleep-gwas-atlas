#!/usr/bin/env python3
"""Join neuronal-line reporter allele measurements to matched-label model tracks."""
import argparse,csv,gzip,pathlib,hashlib,collections
p=argparse.ArgumentParser();p.add_argument('--store',required=True);a=p.parse_args();root=pathlib.Path('outputs');store=pathlib.Path(a.store)
with (root/'stage5_siraj2026_neural_mpra_exact_alleles.csv').open() as f:measure=list(csv.DictReader(f))
with (root/'stage5_variant_multitrait_gwas_matrix.csv').open() as f:gw={r['variant_GRCh38_ref_alt']:r for r in csv.DictReader(f)}
index=collections.defaultdict(list)
for m in measure:index[m['variant_GRCh38']].append(m)
fields=['variant_GRCh38','variant_GRCh37_source','source_file','source_id','A_log2FC','A_logPadj_BF','B_log2FC','B_logPadj_BF','Log2Skew','Skew_logFDR_act','screen_active_and_fdr_lt_0_1','model_output','model_biosample','model_track','model_raw_alt_minus_ref','model_quantile','candidate_gene_context','molecular_pip_max','Jansen_insomnia_p','FinnGen_R13_insomnia_p','Mak_FI_p','Mak_FP_p','Ye_FFS_p','note'];count=collections.Counter()
with (root/'stage5_siraj2026_neural_mpra_model_comparison.csv').open('w',newline='') as out:
 w=csv.DictWriter(out,fieldnames=fields);w.writeheader()
 for key,records in index.items():
  path=store/(hashlib.sha256(key.encode()).hexdigest()[:20]+'.csv.gz')
  with gzip.open(path,'rt',newline='') as f:
   scores=[s for s in csv.DictReader(f) if s['output_type']=='DNASE' and s['biosample_name']=='SK-N-SH']
  for m in records:
   try:
    active=any(abs(float(m[a]))>=1 and float(m[b])>=2 for a,b in [('A_log2FC','A_logPadj_BF'),('B_log2FC','B_logPadj_BF')]);candidate='yes' if active and float(m['Skew_logFDR_act'])>=1 else 'no'
   except (ValueError,TypeError):candidate='not assessable'
   for s in scores:
    g=gw[key];w.writerow({**{k:m[k] for k in ['variant_GRCh38','variant_GRCh37_source','source_file','source_id','A_log2FC','A_logPadj_BF','B_log2FC','B_logPadj_BF','Log2Skew','Skew_logFDR_act']},'screen_active_and_fdr_lt_0_1':candidate,'model_output':s['output_type'],'model_biosample':s['biosample_name'],'model_track':s['track_name'],'model_raw_alt_minus_ref':s['raw_score'],'model_quantile':s['quantile_score'],'candidate_gene_context':g['candidate_gene_context'],'molecular_pip_max':g['molecular_pip_max_across_sets'],**{t+'_p':g[t+'_p'] for t in ['Jansen_insomnia','FinnGen_R13_insomnia','Mak_FI','Mak_FP','Ye_FFS']},'note':'Post hoc source-row screen, not official best-library calls; direction and training independence unresolved; reporter activity and DNase are different measurements; not disease causality'})
    count[candidate]+=1
print('matched score rows',dict(count),'variants',len(index))
