import csv,hashlib,json,re
from pathlib import Path
p=Path(__file__).resolve().parent;parent=p.parent;atlas=next(x for x in p.parents if x.name=='shared-genetics-atlas');provisional=parent.parent/'provisional_allelic_42.tsv';published=parent/'candidate_comparison.tsv';cross=atlas/'biological_mechanisms/cross_assay_prediction/lee_27_crosswalk.tsv';pipeline=Path('/Volumes/Extreme SSD/brain6-crossdisorder-mpra/authors_pipeline.R')
def read(f):
 with f.open(newline='') as h:return list(csv.DictReader(h,delimiter='\t'))
a={x['rsid']:x for x in read(provisional)}; pred={x['rsid']:x for x in read(cross)};rows=[]
for x in read(published):
 if x['allelic_status']!='unique_match':continue
 y=a[x['rsid']];z=pred[x['rsid']]; hg19=re.search(r'_([ACGT]+):([ACGT]+)$',y['source_key']).groups();hg38=re.search(r':([ACGT]+)>([ACGT]+)$',y['grch38']).groups();assert hg19==hg38
 pub=float(x['published_mpra_logfc_unoriented']);reproduced=float(y['v124_logFC']);assert abs(pub-reproduced)<1e-10
 rows.append({'rsid':x['rsid'],'variant38':y['grch38'],'published_logfc':pub,'reconstructed_logfc':reproduced,'absolute_difference':abs(pub-reproduced),'allelic_fdr':x['allelic_fdr'],'element_class':x['element_class'],'dnase_501_alt_vs_ref':z['dnase_501']})
assert len(rows)==27 and len({x['rsid'] for x in rows})==27
s=pipeline.read_text();required=['alt <- c(rep(FALSE,repnum),rep(TRUE,repnum))','design_matrix <- model.matrix(~alt+batch)','colnames(dna_count_matrix) <- c(paste0("Replicate ", 1:repnum, " REF"), paste0("Replicate ", 1:repnum, " ALT"))','mpra_result <- topTable(mpra_lm, coef = 2, number = Inf)'];assert all(x in s for x in required)
with (p/'comparison.tsv').open('w',newline='') as h:w=csv.DictWriter(h,fieldnames=rows[0],delimiter='\t');w.writeheader();w.writerows(rows)
(p/'receipt.json').write_text(json.dumps({'post_hoc':True,'n_matched':len(rows),'maximum_absolute_logfc_difference':max(x['absolute_difference'] for x in rows),'inputs_sha256':{str(x):hashlib.sha256(x.read_bytes()).hexdigest() for x in (provisional,published,cross,pipeline)},'caveat':'Assay labels and pipeline correspondence are supported but original construct and barcode provenance is not independently resolved.'},indent=2)+'\n');print(json.dumps({'matched':len(rows),'max_difference':max(x['absolute_difference'] for x in rows),'examples':[x for x in rows if x['rsid'] in ('rs325485','rs6874138','rs62369151')]},indent=2))
