#!/usr/bin/env python3
"""Make a queue-complete variant index linking phenotype, model and experimental evidence."""
import csv,hashlib,collections
root='outputs/'
def read(n):
 with open(root+n,newline='') as f:return list(csv.DictReader(f))
q=read('stage5_score_queue.csv');matrix=read('stage5_variant_multitrait_gwas_matrix.csv');coverage=read('stage5_reporter_validation_coverage.csv');scores=read('stage5_live_score_coverage_audit.csv');allcell=read('stage5_siraj2026_allcell_reporter_model_comparison.csv')
idx=lambda rows,key:{r[key]:r for r in rows}
m=idx(matrix,'variant_GRCh38_ref_alt');c=idx(coverage,'variant_key');s=idx(scores,'variant_key');ids={r['variant_key'] for r in q}
assert len(ids)==len(q)==len(m)==len(c)==len(s)==488 and ids==set(m)==set(c)==set(s)
exp=collections.defaultdict(list)
for r in allcell:
 if r['emvar_screen']=='TRUE':exp[r['variant_GRCh38']].append(r['cell_type'])
base=list(matrix[0]);extra=['live_score_status','live_returned_output_types','live_missing_output_types','raw_score_zip_member','mpravardb_status','mpravardb_exact_records','ho_primary_tcell_status','ho_primary_tcell_category','ho_jurkat_status','ho_jurkat_category','siraj_grch37_mapping_status','siraj_sknsh_assay_rows','siraj_sknsh_candidate_rows','siraj_selected_any_emvar','siraj_selected_emvar_cell_lines','experimental_interpretation']
with open(root+'stage5_variant_phenotype_experimental_index.csv','w',newline='') as f:
 w=csv.DictWriter(f,fieldnames=base+extra);w.writeheader()
 for qrow in q:
  k=qrow['variant_key'];v=dict(m[k]);cc=c[k];ss=s[k];v.update(live_score_status=ss['status'],live_returned_output_types=ss['returned_output_types'],live_missing_output_types=ss['missing_output_types'],raw_score_zip_member=hashlib.sha256(k.encode()).hexdigest()[:20]+'.csv.gz',mpravardb_status=cc['match_status'],mpravardb_exact_records=cc['exact_allele_records'],ho_primary_tcell_status=cc['primary_tcell_ho_status'],ho_primary_tcell_category=cc['primary_tcell_ho_category'],ho_jurkat_status=cc['jurkat_ho_status'],ho_jurkat_category=cc['jurkat_ho_category'],siraj_grch37_mapping_status=cc['siraj_sknsh_mapping_status'],siraj_sknsh_assay_rows=cc['siraj_sknsh_assay_rows'],siraj_sknsh_candidate_rows=cc['siraj_sknsh_candidate_rows'],siraj_selected_any_emvar=cc['siraj_selected_emvar_screen'],siraj_selected_emvar_cell_lines='|'.join(sorted(set(exp[k]))),experimental_interpretation='Matched reporter effects do not establish target gene, shared disease causality, direct assay concordance, or model-training independence; absent record does not prove zero effect');w.writerow(v)
print('integrated variants',len(q),'selected emvar variants',sum(bool(x) for x in exp.values()))
