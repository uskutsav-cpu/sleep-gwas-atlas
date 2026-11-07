#!/usr/bin/env python3
"""Bind Brain6's exploratory GTEx lead lookups to the six-region AVI pilot."""
import csv,hashlib,json
from pathlib import Path
ROOT=Path(__file__).resolve().parents[2];P=ROOT/'brain6/alphagenome_pilot/results';O=P/'avi_mirror';SRC=ROOT/'brain6/results/exploratory_five_track_v4/lead_eqtl_context';FIN=ROOT/'brain6/paper/final_package_v1'; FUNC=ROOT/'brain6/results/brain6_exploratory_functional_v1'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def read(p):
    with p.open(newline='') as f:return list(csv.DictReader(f,delimiter='\t'))
def write(p,rows):
    with p.open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]),delimiter='\t',lineterminator='\n');w.writeheader();w.writerows(rows)
prov=json.loads((SRC/'provenance.json').read_text())
for name in ['lead_eqtl_associations.tsv','lead_eqtl_query_summary.tsv']:
    if sha(SRC/name)!=prov['output_sha256'][name]:raise ValueError('Upstream eQTL source hash mismatch')
finalprov=json.loads((FIN/'BRAIN6_FINAL_PROVENANCE.json').read_text())
if sha(FIN/'BRAIN6_FINAL_CANDIDATES.tsv')!=finalprov['output_sha256']['BRAIN6_FINAL_CANDIDATES.tsv']:raise ValueError('Final candidate source mismatch')
regions=read(P/'regions.tsv'); ids={r['candidate_locus_id'] for r in regions}
if len(ids)!=6:raise ValueError('Expected six candidates')
avi=read(O/'variant_evidence_avi.tsv'); amap={(r['candidate_locus_id'],r['snp']):r for r in avi}
summary=[r for r in read(SRC/'lead_eqtl_query_summary.tsv') if r['locus_id'] in ids and r['pair_id']=='insomnia__adhd']
if len(summary)!=12 or len({(r['locus_id'],r['dataset_id']) for r in summary})!=12:raise ValueError('Expected 6 x 2 lead-QTL query rows')
lookup={(r['locus_id'],r['lead_snp'],r['dataset_id']):r for r in summary}
qtl=[r for r in read(SRC/'lead_eqtl_associations.tsv') if r['locus_id'] in ids and r['pair_id']=='insomnia__adhd']
for row in qtl:
    key=(row['locus_id'],row['lead_snp'],row['dataset_id']);query=lookup[key]
    if row['allele_status']!='EXACT_UNORDERED' or (row['chr38'],row['bp38'])!=(query['chr38'],query['bp38']):raise ValueError('Unexpected eQTL source identity')
    ar=amap.get((row['locus_id'],row['lead_snp']))
    if ar:
        if ar['chromosome']!='chr'+row['chr38'] or ar['pos_grch38_1based']!=row['bp38'] or {ar['ref_grch38'],ar['alt_grch38']}!={row['qtl_ref'],row['qtl_alt']}:raise ValueError('AVI/eQTL exact allele match failed')
        row.update(avi_match_status='EXACT_LOCUS_RSID_POSITION_UNORDERED_ALLELES',avi_raw_e5=ar['avi_raw_e5'],avi_phred_lower=ar['avi_phred_lower'],avi_phred_upper=ar['avi_phred_upper'],atlas_url=ar['atlas_url'])
    else:row.update(avi_match_status='LEAD_NOT_IN_SCORED_CREDIBLE_SET_UNION',avi_raw_e5='',avi_phred_lower='',avi_phred_upper='',atlas_url='')
for row in summary:
    ar=amap.get((row['locus_id'],row['lead_snp']))
    row['avi_match_status']='LEAD_IN_SCORED_CREDIBLE_SET_UNION' if ar else 'LEAD_NOT_IN_SCORED_CREDIBLE_SET_UNION'
    matches=[x for x in qtl if x['locus_id']==row['locus_id'] and x['lead_snp']==row['lead_snp'] and x['dataset_id']==row['dataset_id']]
    if len(matches)!=int(row['exact_rsid_allele_rows']):raise ValueError('QTL query count mismatch')
if not qtl:raise ValueError('Expected molecular records')
write(O/'lead_eqtl_context.tsv',qtl);write(O/'lead_eqtl_query_status.tsv',summary)
cid='placo_insomnia__adhd_chr5_103447968_104447968';case=[r for r in qtl if r['locus_id']==cid and r['lead_snp']=='rs2431108']
if len(case)!=2 or {r['gene_id'] for r in case}!={'ENSG00000283462'}:raise ValueError('Unexpected chr5 case records')
final={r['candidate_locus_id']:r for r in read(FIN/'BRAIN6_FINAL_CANDIDATES.tsv')}; f=final[cid]
if f['eqtl_coloc_status']!='NOT_ESTIMATED_NO_ELIGIBLE_EQTL_CONTEXT' or f['pair_level_gwas_qtl_coloc_status']!='NOT_ESTIMATED':raise ValueError('Unexpected chr5 final QTL status')
lines=['# Brain6 × AlphaGenome: measured brain eQTL context','','The tables join the six protected candidate IDs to the existing exploratory GTEx v8 exact-lead eQTL lookups. These are observed molecular associations from two bulk brain datasets; they are **not AlphaGenome tissue-specific predictions**. A lead is joined to AVI only when its candidate ID, rsID, GRCh38 coordinate and unordered alleles match. Leads outside the scored credible-set union are retained with an explicit missing AVI status.','','## Chromosome 5: rs2431108','','The existing Brain6 lookups contain two nominal associations for the same molecular trait, ENSG00000283462, at chr5:104612267:T>C (GRCh38). Both have the C allele as the QTL effect allele.','', '| GTEx tissue | Metadata sample size | Nominal eQTL p-value | Beta for C | Standard error |','|---|---:|---:|---:|---:|']
for r in case:lines.append(f"| {r['sample_group']} | {r['sample_size_metadata']} | {r['nominal_pvalue']} | {r['beta_alt']} | {r['se']} |")
lines += ['', 'The signs differ across these two estimates, and their uncertainty is substantial. Neither result establishes a significant molecular effect after multiple-testing correction, a tissue-specific effect, or a causal gene. A non-significant nominal result also does not exclude an effect. These are lead-only observations; they do not substitute for a full-region analysis.','', 'Brain6’s final component-GWAS–eQTL status for this candidate is `NOT_ESTIMATED_NO_ELIGIBLE_EQTL_CONTEXT`, and pair-level GWAS–QTL colocalization is `NOT_ESTIMATED`. The exact-lead rows did not satisfy the separate eligibility gates for regional colocalization. AVI scores, a positional gene and nominal eQTL associations cannot override those gates or the failed primary LAVA family.','', '## Scope and files','',f"The two original query tables have 12 lead-by-tissue status rows and {len(qtl)} exact-lead gene-association rows across the six regions. `lead_eqtl_context.tsv` contains all matching gene rows, while `lead_eqtl_query_status.tsv` keeps the complete query statuses, including no-row results. Sample sizes are source metadata and are not necessarily per-variant analyzed N.",'', 'The eQTL lookup was exploratory and did not apply a gene/tissue multiple-testing correction. A missing lead-level row does not mean no regional eQTL exists. QTL beta is expressed for the listed QTL ALT allele and must not be compared with GWAS beta without checking orientation. Existing model limitations and lack of independent two-trait locus replication are unchanged.','', '[Source project eQTL report](https://github.com/uskutsav-cpu/sleep-gwas-atlas/blob/303c905c229104b7bf623df6de7431abb87e0727/brain6/results/exploratory_five_track_v4/lead_eqtl_context/report.md) · [Brain6 final results](https://github.com/uskutsav-cpu/sleep-gwas-atlas/blob/303c905c229104b7bf623df6de7431abb87e0727/brain6/paper/final_package_v1/BRAIN6_FINAL_RESULTS.md) · [eQTL Catalogue](https://www.ebi.ac.uk/eqtl/)', '']
(O/'MOLECULAR_CONTEXT.md').write_text('\n'.join(lines))
receipt={'upstream_qtl_provenance_sha256':sha(SRC/'provenance.json'),'upstream_qtl_outputs':{n:sha(SRC/n) for n in ['lead_eqtl_associations.tsv','lead_eqtl_query_summary.tsv']},'final_provenance_sha256':sha(FIN/'BRAIN6_FINAL_PROVENANCE.json'),'final_candidate_sha256':sha(FIN/'BRAIN6_FINAL_CANDIDATES.tsv'),'avi_evidence_sha256':sha(O/'variant_evidence_avi.tsv'),'region_sha256':sha(P/'regions.tsv'),'outputs_sha256':{n:sha(O/n) for n in ['lead_eqtl_context.tsv','lead_eqtl_query_status.tsv','MOLECULAR_CONTEXT.md']},'n_candidates':len(ids),'n_status_rows':len(summary),'n_association_rows':len(qtl),'case_status':'EXPLORATORY_LEAD_LEVEL_ONLY_NO_ELIGIBLE_GWAS_QTL_COLOCALIZATION'}
(O/'molecular_context_provenance.json').write_text(json.dumps(receipt,indent=2,sort_keys=True)+'\n')
print(f"{len(ids)} regions, {len(summary)} query states, {len(qtl)} exact-lead gene associations")
