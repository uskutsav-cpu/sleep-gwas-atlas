from openpyxl import load_workbook
from pathlib import Path
import csv,hashlib,json
p=Path('work/ye_supplement.xlsx'); w=load_workbook(p,read_only=True,data_only=True)
def table(s):
 vals=iter(w[s].iter_rows(values_only=True)); header=next(vals);return [dict(zip(header,x)) for x in vals if x[0] is not None]
loci=table(' ST1. COJO independent loci'); links=table('ST3. GWAS Catalog'); assert len(loci)==37
lmap={x['Index']:x for x in loci}; assert len(lmap)==37
url='https://media.springernature.com/original/springer-static/esm/art%3A10.1007%2Fs11357-023-00771-z/MediaObjects/11357_2023_771_MOESM2_ESM.xlsx'
results=[]
for x in loci:
 matches=[y for y in links if y['GenomicLocus']==x['Index'] and any(k in str(y['Trait']).lower() for k in ('insomnia','sleep','chronotype'))]
 results.append({'source_locus_id':x['Index'],'ffs_lead_variant':x['SNP'],'chromosome':x['Chr'],'position_reported':x['bp'],'source_refA_not_yet_author_confirmed_effect_allele':x['refA'],'ffs_marginal_beta_reported':x['b'],'ffs_marginal_se':x['se'],'ffs_marginal_p':x['p'],'ffs_conditional_beta':x['bJ'],'ffs_conditional_se':x['bJ_se'],'ffs_conditional_p':x['pJ'],'ffs_n_reported':x['n'],'nearest_gene_source_annotation_only':x['Nearest Gene'],'annotation':x['Annotation'],'novel_marker_source':x['Novel'] or '', 'prior_sleep_catalog_rows_in_locus':len(matches),'prior_insomnia_catalog_rows_in_locus':sum('insomnia' in str(y['Trait']).lower() for y in matches),'interpretation':'FFS-associated only; prior nearby sleep associations do not establish shared causal variant; not score weights','source_sheet':'ST1. COJO independent loci','source_url':url})
with open('outputs/published_ffs_37_locus_inventory.csv','w',newline='') as f:
 wr=csv.DictWriter(f,fieldnames=results[0]);wr.writeheader();wr.writerows(results)
result2=[]
for y in links:
 if not any(k in str(y['Trait']).lower() for k in ('insomnia','sleep','chronotype')):continue
 x=lmap[y['GenomicLocus']]; distance=abs(int(y['bp'])-int(x['bp']));
 result2.append({'ffs_locus_id':x['Index'],'ffs_lead_variant':x['SNP'],'ffs_lead_position_reported':x['bp'],'ffs_nearest_gene_annotated':x['Nearest Gene'],'sleep_catalog_variant':y['SNP'],'sleep_position_reported':y['bp'],'distance_bp':distance,'sleep_trait':y['Trait'],'reported_sleep_effect_allele':y['Strongest'],'sleep_p_reported':y['P'],'sleep_effect_reported':y['OrBeta'],'reported_gene':y['ReportedGene'],'mapped_gene':y['MappedGene'],'sleep_pmid':y['PMID'],'interpretation':'Same FFS locus grouping in source; variant LD, shared causality, direction and independent replication not established','source_url':url})
assert len(result2)==13 and len({r['ffs_locus_id'] for r in result2})==8
with open('outputs/published_ffs_loci_sleep_catalog_links.csv','w',newline='') as f:
 wr=csv.DictWriter(f,fieldnames=result2[0]);wr.writeheader();wr.writerows(result2)
print(json.dumps({'sha256':hashlib.sha256(p.read_bytes()).hexdigest(),'FFS_loci':len(results),'sleep_catalog_records':len(result2),'loci_with_catalog_sleep':len(set(r['ffs_locus_id'] for r in result2)),'single_insomnia_row':[r for r in result2 if 'insomnia' in r['sleep_trait'].lower()],'marginal_p_above_5e_8':[r['ffs_lead_variant'] for r in results if r['ffs_marginal_p']>=5e-8]},indent=2))
