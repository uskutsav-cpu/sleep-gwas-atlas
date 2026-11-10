from docx import Document
import csv,hashlib,json,collections
p='work/stage3_frailty_twas_supplement.docx';doc=Document(p);assert len(doc.tables)==3
cand=collections.defaultdict(set)
for x in csv.DictReader(open('outputs/stage3_initial_variant_gene_table.csv')):
 if x['Candidate gene']:cand[x['Candidate gene']].add(x['Variant'])
result=[]
for ti,method in [(0,'UTMOST cross-tissue 44 GTEx tissues'),(1,'FUSION GTEx whole blood')]:
 t=doc.tables[ti];header=[x.text.strip() for x in t.rows[0].cells];assert len(set(header))==len(header)
 for row in t.rows[1:]:
  v=[c.text.strip() for c in row.cells];d=dict(zip(header,v));symbol=d['symble' if ti==0 else 'ID'];result.append({'method':method,'gene_symbol':symbol,'p':d['p_value' if ti==0 else 'TWAS.P'],'reported_FDR':d['PFDR' if ti==0 else 'FDR'],'fusion_twas_z':d.get('TWAS.Z',''),'chr':d.get('CHR',''),'best_gwas_id':d.get('BEST.GWAS.ID',''),'n_snps':d.get('NSNP',''),'project_candidate_markers':'|'.join(sorted(cand.get(symbol,[]))),'note':'Published significant-gene table; gene-model association, not marker-specific causation; FI and insomnia GWAS may overlap'})
with open('outputs/stage3_published_fi_twas_supplement_genes.csv','w',newline='') as f:
 w=csv.DictWriter(f,fieldnames=result[0]);w.writeheader();w.writerows(result)
matched=[x for x in result if x['project_candidate_markers']]
json.dump({'source_paper':'https://doi.org/10.3389/fgene.2024.1404456','supplement_url':'https://www.ebi.ac.uk/europepmc/webservices/rest/PMC11272547/supplementaryFiles','supplement_filename':'Table1.DOCX','supplement_sha256':hashlib.sha256(open(p,'rb').read()).hexdigest(),'cross_tissue_rows':sum(x['method'].startswith('UTMOST') for x in result),'whole_blood_rows':sum(x['method'].startswith('FUSION') for x in result),'rows_matching_candidate_symbol':len(matched),'distinct_matching_symbols':len({x['gene_symbol'] for x in matched}),'selection':'source significant supplementary tables, not all tested genes; candidate symbol matching without independent identifier harmonization'},open('outputs/stage3_published_fi_twas_supplement_receipt.json','w'),indent=2)
print('matched',len(matched),len({x['gene_symbol'] for x in matched}));
for x in matched:print(x['method'],x['gene_symbol'],x['p'],x['reported_FDR'],x['project_candidate_markers'])
