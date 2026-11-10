import csv
# Provisional manually transcribed assay-to-gene labels; GDF11/8 maps to both.
panel={'HTRA1','ADAMTSL2','COL28A1','IL1RN','FABP4','LEP','GDF15','TNNT2','FSTL3','ATOX1','RNASE1','FAP','APLP1','WFIKKN2','ANTXR2','LILRA4','VAT1','ENPP5','CACNA2D3','NELL1','GLRX3','GDF11','GDF8','HPGDS','EGFR','NDST1'}
with open('outputs/stage5_variant_gene_tissue_phenotype_integrated.csv', newline='') as f:
 genes={r['gene_symbol'] for r in csv.DictReader(f)}
print('panel_symbols',len(panel),'inventory_gene_symbols',len(genes),'exact_overlap',sorted(panel & genes))
