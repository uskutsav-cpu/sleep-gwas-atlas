import pathlib,csv,json,collections
p=pathlib.Path('outputs/shared-genetics-atlas/evidence_graph')
def read(n):
 with (p/n).open() as f:return list(csv.DictReader(f,delimiter='\t'))
n=read('nodes.tsv');e=read('edges.tsv');receipt=json.loads((p/'receipt.json').read_text());ids={x['id'] for x in n};assert len(ids)==len(n) and len(e)==1167 and all(x['source'] in ids and x['target'] in ids for x in e);assert len({tuple(x[k] for k in ('source','target','relation','detail')) for x in e})==len(e)
assert dict(collections.Counter(x['relation'] for x in e))==receipt['counts_edges_by_type'];v=[x for x in e if x['source']=='variant:rs62264764'];assert any(x['relation']=='positional_gene_overlap' and x['target']=='gene:ENSG00000185565' for x in v);assert any(x['relation']=='gtex_v8_significant_eqtl' and x['target']=='gene:ENSG00000185565' and 'tissue=Adrenal_Gland' in x['detail'] for x in v);assert not any('Brain_Cortex' in x['detail'] and x['relation']=='gtex_v8_significant_eqtl' for x in v)
assert any(x['source']=='variant:rs2155292' and x['target']=='trait:FFS_Ye' and x['value']=='0.54' and x['relation']=='postselected_marginal_association' for x in e)
assert any(x['source']=='variant:rs12536335' and x['target']=='trait:HFRS_FinnGen' for x in e)
assert sum(x['relation']=='ldsc_genetic_correlation' for x in e)==100
assert any(x['source']=='trait:insomnia' and x['target']=='trait:adhd' and x['relation']=='ldsc_genetic_correlation' and abs(float(x['value'])-.4033)<1e-9 for x in e)
assert any(x['source']=='trait:FI_Atkins' and x['target']=='trait:FFS_Ye' and x['relation']=='ldsc_genetic_correlation' and abs(float(x['value'])-.753411068344601)<1e-12 for x in e)
print('PASS: references, uniqueness, 100 genome-wide edges and source-specific variant evidence')
