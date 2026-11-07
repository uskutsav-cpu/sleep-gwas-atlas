#!/usr/bin/env python3
"""Report direct typed evidence and paths with a single reference-LD edge."""
import argparse,csv,json,pathlib,sys
p=pathlib.Path(__file__).resolve().parent
def read(name):
 with (p/name).open() as f:return list(csv.DictReader(f,delimiter='\t'))
ap=argparse.ArgumentParser();ap.add_argument('rsid');args=ap.parse_args();nodes={x['id']:x for x in read('nodes.tsv')};edges=read('edges.tsv');start='variant:'+args.rsid
if start not in nodes:sys.exit('Variant absent from this selected graph; absence is not evidence of no association.')
def decorate(x):
 return dict(x,source_label=nodes[x['source']]['label'],target_label=nodes[x['target']]['label'])
direct=[decorate(x) for x in edges if x['source']==start or (x['target']==start and x['relation']=='reference_ld_r2_undirected')];paths=[]
for e in direct:
 if e['relation']!='reference_ld_r2_undirected':continue
 neighbor=e['target'] if e['source']==start else e['source'];
 for q in edges:
  if q['source']==neighbor and q['relation'] in ('postselected_marginal_association','positional_gene_overlap','nearest_positional_gene','gtex_v8_significant_eqtl','conditional_abf_credible_set_member'):
   paths.append({'ld_edge':e,'neighbor_evidence':decorate(q)})
print(json.dumps({'query':nodes[start],'direct_edges':direct,'one_ld_hop_paths':paths,'interpretation':'Descriptive, selected evidence only. Paths are not proof of shared causality, mediation, replication, or a direct effect on the query variant.'},indent=2))
