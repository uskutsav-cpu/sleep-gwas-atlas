#!/usr/bin/env python3
"""Build source-tracked MEIS1 evidence as an explicitly separate graph extension."""
import hashlib,json
from pathlib import Path
p=Path('outputs')
base_nodes={json.loads(x)['id'] for x in (p/'stage6_evidence_graph_nodes.jsonl').read_text().splitlines()}
base_edges=[json.loads(x) for x in (p/'stage6_evidence_graph_edges.jsonl').read_text().splitlines()]
proxy='variant:GRCh37:2:66660403:T:G'; sentinel='variant:GRCh37:2:66677816:C:G'
assert proxy not in base_nodes and sentinel not in base_nodes and 'gene_symbol:MEIS1' in base_nodes
nodes=[{'id':proxy,'type':'variant','rsid':'rs13033745','scope':'external proxy, not a frozen score variant'},{'id':sentinel,'type':'variant','rsid':'rs1519102','scope':'external sentinel, not a frozen score variant'},{'id':'assay:Palermo2023:neural_progenitor:MEIS1','type':'molecular_assay','species':'human','context':'iPSC-derived neural progenitors','url':'https://pmc.ncbi.nlm.nih.gov/articles/PMC9821868/'},{'id':'assay:Zimmerman2026:meis1b_zebrafish','type':'gene_perturbation_assay','species':'zebrafish','gene':'meis1b','url':'https://pubmed.ncbi.nlm.nih.gov/42532834/'}]
edges=[]
def add(source,target,relation,**kw): edges.append(dict(source=source,target=target,relation=relation,**kw))
add(proxy,'assay:Palermo2023:neural_progenitor:MEIS1','variant_in_contacting_open_region',caveat='Region-level chromatin contact; not an allele edit or proof of causal variant.')
add('assay:Palermo2023:neural_progenitor:MEIS1','gene_symbol:MEIS1','promoter_contact',caveat='Physical contact in human neural progenitors; not proof of functional or disease mediation.')
add('assay:Zimmerman2026:meis1b_zebrafish','gene_symbol:MEIS1','ortholog_gene_perturbation_context',caveat='Zebrafish gene perturbation; not human variant perturbation or frailty evidence.')
ld=p/'stage6_meis1_npc_proxy_reference_ld.json';data=json.loads(ld.read_text()); sha=hashlib.sha256(ld.read_bytes()).hexdigest()
lookup={'rs13033745':proxy,'rs1519102':sentinel,'rs113851554':'variant:GRCh37:2:66750564:T:G','rs78206187':'variant:GRCh37:2:67022234:G:A'}
for row in data['pairs']:
 a,b=row['variants']
 if row['status']=='available' and proxy in (lookup[a],lookup[b]):
  add(lookup[a],lookup[b],'reference_haplotype_ld',r2=row['r2'],n_individuals=row['n_individuals'],source_sha256=sha,source_file=ld.name,caveat='EUR reference LD only; not GWAS association or causal variant.')
for tag,target in [('Atkins_FI','study:FI_meta'),('Ye_FFS','study:Ye_FFS'),('Jansen_insomnia','study:Jansen_UKB_insomnia'),('MVP_insomnia','study:MVP_insomnia'),('FinnGen_insomnia','study:FinnGen_insomnia'),('FinnGen_HFRS','study:FinnGen_HFRS')]:
 f=p/f'stage6_meis1_proxy_{tag}_lookup.json';j=json.loads(f.read_text());h=hashlib.sha256(f.read_bytes()).hexdigest()
 for row in j['hits']:
  rs=row.get('variant_id',row.get('SNP',row.get('rsid',row.get('rs_id',row.get('rsids')))))
  if rs not in lookup:continue
  if tag=='Atkins_FI': allele=row['effect_allele'];other=row['other_allele'];value=float(row['beta']);scale='beta';pv=float(row['p_value'])
  elif tag=='Ye_FFS':allele=row['A1'];other=row['A2'];value=float(row['BETA']);scale='beta';pv=float(row['P'])
  elif tag=='Jansen_insomnia':allele=row['A1'];other=row['A2'];value=float(row['OR']);scale='odds_ratio';pv=float(row['P'])
  elif tag=='MVP_insomnia':allele=row['effect_allele'];other=row['other_allele'];value=float(row['odds_ratio']);scale='odds_ratio';pv=float(row['p_value'])
  elif tag=='FinnGen_insomnia':allele=row['alt'];other=row['ref'];value=float(row['beta']);scale='beta';pv=float(row['pval'])
  else:allele=row['effect_allele'];other=row['other_allele'];value=float(row['beta']);scale='beta';pv=float(row['p_value'])
  add(lookup[rs],target,'gwas_marginal_association',effect_allele=allele,other_allele=other,effect=value,effect_scale=scale,p_value=pv,input_file=f.name,input_sha256=h,caveat='Reported source orientation; postselected; study overlap and multiallelic source coding require review; not score validation.')
ids=base_nodes|{n['id'] for n in nodes};assert len(ids)==len(base_nodes)+len(nodes)
assert all(e['source'] in ids and e['target'] in ids for e in edges)
assert all((e.get('p_value') is None or 0<=e['p_value']<=1) for e in edges)
for name,vals in [('nodes',nodes),('edges',edges)]:
 (p/f'stage6_meis1_graph_extension_{name}.jsonl').write_text(''.join(json.dumps(x,sort_keys=True)+'\n' for x in vals))
print(json.dumps({'new_nodes':len(nodes),'new_edges':len(edges),'base_nodes':len(base_nodes),'base_edges':len(base_edges),'note':'separate extension; candidate score membership unchanged'},indent=2))
