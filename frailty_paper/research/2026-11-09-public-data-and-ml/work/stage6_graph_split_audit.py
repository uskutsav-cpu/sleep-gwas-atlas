import json, hashlib, collections
from pathlib import Path
root=Path('outputs'); paths=[root/'stage6_evidence_graph_nodes.jsonl',root/'stage6_evidence_graph_edges.jsonl']
nodes,edges=[[json.loads(l) for l in p.open()] for p in paths]
v={n['id']:n for n in nodes if n['type']=='variant'}
chrom=collections.Counter(n['chromosome'] for n in v.values())
# Descriptive single-linkage blocks only: one megabase is not an LD guarantee.
groups=[]
for c in sorted(chrom,key=int):
 rows=sorted((n for n in v.values() if n['chromosome']==c),key=lambda n:n['position'])
 for n in rows:
  if not groups or groups[-1]['chromosome']!=c or n['position']-groups[-1]['last_position']>1000000:
   groups.append({'chromosome':c,'last_position':n['position'],'variants':[]})
  groups[-1]['last_position']=n['position'];groups[-1]['variants'].append(n['rsid'])
target_chrom=collections.defaultdict(set); target_rels=collections.defaultdict(set)
for e in edges:
 if e['source'] in v:
  target_chrom[e['target']].add(v[e['source']]['chromosome']); target_rels[e['target']].add(e['relation'])
shared=[{'node':k,'chromosomes':sorted(s,key=int),'relations':sorted(target_rels[k])} for k,s in target_chrom.items() if len(s)>1]
result={'input_sha256':{p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in paths},'variants':len(v),'chromosome_counts':dict(sorted(chrom.items(),key=lambda kv:int(kv[0]))),'one_megabase_single_linkage_blocks':groups,'shared_cross_chromosome_targets':sorted(shared,key=lambda x:x['node']),'note':'Descriptive blocks are not LD blocks, independent labels, a cross-validation experiment, or proof against leakage.'}
(root/'stage6_graph_split_audit.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps({'variants':len(v),'chromosomes':len(chrom),'blocks':len(groups),'largest_chromosome':max(chrom.values()),'multi_variant_blocks':[g for g in groups if len(g['variants'])>1],'shared_targets':shared},indent=2))
