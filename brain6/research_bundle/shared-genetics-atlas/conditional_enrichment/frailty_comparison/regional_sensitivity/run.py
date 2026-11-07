import csv,collections,hashlib,json,pathlib
p=pathlib.Path(__file__).resolve().parent;source=p.parent/'block_counts.tsv'
with source.open() as f:rows=[r for r in csv.DictReader(f,delimiter='\t') if float(r['block_bp'])==5000000]
assert len(rows)==562 and len({r['block'] for r in rows})==len(rows)
keys=['n','n_insomnia','FI_cond','FI_joint','FFS_cond','FFS_joint','HFRS_cond','HFRS_joint']; totals={k:sum(int(r[k]) for r in rows) for k in keys};assert totals['n']==1137927
out=[]
for label,exclude in [('none',lambda b:False)]+[(f'chr{c}',lambda b,c=c:b.split(':')[0]==str(c)) for c in range(1,23)]+[('chr19_45_50Mb',lambda b:b=='19:9')]:
 removed={k:sum(int(r[k]) for r in rows if exclude(r['block'])) for k in keys};a={k:totals[k]-removed[k] for k in keys};assert a['n']>0 and all(a[k+'_cond']>0 for k in ['FI','FFS','HFRS'])
 result={'excluded':label,'n_remaining':a['n'],'n_removed':removed['n'],'baseline_insomnia_fraction':a['n_insomnia']/a['n']}
 for trait in ['FI','FFS','HFRS']:result.update({trait+'_selected':a[trait+'_cond'],trait+'_joint':a[trait+'_joint'],trait+'_fraction':a[trait+'_joint']/a[trait+'_cond'],trait+'_removed_selected':removed[trait+'_cond'],trait+'_removed_joint':removed[trait+'_joint']})
 out.append(result)
with (p/'exclusion_results.tsv').open('w',newline='') as f:w=csv.DictWriter(f,list(out[0]),delimiter='\t');w.writeheader();w.writerows(out)
ranked=sorted(rows,key=lambda r:int(r['HFRS_joint']),reverse=True)
with (p/'hfrs_block_counts_ranked.tsv').open('w',newline='') as f:w=csv.DictWriter(f,list(ranked[0]),delimiter='\t');w.writeheader();w.writerows(ranked)
receipt={'source_sha256':hashlib.sha256(source.read_bytes()).hexdigest(),'source_blocks':len(rows),'total_support':totals['n'],'chr19_window_definition':'1-based inclusive 45000001–50000000 GRCh37','script_sha256':hashlib.sha256((p/'run.py').read_bytes()).hexdigest(),'caveat':'post-hoc physical exclusions; no genetic or causal inference'}
(p/'receipt.json').write_text(json.dumps(receipt,indent=2)+'\n')
for x in out:
 if x['excluded'] in ('none','chr19','chr19_45_50Mb'):print(x)
print('top HFRS joint',[(r['block'],r['HFRS_joint'],r['HFRS_cond']) for r in ranked[:5]])
