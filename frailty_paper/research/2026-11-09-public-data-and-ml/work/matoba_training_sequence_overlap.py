import runpy,re,json
m=runpy.run_path('work/matoba_fixed_baseline_transfer.py');b=m['b'];trained=[]
for ident,rs in b['groups'].items():
 if len(rs)!=2 or {r['type'] for r in rs}!={'ref','alt'}:continue
 labels={r['emVar.group'] for r in rs}
 if len(labels)!=1 or next(iter(labels)) not in ('fdr <= 0.05','n.s'):continue
 cm=re.fullmatch(r'chr(\d+):rs\d+:\d+',ident)
 if not cm or not 1<=int(cm[1])<=22 or int(cm[1]) in (1,2):continue
 a,c=(next(r for r in rs if r['type']==t) for t in ('ref','alt'))
 sa=b['seq'].get(a['CRS']);sb=b['seq'].get(c['CRS'])
 if not sa or not sb or len(sa)!=199 or len(sb)!=199 or set(sa+sb)-set('ACGT'):continue
 if sum(x!=y for x,y in zip(sa,sb))!=1:continue
 trained.append((sa,sb))
if len(trained)!=801:raise ValueError('training reproduction mismatch')
comp=str.maketrans('ACGT','TGCA');train_strings={s for p in trained for s in p};matches=[]
for a,c,y in m['pairs']:
 test=[a,c,a.translate(comp)[::-1],c.translate(comp)[::-1]]
 matches.append(any(t in source for t in test for source in train_strings))
r={'training_pairs':len(trained),'test_pairs':len(matches),'test_pairs_with_150_base_exact_substring_in_either_orientation':sum(matches),'method':'each test allele and its reverse complement compared against both 199-base training sequences','warning':'Only exact full-length substring matches; near matches, homologous regions, reverse containment, donor and pretraining overlap remain unassessed'}
open('outputs/stage6_matoba_training_sequence_overlap.json','w').write(json.dumps(r,indent=2)+'\n');print(json.dumps(r,indent=2))
