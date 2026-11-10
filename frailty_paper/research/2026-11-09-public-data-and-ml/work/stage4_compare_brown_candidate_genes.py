import csv,collections
m={}
for x in csv.DictReader(open('outputs/stage3_initial_variant_gene_table.csv')):
 g=x['gene ID'].split('.')[0]
 if g.startswith('ENSG'):
  z=m.setdefault(g,{'markers':set(),'symbols':set()});z['markers'].add(x['Variant']);z['symbols'].add(x['Candidate gene'])
traits=['atkins_fi_brown','mak_fi_brown','ye_ffs_brown','jansen_insomnia_brown'];data={t:{x['GENE']:x for x in csv.DictReader(open('outputs/stage4_'+t+'_magma_gene_results.csv'))} for t in traits}
rows=[]
for g,v in sorted(m.items()):
 row={'gene_id':g,'candidate_symbols':'|'.join(sorted(v['symbols'])),'candidate_markers':'|'.join(sorted(v['markers']))}
 for t in traits:
  x=data[t].get(g,{});row[t+'_p']=x.get('P','');row[t+'_BH_q']=x.get('bh_q_all_genes','')
 row['interpretation']='separate approximate gene tests; no joint test, shared causal variant, marker causality, or mediation inferred';rows.append(row)
with open('outputs/stage4_candidate_brown_four_trait_gene_comparison.csv','w',newline='') as f:w=csv.DictWriter(f,fieldnames=rows[0]);w.writeheader();w.writerows(rows)
for t in traits: print(t,'tested',sum(bool(x[t+'_p']) for x in rows),'q<.05',sum(bool(x[t+'_p']) and float(x[t+'_BH_q'])<.05 for x in rows))
trip=['atkins_fi_brown','ye_ffs_brown','jansen_insomnia_brown'];common=[x for x in rows if all(x[t+'_p'] for t in trip)];hits=[x for x in common if all(float(x[t+'_BH_q'])<.05 for t in trip)];print('three tested',len(common),'three significant',len(hits))
for x in hits:print(x['gene_id'],x['candidate_symbols'],x['candidate_markers'],*[x[t+'_BH_q'] for t in trip])
