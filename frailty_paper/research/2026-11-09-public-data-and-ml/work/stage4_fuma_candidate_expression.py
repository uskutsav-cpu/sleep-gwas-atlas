import csv,collections,pathlib,json,math
q=list(csv.DictReader(open('outputs/stage3_initial_variant_gene_table.csv')));genes=collections.defaultdict(set)
for r in q:
 g=r['gene ID'].split('.')[0]
 if g.startswith('ENSG'):genes[g].add(r['Candidate gene'])
datasets={'GSE67835_Human_Cortex_woFetal':('brain','human'),'GSE87544_Mouse_Hypothalamus':('brain','mouse, mapped to human gene IDs'),'PBMC_10x_68k':('immune','human'),'GSE98816_Mouse_Brain_Vascular':('vascular','mouse, mapped to human gene IDs'),'MouseCellAtlas_Muscle':('muscle','mouse, mapped to human gene IDs'),'TabulaMuris_FACS_Limb_Muscle':('muscle','mouse, mapped to human gene IDs'),'GSE168408_Human_Prefrontal_Cortex_level2_Adult':('brain','human')}
out=[];coverage={}
for name,(domain,species) in datasets.items():
 p=pathlib.Path('work/stage4_magma/matrices')/(name+'.txt');seen=collections.Counter()
 with p.open() as f:
  reader=csv.DictReader(f,delimiter='\t');fields=reader.fieldnames; assert fields[0]=='GENE' and fields[-1]=='Average'; types=fields[1:-1]
  for r in reader:
   g=r['GENE']
   if g not in genes:continue
   seen[g]+=1
   for ct in types:
    v=float(r[ct]);assert math.isfinite(v) and v>=0
    out.append({'gene_id':g,'gene_symbols':'|'.join(sorted(genes[g])),'dataset':name,'domain':domain,'species':species,'cell_type':ct,'processed_expression':v,'evidence':'source processed expression; not genetic association, raw-cell fraction, donor-aware test or causality'})
 assert not [g for g,n in seen.items() if n>1],(name,seen)
 coverage[name]={'cell_types':len(types),'candidate_genes_present':len(seen),'candidate_gene_ids_missing':sorted(genes.keys()-seen.keys()),'species':species,'note':'mouse expression values mapped by source to human IDs; exact homolog mapping requires source review' if species!='human' else 'human'}
with open('outputs/stage4_fuma_candidate_gene_expression.csv','w',newline='') as f:
 w=csv.DictWriter(f,fieldnames=out[0]);w.writeheader();w.writerows(out)
pathlib.Path('outputs/stage4_fuma_candidate_expression_coverage.json').write_text(json.dumps(coverage,indent=2)+'\n')
print('candidate genes',len(genes),'rows',len(out),'positive',sum(x['processed_expression']>0 for x in out))
for k,v in coverage.items():print(k,v['cell_types'],v['candidate_genes_present'])
