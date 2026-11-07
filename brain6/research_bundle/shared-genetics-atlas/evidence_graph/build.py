#!/usr/bin/env python3
import csv,pathlib,json,hashlib,collections,re
root=pathlib.Path('outputs');out=root/'shared-genetics-atlas/evidence_graph';paths={'variants':root/'brain6-alphagenome-pilot/results/variant_evidence.tsv','context':root/'shared-genetics-atlas/variant_context/selected_variant_gene_context.tsv','gtex':root/'shared-genetics-atlas/variant_context/gtex_v8/significant_associations.tsv','ld':root/'shared-genetics-atlas/cross_project_ld/pairwise_ld.tsv','frailty_candidates':root/'shared-genetics-atlas/locus_atlas/frailty_six_candidate_atlas.tsv','frailty_joint':root/'shared-genetics-atlas/evidence/song124_mvp_hfrs_joint_screen.csv','rg_matrix':root/'shared-genetics-atlas/multiscale_embedding/source_matrix.tsv','hfrs_rg':root/'shared-genetics-atlas/evidence/NEW_VALIDATION_FAMILY_RESULTS.tsv','FI_R13':root/'shared-genetics-atlas/frailty_sensitivity/FI_R13_ukbcommon.exact.json','FFS_R13':root/'shared-genetics-atlas/frailty_sensitivity/FFS_R13_ukbcommon.exact.json','FI_FFS':root/'shared-genetics-atlas/global_audit/fi_ffs_native_exact.json'}
def read(k):
 with paths[k].open() as f:return list(csv.DictReader(f,delimiter=',' if k=='frailty_joint' else '\t'))
v=read('variants');c=read('context');q=read('gtex');l=read('ld');frailty=read('frailty_candidates');joint={x['variant']:x for x in read('frailty_joint')};assert len(v)==len({r['snp'] for r in v})==425
nodes={};edges=[];variants={x['snp']:x for x in v}
def node(id,kind,label):
 if id in nodes:assert nodes[id]['type']==kind;return
 nodes[id]={'id':id,'type':kind,'label':label}
def edge(src,dst,typ,source,value='',detail=''):edges.append({'source':src,'target':dst,'relation':typ,'value':value,'detail':detail,'provenance':source})
for x in v:
 snp='variant:'+x['snp'];reg='region:'+x['candidate_locus_id'];node(snp,'variant',x['snp']+' '+x['atlas_variant']);node(reg,'region',x['candidate_locus_id']);edge(snp,reg,'in_selected_region','variants',detail='GRCh37 region; variant representation GRCh38')
 for t in ('insomnia','adhd'):
  if x[t+'_credible_set_member']=='YES':
   trait='trait:'+t;node(trait,'trait',t+' (Brain6 GWAS source)');edge(snp,trait,'conditional_abf_credible_set_member','variants',x[t+'_conditional_abf_pip'],'PIP conditional on single-signal model; source A1='+x['source_a1'])
for x in c:
 if x['status']!='ok' or not x['overlapping_genes']:continue
 for item in x['overlapping_genes'].split(';'):
  name,gid=item.split('|',1);gid=gid.split('.')[0];gene='gene:'+gid;snp='variant:'+x['snp'];assert snp in nodes;node(gene,'gene',name+' '+gid);edge(snp,gene,'positional_gene_overlap','context',detail='Ensembl GRCh38 overlap; not causal gene assignment')
qtl_count=0
for x in q:
 if x['type']!='Eqtl' or x['snp'] not in variants or not x['gencode_id']:continue
 match=re.fullmatch(r'(chr[0-9XY]+)_([0-9]+)_([ACGT])_([ACGT])_b38',x['variant_id']);assert match;chrom,pos,ref,alt=match.groups();expected=f'{chrom}:{pos}:{ref}>{alt}';assert variants[x['snp']]['atlas_variant']==expected
 gid=x['gencode_id'].split('.')[0];gene='gene:'+gid;snp='variant:'+x['snp'];node(gene,'gene',x['gene_symbol']+' '+gid);edge(snp,gene,'gtex_v8_significant_eqtl','gtex',x['p_value'],'tissue='+x['tissue']+'; GTEx_gene_symbol='+x['gene_symbol']+'; GTEx_gencode_id='+x['gencode_id']+'; NES_ALT='+x['normalized_effect_size']+'; GTEx significant-pair selection; query_sha256='+x['query_sha256']);qtl_count+=1
for x in l:
 a='variant:'+x['variant_1'];b='variant:'+x['variant_2'];assert a in nodes and b in nodes;edge(a,b,'reference_ld_r2_undirected','ld',x['r2'],'1000 Genomes EUR n='+x['individuals']+'; allele1='+x['counted_allele_1']+'; allele2='+x['counted_allele_2']+'; r='+x['r'])
assert len(frailty)==6 and len({x['variant'] for x in frailty})==6
for x in frailty:
 rs=x['variant'];snp='variant:'+rs;node(snp,'variant',rs+'; source-reported GRCh37 build assumed');y=joint[rs];assert y['eligible_for_joint_screen']=='True' and float(x['MVP_insomnia_p'])==float(y['MVP_p']) and float(x['FinnGen_HFRS_p'])==float(y['HFRS_p'])
 for ident,label,pcol,eallele,oallele in [('trait:insomnia_MVP','MVP insomnia','MVP_p','MVP_effect_allele','MVP_other_allele'),('trait:HFRS_FinnGen','FinnGen HFRS','HFRS_p','HFRS_effect_allele','HFRS_other_allele')]:
  node(ident,'trait',label);edge(snp,ident,'postselected_marginal_association','frailty_joint',y[pcol], 'source_effect_allele='+y[eallele]+'; source_other_allele='+y[oallele]+'; GRCh37 mapping assumed; multiallelic representation and source overlap require review')
 node('trait:FFS_Ye','trait','Ye Fried frailty score');edge(snp,'trait:FFS_Ye','postselected_marginal_association','frailty_candidates',x['FFS_marginal_p'],'effect allele convention inferred; not independent replication; source p is not multiplicity-corrected')
 gene='gene:'+x['nearest_gene_id'].split('.')[0];node(gene,'gene',x['nearest_gene_positional_only']+' '+x['nearest_gene_id']);edge(snp,gene,'nearest_positional_gene','frailty_candidates',x['gene_distance_bp'],'distance in bp; source gene label='+x['nearest_gene_positional_only']+'; selection is post hoc; not gene assignment')
matrix=read('rg_matrix');assert len(matrix)==96 and len({(x['sleep_trait'],x['outcome']) for x in matrix})==96
for x in matrix:
 st='trait:'+x['sleep_trait'];ot={'FI':'trait:FI_Atkins','FFS':'trait:FFS_Ye'}.get(x['outcome'],'trait:'+x['outcome']);node(st,'trait',x['sleep_trait']+' (source as recorded in matrix)');node(ot,'trait',{'FI':'Atkins frailty index','FFS':'Ye Fried frailty score'}.get(x['outcome'],x['outcome']+' (Brain6 GWAS source)'));edge(st,ot,'ldsc_genetic_correlation','rg_matrix',x['rg'],'SE='+x['se']+'; '+x['source']+'; no causal or individual-level interpretation')
h=[x for x in read('hfrs_rg') if x['test_id']=='insomnia_x_HFRS'];assert len(h)==1 and h[0]['status']=='ADMITTED_EXPLORATORY_GLOBAL_RESULT';x=h[0];edge('trait:insomnia','trait:HFRS_FinnGen','ldsc_genetic_correlation','hfrs_rg',x['rg_or_contrast'],'SE='+x['SE']+'; source1='+x['source1']+'; source2='+x['source2']+'; exploratory 16-test family')
node('trait:insomnia_R13','trait','FinnGen R13 registry insomnia')
for key,target in [('FI_R13','trait:FI_Atkins'),('FFS_R13','trait:FFS_Ye'),('FI_FFS','trait:FFS_Ye')]:
 z=json.loads(paths[key].read_text());z=z[0] if isinstance(z,list) else z;start='trait:FI_Atkins' if key=='FI_FFS' else 'trait:insomnia_R13';edge(start,target,'ldsc_genetic_correlation',key,z['rg'],'SE='+str(z['rg_se'])+'; SNPs='+str(z['n_snps'])+'; archived LDSC; cohort overlap/reference uncertainty')
with (out/'nodes.tsv').open('w',newline='') as f:w=csv.DictWriter(f,fieldnames=['id','type','label'],delimiter='\t');w.writeheader();w.writerows(nodes[k] for k in sorted(nodes))
with (out/'edges.tsv').open('w',newline='') as f:w=csv.DictWriter(f,fieldnames=['source','target','relation','value','detail','provenance'],delimiter='\t');w.writeheader();w.writerows(edges)
assert all(x['source'] in nodes and x['target'] in nodes for x in edges);assert qtl_count>0
receipt={'counts_nodes_by_type':dict(collections.Counter(x['type'] for x in nodes.values())),'counts_edges_by_type':dict(collections.Counter(x['relation'] for x in edges)),'source_paths_and_sha256':{k:{'path':str(p),'sha256':hashlib.sha256(p.read_bytes()).hexdigest()} for k,p in paths.items()},'warning':'GTEx and context targets were selected in part using AlphaGenome scores; do not assume graph permitted for model training; review governing terms.'}
(out/'receipt.json').write_text(json.dumps(receipt,indent=2)+'\n');print(json.dumps(receipt,indent=2))
