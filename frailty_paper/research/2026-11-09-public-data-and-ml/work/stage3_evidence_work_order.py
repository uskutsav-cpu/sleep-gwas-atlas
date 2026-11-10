import csv,json,collections
loc=list(csv.DictReader(open('outputs/stage3_provisional_locus_triage.csv')));selected={}
def add(c,gene,tissue,source,h,cs_count,cover):selected[c]={'gene':gene,'tissue':tissue,'source':source,'h4':h,'cs_count':cs_count,'coverage':cover}
def baseline(path, key, val):
 return [x for x in csv.DictReader(open(path)) if x[key]==val and float(x['prior_scale'])==1 and float(x['p12'])==1e-5]
def covtext(d):
 return ';'.join(t+':'+str(z['molecular_credible_variants_matched'])+'/'+str(z['molecular_credible_set_size']) for t,z in d.items())
a=json.load(open('outputs/stage3_original_gwas_coloc_coverage.json'))
for c,an,g,t,src,n in [('rs323509','rs323509_testis','ENSG00000251574','testis','Catalogue QTD000340 splicing',1),('rs34290943','inka1','INKA1','skeletal muscle','Catalogue QTD000281 expression',1),('rs12536335','foxp2_pituitary','FOXP2','pituitary','Catalogue QTD000301 expression',1),('rs2838788','adarb1_muscle_splice','ADARB1','skeletal muscle','Catalogue QTD000285 splicing',1)]:
 r=baseline('outputs/stage3_original_gwas_coloc_sensitivity.csv','analysis',an);h={x['trait']:float(x['H4']) for x in r};add(c,g,t,src,h,n,covtext(a[an]))
 bfile={'rs323509':'rs323509_testis','rs34290943':'inka1_muscle','rs12536335':'foxp2_pituitary','rs2838788':'adarb1_muscle_splice'}[c]
 path='outputs/stage3_'+bfile+'_coloc_sensitivity.csv';rr=[x for x in csv.DictReader(open(path)) if x['trait']=='Ye_FFS' and float(x['prior_scale'])==1 and float(x['p12'])==1e-5];assert len(rr)==1;selected[c]['h4']['Ye_FFS']=float(rr[0]['H4'])
r=[x for x in csv.DictReader(open('outputs/stage3_remaining_coloc_baseline_summary.csv'))];summary={x['candidate']:x for x in csv.DictReader(open('outputs/stage3_remaining_molecular_finemapping_summary.csv'))};cov=json.load(open('outputs/stage3_remaining_coloc_coverage.json'))['results']
for c in sorted({x['candidate'] for x in r}):
 rows=[x for x in r if x['candidate']==c];z=rows[0];h={x['trait']:float(x['H4']) for x in rows};sets=cov[c]['Atkins_FI']['credible_sets'];coverage=';'.join(t+':'+','.join(str(v['matched'])+'/'+str(v['size']) for v in cov[c][t]['credible_sets'].values()) for t in ('Atkins_FI','Jansen_insomnia')) if sets else 'no reported molecular credible set';add(c,z['gene_symbol'],z['tissue'],'Catalogue '+z['dataset']+' expression',h,int(summary[c]['credible_set_count']),coverage)
for prefix,info in [('negr1',{'rs2815764':('NEGR1','brain caudate')}),('positional_fallback',{'rs9845387':('LSAMP','brain cortex'),'rs4886140':('DIAPH3','brain cortex')})]:
 rows=[x for x in csv.DictReader(open('outputs/stage3_'+prefix+'_coloc_sensitivity.csv')) if float(x['prior_scale'])==1 and float(x['p12'])==1e-5]
 for c,(g,t) in info.items():add(c,g,t,'Catalogue positional/fallback expression',{x['trait']:float(x['H4']) for x in rows if x['candidate']==c},0,'no reported molecular credible set')
rows=baseline('outputs/stage3_smim32_original_gwas_coloc_sensitivity.csv','analysis','smim32_muscle');add('rs7734045','SMIM32','skeletal muscle','Catalogue QTD000281 expression',{x['trait']:float(x['H4']) for x in rows},1,covtext(json.load(open('outputs/stage3_smim32_original_gwas_coloc_coverage.json'))['smim32_muscle']))
out=[]
for x in loc:
 cs=[c for c in x['markers'].split(';') if c in selected];assert len(cs)==1,(x['markers'],cs);c=cs[0];z=selected[c];h=z['h4'];fi=h.get('Atkins_FI');ins=h.get('Jansen_insomnia');ffs=h.get('Ye_FFS');assert fi is not None and ins is not None
 if c=='rs323509':tier=1;reason='Both original pairwise H4 >=0.7 and 4/4 molecular variants matched; FFS H4 low and FI prior sensitive; testis context unvalidated'
 elif c=='rs1909650':tier=2;reason='Both original pairwise H4 >=0.5 but no molecular credible set reported; signal multiplicity and coverage unknown'
 elif c in ('rs9845387','rs4886140'):tier=4;reason='No GTEx v10 exact-marker QTL; positional fallback marker nominal p not significant; no molecular credible set'
 else:tier=3;reason='Other selected molecular evidence; no basis for strong shared-trait conclusion; inspect coverage, multiplicity and other genes'
 out.append({'priority_tier':tier,'previous_triage_rank':int(x['triage_rank']),'markers':x['markers'],'representative_gwas_marker':x['triage_marker'],'molecular_anchor_marker':c,'atkins_fi_p_representative':x['atkins_FI_p'],'ye_ffs_p_representative':x['ye_FFS_p'],'jansen_insomnia_p_representative':x['ukb_insomnia_p'],'finngen_insomnia_p_representative':x['finngen_insomnia_p'],'selected_molecular_gene':z['gene'],'selected_molecular_tissue':z['tissue'],'qtl_source':z['source'],'reported_molecular_cs_count':z['cs_count'],'credible_set_coverage_original_gwas':z['coverage'],'atkins_fi_H4_baseline':fi,'jansen_insomnia_H4_baseline':ins,'ye_ffs_H4_baseline':ffs if ffs is not None else 'not calculated for selected QTL','priority_rationale':reason,'warning':'post hoc work priority only; pairwise H4 is model- and prior-dependent, cohorts/phenotypes overlap, and neither signal causality nor risk is established'})
out.sort(key=lambda x:(x['priority_tier'],x['previous_triage_rank']))
for i,x in enumerate(out,1):x['work_order']=i
fields=['work_order']+list(out[0].keys())[:-1]
with open('outputs/stage3_evidence_informed_work_order.csv','w',newline='') as f:w=csv.DictWriter(f,fieldnames=fields);w.writeheader();w.writerows(out)
for x in out[:5]:print(x['work_order'],x['markers'],x['selected_molecular_gene'],x['atkins_fi_H4_baseline'],x['jansen_insomnia_H4_baseline'],x['ye_ffs_H4_baseline'])
