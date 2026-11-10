import csv,pathlib
l=list(csv.DictReader(open('outputs/stage3_provisional_locus_triage.csv')));q=list(csv.DictReader(open('outputs/stage3_gtex_v10_candidate_qtls.csv')));meta=list(csv.DictReader(open('work/eqtl_catalogue_r7_metadata.tsv'),delimiter='\t'));m={(x['sample_group'],x['quant_method']):x for x in meta if x['study_id']=='QTS000015'}
def tissue(t):
 t=t.lower();aliases={'brain_anterior_cingulate_cortex_ba24':'brain_anterior_cingulate_cortex','brain_frontal_cortex_ba9':'brain_frontal_cortex','brain_nucleus_accumbens_basal_ganglia':'brain_nucleus_accumbens','brain_caudate_basal_ganglia':'brain_caudate','brain_putamen_basal_ganglia':'brain_putamen'};return aliases.get(t,t)
def tier(t):
 t=t.lower();return 0 if t.startswith('brain_') or t in ('muscle_skeletal','whole_blood','artery_tibial','artery_aorta','artery_coronary') else 1 if t=='nerve_tibial' else 2 if t=='pituitary' else 3
out=[]
for x in l:
 if int(x['triage_rank']) in (1,2,3,4,6):continue
 a=[z for z in q if z['candidate'] in x['markers'].split(';') and z['qtl_type']=='Eqtl' and (tissue(z['tissue']),'ge') in m];a.sort(key=lambda z:(tier(z['tissue']),float(z['p_value']),z['candidate'],z['gencode_id']))
 if not a:out.append({'rank':x['triage_rank'],'markers':x['markers'],'status':'no eligible exact-marker GTEx v10 eQTL in mapped tissue; splicing and proxy mapping remain','candidate':'','gene_symbol':'','gene_id':'','tissue_v10':'','dataset_r7':'','n_samples':'','marker_grch38':'','candidate_ref':'','candidate_alt':'','p_v10':''});continue
 z=a[0];v=z['gtex_variant_id'].split('_');out.append({'rank':x['triage_rank'],'markers':x['markers'],'status':'selected for exploratory release-7 regional review','candidate':z['candidate'],'gene_symbol':z['gene_symbol'],'gene_id':z['gencode_id'].split('.')[0],'tissue_v10':z['tissue'],'dataset_r7':m[(tissue(z['tissue']),'ge')]['dataset_id'],'n_samples':m[(tissue(z['tissue']),'ge')]['sample_size'],'marker_grch38':v[0]+':'+v[1],'candidate_ref':v[2],'candidate_alt':v[3],'p_v10':z['p_value']})
with open('outputs/stage3_remaining_molecular_worklist.csv','w',newline='') as f:w=csv.DictWriter(f,fieldnames=out[0]);w.writeheader();w.writerows(out)
print('remaining groups',len(out),'selected',sum(bool(x['candidate']) for x in out),'unselected',sum(not x['candidate'] for x in out))
