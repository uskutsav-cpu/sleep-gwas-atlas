import csv,hashlib,json,math,pathlib
p=pathlib.Path(__file__).resolve().parent;root=p.parents[1]
def read(q):
 with q.open(newline='') as f:return list(csv.DictReader(f,delimiter='\t'))
mvp=read(p/'mvp_result.tsv');receipt=json.loads((p/'query_receipt.json').read_text());ref=read(p/'pairwise_ld.tsv');alpha=root.parent/'brain6-alphagenome-pilot/results/jbrowse_expanded/expanded_variant_evidence.tsv';a=[x for x in read(alpha) if x['snp'] in ('rs62264764','rs62264767')];assert len(mvp)==len(ref)==1 and len(a)==2 and len({x['snp'] for x in a})==2
q=receipt['exact_records'];assert len(q)==1 and q[0]['rsid']=='rs62264764' and mvp[0]['status']=='MATCHED';v=q[0];aligned=math.log(float(v['odds_ratio']))*(1 if v['effect_allele']=='A' else -1);assert math.isclose(aligned,float(mvp[0]['mvp_logodds_aligned']),abs_tol=1e-14) and float(v['p_value'])==float(mvp[0]['mvp_p'])
r=ref[0];assert receipt['reference_haplotype_discordance']==round(int(r['haplotypes'])*(float(r['allele1_frequency'])+float(r['allele2_frequency'])-2*float(r['p11']))) and receipt['reference_haplotype_discordance']==0
with (p/'alpha_comparison.tsv').open('w',newline='') as f:
 w=csv.DictWriter(f,fieldnames=['snp','atlas_variant','dnase_log2_ratio_501bp','prediction_metadata_sha256','prediction_window_sha256'],delimiter='\t');w.writeheader();w.writerows({k:z[k] for k in w.fieldnames} for z in a)
check={'mvp_status':mvp[0]['status'],'eur_panel_haplotypes':receipt['reference_haplotypes'],'discordant_haplotypes':receipt['reference_haplotype_discordance'],'prediction_input_sha256':hashlib.sha256(alpha.read_bytes()).hexdigest(),'comparison_sha256':hashlib.sha256((p/'alpha_comparison.tsv').read_bytes()).hexdigest(),'note':'Allele, arithmetic and data-binding checks; not scientific or causal validation.'};(p/'verification.json').write_text(json.dumps(check,indent=2)+'\n');print(json.dumps(check,indent=2))
