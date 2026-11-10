import pathlib,json,math,hashlib,csv,statistics
base=pathlib.Path('/Users/swethasunilkumar/Documents/Codex/2026-10-08/g/work/sleep-gwas-atlas/frailty_paper/nature_communications_continuation_v6/repro_inputs/native_exact')
expected={'aux_R13_common_FI_FFS.exact.json':'8a48e11b9e93b3c49233422e28784943a51224d374aa928f5269c5d00017107c','aux_MVP_common_FI_FFS.exact.json':'19cd1edaff3e16ba0ba0d6982dc1510bde9cfdda65a9d38b90c4e465db752f4d'}
rows=[]
for name,sha in expected.items():
 f=base/name; data=f.read_bytes(); assert hashlib.sha256(data).hexdigest()==sha
 d=json.loads(data)[0]; a=d['hsq1']['total'];b=d['hsq2']['total']; rg=d['gencov']['total']/math.sqrt(a*b);assert abs(rg-d['rg'])<1e-12
 ds=[z[0] for z in d['delete_rg']]; mean=statistics.mean(ds); se=math.sqrt((len(ds)-1)/len(ds)*sum((v-mean)**2 for v in ds)); assert len(ds)==200 and abs(se-d['rg_se'])<1e-12
 rows.append({'support': 'FinnGen-insomnia matched support' if 'R13' in name else 'MVP-insomnia matched support', 'trait_1':'FI Atkins 2021 UKB + TwinGene','trait_2':'FFS Ye 2023 UKB', 'rg':rg,'jackknife_se':se,'wald_ci_lower':rg-1.959963984540054*se,'wald_ci_upper':rg+1.959963984540054*se,'n_snps':d['n_snps'],'n_jackknife_blocks':len(ds),'ordered_snp_sha256':d['row_order_SNP_sha256'],'exact_json_sha256':sha,'status':'Auxiliary matched-support LDSC; same two frailty sources; not independent replication','source_file':str(f)})
with open('outputs/direct_fi_ffs_genetic_correlation.csv','w',newline='') as f:
 writer=csv.DictWriter(f,fieldnames=rows[0]);writer.writeheader();writer.writerows(rows)
print(json.dumps(rows,indent=2))
