import pathlib,csv,json,hashlib,math
p=pathlib.Path(__file__).resolve().parent;parent=p.parent
def read(path):
 with path.open() as f:return list(csv.DictReader(f,delimiter='\t'))
def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()
variants=read(parent/'variants.tsv');diag=read(parent/'kriging_diagnostics.tsv');retained=read(p/'retained_variants.tsv');flag={x['variant'] for x in diag if x['fit'] in ('qtl','insomnia','adhd') and abs(float(x['z_std_diff']))>5};rset={x['variant'] for x in retained};assert len(flag)==10 and len(rset)==1538 and rset=={x['variant'] for x in variants}-flag
removed=read(p/'removed_diagnostics.tsv');assert {x['variant'] for x in removed}==flag
cs0={x['fit']:set(x['variants'].split(',')) for x in read(parent/'credible_sets.tsv')};cs1={x['fit']:set(x['variants'].split(',')) for x in read(p/'credible_sets.tsv')};overlap={}
for t in cs0.keys()&cs1.keys():overlap[t]={'before':len(cs0[t]),'after':len(cs1[t]),'intersection':len(cs0[t]&cs1[t]),'jaccard':len(cs0[t]&cs1[t])/len(cs0[t]|cs1[t])}
post={}
for prior in ('1e-05','1e-06'):
 a=read(parent/('coloc_insomnia_'+prior+'.tsv'));b=read(p/('coloc_insomnia_'+prior+'.tsv'));assert len(a)==len(b)==1
 for x in (a[0],b[0]):assert math.isclose(sum(float(x['PP.H'+str(i)+'.abf']) for i in range(5)),1,abs_tol=1e-10)
 post[prior]={'before_H4':float(a[0]['PP.H4.abf']),'after_H4':float(b[0]['PP.H4.abf'])}
result={'removed_n':len(flag),'retained_n':len(rset),'candidate_removed':'chr3_117920728_G_A' in flag,'eqtl_lead_removed':'chr3_117873501_C_A' in flag,'credible_set_overlap':overlap,'h4':post,'parent_variants_sha256':sha(parent/'variants.tsv'),'parent_kriging_sha256':sha(parent/'kriging_diagnostics.tsv'),'parent_ld_sha256':sha(parent/'ld_matrix.tsv.gz'),'retained_sha256':sha(p/'retained_variants.tsv'),'note':'subset and posterior arithmetic validation; no scientific calibration'}
(p/'verification.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))
