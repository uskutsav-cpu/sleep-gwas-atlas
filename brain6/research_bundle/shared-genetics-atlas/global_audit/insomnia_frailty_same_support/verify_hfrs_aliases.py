import csv,json,pathlib,subprocess,math,hashlib
p=pathlib.Path(__file__).resolve().parent;receipt=json.loads((p/'hfrs_raw_alignment_receipt.json').read_text());source=pathlib.Path(receipt['raw_path']);items=receipt['missing_filtered_rows'];assert len(items)==4 and all(not x['in_final_support'] for x in items); ids={x['rsid']:x for x in items};argv=['rg','-F']
for rs in ids:argv+=['-e',rs]
argv.append(str(source));text=subprocess.check_output(argv,text=True);matches={};comp=str.maketrans('ACGT','TGCA')
for line in text.splitlines():
 c=line.split('\t');assert len(c)==9
 for rs in set(c[8].split(','))&ids.keys():
  t=ids[rs]['filtered'];a,b,beta,se,pval,n,eaf,ch=t;assert c[0]==ch
  if (c[2],c[3])==(a,b):flip=1
  elif (c[3],c[2])==(a,b):flip=-1
  elif (c[2].translate(comp),c[3].translate(comp))==(a,b):flip=1
  elif (c[3].translate(comp),c[2].translate(comp))==(a,b):flip=-1
  else:continue
  expectedf=float(c[6]) if flip==1 else 1-float(c[6]);assert all(math.isclose(float(x),float(y),rel_tol=1e-10,abs_tol=1e-12) for x,y in [(beta,float(c[4])*flip),(se,c[5]),(pval,c[7]),(eaf,expectedf)]) and float(n)==500737
  assert rs not in matches;matches[rs]={'source_row':c,'effect_flip':flip,'in_final_support':False}
assert set(matches)==set(ids)
(p/'hfrs_alias_receipt.json').write_text(json.dumps({'input_raw_sha256':receipt['raw_sha256'],'alias_matches':matches,'status':'PASS','limitation':'Source rs_id lists comma-separated aliases; matched by token, chromosome and alleles, with no independent coordinate-liftover verification.'},indent=2)+'\n');print('all four excluded markers match a comma-separated source alias, chromosome, alleles and values')
