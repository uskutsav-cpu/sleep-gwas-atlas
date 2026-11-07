import pathlib,hashlib,zlib,csv,io,collections,json
p=pathlib.Path(__file__).resolve().parent
with (p/'selected_rows.tsv').open() as f:selected=list(csv.DictReader(f,delimiter='\t'))
with (p/'build_matches.tsv').open() as f:build={x['rsid']:x for x in csv.DictReader(f,delimiter='\t')}
assert len(selected)==20
names=('Low_Grip_strength','Low_Physical_activity','Slow_Walking_speed','Weight_loss');out=[];meta={}
for name in names:
 raw=(p/(name+'_first_131072_bytes.bin')).read_bytes();assert len(raw)==131072
 z=zlib.decompressobj(31);s=z.decompress(raw).decode();rows=list(csv.DictReader(io.StringIO(s[:s.rfind('\n')]),delimiter='\t'));index=collections.defaultdict(list)
 for x in rows:index[x['ID']].append(x)
 for x in selected:
  ys=index[x['ID']];matching=[y for y in ys if y['CHROM']==x['CHROM'] and y['GENPOS']==x['GENPOS'] and {y['ALLELE0'],y['ALLELE1']}=={x['ALLELE0'],x['ALLELE1']}]
  out.append({'component':name,'rsid':x['ID'],'id_rows':len(ys),'matching_rows':len(matching),'status':'UNIQUE_EXACT' if len(ys)==len(matching)==1 else 'ABSENT' if not ys else 'NONUNIQUE_OR_DIFFERENT','exhaustion_direct_reference_matches':build[x['ID']]['direct_exact_matches']})
 meta[name]={'prefix_sha256':hashlib.sha256(raw).hexdigest(),'complete_rows':len(rows),'unique_exact':sum(q['component']==name and q['status']=='UNIQUE_EXACT' for q in out)}
with (p/'component_prefix_matches.tsv').open('w',newline='') as f:w=csv.DictWriter(f,list(out[0]),delimiter='\t');w.writeheader();w.writerows(out)
(p/'other_prefix_receipt.json').write_text(json.dumps({'components':meta,'selection_sha256':hashlib.sha256((p/'selected_rows.tsv').read_bytes()).hexdigest(),'results_sha256':hashlib.sha256((p/'component_prefix_matches.tsv').read_bytes()).hexdigest(),'caveat':'partial files; does not establish provenance of full archives or effect coding'},indent=2)+'\n');print(json.dumps(meta,indent=2))
