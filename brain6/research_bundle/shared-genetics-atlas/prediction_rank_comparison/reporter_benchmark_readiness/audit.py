import collections,csv,hashlib,json,pathlib
p=pathlib.Path(__file__).resolve().parent;atlas=p.parent.parent
assay=atlas/'variant_context/siraj_2026_mpra/sknsh_results.tsv';ranking=p.parent/'variant_ranks.tsv'
read=lambda f:list(csv.DictReader(f.open(),delimiter='\t'))
a=read(assay);r=read(ranking);by=collections.defaultdict(list)
for x in r:by[x['snp']].append(x)
out=[]
for x in a:
 candidates=by[x['rsid']];assert len(candidates)==1,(x['rsid'],len(candidates));y=candidates[0]
 f=float(x['source_fdr']);a_p=float(x['a_activity_bonferroni_p']);b_p=float(x['b_activity_bonferroni_p'])
 assert all(0<=v<=1 for v in (f,a_p,b_p))
 out.append(dict(rsid=x['rsid'],library=x['library'],region=y['region'],dnase=y['dnase_alt_ref_log2_501bp'],insomnia_pip=y['insomnia_pip'],adhd_pip=y['adhd_pip'],source_fdr=f,activity_a_bonf=a_p,activity_b_bonf=b_p,skew_fdr_pass=f<.05,proxy_both=f<.05 and min(a_p,b_p)<.05))
with (p/'joined_rows.tsv').open('w',newline='') as f:w=csv.DictWriter(f,fieldnames=list(out[0]),delimiter='\t');w.writeheader();w.writerows(out)
ids={x['rsid'] for x in out};libraries=collections.Counter(x['rsid'] for x in out)
summary=dict(assay_rows=len(out),distinct_variants=len(ids),unique_variant_region_counts=dict(collections.Counter(next(x['region'] for x in out if x['rsid']==i) for i in ids)),variant_multiple_library=sum(n>1 for n in libraries.values()),all_rows_have_dnase=sum(bool(x['dnase']) for x in out),all_rows_have_both_pips=sum(bool(x['insomnia_pip']) and bool(x['adhd_pip']) for x in out),skew_fdr_rows=sum(x['skew_fdr_pass'] for x in out),row_proxy_positive=sum(x['proxy_both'] for x in out),proxy_positive_variants=len({x['rsid'] for x in out if x['proxy_both']}),validated_negative_labels=0)
assert summary['row_proxy_positive']==0 and summary['skew_fdr_rows']==2
sha=lambda f:hashlib.sha256(f.read_bytes()).hexdigest()
receipt=dict(summary=summary,inputs={str(assay):sha(assay),str(ranking):sha(ranking)},joined_sha256=sha(p/'joined_rows.tsv'),note='row proxy is not the complete author classification; missing significant effects are not true negatives; no model fitted')
(p/'receipt.json').write_text(json.dumps(receipt,indent=2)+'\n');print(json.dumps(summary,indent=2))
