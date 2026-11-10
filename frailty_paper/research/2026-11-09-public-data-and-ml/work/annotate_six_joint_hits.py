import csv,json,hashlib,concurrent.futures,subprocess
rs=list(csv.DictReader(open('outputs/song124_mvp_hfrs_joint_screen.csv')));hits=[x for x in rs if x.get('joint_bonferroni_124') and float(x['joint_bonferroni_124'])<0.05]
def one(r):
 p=int(r['song_position_reported']); c=r['song_chr'];urls=[f'https://grch37.rest.ensembl.org/overlap/region/human/{c}:{p-100000}-{p+100000}?feature=gene;content-type=application/json',f'https://grch37.rest.ensembl.org/variation/human/{r["variant"]}?content-type=application/json']; vals=[]
 for url in urls:
  z=subprocess.run(['curl','-fsSL','--max-time','20',url],capture_output=True); vals.append((json.loads(z.stdout),hashlib.sha256(z.stdout).hexdigest()) if z.returncode==0 else (None,''))
 genes=vals[0][0];v=vals[1][0];out={'variant':r['variant'],'chr_GRCh37_assumed':c,'position_GRCh37_assumed':p,'gene_search_url':urls[0],'gene_response_sha256':vals[0][1],'variant_search_url':urls[1],'variant_response_sha256':vals[1][1]}
 if genes is not None:
  valid=[(max(x['start']-p,0,p-x['end']),x) for x in genes if 'start'in x and 'end'in x];valid.sort(key=lambda x:x[0]);near=[x for d,x in valid if d==valid[0][0]] if valid else []
  out.update(nearest_names=';'.join(x.get('external_name','') for x in near),nearest_ids=';'.join(x.get('id','') for x in near),nearest_biotypes=';'.join(x.get('biotype','') for x in near),distance_bp=valid[0][0] if valid else '',gene_status='found' if near else 'none_in_100kb')
 else:out['gene_status']='query_error'
 if v is not None:
  map_=[x for x in v.get('mappings',[]) if x.get('seq_region_name')==c and x.get('start')==p];out['reference_alleles']=';'.join(x.get('allele_string','') for x in map_);out['matching_reference_mappings']=len(map_)
 return out
with concurrent.futures.ThreadPoolExecutor(max_workers=4) as ex: out=list(ex.map(one,hits))
p='outputs/song124_six_hit_positional_genes.csv'; keys=list(dict.fromkeys(k for r in out for k in r));
with open(p,'w',newline='') as f:w=csv.DictWriter(f,fieldnames=keys);w.writeheader();w.writerows(out)
for x in out:print(x['variant'],x.get('nearest_names'),x.get('nearest_biotypes'),x.get('distance_bp'),x.get('reference_alleles'),x.get('gene_status'))
