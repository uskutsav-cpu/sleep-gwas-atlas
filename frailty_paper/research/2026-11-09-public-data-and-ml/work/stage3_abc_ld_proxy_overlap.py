import csv,gzip,collections,bisect,json,pathlib,time
p=pathlib.Path('work/stage3_region_cache/ABC_Nasser2021_predictions_GRCh37.txt.gz')
a=json.load(open('outputs/stage3_abc_candidate_regulatory_receipt.json'))
proxy=list(csv.DictReader(open('outputs/stage3_candidate_reference_ld_proxies.csv')))
by=collections.defaultdict(lambda:collections.defaultdict(list))
for x in proxy:
 if x['proxy_pos_GRCh37'] != x['candidate_pos_GRCh37']:
  by['chr'+x['candidate_chr_GRCh37']][int(x['proxy_pos_GRCh37'])-1].append(x)
positions={k:sorted(v) for k,v in by.items()};out=[];n=0;start=time.time()
with gzip.open(p,'rt') as f:
 header=f.readline().rstrip('\n').split('\t');assert header[:3]==['chr','start','end']
 for line in f:
  n+=1;v=line.rstrip('\n').split('\t');ch=v[0]
  if ch not in by:continue
  lo,hi=int(v[1]),int(v[2]);ps=positions[ch]; i=bisect.bisect_left(ps,lo);j=bisect.bisect_left(ps,hi)
  if i==j:continue
  for pos in ps[i:j]:
   for x in by[ch][pos]:
    out.append({'candidate':x['candidate'],'proxy_reference_id':x['proxy_reference_id'],'proxy_chr_GRCh37':x['candidate_chr_GRCh37'],'proxy_pos_GRCh37':x['proxy_pos_GRCh37'],'proxy_bim_A1_A2':x['proxy_bim_A1_A2'],'proxy_r2':x['r2_genotype_dosage'],**dict(zip(header,v))})
fields=['candidate','proxy_reference_id','proxy_chr_GRCh37','proxy_pos_GRCh37','proxy_bim_A1_A2','proxy_r2']+header
with open('outputs/stage3_abc_ld_proxy_regulatory_links.csv','w',newline='') as f:
 w=csv.DictWriter(f,fieldnames=fields);w.writeheader();w.writerows(out)
count=collections.Counter(x['candidate'] for x in out)
json.dump({'source':a['source'],'source_sha256_prior_verified':a['local_source_sha256'],'source_rows_scanned':n,'gzip_read_to_end':True,'candidate_proxy_pairs_input':len(proxy),'nonself_candidate_proxy_pairs_input':sum(x['proxy_pos_GRCh37']!=x['candidate_pos_GRCh37'] for x in proxy),'regulatory_prediction_rows':len(out),'candidate_counts':count,'interval_rule':'1-based proxy position minus one within source 0-based half-open interval; source build is GRCh37 per prior documentation','note':'reference LD plus predicted regulatory links, not GWAS finemapping or causal evidence','seconds':round(time.time()-start,1)},open('outputs/stage3_abc_ld_proxy_overlap_receipt.json','w'),indent=2)
print('source rows',n,'prediction rows',len(out),'candidates',dict(count),flush=True)
