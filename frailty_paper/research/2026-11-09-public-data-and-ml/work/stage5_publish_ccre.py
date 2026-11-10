import csv,json,pathlib,zipfile,collections,hashlib,datetime
queue=list(csv.DictReader(open('outputs/stage5_score_queue.csv')));raw=pathlib.Path('work/stage5_ccre_raw');index={}
for f in raw.glob('*.json'):
 obj=json.loads(f.read_text());d=obj['response'];index[(d['chrom'],d['start'],d['end'])]=(obj['url'],d['cCREregistry'])
rows=[];hits=[]
for v in queue:
 chrom=v['chr'];a=int(v['position_GRCh38'])-1;b=a+len(v['reference_allele']);cover=[x for k,x in index.items() if k[0]==chrom and k[1]<=a and k[2]>=b];assert len(cover)==1,(v['variant_key'],len(cover));url,features=cover[0]
 match=[x for x in features if x['chromStart']<b and x['chromEnd']>a];rows.append({'variant_GRCh38':v['variant_key'],'input_classes':v['input_classes'],'reference_span_start_0based':a,'reference_span_end_0based_exclusive':b,'overlap_count':len(match),'overlapping_ccre_ids':'|'.join(x['name'] for x in match),'overlapping_ccre_classes':'|'.join(x['cCRE_class'] for x in match),'query_status':'complete','source_query':url,'interpretation':'Reference-span overlap only; insertion uses anchor base; combined registry, tissue unknown; not allele-specific proof'})
 for x in match:hits.append({'variant_GRCh38':v['variant_key'],'ccre_id':x['name'],'ccre_start_0based':x['chromStart'],'ccre_end_0based_exclusive':x['chromEnd'],'ccre_class':x['cCRE_class'],'DNase_maxZ':x['DNase_maxZ'],'H3K4me3_maxZ':x['H3K4me3_maxZ'],'H3K27ac_maxZ':x['H3K27ac_maxZ'],'CTCF_maxZ':x['CTCF_maxZ'],'source_query':url})
def write(path,data,fields):
 with open(path,'w',newline='') as f:
  w=csv.DictWriter(f,fieldnames=fields);w.writeheader();w.writerows(data)
write('outputs/stage5_encode_ccre_variant_overlap.csv',rows,list(rows[0]));write('outputs/stage5_encode_ccre_overlap_detail.csv',hits,['variant_GRCh38','ccre_id','ccre_start_0based','ccre_end_0based_exclusive','ccre_class','DNase_maxZ','H3K4me3_maxZ','H3K27ac_maxZ','CTCF_maxZ','source_query'])
with zipfile.ZipFile('outputs/stage5_encode_ccre_query_responses.zip','w',zipfile.ZIP_DEFLATED) as z:
 for f in sorted(raw.glob('*.json')):z.write(f,f.name)
with zipfile.ZipFile('outputs/stage5_encode_ccre_query_responses.zip') as z:assert z.testzip() is None
stats={'variants':len(rows),'successful_query_groups':len(index),'variants_with_overlap':sum(int(x['overlap_count'])>0 for x in rows),'overlap_pairs':len(hits),'classes':dict(collections.Counter(x['ccre_class'] for x in hits)),'build':'hg38','source':'UCSC ENCODE4 cCREregistry combined registry','access_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'sha256_overlap_csv':hashlib.sha256(pathlib.Path('outputs/stage5_encode_ccre_variant_overlap.csv').read_bytes()).hexdigest()}
pathlib.Path('outputs/stage5_encode_ccre_overlap_receipt.json').write_text(json.dumps(stats,indent=2)+'\n');print(json.dumps(stats,indent=2))
