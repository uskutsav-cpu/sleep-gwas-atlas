#!/usr/bin/env python3
"""Extract a 2,001-bp window from a ready JBrowse AlphaGenome variant prediction.

The input is a live status JSON with short-lived signed URLs. Do not publish it.
Outputs omit URLs and credentials. This script calculates local summaries, not
AlphaGenome gene-level RNA scorers or causal inference.
"""
import argparse,csv,hashlib,json,math,ssl,struct,urllib.request
from pathlib import Path

def digest(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def main():
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--status',required=True,type=Path);p.add_argument('--out',required=True,type=Path);p.add_argument('--expected-variant',default='chr5:104612267:T>C');a=p.parse_args()
 source=json.loads(a.status.read_text());m=source['manifest'];v=m['variant'];interval=m['interval']
 if source['status']!='ready' or m['mode']!='variant' or m['organism']!='human':raise ValueError('Not a ready human variant prediction')
 if f"{v['chromosome']}:{v['position']}:{v['reference_bases']}>{v['alternate_bases']}"!=a.expected_variant:raise ValueError('Unexpected variant')
 if m['sequence_length']!=1048576 or m['ontology_terms']!=['UBERON:0001870']:raise ValueError('Unexpected context or tissue')
 p0=v['position']-1; start=p0-1000;end=p0+1001; offset=start-interval['start']; length=end-start
 if length!=2001 or offset<0 or end>interval['end']:raise ValueError('Invalid extraction window')
 a.out.mkdir(parents=True,exist_ok=False);ctx=ssl.create_default_context(cafile='/etc/ssl/cert.pem');tracks=[];raws=[]
 for output in ('dnase','rna_seq'):
  ref=m['outputs'][output]['reference'];alt=m['outputs'][output]['alternate']
  if ref['resolution']!=1 or alt['resolution']!=1 or ref['num_bins']!=m['sequence_length'] or ref['tracks']!=alt['tracks']:raise ValueError('Output metadata mismatch')
  for i,t in enumerate(ref['tracks']):
   track={k:v for k,v in t.items() if k not in ('url','key')};track['output_type']=output;track['track_ordinal']=i; arrays=[]
   for variant,entry in [('reference',ref),('alternate',alt)]:
    begin=i*entry['bytes_per_track']+offset*4;last=begin+length*4-1
    req=urllib.request.Request(entry['url'],headers={'Range':f'bytes={begin}-{last}'})
    with urllib.request.urlopen(req,context=ctx,timeout=45) as response:
     body=response.read();header=response.headers.get('Content-Range','');status=response.status
    if status!=206 or not header.startswith(f'bytes {begin}-{last}/') or len(body)!=length*4:raise ValueError(f'Incorrect HTTP range response {output} {variant} track {i}')
    vals=struct.unpack('<'+str(length)+'f',body)
    if not all(math.isfinite(x) for x in vals):raise ValueError('Non-finite predicted value')
    arrays.append(vals);raws.append({'output_type':output,'track_ordinal':i,'variant_type':variant,'range_bytes':f'{begin}-{last}','sha256':hashlib.sha256(body).hexdigest()})
   tracks.append((track,*arrays))
 with (a.out/'track_window.tsv').open('w',newline='') as f:
  fields=['chromosome','position_grch38_1based','offset_from_variant_bp','output_type','track_ordinal','track_name','biosample_name','life_stage','data_source','strand','reference','alternate','delta_alt_minus_ref']
  w=csv.DictWriter(f,fieldnames=fields,delimiter='\t',lineterminator='\n');w.writeheader()
  for t,r,b in tracks:
   for i,(ref,alt) in enumerate(zip(r,b)):
    w.writerow({'chromosome':v['chromosome'],'position_grch38_1based':start+i+1,'offset_from_variant_bp':i-1000,'output_type':t['output_type'],'track_ordinal':t['track_ordinal'],'track_name':t.get('name'),'biosample_name':t.get('biosample_name'),'life_stage':t.get('biosample_life_stage'),'data_source':t.get('data_source'),'strand':t.get('strand'),'reference':repr(ref),'alternate':repr(alt),'delta_alt_minus_ref':repr(alt-ref)})
 summary=[]
 for t,r,b in tracks:
  refsum=sum(r[750:1251]);altsum=sum(b[750:1251]);ratio=math.log2((altsum+1)/(refsum+1)) if refsum+1>0 and altsum+1>0 else None
  summary.append({'output_type':t['output_type'],'track_ordinal':t['track_ordinal'],'track_name':t.get('name'),'life_stage':t.get('biosample_life_stage'),'data_source':t.get('data_source'),'strand':t.get('strand'),'ref_sum_501bp':refsum,'alt_sum_501bp':altsum,'delta_sum_501bp':altsum-refsum,'local_log2_ratio_plus_one':ratio,'max_abs_delta_2001bp':max(abs(x-y) for x,y in zip(r,b)),'note':'local signal summary, not a gene-level score' if t['output_type']=='rna_seq' else '501-bp sum ratio; predicted track only'})
 with (a.out/'track_summary.tsv').open('w',newline='') as f:
  w=csv.DictWriter(f,fieldnames=list(summary[0]),delimiter='\t',lineterminator='\n');w.writeheader();w.writerows(summary)
 receipt={'source':'JBrowse public AlphaGenome prediction service','source_documentation':'https://www.jbrowse.org/jb2/docs/tutorials/alphagenome/','alpha_genome_model_release':'unknown; not provided by proxy','request':source['request'],'manifest_created_epoch':m['created'],'tracks':[t for t,_,_ in tracks],'range_sha256':raws,'output_sha256':{k:digest(a.out/k) for k in ['track_window.tsv','track_summary.tsv']},'coordinates':'model interval 0-based half-open; output positions 1-based GRCh38','window':'variant ±1000bp; 501-bp summary ±250bp','limitations':'predicted signals not measurements, not evidence of shared trait causality; API outputs for noncommercial research; no model training'}
 (a.out/'provenance.json').write_text(json.dumps(receipt,indent=2,sort_keys=True)+'\n')
 print(json.dumps(summary,indent=2))
if __name__=='__main__':main()
