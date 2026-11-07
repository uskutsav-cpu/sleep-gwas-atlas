import csv,gzip,hashlib,json,math,struct
from pathlib import Path
p=Path(__file__).resolve().parent;root=p.parents[2].parent/'brain6-alphagenome-pilot/results/jbrowse_expanded'
# root is outputs/brain6-alphagenome-pilot/results/jbrowse_expanded.
variant='chr3:117920728:G>A';folder=root/'variants'/hashlib.sha256(variant.encode()).hexdigest()[:20]
m=json.loads((folder/'provenance.json').read_text());raw=(folder/'track_window.f32.gz').read_bytes();assert hashlib.sha256(raw).hexdigest()==m['track_window_sha256'];a=struct.unpack('<16008f',gzip.decompress(raw));tracks=[i for i,t in enumerate(m['tracks']) if t['output_type']=='dnase' and t['metadata'].get('biosample_life_stage')=='adult'];assert len(tracks)==1;i=tracks[0];ref=a[i*4002:i*4002+2001];alt=a[i*4002+2001:(i+1)*4002]
start=m['interval']['start'];entry=next(x for x in m['source_ranges'] if x['output_type']=='dnase' and x['variant_type']=='reference');byte0=int(entry['range_bytes'].split('-')[0]);assert byte0%4==0
origin=start+byte0//4;center=117920728-1;assert origin+1000==center
win=(center-250,center+251);peak=(117920181,117920630);assert win==(117920477,117920978) and origin<=peak[0]<peak[1]<=origin+2001
with (p.parents[2]/'alphagenome_window_sensitivity/window_scores.tsv').open() as f:prior=next(r for r in csv.DictReader(f,delimiter='\t') if r['snp']=='rs62264764')
def calc(name,left,right):
 assert origin<=left<right<=origin+2001;v=ref[left-origin:right-origin];w=alt[left-origin:right-origin];assert min(v)>=0 and min(w)>=0
 x=sum(v);y=sum(w);return {'segment':name,'start0':left,'end0':right,'bases':right-left,'ref_sum':x,'alt_sum':y,'alt_minus_ref':y-x,'log2ratio_offset1':math.log2((y+1)/(x+1))}
rows=[calc('window501',*win),calc('peak_within_window',max(win[0],peak[0]),min(win[1],peak[1])),calc('window_outside_peak',peak[1],win[1]),calc('entire_peak',*peak)]
assert abs(rows[0]['log2ratio_offset1']-float(prior['log2ratio_501_pc1']))<1e-10
assert abs(rows[0]['alt_minus_ref']-rows[1]['alt_minus_ref']-rows[2]['alt_minus_ref'])<1e-12
for shift in (-1,1):rows.append(calc('peak_in_window_shift_'+str(shift),max(win[0],peak[0]+shift),min(win[1],peak[1]+shift)))
with (p/'segment_scores.tsv').open('w',newline='') as f:w=csv.DictWriter(f,fieldnames=list(rows[0]),delimiter='\t');w.writeheader();w.writerows(rows)
receipt={'variant':variant,'track_sha256':m['track_window_sha256'],'source_model_release':m['model_release'],'source_track_type':'adult frontal cortex DNase','input_interval_start':start,'reference_byte_offset':byte0,'track_window_start_0based':origin,'center_0based':center,'peak_0based_assumption':list(peak),'prior_501_log2ratio_match':True,'delta_decomposition_exact':True,'coordinate_caveat':'Peak is from developmental cells and assumed BED format; predicted track is adult frontal cortex. No genotype accessibility effect is measured.'}
(p/'receipt.json').write_text(json.dumps(receipt,indent=2)+'\n');print(json.dumps(rows,indent=2))
