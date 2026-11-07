import collections,csv,gzip,hashlib,json,math
from pathlib import Path
p=Path(__file__).resolve().parent;receipt=json.loads((p/'retrieval_receipt.json').read_text());centers={'rs62264764':117920728,'rs62264767':117923158};by={}
for src in receipt['members']:
 peaks={}
 with gzip.open(src['file'],'rt') as f:
  next(f)
  for line in f:
   var,beta,pv,peak,effect=line.split();parts=peak.split('_');start,end=int(parts[1]),int(parts[2]);eligible=[k for k,v in centers.items() if start < v+250 and end > v-251]
   if not eligible:continue
   key=peak;z=float(pv)
   if key not in peaks:peaks[key]={'cell':src['cell'],'peak':peak,'start':start,'end':end,'overlapping_windows':','.join(eligible),'source_rows':0,'smallest_nominal_p':1.,'best_variant':''}
   peaks[key]['source_rows']+=1
   if math.isfinite(z) and z<peaks[key]['smallest_nominal_p']:peaks[key]['smallest_nominal_p']=z;peaks[key]['best_variant']=var
 by[src['cell']]=list(peaks.values())
rows=[x for a in by.values() for x in a];fields=['cell','peak','start','end','overlapping_windows','source_rows','smallest_nominal_p','best_variant']
with (p/'overlap_peaks.tsv').open('w',newline='') as f:w=csv.DictWriter(f,fieldnames=fields,delimiter='\t');w.writeheader();w.writerows(rows)
out={'counts':{k:len(v) for k,v in by.items()},'coordinate_assumption':'peak identifiers are 0-based half-open; this has not been independently verified','interpretation':'Peak detected in filtered source is not a caQTL at either candidate variant; minima are post-selected uncorrected over many tests. Absence does not establish no accessibility.','file_sha256':hashlib.sha256((p/'overlap_peaks.tsv').read_bytes()).hexdigest()}
(p/'peak_receipt.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out,indent=2));print(rows)
