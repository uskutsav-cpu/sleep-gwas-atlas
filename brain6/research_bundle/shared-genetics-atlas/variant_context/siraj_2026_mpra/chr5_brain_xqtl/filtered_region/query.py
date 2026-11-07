import urllib.request,csv,io,hashlib,json,collections
from pathlib import Path
p=Path(__file__).resolve().parent
url='https://xqtl.niagads.org/download/region/chr5%3A104400000-104700000';req=urllib.request.Request(url,headers={'Accept-Encoding':'identity'})
h=hashlib.sha256();n=0;size=0;out=[];counts=collections.Counter()
with urllib.request.urlopen(req,timeout=90) as response:
 lines=iter(response)
 header=next(lines);h.update(header);size+=len(header);fields=header.decode().rstrip('\r\n').split('\t');assert 'rsID' in fields and 'target' in fields
 for line in lines:
  h.update(line);size+=len(line);n+=1
  vals=line.decode().rstrip('\r\n').split('\t');assert len(vals)==len(fields)
  r=dict(zip(fields,vals));counts[r['xQTLtype']]+=1
  if r['xQTLtype']=='mQTL' and r['context']=='PC' and r['study']=='Knight-ADRC' and r['target']=='chr5:104484889-104484891:cg24905254':out.append(r)
with (p/'target_filtered.tsv').open('w',newline='') as f:
 w=csv.DictWriter(f,fieldnames=fields,delimiter='\t');w.writeheader();w.writerows(out)
receipt={'url':url,'source_rows':n,'source_bytes':size,'source_stream_sha256':h.hexdigest(),'source_not_saved_in_full':True,'counts_by_type':dict(counts),'target_rows':len(out),'output_sha256':hashlib.sha256((p/'target_filtered.tsv').read_bytes()).hexdigest(),'filter':'HMT-significant portal results only; not all tested SNPs and not valid input to dense colocalization'}
(p/'receipt.json').write_text(json.dumps(receipt,indent=2)+'\n');print(json.dumps(receipt,indent=2))
