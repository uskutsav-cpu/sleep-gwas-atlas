#!/usr/bin/env python3
"""Stream and checksum the public CHARGE/UKB low-grip GWAS; retain selected rows."""
import requests,zlib,hashlib,csv,io,json
from pathlib import Path
import pandas as pd
O=Path('outputs');s15=pd.read_csv(O/'unvalidated_15_variant_score_prototype.csv');s22=pd.read_csv(O/'unvalidated_22_shared_direction_score.csv');ids=set(s15.variant)|set(s22.variant);url='https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90007001-GCST90008000/GCST90007526/GCST90007526_buildGRCh37.tsv.gz';expected='10ed84976ef644b1ac7ceb06cd34e2dc';r=requests.get(url,stream=True,timeout=(20,90));r.raise_for_status();z=zlib.decompressobj(31);h=hashlib.md5();buffer=b'';header=None;matches=[];n=0;total=0;nextprogress=64000000
for chunk in r.iter_content(1048576):
 if not chunk:continue
 h.update(chunk);total+=len(chunk);buffer+=z.decompress(chunk);lines=buffer.split(b'\n');buffer=lines.pop()
 for line in lines:
  if header is None:header=line.decode().split('\t');continue
  n+=1
  if line.split(b'\t',1)[0].decode() in ids:matches.append(dict(zip(header,line.decode().split('\t'))))
 if total>=nextprogress:print('read bytes',total,flush=True);nextprogress+=64000000
buffer+=z.flush();assert z.eof and not buffer.strip() and h.hexdigest()==expected,(z.eof,h.hexdigest());out=O/'stage6_grip_gwas_raw_candidates.csv';pd.DataFrame(matches,columns=header).to_csv(out,index=False);receipt={'source':url,'expected_md5':expected,'observed_md5':h.hexdigest(),'source_checksum_matches':True,'bytes':total,'data_rows':n,'selected_ids':len(ids),'selected_source_rows':len(matches),'warning':'Low grip strength is one physical frailty component, not FI or FFS. Source includes UK Biobank and is not independent of discovery; exact allele and coordinates require additional checks.'};(O/'stage6_grip_gwas_extract_receipt.json').write_text(json.dumps(receipt,indent=2)+'\n');print(json.dumps(receipt,indent=2),flush=True)
