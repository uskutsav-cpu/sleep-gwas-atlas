#!/usr/bin/env python3
import hashlib,json,zipfile
from pathlib import Path
from remote_zip import Remote
p=Path(__file__).resolve().parent;meta=json.loads((p/'archive_inventory.json').read_text());remote=Remote(meta['url'],meta['archive_bytes'])
rows=0;selected=[]
with zipfile.ZipFile(remote) as z:
 info=z.getinfo('ldm13M/snp.info')
 with z.open(info) as f:
  header=next(f).decode().split()
  assert header==['Chrom','ID','Index','GenPos','PhysPos','A1','A2','A1Freq','N','Block']
  for line in f:
   rows+=1
   if 2682830<=rows-1<=2696109:selected.append(line.decode())
 assert len(selected)==13280
 for i,line in enumerate(selected):
  parts=line.split();assert int(parts[2])==2682830+i and int(parts[-1])==341
 text='\t'.join(header)+'\n'+''.join('\t'.join(row.split())+'\n' for row in selected)
 (p/'block341_snp.tsv').write_text(text)
 result={'url':meta['url'],'archive_bytes':meta['archive_bytes'],'member':info.filename,'member_compressed_bytes':info.compress_size,'member_uncompressed_bytes':info.file_size,'member_crc32':f'{info.CRC:08x}','zip_member_crc_passed':True,'source_rows':rows,'selected_block_rows':len(selected),'selected_sha256':hashlib.sha256(text.encode()).hexdigest(),'http_range_bytes_transferred':remote.transferred,'limitations':'Full SNP member CRC checked, full archive not downloaded or hashed.'}
 (p/'snp_receipt.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))
