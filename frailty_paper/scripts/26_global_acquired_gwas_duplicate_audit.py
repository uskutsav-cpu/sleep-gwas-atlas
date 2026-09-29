#!/usr/bin/env python3
"""Global coordinate/allele duplicate audit of acquired component/context GWAS.

Coordinates are checked in a streaming pass. Ordered inputs can be deduplicated
within each coordinate in constant memory; unordered inputs use external sort.
Source files are read-only and are never filtered or rewritten.
"""
from __future__ import annotations
import argparse,csv,os,shutil,subprocess,tempfile
from pathlib import Path

IDS={
 'zenodo_14011550_weight_loss','zenodo_14011550_exhaustion',
 'zenodo_14011550_low_physical_activity','zenodo_14011550_slow_walking_speed',
 'zenodo_14011550_low_grip_strength','zenodo_1302861_healthspan',
 'nalls_2019_parkinson_public_proxy',
}
FIELDS=['resource_id','file','bytes','sha256','rows','coordinate_order','coordinate_order_descents','duplicate_key_groups','duplicate_excess_rows','global_key_uniqueness','key_definition','method']
ALIASES={'chromosome':('chromosome','chr','chrom'),'position':('base_pair_location','position','pos','bp','genpos'),'effect_allele':('effect_allele','ea','a1','allele1'),'other_allele':('other_allele','ra','nea','oa','a0','a2','allele0')}

def columns(path:Path, resource_id:str):
 delim=',' if path.name.lower().endswith(('.csv','.csv.gz')) else '\t'
 opener=__import__('gzip').open if path.suffix.lower()=='.gz' else open
 with opener(path,'rt',encoding='utf-8-sig',newline='') as f: header=next(csv.reader(f,delimiter=delim))
 lookup={x.strip().lower():i+1 for i,x in enumerate(header)}; out=[]
 for key in ('chromosome','position','effect_allele','other_allele'):
  ix=next((lookup[x] for x in ALIASES[key] if x in lookup),None)
  if ix is None:raise ValueError(f'{resource_id}: missing {key} column; header={header}')
  out.append(ix)
 return delim,out

def decompressor(path:Path): return ['gzip','-cd',str(path)] if path.suffix.lower()=='.gz' else ['cat',str(path)]

def stream_scan(path:Path, delim:str, ix:list[int]):
 c,p,ea,oa=ix; fs=r'\t' if delim=='\t' else ','
 program=(r'BEGIN { FS="__FS__"; OFS="\t" } '
  r'function norm(s,u) { u=toupper(s); sub(/^CHR/,"",u); if (u ~ /^[0-9]+$/) return sprintf("%d",u+0); return u } '
  r'function rank(s,u) { u=norm(s); if (u ~ /^[0-9]+$/) return u+0; if (u=="X") return 23; if (u=="Y") return 24; if (u=="M" || u=="MT") return 25; return 1000 } '
  r'NR==1 { next } { rows++; chr=norm($CH); cr=rank($CH); pos=$POS+0; '
  r'if (rows>1) { if (cr<prevcr || (cr==prevcr && cr<1000 && pos<prevpos) || (cr==prevcr && cr>=1000 && (chr<prevchr || (chr==prevchr && pos<prevpos)))) desc++ } '
  r'site=chr SUBSEP sprintf("%.0f",pos); if (site!=prevsite) { delete seen; prevsite=site } '
  r'key=toupper($EA) SUBSEP toupper($OA); if (seen[key]++) { excess++; if (seen[key]==2) groups++ } '
  r'prevcr=cr; prevchr=chr; prevpos=pos } END { printf "%d\t%d\t%d\t%d\n", rows,desc,groups,excess }')
 program=(program.replace('__FS__',fs).replace('$CH',f'${c}').replace('$POS',f'${p}').replace('$EA',f'${ea}').replace('$OA',f'${oa}'))
 dec=subprocess.Popen(decompressor(path),stdout=subprocess.PIPE,stderr=subprocess.PIPE)
 awk=subprocess.Popen(['awk',program],stdin=dec.stdout,stdout=subprocess.PIPE,stderr=subprocess.PIPE,env={**os.environ,'LC_ALL':'C'})
 assert dec.stdout is not None; dec.stdout.close()
 output,awk_err=awk.communicate(); dec_err=dec.stderr.read() if dec.stderr else b''; dec_code=dec.wait()
 if awk.returncode or dec_code:raise RuntimeError(f'awk={awk.returncode} gzip={dec_code}: {awk_err.decode(errors="replace")} {dec_err.decode(errors="replace")}')
 vals=[int(v) for v in output.decode().strip().split('\t')]
 return vals

def sorted_fallback(path:Path, delim:str, ix:list[int], external_root:Path):
 c,p,ea,oa=ix;fs=r'\t' if delim=='\t' else ','
 program=(r'BEGIN { FS="__FS__"; OFS="\t" } NR>1 { print norm($CH), sprintf("%.0f",$POS+0), toupper($EA), toupper($OA) } '
  r'function norm(s,u) { u=toupper(s); sub(/^CHR/,"",u); if (u ~ /^[0-9]+$/) return sprintf("%d",u+0); return u }'
  .replace('__FS__',fs).replace('$CH',f'${c}').replace('$POS',f'${p}').replace('$EA',f'${ea}').replace('$OA',f'${oa}'))
 tmp=Path(tempfile.mkdtemp(prefix='acquired-gwas-sort-',dir=external_root))
 dec=awk=sorter=None
 try:
  dec=subprocess.Popen(decompressor(path),stdout=subprocess.PIPE,stderr=subprocess.PIPE)
  awk=subprocess.Popen(['awk',program],stdin=dec.stdout,stdout=subprocess.PIPE,stderr=subprocess.PIPE,env={**os.environ,'LC_ALL':'C'});assert dec.stdout is not None;dec.stdout.close()
  sorter=subprocess.Popen(['sort','-S','1G','-T',str(tmp),'-t','\t','-k1,1','-k2,2n','-k3,3','-k4,4'],stdin=awk.stdout,stdout=subprocess.PIPE,stderr=subprocess.PIPE,env={**os.environ,'LC_ALL':'C'});assert awk.stdout is not None;awk.stdout.close();assert sorter.stdout is not None
  rows=groups=excess=0;prev=prev_group=None
  for raw in sorter.stdout:
   key=raw.rstrip(b'\n');rows+=1
   if key==prev:
    excess+=1
    if prev_group!=key:groups+=1
   prev_group=key;prev=key
  sorter.stdout.close(); sort_err=sorter.stderr.read() if sorter.stderr else b'';awk_err=awk.stderr.read() if awk.stderr else b'';dec_err=dec.stderr.read() if dec.stderr else b''
  sc=sorter.wait();ac=awk.wait();dc=dec.wait()
  if sc or ac or dc:raise RuntimeError(f'awk={ac} gzip={dc} sort={sc}: {awk_err.decode(errors="replace")} {dec_err.decode(errors="replace")} {sort_err.decode(errors="replace")}')
  return rows,groups,excess
 finally:
  for proc in (sorter,awk,dec):
   if proc is not None and proc.poll() is None:proc.terminate();proc.wait()
  shutil.rmtree(tmp,ignore_errors=True)

def main():
 ap=argparse.ArgumentParser();ap.add_argument('--repo',type=Path,default=Path(__file__).resolve().parents[2]);ap.add_argument('--external-storage-root',type=Path,required=True);ap.add_argument('--output',type=Path,default=Path('frailty_paper/manifests/gwas_component_duplicate_audit.tsv'));ap.add_argument('--resource-ids',nargs='+',choices=sorted(IDS),default=sorted(IDS));a=ap.parse_args();repo=a.repo.resolve();ext=a.external_storage_root.resolve();selected=set(a.resource_ids)
 if shutil.disk_usage(ext).free<20*1024**3:raise SystemExit('external storage has less than 20 GiB free; no audit started')
 items={r['resource_id']:r for r in csv.DictReader((repo/'frailty_paper/manifests/all_acquired_resources.tsv').open(),delimiter='\t') if r['resource_type']=='GWAS summary statistics' and r['resource_id'] in IDS}
 if set(items)!=IDS:raise SystemExit(f'manifest resource set mismatch: missing={sorted(IDS-set(items))}')
 dest=a.output if a.output.is_absolute() else repo/a.output;dest.parent.mkdir(parents=True,exist_ok=True)
 retained={}
 if dest.exists():
  with dest.open(newline='') as prior:
   retained={r['resource_id']:r for r in csv.DictReader(prior,delimiter='\t') if r.get('resource_id') not in selected}
 with dest.open('w',newline='') as f:
  writer=csv.DictWriter(f,fieldnames=FIELDS,delimiter='\t',lineterminator='\n');writer.writeheader()
  for row in retained.values():writer.writerow(row)
  f.flush()
  for rid in sorted(selected):
   item=items[rid];path=repo/item['file']
   if not path.is_file():raise SystemExit(f'missing input: {path}')
   delim,ix=columns(path,rid);rows,descents,groups,excess=stream_scan(path,delim,ix);method='streaming coordinate-group key count'
   if descents:
    sorted_rows,groups,excess=sorted_fallback(path,delim,ix,ext);rows=sorted_rows;method='input-order descents; global external sort key count'
   result={'resource_id':rid,'file':item['file'],'bytes':item['bytes'],'sha256':item['sha256'],'rows':rows,'coordinate_order':'NONDECREASING' if not descents else 'DESCENTS_REPORTED','coordinate_order_descents':descents,'duplicate_key_groups':groups,'duplicate_excess_rows':excess,'global_key_uniqueness':'UNIQUE' if not excess else 'DUPLICATE_KEYS_PRESENT','key_definition':'chromosome + numeric position + uppercase effect allele + uppercase other allele','method':method+'; no filtering'}
   writer.writerow(result);f.flush();print(f"{rid} rows={rows} descents={descents} duplicate_groups={groups} excess={excess}",flush=True)
 print(f'ACQUIRED_GWAS_DUPLICATE_AUDIT_WRITTEN rows={len(selected)} output={dest}')
if __name__=='__main__':main()
