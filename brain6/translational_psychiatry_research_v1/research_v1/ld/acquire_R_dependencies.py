#!/usr/bin/env python3
import csv,hashlib,json,pathlib,re,subprocess
HERE=pathlib.Path(__file__).resolve().parent
RAW=HERE.parents[3]/'work'/'R_source_packages_research_v1'
RAW.mkdir(parents=True,exist_ok=True)
rows=list(csv.DictReader((HERE/'R_dependency_manifest.tsv').open(),delimiter='\t'))
manifest={'before_native_fits':True,'source_index_sha256':hashlib.sha256((HERE.parents[3]/'work/ld_genotypes_research_v1/CRAN_SOURCE_PACKAGES.gz').read_bytes()).hexdigest(),'combined_compressed_cap':50000000,'installed_library_disk_cap_bytes':150*2**20,'sources':[{'package':r['Package'],'version':r['Version'],'source_md5':r['MD5sum'],'url':f'https://cran.r-project.org/src/contrib/{r["Package"]}_{r["Version"]}.tar.gz','per_object_cap':12000000} for r in rows]}
(HERE/'R_full_dependency_acquisition_manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
receipts=[]
for obj in manifest['sources']:
 path=RAW/(obj['package']+'_'+obj['version']+'.tar.gz')
 result=subprocess.run(['curl','--fail','--silent','--show-error','--location','--max-time','40','--max-filesize',str(obj['per_object_cap']),'-o',str(path),obj['url']],capture_output=True,text=True)
 rec={**obj,'path':str(path),'returncode':result.returncode,'stderr':result.stderr}
 if path.exists():
  body=path.read_bytes();rec.update(bytes=len(body),sha256=hashlib.sha256(body).hexdigest(),actual_md5=hashlib.md5(body).hexdigest())
 receipts.append(rec);(HERE/'R_full_dependency_acquisition_receipts.json').write_text(json.dumps(receipts,indent=2)+'\n')
 if result.returncode or rec['actual_md5']!=obj['source_md5']:raise RuntimeError('Exact package acquisition failed '+json.dumps(rec))
 if sum(x.get('bytes',0) for x in receipts)>manifest['combined_compressed_cap']:raise RuntimeError('Source package cap exceeded')
 print(obj['package'],rec['bytes'],flush=True)
# Topological source-install order, excluding existing base/recommended dependencies.
by={r['Package']:r for r in rows};ordered=[];visiting=set()
def visit(name):
 if name in ordered:return
 if name in visiting:raise ValueError('Dependency cycle')
 visiting.add(name)
 for col in ['Depends','Imports','LinkingTo']:
  for x in by[name][col].split(','):
   d=re.sub(r'\s*\(.*','',x).strip()
   if d in by:visit(d)
 visiting.remove(name);ordered.append(name)
for n in by:visit(n)
(HERE/'R_install_source_order.txt').write_text('\n'.join(str(RAW/(n+'_'+by[n]['Version']+'.tar.gz')) for n in ordered)+'\n')
print(json.dumps({'packages':len(receipts),'source_bytes':sum(x['bytes'] for x in receipts)}),flush=True)
