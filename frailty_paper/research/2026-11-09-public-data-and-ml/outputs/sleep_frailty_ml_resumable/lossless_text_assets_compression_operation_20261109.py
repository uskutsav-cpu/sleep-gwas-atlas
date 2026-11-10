import pathlib,subprocess,hashlib,json,datetime,shutil,mmap,os
root=pathlib.Path.cwd();out=root/'outputs/sleep_frailty_ml_resumable/lossless_text_assets_filesystem_compression_20261109.json'
assert not out.exists(),'Preserve prior operation receipt'
raw=subprocess.run(['ps','-axo','command'],capture_output=True,text=True,check=True).stdout
assert not any('python' in row and ('readonly_native_runtime' in row or 'train_permitted_AVI_development.py --' in row or 'run_resumable_comparison.py' in row) for row in raw.splitlines()),'Scientific worker or controller active'
allpaths=subprocess.run(['rg','--files','work','outputs','data'],capture_output=True,text=True).stdout.splitlines()
paths=[]
for name in allpaths:
 p=root/name
 if p.suffix.lower() not in {'.csv','.tsv','.json','.jsonl','.txt','.bed','.vcf'} or p.is_symlink():continue
 s=p.stat()
 if s.st_size>=4*1024**2 and s.st_blocks*512>s.st_size*.8 and s.st_nlink==1:paths.append(p)
paths.sort(key=lambda p:p.stat().st_size)
def sha(p):
 h=hashlib.sha256()
 with p.open('rb') as f:
  for b in iter(lambda:f.read(1024**2),b''):h.update(b)
 return h.hexdigest()
result={'time_UTC':datetime.datetime.now(datetime.timezone.utc).isoformat(),'operation':'ditto transparent lossless filesystem compression of completed project text assets','free_disk_before':shutil.disk_usage(root).free,'scientific_training_active_during_operation':False,'files':[],'data_paths_and_logical_contents_preserved':True,'user_documents_or_models_deleted':False}
def write():out.write_text(json.dumps(result,indent=2)+'\n')
write()
for p in paths:
 s=p.stat();tmp=p.with_name(p.name+'.filesystem_compression_pending')
 assert not tmp.exists(),str(tmp)
 opened=subprocess.run(['/usr/sbin/lsof','-F','p','--',str(p)],capture_output=True,text=True)
 if opened.returncode!=1 or opened.stdout.strip():
  result['files'].append({'path':str(p.relative_to(root)),'skipped':'File open or lsof unable to establish closed status'});write();continue
 if shutil.disk_usage(root).free<s.st_size+300*1024**2:
  result['files'].append({'path':str(p.relative_to(root)),'skipped':'Insufficient full temporary-copy reserve'});write();continue
 digest=sha(p)
 subprocess.run(['/usr/bin/ditto','--hfsCompression','--noclone','--nocache',str(p),str(tmp)],check=True,capture_output=True,text=True)
 assert tmp.stat().st_size==s.st_size and sha(tmp)==digest
 with p.open('rb') as a,tmp.open('rb') as b,mmap.mmap(a.fileno(),0,access=mmap.ACCESS_READ) as x,mmap.mmap(b.fileno(),0,access=mmap.ACCESS_READ) as y:
  for off in [0,max(0,s.st_size//2-2048),max(0,s.st_size-4096)]:assert x[off:off+4096]==y[off:off+4096]
 now=p.stat();assert (now.st_ino,now.st_size,now.st_mtime_ns)==(s.st_ino,s.st_size,s.st_mtime_ns) and sha(p)==digest
 allocated_after=tmp.stat().st_blocks*512;applied=allocated_after<s.st_blocks*512
 if applied:
  os.utime(tmp,ns=(s.st_atime_ns,s.st_mtime_ns));os.replace(tmp,p)
 else:tmp.unlink()
 row={'path':str(p.relative_to(root)),'logical_bytes':s.st_size,'SHA256_unchanged':digest,'read_only_mmap_samples_match':True,'compression_applied':applied,'allocated_before':s.st_blocks*512,'allocated_after':p.stat().st_blocks*512,'same_path_format_and_mtime_preserved':p.stat().st_mtime_ns==s.st_mtime_ns}
 assert row['same_path_format_and_mtime_preserved'];result['files'].append(row);result['free_disk_current']=shutil.disk_usage(root).free;write()
 print(json.dumps({'path':row['path'],'applied':applied,'free_bytes':result['free_disk_current']}),flush=True)
result['completed']=True;result['completed_time_UTC']=datetime.datetime.now(datetime.timezone.utc).isoformat();result['free_disk_after']=shutil.disk_usage(root).free;write()
print(json.dumps({'completed':True,'files':len(result['files']),'free_before':result['free_disk_before'],'free_after':result['free_disk_after']}),flush=True)
