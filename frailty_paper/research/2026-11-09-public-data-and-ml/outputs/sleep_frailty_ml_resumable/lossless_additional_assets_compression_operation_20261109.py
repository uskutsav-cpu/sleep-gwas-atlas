import pathlib,subprocess,hashlib,json,datetime,shutil,mmap,os,shlex
root=pathlib.Path.cwd();out=root/'outputs/sleep_frailty_ml_resumable/lossless_additional_assets_filesystem_compression_20261109.json';assert not out.exists()
raw=subprocess.run(['ps','-axo','pid,command'],capture_output=True,text=True,check=True).stdout
names={'run_resumable_comparison.py','readonly_native_runtime.py','readonly_native_runtime_768.py','train_permitted_AVI_development.py','bench_permitted_AVI_resources.py','audit_completed_readonly_fit.py'}
for line in raw.splitlines()[1:]:
 f=line.strip().split(None,1)
 if len(f)!=2:continue
 try:a=shlex.split(f[1])
 except ValueError:continue
 if a and pathlib.Path(a[0]).name.startswith('python'):assert not any(pathlib.Path(x).name in names for x in a[1:]),line
explicit=['work/stage4_magma/ensembl_v92_1kb.genes.annot','work/gepitel/GEPITEL.bim','work/e2g_portal_bundle.js','work/PsychAD_snTWAS_EUR_Class_Mural.db','work/PsychAD_snTWAS_EUR_Class_Endo.db','work/stage6_tenk10k_susie_expanded/ld_f64.bin','work/niagads_smad5_susie/ld_f64.bin','work/stage5_susie_diagnostic/eur_ld_float64.bin','work/stage6_tenk10k_susie_diagnostic/ld_f64.bin']
explicit += [str(p) for p in pathlib.Path('work/stage4_magma').glob('*.genes.raw')]
paths=[root/p for p in explicit]
for name in subprocess.run(['rg','--files','work','outputs'],capture_output=True,text=True).stdout.splitlines():
 p=root/name
 if p.suffix.lower() not in {'.csv','.tsv','.json','.jsonl','.txt','.bed'} or p.is_symlink() or any(x in p.parts for x in ['site-packages','lib','node_modules']):continue
 if any(x in name.lower() for x in ['reserved','ulirsch','external_outcome']):continue
 s=p.stat()
 if 1024**2<=s.st_size<4*1024**2 and s.st_blocks*512>s.st_size*.8 and s.st_nlink==1:paths.append(p)
paths=sorted(set(paths),key=lambda p:p.stat().st_size)
result={'time_UTC':datetime.datetime.now(datetime.timezone.utc).isoformat(),'operation':'Lossless transparent compression of closed completed project assets; identical paths and bytes.','scientific_training_active_during_operation':False,'compression_operation_minimum_reserve_bytes':128*1024**2,'compression_reserve_justification':'No fitting or model loading is active. Full uncompressed temporary-file size plus 128 MiB must fit before each filesystem copy. Scientific 768 MiB admission and 300 MiB live-save guards are unchanged.','free_disk_before':shutil.disk_usage(root).free,'files':[],'data_paths_and_logical_contents_preserved':True,'installed_libraries_models_or_user_documents_deleted':False}
def sha(p):
 h=hashlib.sha256()
 with p.open('rb') as f:
  for b in iter(lambda:f.read(1024**2),b''):h.update(b)
 return h.hexdigest()
def write():out.write_text(json.dumps(result,indent=2)+'\n')
write()
for p in paths:
 s=p.stat();tmp=p.with_name(p.name+'.filesystem_compression_pending');assert not tmp.exists()
 if p.is_symlink() or s.st_nlink!=1:continue
 opened=subprocess.run(['/usr/sbin/lsof','-F','p','--',str(p)],capture_output=True,text=True)
 if opened.returncode!=1 or opened.stdout.strip() or any(p.with_name(p.name+x).exists() for x in ['-wal','-journal']):
  result['files'].append({'path':str(p.relative_to(root)),'skipped':'Open file or SQLite journal present'});write();continue
 if shutil.disk_usage(root).free<s.st_size+128*1024**2:
  result['files'].append({'path':str(p.relative_to(root)),'skipped':'Full temporary-copy reserve not available'});write();continue
 digest=sha(p);subprocess.run(['/usr/bin/ditto','--hfsCompression','--noclone','--nocache',str(p),str(tmp)],check=True,capture_output=True,text=True)
 assert tmp.stat().st_size==s.st_size and sha(tmp)==digest
 with p.open('rb') as a,tmp.open('rb') as b,mmap.mmap(a.fileno(),0,access=mmap.ACCESS_READ) as x,mmap.mmap(b.fileno(),0,access=mmap.ACCESS_READ) as y:
  for off in [0,max(0,s.st_size//2-2048),max(0,s.st_size-4096)]:assert x[off:off+4096]==y[off:off+4096]
 now=p.stat();assert (now.st_ino,now.st_size,now.st_mtime_ns)==(s.st_ino,s.st_size,s.st_mtime_ns) and sha(p)==digest
 applied=tmp.stat().st_blocks*512<s.st_blocks*512
 if applied:os.utime(tmp,ns=(s.st_atime_ns,s.st_mtime_ns));os.replace(tmp,p)
 else:tmp.unlink()
 row={'path':str(p.relative_to(root)),'logical_bytes_unchanged':s.st_size,'SHA256_unchanged':digest,'mmap_samples_match':True,'compression_applied':applied,'allocated_before':s.st_blocks*512,'allocated_after':p.stat().st_blocks*512,'same_name_format_mtime_preserved':p.stat().st_mtime_ns==s.st_mtime_ns};assert row['same_name_format_mtime_preserved'];result['files'].append(row);result['free_disk_current']=shutil.disk_usage(root).free;write();print(json.dumps({'path':row['path'],'applied':applied,'free_bytes':result['free_disk_current']}),flush=True)
result['completed']=True;result['free_disk_after']=shutil.disk_usage(root).free;result['completed_time_UTC']=datetime.datetime.now(datetime.timezone.utc).isoformat();write();print(json.dumps({'completed':True,'file_count':len(result['files']),'free_before':result['free_disk_before'],'free_after':result['free_disk_after']}),flush=True)
