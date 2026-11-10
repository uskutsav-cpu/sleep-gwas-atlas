import tarfile,pathlib,gzip,hashlib,json
archive='/Volumes/Extreme SSD/Utsav-Research-Archive/Sleep-GWAS/FAILED-PARTIAL-2026-09-01-sleep_gwas_atlas/work/track_b_completion/local_dependency_copies/ref/interpretation/cell_type/fuma_scrna/FUMA_scRNA_data-dd526163ea80af1a80a6cdc80db167144500694b.tar.gz'
names=['GSE67835_Human_Cortex_woFetal','PBMC_10x_68k','GSE98816_Mouse_Brain_Vascular','MouseCellAtlas_Muscle','TabulaMuris_FACS_Limb_Muscle','GSE87544_Mouse_Hypothalamus']
out=pathlib.Path('work/stage4_magma/matrices');out.mkdir(parents=True,exist_ok=True);found={}
with tarfile.open(archive,'r|gz') as tar:
 for m in tar:
  for name in names:
   if m.name.endswith('/processed_data/'+name+'.txt.gz'):
    data=tar.extractfile(m).read(); raw=gzip.decompress(data);(out/(name+'.txt')).write_bytes(raw)
    found[name]={'archive_member':m.name,'compressed_sha256':hashlib.sha256(data).hexdigest(),'text_sha256':hashlib.sha256(raw).hexdigest(),'header':raw.splitlines()[0].decode(),'rows':len(raw.splitlines())-1}
assert len(found)==len(names),set(names)-found.keys()
pathlib.Path('outputs/stage4_fuma_matrix_receipt.json').write_text(json.dumps({'archive':archive,'selection':'exploratory after source inventory; different species and datasets','matrices':found},indent=2)+'\n')
for k,v in found.items():print(k,v['rows'],v['header'][:300],flush=True)
