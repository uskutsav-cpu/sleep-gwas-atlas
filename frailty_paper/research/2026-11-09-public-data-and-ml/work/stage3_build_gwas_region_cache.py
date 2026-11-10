import csv,gzip,pathlib,collections,json,hashlib,time
rows=list(csv.DictReader(open('outputs/shared22_eight_source_evidence_matrix.csv')));extra=next(x for x in csv.DictReader(open('outputs/integrated_15_candidate_review.csv')) if x['variant']=='rs323509')
coords=[(x['chr_GRCh37_assumed'],int(x['position_GRCh37_assumed']),x['chr_GRCh38'],int(x['position_GRCh38']),x['variant']) for x in rows]+[(extra['GRCh37_chr'],int(extra['GRCh37_pos']),extra['GRCh37_chr'],int(extra['GRCh38_pos']),'rs323509')]
assert len(coords)==23
windows={37:collections.defaultdict(list),38:collections.defaultdict(list)}
for c37,p37,c38,p38,name in coords:
 windows[37][c37].append((max(1,p37-2500000),p37+2500000,name));windows[38][c38].append((max(1,p38-2500000),p38+2500000,name))
root=pathlib.Path('/Volumes/Extreme SSD/sleep-gwas-atlas-frailty-v1/analysis-workspace/frailty_v1/harmonized');mak=pathlib.Path('/Volumes/Extreme SSD/Codex Archive/mak2023_same_cohort');base=pathlib.Path('/Volumes/Extreme SSD/sleep-frailty-full-autonomous-v5-20261009/frailty_paper/full_autonomous_research_v5/genetics/new_public_validation/raw');sources=[('Atkins_FI',root/'frailty.harmonized.tsv.gz','CHR','BP',37),('Jansen_insomnia',root/'insomnia.harmonized.tsv.gz','CHR','BP',37),('Mak_FI',mak/'fi.fastGWA','CHR','POS',37),('Mak_FP',mak/'fp.fastGWA','CHR','POS',37),('Ye_FFS',pathlib.Path('/Volumes/Extreme SSD/sleep-gwas-atlas-frailty-v1/frailty_paper/data/gwas/fried_frailty_score/Frailty_GWAS_maf0.01_info0.9_hwe1e-6_FUMA.tsv'),'CHROM','POS',37),('FinnGen_insomnia',base/'finngen_R13_F5_INSOMNIA.gz','#chrom','pos',38)];receipt={}
for name,path,chrom_col,pos_col,build in sources:
 start=time.time();out=pathlib.Path('work/stage3_region_cache')/(name+'.regional.tsv.gz');n=kept=0
 with (gzip.open(path,'rt') if path.suffix=='.gz' else path.open()) as f,gzip.open(out,'wt',compresslevel=1) as o:
  header=f.readline();o.write(header);cols=header.rstrip('\n').split('\t');ci=cols.index(chrom_col);pi=cols.index(pos_col)
  for line in f:
   n+=1;v=line.split('\t');chrom=v[ci];regions=windows[build].get(chrom)
   if not regions:continue
   try:p=int(v[pi])
   except ValueError:continue
   if any(lo<=p<=hi for lo,hi,_ in regions):o.write(line);kept+=1
 h=hashlib.sha256()
 with out.open('rb') as f:
  while b:=f.read(8*1024*1024):h.update(b)
 receipt[name]={'source':str(path),'source_build_assumed':build,'source_rows':n,'retained_rows':kept,'output':str(out),'output_sha256':h.hexdigest(),'seconds':round(time.time()-start,1)}
 print(name,receipt[name],flush=True)
pathlib.Path('outputs/stage3_gwas_regional_cache_receipt.json').write_text(json.dumps({'radius_bp':2500000,'marker_windows':coords,'traits':receipt,'limitations':'candidate coordinates and source builds assumed; this cache preserves all source rows within any window without allele QC, not independently harmonized or a discovery test'},indent=2)+'\n')
