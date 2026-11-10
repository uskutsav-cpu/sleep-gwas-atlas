#!/usr/bin/env python3
"""Join GRCh37 source reporter variants to prevalidated reciprocal mappings."""
import argparse,csv,gzip,pathlib,collections
p=argparse.ArgumentParser();p.add_argument('--directory',required=True);a=p.parse_args();root=pathlib.Path('outputs');source=pathlib.Path(a.directory)
with (root/'stage5_reciprocal_snv_hg19_reference_check.csv').open() as f:
 mapping={f"{r['chr37']}:{r['pos37']}:{r['ref38']}:{r['alt38']}":r['variant_GRCh38'] for r in csv.DictReader(f) if r['reference_status']=='reference_same'}
with (root/'stage5_indel_hg19_mapping.csv').open() as f:
 for r in csv.DictReader(f):
  if r['status']=='mapped_reference_checked_normalized':mapping[f"{r['hg19_chrom']}:{r['hg19_position']}:{r['hg19_ref']}:{r['hg19_alt']}"]=r['variant_GRCh38']
files=sorted(source.glob('*_SKNSH_emVAR_glm_paired_20221220.out.gz'));files=[x for x in files if not x.name.startswith('._')];assert len(files)==9,len(files)
fields=['variant_GRCh38','variant_GRCh37_source','source_file','source_id','cell_line','A_log2FC','A_logPadj_BF','B_log2FC','B_logPadj_BF','Log2Skew','Skew_SE','Skew_logP','Skew_logFDR','Skew_logFDR_act','note'];matches=[];counts={}
for file in files:
 n=0;hits=0
 with gzip.open(file,'rt',newline='') as f:
  for x in csv.DictReader(f,delimiter='\t'):
   n+=1;key=mapping.get(x['SNP'])
   if not key:continue
   hits+=1;matches.append({'variant_GRCh38':key,'variant_GRCh37_source':x['SNP'],'source_file':file.name,'source_id':x['ID'],'cell_line':'SK-N-SH','A_log2FC':x['A_log2FC'],'A_logPadj_BF':x['A_logPadj_BF'],'B_log2FC':x['B_log2FC'],'B_logPadj_BF':x['B_logPadj_BF'],'Log2Skew':x['Log2Skew'],'Skew_SE':x['Skew_SE'],'Skew_logP':x['Skew_logP'],'Skew_logFDR':x['Skew_logFDR'],'Skew_logFDR_act':x['Skew_logFDR_act'],'note':'GRCh37 exact allele labels mapped via reciprocal/reference-checked project table; verify allele direction and source significance conventions before comparison'})
 counts[file.name]=(n,hits)
with (root/'stage5_siraj2026_neural_mpra_exact_alleles.csv').open('w',newline='') as f:
 w=csv.DictWriter(f,fieldnames=fields);w.writeheader();w.writerows(matches)
print('files',len(files),'reference-mapped variants',len(mapping),'match rows',len(matches),'distinct variants',len({x['variant_GRCh38'] for x in matches}));print('file counts',counts)
