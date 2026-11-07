import hashlib,json,subprocess,zlib
from pathlib import Path
import pysam
p=Path(__file__).resolve().parent
url='https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90565001-GCST90566000/GCST90565837/harmonised/GCST90565837.h.tsv.gz'
chunk=subprocess.check_output(['curl','-fsSL','--range','0-65535','--max-time','45',url])
head=zlib.decompressobj(31).decompress(chunk).split(b'\n',1)[0].decode().rstrip('\r').split('\t')
expected=['chromosome','base_pair_location','effect_allele','other_allele','beta','standard_error','effect_allele_frequency','p_value','rsid','n','hm_coordinate_conversion','hm_code','variant_id']
assert head==expected
idx=p/'GCST90565837.h.tsv.gz.tbi'
f=pysam.TabixFile(url,index=str(idx))
anchors=[('rs62264764',117920728,{'G','A'}),('rs62264767',117923158,{'A','C'})]
output=[]
for rs,pos,alleles in anchors:
 rows=[dict(zip(head,line.split('\t'))) for line in f.fetch('3',pos-1,pos)]
 matching=[r for r in rows if r['rsid']==rs and {r['effect_allele'],r['other_allele']}==alleles]
 output.append({'expected_rsid':rs,'position_grch38':pos,'expected_alleles':sorted(alleles),'all_rows_at_position':rows,'exact_matching_count':len(matching)})
result={'url':url,'publisher_declared_compressed_md5':'204e0d2198afb72af7b42f81317bc0b0','scope':'Partial indexed query only; compressed file MD5 not verified by this query.','index_sha256':hashlib.sha256(idx.read_bytes()).hexdigest(),'pysam_version':pysam.__version__,'rows':output}
(p/'receipt.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps(output,indent=2))
