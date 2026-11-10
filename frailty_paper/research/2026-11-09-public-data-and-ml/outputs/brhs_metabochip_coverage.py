#!/usr/bin/env python3
"""Compare prototype coordinates with the public Illumina mixed-build manifest.
Usage: python3 outputs/brhs_metabochip_coverage.py MANIFEST.csv hg19ToHg18.over.chain.gz
Requires pyliftover. The output identifies array loci, not cohort genotypes.
"""
import csv, gzip, hashlib, json, sys
from pathlib import Path
from pyliftover import LiftOver

manifest, chain = map(Path, sys.argv[1:3])
def sha(path):
    h = hashlib.sha256()
    with path.open('rb') as f:
        for chunk in iter(lambda: f.read(1024*1024), b''): h.update(chunk)
    return h.hexdigest()
lo=LiftOver(str(chain))
scores=[('15','outputs/unvalidated_15_variant_score_prototype.csv','chr','position_GRCh37_inferred'),('22','outputs/unvalidated_22_shared_direction_score.csv','chr_GRCh37_assumed','position_GRCh37_assumed')]
variants=[]
for score,path,chrom,pos in scores:
    with open(path,newline='') as f:
        for x in csv.DictReader(f):
            conversions=lo.convert_coordinate('chr'+x[chrom],int(x[pos])-1)
            variants.append(dict(score=score,variant=x['variant'],chr37=x[chrom],pos37=x[pos],mappings=conversions))
with manifest.open(newline='',encoding='utf-8-sig') as f:
    for line in f:
        if line.strip()=='[Assay]': break
    reader=csv.DictReader(f)
    bycoord={}; byname={}; builds={}
    for x in reader:
        if not x.get('Name'): continue
        build=x.get('GenomeBuild',''); builds[build]=builds.get(build,0)+1
        
        if build in ('36','36.2','18'): key='hg18'
        elif build in ('37','37.1'): key='hg19'
        else: key='other'
        bycoord.setdefault((key,str(x.get('Chr')),str(x.get('MapInfo'))),[]).append(x['Name'])
        byname.setdefault(x['Name'],[]).append((x.get('Chr'),x.get('MapInfo')))
rows=[]
for x in variants:
    mappings=x['mappings']; coordnames=[]
    if len(mappings)==1:
        chrom,p0,strand,quality=mappings[0]
        coordnames=bycoord.get(('hg18',chrom.removeprefix('chr'),str(p0+1)),[])
        coord=f'{chrom}:{p0+1}:{strand}'
    else: coord=''
    coordnames37=bycoord.get(('hg19',x['chr37'],x['pos37']),[])
    coordnames=list(dict.fromkeys(coordnames+coordnames37))
    rows.append(dict(score=x['score'],variant=x['variant'],chr37=x['chr37'],pos37=x['pos37'],hg18_liftover=coord,mapping_count=len(mappings),exact_name_match=bool(byname.get(x['variant'])),build_matched_coordinate_match=bool(coordnames),manifest_names=';'.join(coordnames)))
out=Path('outputs/brhs_metabochip_score_marker_screen.csv')
with out.open('w',newline='') as f:
    w=csv.DictWriter(f,fieldnames=rows[0].keys()); w.writeheader(); w.writerows(rows)
receipt={'manifest_sha256':sha(manifest),'chain_sha256':sha(chain),'manifest_build_counts':builds,'manifest_source':'https://webdata.illumina.com/downloads/productfiles/cardio-metabochip/cardio-metabo_chip_11395247_c.csv','chain_source':'https://hgdownload.soe.ucsc.edu/goldenPath/hg19/liftOver/hg19ToHg18.over.chain.gz','scores':{score:{'rows':sum(x['score']==score for x in rows),'unique_liftover':sum(x['score']==score and x['mapping_count']==1 for x in rows),'coordinate_matches':sum(x['score']==score and x['build_matched_coordinate_match'] for x in rows),'exact_name_matches':sum(x['score']==score and x['exact_name_match'] for x in rows)} for score in ['15','22']}}
Path('outputs/brhs_metabochip_score_marker_screen.json').write_text(json.dumps(receipt,indent=2)+'\n')
print(json.dumps(receipt,indent=2))
