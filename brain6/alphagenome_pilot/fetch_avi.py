#!/usr/bin/env python3
"""Join a pinned public AVI mirror to the prepared Brain6 variants by exact allele.

Requires duckdb with httpfs. This queries the permissive static AVI artifact's
third-party re-encoding; it does not use Atlas API output or predict tissues.
"""
import argparse,csv,hashlib,json,decimal,urllib.request,ssl
from pathlib import Path
import duckdb
REV='04820dafc5af625d1d941e080b5ae7d2b17dae87'
BASE=f'https://huggingface.co/datasets/just-dna-seq/alphagenome_avi/resolve/{REV}/'
ROOT=Path(__file__).resolve().parent

def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()
def exact(v):return str(decimal.Decimal(int(v))/decimal.Decimal(100000))
def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--evidence',type=Path,default=ROOT/'results/variant_evidence.tsv');p.add_argument('--out',type=Path,required=True);p.add_argument('--extension-dir',type=Path,required=True);a=p.parse_args()
    with a.evidence.open() as f: rows=list(csv.DictReader(f,delimiter='\t'))
    if not rows:raise ValueError('Empty evidence')
    a.out.mkdir(parents=True,exist_ok=False)
    c=duckdb.connect();c.execute('SET extension_directory = ?', [str(a.extension_dir.resolve())]);c.execute('LOAD httpfs')
    # On the 2026-09-10 mirror, contrary to its README, each row is wide:
    # alt0/1/2 are inferred to be the ordered A,C,G,T alternatives excluding REF.
    # The pre-outcome reference examples and independent lookups below check this.
    def parquet(chrom):return BASE+f'data/alphagenome_avi-{chrom}.parquet'
    def at(chrom,pos):return c.execute('SELECT chrom,pos,ref,alt0,alt1,alt2 FROM read_parquet(?) WHERE pos=?',[parquet(chrom),pos]).fetchall()
    anchors=[('chr1',10001,'T',(-3868,-3200,-3720)),('chr1',10002,'A',(-3583,-3301,-3668)),('chr5',104565723,'C',(-5408,-5977,-5709)),('chr5',104628884,'G',(-41460,-42840,-42260))]
    for chrom,pos,ref,values in anchors:
        got=at(chrom,pos)
        if len(got)!=1 or got[0][2]!=ref or tuple(got[0][3:])!=values:raise ValueError(f'Mirror encoding anchor failed: {chrom}:{pos}')
    knots={int(k):(int(n),float(lo),float(hi)) for k,n,lo,hi in c.execute('SELECT raw_score_e5,n,phred_lo,phred_hi FROM read_parquet(?)',[BASE+'avi_knots.parquet']).fetchall()}
    if not 24.65 < knots[106100][1] < 24.67: raise ValueError('Independent published G>A PHRED anchor mismatch')
    grouped={}
    for row in rows:
        if row['mapping_status']=='READY_UNIQUE_REFERENCE_VERIFIED':grouped.setdefault((row['candidate_locus_id'],row['chromosome']),[]).append(row)
    mapped={}
    for (cid,chrom),group in grouped.items():
        positions=sorted({int(r['pos_grch38_1based']) for r in group}); placeholders=','.join('?' for _ in positions)
        sql=f'SELECT chrom,pos,ref,alt0,alt1,alt2 FROM read_parquet(?) WHERE pos BETWEEN ? AND ? AND pos IN ({placeholders})'
        records=c.execute(sql,[parquet(chrom),min(positions),max(positions),*positions]).fetchall()
        found={}
        for item in records:
            if item[1] in found:raise ValueError(f'Duplicate source position {chrom}:{item[1]}')
            found[item[1]]=item
        for row in group:
            key=(cid,row['snp']);item=found.get(int(row['pos_grch38_1based']))
            if item is None: mapped[key]={'avi_status':'NOT_FOUND_IN_MIRROR'};continue
            ref=row['ref_grch38'];alt=row['alt_grch38']
            if item[0]!=chrom or item[2]!=ref: mapped[key]={'avi_status':'REFERENCE_MISMATCH_HOLD'};continue
            alternatives=[b for b in 'ACGT' if b!=ref]
            if ref not in 'ACGT' or alt not in alternatives:raise ValueError(f'Invalid allele {key}')
            raw=item[3+alternatives.index(alt)]
            if raw is None: mapped[key]={'avi_status':'MISSING_SCORE'};continue
            if int(raw) not in knots: raise ValueError(f'Missing PHRED knot for {raw}')
            _,lo,hi=knots[int(raw)]
            mapped[key]={'avi_status':'SCORED_MIRROR_EXPLORATORY','avi_raw_e5':int(raw),'avi_raw_score_decimal':exact(raw),'avi_phred_lower':lo,'avi_phred_upper':hi}
        print(f'{cid}: {len(group)} variants, {len(found)} positions found',flush=True)
    extra=['avi_status','avi_raw_e5','avi_raw_score_decimal','avi_phred_lower','avi_phred_upper']
    dest=a.out/'variant_evidence_avi.tsv'
    with dest.open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0])+extra,delimiter='\t',lineterminator='\n');w.writeheader()
        for row in rows:w.writerow({**row,**mapped.get((row['candidate_locus_id'],row['snp']),{'avi_status':'MAPPING_HOLD'})})
    ctx=ssl.create_default_context(cafile='/etc/ssl/cert.pem')
    for name in ('release.json','README.md','LICENSE.txt'):
        with urllib.request.urlopen(BASE+name,context=ctx,timeout=30) as response:(a.out/('mirror_'+name)).write_bytes(response.read())
    status={}
    for m in mapped.values():status[m['avi_status']]=status.get(m['avi_status'],0)+1
    receipt={'mirror_revision':REV,'source_url':BASE,'source_dataset':'third-party re-encoding of AlphaGenome AVI SNV static download','source_artifact_sha256_claimed_by_mirror':json.loads((a.out/'mirror_release.json').read_text())['source_sha256'],'source_artifact_not_independently_downloaded':True,'evidence_sha256':sha(a.evidence),'output_sha256':sha(dest),'n_rows':len(rows),'counts':status,'duckdb_version':duckdb.__version__,'layout_discrepancy':'Mirror README describes long-format alt and raw_score_e5; pinned parquet contains alt0/alt1/alt2. Ordering inferred as A,C,G,T excluding reference and checked against cited examples for all four reference bases.','validation_anchors':anchors,'independent_anchor_sources':['https://switzerlandomics.ch/blog/2026-09-09-alphagenome-atlas/','https://github.com/broadinstitute/SpliceAI-lookup/tree/main/google_cloud_run_services/score_lookup','https://www.medrxiv.org/content/10.64898/2026.09.16.26363192v1.full'],'limitations':'Unverified mirror, row-group responses not source-file hash verified; PHRED is interval reconstructed from knots; scores not tissue-specific or causal; do not train genetic variant effect models.'}
    (a.out/'avi_provenance.json').write_text(json.dumps(receipt,indent=2)+'\n')
    print(json.dumps(status))
if __name__=='__main__':main()
