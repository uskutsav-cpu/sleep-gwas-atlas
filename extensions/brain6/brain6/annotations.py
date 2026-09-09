"""PIP-weighted regulatory overlap with reviewed BED0 cell annotations.

This produces descriptive independent-locus scores for the competitive test.
A weighted overlap is NOT evidence that a cell type mediates disease risk.
"""
from __future__ import annotations
import bisect
import collections
import math
from pathlib import Path
from .artifacts import transaction
from .io import check_hash, open_text, read_json, read_tsv, require, write_json, write_tsv


class IntervalIndex:
    """Merge overlapping BED intervals so SNPs are never counted twice."""
    def __init__(self,path):
        intervals=collections.defaultdict(list)
        self.skipped_non_autosomal=0
        with open_text(path) as f:
            for line in f:
                if not line.strip() or line.startswith(('#','track','browser')):
                    continue
                row=line.rstrip('\n').split('\t')
                require(len(row)>=3,'BED requires chromosome/start/end')
                contig=row[0].lower().removeprefix('chr')
                if not contig.isdigit() or int(contig) not in range(1,23):
                    self.skipped_non_autosomal+=1;continue
                ch=int(contig);start,stop=int(row[1]),int(row[2])
                require(0<=start<stop,'Invalid BED0 half-open interval')
                if ch in range(1,23):intervals[ch].append((start,stop))
        self.intervals={};self.starts={}
        for ch,values in intervals.items():
            merged=[]
            for start,stop in sorted(values):
                if merged and start<=merged[-1][1]:
                    merged[-1]=(merged[-1][0],max(stop,merged[-1][1]))
                else:merged.append((start,stop))
            self.intervals[ch]=merged;self.starts[ch]=[r[0] for r in merged]
        require(bool(self.intervals),'No usable autosomal BED annotations')

    def contains(self,ch,bp1):
        require(bp1>=1,'GWAS positions are one-based')
        bp0=bp1-1
        i=bisect.bisect_right(self.starts.get(ch,[]),bp0)-1
        return i>=0 and bp0<self.intervals[ch][i][1]


def score(variants,annotation_manifest,root,name,*,synthetic=False):
    c=read_json(annotation_manifest)
    require(c.get('reviewed') is True and c.get('synthetic',False)==synthetic,'Unreviewed/mixed annotation manifest')
    require(c.get('genome_build') in {'GRCh37','GRCh38'},'Declare locus/annotation build')
    require(c.get('independent_loci_reviewed') is True,'Do not permute individual SNPs as independent loci')
    require(c.get('matched_controls_reviewed') is True,'Control locus matching needs review')
    require(c.get('annotations'),'No cell annotations')
    cells={};inputs=[variants,annotation_manifest]
    for a in c['annotations']:
        require(a.get('cell_type') and a['cell_type'] not in cells,'Duplicate/empty cell label')
        require(a.get('source_uri') and a.get('citation'),'Unsourced cell annotations')
        require(a.get('genome_build')==c['genome_build'] and a.get('coordinate_system')=='BED0',
                'No implicit genome-build or coordinate conversion')
        check_hash(a['path'],a['sha256']);inputs.append(a['path'])
        cells[a['cell_type']]=IntervalIndex(a['path'])
    units={};seen=set()
    for r in read_tsv(variants,['pair_id','locus_id','SNP','CHR','BP','PIP','selected','stratum']):
        key=r['pair_id'],r['locus_id'];snp_key=key+(r['SNP'],)
        require(snp_key not in seen,'Duplicate variant within locus score')
        seen.add(snp_key)
        require(r['selected'] in {'0','1'} and r['stratum'],'Invalid selected flag/matching stratum')
        pip=float(r['PIP']);ch,bp=int(r['CHR']),int(r['BP'])
        require(math.isfinite(pip) and 0<=pip<=1 and ch in range(1,23) and bp>=1,'Invalid variant/PIP')
        if key not in units:
            units[key]={'selected':r['selected'],'stratum':r['stratum'],'total':0.,'variants':0,
                        'mass':{cell:0. for cell in cells}}
        u=units[key]
        require(u['selected']==r['selected'] and u['stratum']==r['stratum'],'Inconsistent locus membership')
        u['total']+=pip;u['variants']+=1
        for cell,index in cells.items():
            if index.contains(ch,bp):u['mass'][cell]+=pip
    require(units,'No locus variants')
    rows=[]
    for (pair,locus),u in sorted(units.items()):
        require(u['total']>0,'Zero posterior mass is unresolved, not a zero-overlap result')
        for cell in sorted(cells):
            rows.append({'pair_id':pair,'cell_type':cell,'locus_id':locus,
                         'score':u['mass'][cell]/u['total'],'selected':u['selected'],
                         'stratum':u['stratum'],'n_variants':u['variants'],'pip_sum':u['total']})
    with transaction(root,name,stage='cell_annotation_scores',inputs=inputs,
                     parameters=c,synthetic=synthetic) as (work,meta):
        write_tsv(work/'scores.tsv',list(rows[0]),rows)
        meta['scientific_status']='PASS'
        write_json(work/'status.json',{'status':'PASS','loci':len(units),'cell_types':len(cells),
                   'skipped_non_autosomal_intervals':{k:v.skipped_non_autosomal for k,v in cells.items()},
                   'note':'Descriptive fraction of PIP mass overlapping annotations; not a cell-specific causal probability. Competitive enrichment requires reviewed matched control loci.'})
    return Path(root)/name
