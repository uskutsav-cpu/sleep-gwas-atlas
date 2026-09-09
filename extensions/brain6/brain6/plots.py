"""Plot measured outputs without altering inference or silently dropping failures."""
from __future__ import annotations
import math
import sqlite3
from pathlib import Path
import numpy as np
from .io import read_tsv,require


def _plt():
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    return plt


def heatmap(ranking,out):
    from .pairs import BRAIN,SLEEP
    rows=list(read_tsv(ranking,['sleep_trait','disease_trait','rg','fdr','analysis_status']))
    require(len(rows)==72,'Need the complete 12 x 6 table, not only significant rows')
    values=np.full((12,6),np.nan);stars=[]
    seen=set()
    for r in rows:
        i=SLEEP.index(r['sleep_trait']);j=BRAIN.index(r['disease_trait'])
        require((i,j) not in seen,'Duplicate heatmap cell');seen.add((i,j))
        try:rg=float(r['rg']);q=float(r['fdr'])
        except ValueError:continue
        if math.isfinite(rg) and abs(rg)<=1:
            values[i,j]=rg
            if r.get('eligible') in {'True','true','1'} and q<.05:stars.append((i,j))
    plt=_plt();fig,ax=plt.subplots(figsize=(9,8))
    image=ax.imshow(np.ma.masked_invalid(values),vmin=-1,vmax=1,aspect='auto')
    ax.set_xticks(range(6),['ADHD','MDD','Schizophrenia','Bipolar','Alzheimer’s','Parkinson’s'],rotation=35,ha='right')
    ax.set_yticks(range(12),[s.replace('_',' ') for s in SLEEP])
    for i,j in stars:ax.text(j,i,'*',ha='center',va='center')
    fig.colorbar(image,ax=ax,label='Genetic correlation')
    ax.set_title(('SYNTHETIC — ' if 'SYNTHETIC' in str(ranking).upper() else '')+'Measured sleep–brain atlas\n* original-family FDR < 0.05 and eligible QC')
    fig.tight_layout();fig.savefig(out,dpi=180);plt.close(fig)


def manhattan(results,out,*,p_column='P_PLACO',threshold=None,bin_bp=100000):
    require(bin_bp>=1,'Invalid plotting bin')
    maxima={};lengths={};counts={'input':0,'non_tested':0,'zero_p':0}
    for r in read_tsv(results,['CHR','BP',p_column]):
        counts['input']+=1
        if r.get('status','TESTED')!='TESTED' or r[p_column] in {'NA','','nan'}:
            counts['non_tested']+=1;continue
        ch=int(r['CHR']);bp=int(r['BP']);p=float(r[p_column])
        require(1<=ch<=22 and bp>=1 and math.isfinite(p) and 0<=p<=1,'Invalid Manhattan row')
        if p==0:counts['zero_p']+=1
        score=-math.log10(max(p,float(np.nextafter(0.,1.))))
        key=ch,bp//bin_bp
        if key not in maxima or score>maxima[key][1]:maxima[key]=bp,score
        lengths[ch]=max(lengths.get(ch,0),bp)
    require(maxima,'No tested variants to plot')
    offset={};total=0;ticks=[];labels=[]
    for ch,length in sorted(lengths.items()):
        offset[ch]=total;ticks.append(total+length/2);labels.append(str(ch));total+=length+1000000
    x=[];y=[]
    for (ch,_),(bp,score) in sorted(maxima.items()):x.append(offset[ch]+bp);y.append(score)
    plt=_plt();fig,ax=plt.subplots(figsize=(13,4.5))
    ax.scatter(x,y,s=5);ax.set_xticks(ticks,labels)
    if threshold is not None:
        require(0<threshold<1,'Invalid threshold');ax.axhline(-math.log10(threshold),linestyle='--')
    ax.set_xlabel('Autosome');ax.set_ylabel('-log10(P)')
    ax.set_title(f'Pleiotropy results — plotting maxima per {bin_bp:,} bp bin; inference uses all variants')
    fig.tight_layout();fig.savefig(out,dpi=180);plt.close(fig)
    return counts


def qq(sqlite_path,out,*,max_points=5000):
    require(max_points>=10,'Too few QQ display points')
    db=sqlite3.connect(f'{Path(sqlite_path).resolve().as_uri()}?mode=ro',uri=True)
    n=db.execute("SELECT count(*) FROM results WHERE status='TESTED'").fetchone()[0]
    require(n>0,'No tested P values')
    ranks=set(np.unique(np.geomspace(1,n,min(n,max_points)).astype(int)))
    expected=[];observed=[]
    for i,(p,) in enumerate(db.execute("SELECT p FROM results WHERE status='TESTED' ORDER BY p"),1):
        if i in ranks:
            expected.append(-math.log10((i-.5)/n));observed.append(-math.log10(max(p,float(np.nextafter(0.,1.)))))
    db.close()
    plt=_plt();fig,ax=plt.subplots(figsize=(5.5,5.5))
    ax.scatter(expected,observed,s=6)
    top=max(max(expected),max(observed));ax.plot([0,top],[0,top],linestyle='--')
    ax.set_xlabel('Expected -log10(P)');ax.set_ylabel('Observed -log10(P)')
    ax.set_title('Pleiotropy QQ — tested variants only\nReview exclusions and failure denominator separately')
    fig.tight_layout();fig.savefig(out,dpi=180);plt.close(fig)
