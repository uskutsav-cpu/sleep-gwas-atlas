#!/usr/bin/env python3
"""Static historical-evidence figures, traceable to sealed numerical tables."""
import csv
from collections import Counter
from datetime import datetime,timezone
import hashlib
import json
from pathlib import Path
import textwrap
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
import numpy as np

PACKAGE=Path(__file__).resolve().parents[1]
FIG=PACKAGE/'figures'
INPUTS={};ARTIFACTS=[]
def read(relative):
    p=PACKAGE/relative;INPUTS[relative]=hashlib.sha256(p.read_bytes()).hexdigest()
    with p.open(newline='') as f:return list(csv.DictReader(f,delimiter='\t'))
def save(fig,name,tables,claims):
    paths=[]
    for suffix in ['pdf','svg','png']:
        p=FIG/(name+'.'+suffix);fig.savefig(p,dpi=220,bbox_inches='tight',facecolor='white');paths.append(str(p.relative_to(PACKAGE)))
    plt.close(fig)
    ARTIFACTS.append({'figure':name,'outputs':paths,'source_tables':tables,'claims':claims,'evidence_level':'HISTORICAL_TABLES_AND_CURRENT_EXECUTION_STATUS; NOT_FULL_NATIVE_REPRODUCTION'})
def style():
    plt.rcParams.update({'font.family':'DejaVu Sans','font.size':9,'axes.spines.top':False,'axes.spines.right':False,'pdf.fonttype':42,'svg.fonttype':'none','axes.titlesize':12})
def main():
    style();FIG.mkdir(exist_ok=True)
    core=read('tables/original_core_396.tsv');ext=read('tables/original_extension_1200.tsv');rep=read('REPLICATION_RESULTS.tsv')
    counts=read('tables/study_design_counts.tsv');meta=read('tables/phenotype_and_source_metadata.tsv')
    labels={r['trait_id']:r['phenotype'] for r in meta}
    sleep=list(dict.fromkeys(r['sleep_trait'] for r in core));disease=list(dict.fromkeys(r['disease_trait'] for r in core))
    assert len(sleep)==12 and len(disease)==33 and len(core)==396 and len(ext)==1200
    blue='#215a88';gray='#68737d';red='#a52c35'
    fig,ax=plt.subplots(figsize=(12,5));ax.set(xlim=(0,12),ylim=(0,5));ax.axis('off')
    ax.text(6,4.7,'Unified sleep genetics evidence · current status',ha='center',fontsize=16,weight='bold')
    ax.text(6,4.28,'12 original sleep traits; three preserved testing families',ha='center',color=gray)
    descriptions=[('Original atlas','12 × 33 = 396 tests','153 primary FDR positives\n161 including 8 sensitivity rows','Native replay: 1 / 396 pairs'),('Phenome extension','12 × 100 = 1,200 tests','603 historical FDR positives','Native replay: 0 / 1,200 pairs'),('External validation','217 locked candidates','23 qualified outcome-side positives\n18 concordant; 17 QC excluded; 159 absent','Native replay: 0 / 41 estimates')]
    for x,(title,n,positive,native) in zip([.15,4.15,8.15],descriptions):
        ax.add_patch(Rectangle((x,1.45),3.7,2.35,facecolor='#f0f5f8',edgecolor=blue,lw=1.2))
        for y,text,size in [(3.42,title,13),(2.95,n,11),(2.36,positive,9),(1.76,native,10)]:
            ax.text(x+1.85,y,text,ha='center',va='center',fontsize=size,color=blue if y==1.76 else '#182735')
    ax.text(6,.95,'18 / 18 checkpoint files recovered exactly · 100 extension + 13 validation inputs receipt-matched',ha='center')
    ax.text(6,.55,'0 completed independent both-trait validations · no new primary biological question admitted',ha='center',color=red,weight='bold')
    ax.text(6,.15,'Exact archival recovery and printed-log concordance do not complete native source-to-result reproduction.',ha='center',fontsize=8,color=gray)
    save(fig,'01_unified_study_design',['tables/study_design_counts.tsv','REPLICATION_RESULTS.tsv'],['Separate families396/1200/217','Native counts1/0/0','Independent completed count0'])

    def atlas(rows,traits,key,qcol,name,size,primary=False):
        lookup={(r['sleep_trait'],r[key]):r for r in rows}
        matrix=np.array([[float(lookup[(s,t)]['rg']) for s in sleep] for t in traits])
        fig,ax=plt.subplots(figsize=size)
        im=ax.imshow(matrix,cmap='RdBu_r',vmin=-1,vmax=1,aspect='auto',interpolation='none')
        for i,t in enumerate(traits):
            is_sensitivity=primary and lookup[(sleep[0],t)]['analysis_tier']!='PRIMARY_PHASE1'
            if is_sensitivity:ax.add_patch(Rectangle((-.5,i-.5),12,1,fill=False,hatch='////',edgecolor='#777',lw=0))
            for j,s in enumerate(sleep):
                r=lookup[(s,t)]
                if float(r[qcol])<.05:ax.plot(j,i,'o',ms=2.8 if not is_sensitivity else 2,color='black' if not is_sensitivity else '#777')
        ax.set_xticks(range(12),[labels[s] for s in sleep],rotation=55,ha='right',fontsize=8)
        ax.set_yticks(range(len(traits)),[textwrap.fill(labels[t],53) for t in traits],fontsize=7 if len(traits)>40 else 9)
        ax.set_title('Original core · frozen 396-test family' if primary else 'Original extension · frozen 1,200-test family',loc='left',pad=15)
        bar=fig.colorbar(im,ax=ax,pad=.025,fraction=.018);bar.set_label('Global genetic correlation (rg)')
        ax.set_xlabel('Historical estimates; dots: within-family BH q < 0.05. '+('Hatching: sensitivity only.' if primary else 'Native replay incomplete.'),labelpad=12,fontsize=8)
        ax.tick_params(length=0)
        save(fig,name,['tables/original_core_396.tsv' if primary else 'tables/original_extension_1200.tsv','tables/phenotype_and_source_metadata.tsv'],['Historical global rg only','Separate frozen BH correction','No objective/self-report significance-count comparison'])
    atlas(core,disease,'disease_trait','fdr','02_global_core_atlas',(10,11),True)
    etraits=list(dict.fromkeys(r['extension_trait_id'] for r in ext))
    atlas(ext,etraits,'extension_trait_id','extension_fdr','03_global_extension_atlas',(11,27))

    classes=Counter(r['historical_replication_class'] for r in rep)
    fig,axs=plt.subplots(1,2,figsize=(11,4.5),gridspec_kw={'width_ratios':[1.5,1]})
    names=['Qualified outcome-side positive','Concordant below threshold','QC ineligible','No eligible external source']
    vals=[classes[k] for k in ['REPLICATED','DIRECTIONALLY_CONCORDANT','UNDERPOWERED','NO_INDEPENDENT_DATASET']]
    bars=axs[0].barh(range(4),vals,color=[blue,'#7ba6c3','#d9a441','#a8afb5'])
    axs[0].invert_yaxis();axs[0].set_yticks(range(4),names);axs[0].set_xlabel('Locked candidates (total 217)');axs[0].set_xlim(0,183)
    for bar,value in zip(bars,vals):axs[0].text(value+2,bar.get_y()+bar.get_height()/2,str(value),va='center')
    axs[0].set_title('Historical classification',loc='left')
    axs[1].axis('off');axs[1].text(.5,.83,'0',ha='center',fontsize=54,color=red)
    axs[1].text(.5,.58,'Completed independent\nboth-trait validations',ha='center',fontsize=13)
    axs[1].text(.5,.32,'All 41 estimated historical pairs\nreuse the discovery sleep source.',ha='center',fontsize=10)
    axs[1].text(.5,.1,'New MVP clinical-insomnia source:\n10 conditional transport candidates;\nno new rg estimated.',ha='center',fontsize=9,color=gray)
    fig.suptitle('Validation evidence and independence boundary',fontsize=14);fig.tight_layout()
    save(fig,'04_validation_evidence',['REPLICATION_RESULTS.tsv','tables/study_design_counts.tsv'],['23positive18concordant17QC159absent','0 completed independent validations','MVP proposal metadata only'])

    hyp=read('HYPOTHESIS_FEASIBILITY_MATRIX.tsv')
    keys=['novelty_0_4','rationale_0_4','data_0_4','statistical_validity_0_4','power_0_4','clinical_relevance_0_4','independent_test_0_4','feasibility_0_4','stopping_rule_readiness_0_4']
    fig,ax=plt.subplots(figsize=(10,4.5));h=np.array([[int(r[k]) for k in keys] for r in hyp])
    ax.imshow(h,cmap='Blues',vmin=0,vmax=4,aspect='auto')
    ax.set_yticks(range(5),[r['question_id']+' · '+textwrap.shorten(r['candidate_estimand'],width=55,placeholder='…') for r in hyp],fontsize=8)
    ax.set_xticks(range(9),['Novelty','Rationale','Data','Statistical\nvalidity','Power','Clinical\nrelevance','Independent\ntest','Feasibility','Stopping\nrule'],rotation=30,ha='right',fontsize=8)
    for i in range(5):
        for j in range(9):ax.text(j,i,str(h[i,j]),ha='center',va='center',color='white' if h[i,j]>=3 else '#192632')
    ax.set_title('Prospective question audit · no primary question admitted',loc='left',pad=15)
    ax.set_xlabel('Reviewer ordinal scores (0–4), not biological results. All five fail at least one hard admission gate.',labelpad=14,fontsize=8);ax.tick_params(length=0)
    save(fig,'05_question_admission_audit',['HYPOTHESIS_FEASIBILITY_MATRIX.tsv'],['All five prospective hypotheses currentlyNO_GO','Scores ordinal not scientific estimates'])

    estimated=[r for r in rep if r['replication_rg'] not in ['NA','']]
    fig,axs=plt.subplots(1,2,figsize=(11,5))
    for positive,col in [(False,'#9aa4ad'),(True,blue)]:
        rows=[r for r in estimated if (r['historical_replication_class']=='REPLICATED')==positive]
        axs[0].errorbar([float(r['discovery_rg']) for r in rows],[float(r['replication_rg']) for r in rows],
                       xerr=[1.96*float(r['discovery_se']) for r in rows],yerr=[1.96*float(r['replication_se']) for r in rows],
                       fmt='o',ms=4,elinewidth=.6,capsize=0,color=col,alpha=.8,label='Outcome-side positive' if positive else 'Concordant below threshold')
    axs[0].plot([-.7,.9],[-.7,.9],ls='--',lw=.8,color=gray);axs[0].set(xlabel='Discovery rg (95% marginal interval)',ylabel='Outcome-side validation rg (95% marginal interval)',title='41 historical estimated pairs');axs[0].legend(fontsize=7)
    z=[float(r['heterogeneity_z']) for r in estimated];axs[1].hist(z,bins=np.arange(-5,5.5,.5),color='#8caec4',edgecolor='white')
    for v in [-1.959963984540054,1.959963984540054]:axs[1].axvline(v,color=red,ls='--',lw=1)
    axs[1].set(xlabel='Effect-difference Z under Cov = 0',ylabel='Historical pair count',title='Uncalibrated heterogeneity diagnostic')
    allnom=sum(float(r['heterogeneity_p'])<.05 for r in estimated);posnom=sum(float(r['heterogeneity_p'])<.05 and r['historical_replication_class']=='REPLICATED' for r in estimated)
    assert (allnom,posnom)==(16,7)
    axs[1].text(.97,.95,'16 / 41 nominal flags\n7 / 23 positive-subset flags',transform=axs[1].transAxes,ha='right',va='top',fontsize=9)
    fig.text(.5,.01,'Shared sleep estimates, unestimated covariance and discovery selection prevent confirmed effect-difference inference.',ha='center',fontsize=8,color=red);fig.tight_layout(rect=(0,.05,1,1))
    save(fig,'06_covariance_heterogeneity_diagnostic',['REPLICATION_RESULTS.tsv'],['Covariancezero only','16 nominal all41;7positive23','Marginal intervals not joint covariance-aware contrasts'])

    standalone=read('sources/recovered/results/tables/h2_summary.tsv');hm={r['trait']:r for r in standalone}
    qc=[]
    for t in disease:
        rows=[r for r in core if r['disease_trait']==t];raw=hm[t]
        intkey='intercept' if 'intercept' in raw else 'h2_int'
        qc.append({'trait_id':t,'label':labels[t],'standalone_intercept':float(raw[intkey]),'pairwise_intercept_min':min(float(r['h2_int']) for r in rows),'pairwise_intercept_max':max(float(r['h2_int']) for r in rows),'historical_tier':rows[0]['analysis_tier']})
    p=PACKAGE/'tables/core_intercept_stage_diagnostic.tsv'
    with p.open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(qc[0]),delimiter='\t');w.writeheader();w.writerows(qc)
    INPUTS[str(p.relative_to(PACKAGE))]=hashlib.sha256(p.read_bytes()).hexdigest()
    fig,ax=plt.subplots(figsize=(9,11))
    for y,r in enumerate(qc):
        col=red if r['trait_id'] in ['hdl','ldl','triglycerides'] else blue
        ax.plot([r['pairwise_intercept_min'],r['pairwise_intercept_max']],[y,y],color=col,lw=3)
        ax.plot(r['standalone_intercept'],y,'o',mfc='white',mec=col,ms=5)
    ax.axvline(1.2,color=gray,ls='--',lw=1);ax.set_yticks(range(33),[r['label'] for r in qc]);ax.invert_yaxis()
    ax.set(xlabel='Outcome LDSC intercept (circle: standalone; line: pairwise range)',title='Core QC diagnostic · estimator-stage differences')
    ax.text(.01,-.08,'All 36 sleep–lipid pairwise intercepts exceed 1.2. Historical standalone-based classifications are preserved.\nResolve estimator/SNP-set differences before new lipid-conditioned inference.',transform=ax.transAxes,fontsize=8,color=red)
    save(fig,'07_core_QC_intercept_diagnostic',['tables/core_intercept_stage_diagnostic.tsv','tables/original_core_396.tsv','sources/recovered/results/tables/h2_summary.tsv'],['All36lipidpairwiseinterceptsabove1.2','Historicalclassificationpreserved'])

    manifest={'created_utc':datetime.now(timezone.utc).isoformat(),'matplotlib_version':matplotlib.__version__,'numpy_version':np.__version__,
              'input_sha256':INPUTS,'figures':ARTIFACTS,'output_sha256':{p:hashlib.sha256((PACKAGE/p).read_bytes()).hexdigest() for a in ARTIFACTS for p in a['outputs']},
              'new_biological_claims':False,'mechanism_figure_admitted':False}
    (PACKAGE/'manifests/figure_manifest_v1.json').write_text(json.dumps(manifest,indent=2)+'\n')
    print(json.dumps({'figures':len(ARTIFACTS),'outputs':sum(len(a['outputs']) for a in ARTIFACTS)},indent=2))
if __name__=='__main__':main()
