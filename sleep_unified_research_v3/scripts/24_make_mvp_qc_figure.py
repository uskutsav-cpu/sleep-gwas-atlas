#!/usr/bin/env python3
"""Source-only QC figure, directly traceable to diagnostic numerical tables."""
import csv
import hashlib
import json
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

PACKAGE=Path(__file__).resolve().parents[1]

def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()

def write_table(path,rows,fields):
    with path.open('x',newline='') as f:
        w=csv.DictWriter(f,delimiter='\t',fieldnames=fields,lineterminator='\n');w.writeheader();w.writerows(rows)

def main():
    receipt_path=PACKAGE/'logs/mvp_preprocessing_worker_receipt_v3.json'
    d=json.loads(receipt_path.read_text())
    if d['status']!='PREPROCESSING_ONLY_PASS':raise RuntimeError('NO_SUCCESSFUL_PREPROCESSING_RECEIPT')
    tables=PACKAGE/'tables';figures=PACKAGE/'figures';figures.mkdir(exist_ok=True)
    counts=[{'metric':k,'count':v,'denominator':d['counts']['source_rows'],'scope':'mutually_exclusive_source_filter'} for k,v in d['exclusive_exclusions'].items()]
    counts.append({'metric':'retained_HM3','count':d['counts']['output_rows'],'denominator':d['counts']['source_rows'],'scope':'retained_output'})
    count_path=tables/'MVP_PREPROCESSING_COUNTS.tsv';write_table(count_path,counts,['metric','count','denominator','scope'])
    order=['0.001','0.01','0.05','0.1','1','10','gt10']
    buckets=[{'upper_bound_relative_difference':k,'count':d['retained_P_CI_relative_difference_buckets'].get(k,0),
              'denominator_retained_rows':d['counts']['output_rows'],'formula':'abs(P_CI_Z-P_source)/max(P_source,1e-300)',
              'is_filter':'False'} for k in order]
    p_path=tables/'MVP_RETAINED_P_CI_DIAGNOSTICS.tsv';write_table(p_path,buckets,list(buckets[0]))
    center=[{'upper_bound_SE_units':k,'count':d['whole_file_log_CI_centering_SE_unit_buckets'].get(k,0),
             'denominator_source_rows':d['counts']['source_rows'],'is_filter':'False'} for k in order]
    c_path=tables/'MVP_WHOLE_SOURCE_CI_CENTERING.tsv';write_table(c_path,center,list(center[0]))
    with count_path.open(newline='') as f:counts=list(csv.DictReader(f,delimiter='\t'))
    with p_path.open(newline='') as f:buckets=list(csv.DictReader(f,delimiter='\t'))
    plt.rcParams.update({'font.family':'DejaVu Sans','font.size':9,'pdf.fonttype':42,'svg.fonttype':'none',
                         'axes.spines.top':False,'axes.spines.right':False})
    fig,axes=plt.subplots(1,2,figsize=(12.4,4.8),gridspec_kw={'width_ratios':[1.3,1]})
    labels={'not_nonambiguous_HM3':'Outside nonambiguous HM3','invalid_OR_P_CI_or_R2_le_0_9':'OR/P/CI/R2 exclusion',
            'invalid_AF_or_MAF_le_0_01':'AF/MAF exclusion','duplicate_eligible_HM3_rsID':'Duplicate eligible rsID',
            'allele_mismatch_or_ambiguous':'Allele mismatch/ambiguous','chain_unmapped':'Unmapped coordinate','retained_HM3':'Retained HM3'}
    values=[int(r['count']) for r in counts];names=[labels[r['metric']] for r in counts]
    axes[0].barh(names,values,color=['#215a88' if r['metric']=='retained_HM3' else '#8093a1' for r in counts]);axes[0].invert_yaxis()
    axes[0].set_xscale('log');axes[0].set_xlim(1,1e8);axes[0].set_xlabel('Rows (log scale; mutually exclusive categories)')
    for i,n in enumerate(values):axes[0].text(n*1.15,i,format(n,','),va='center',fontsize=8)
    axes[0].set_title('A  Source retention and exclusions',loc='left',fontweight='bold',pad=13)
    pct=[100*int(r['count'])/int(r['denominator_retained_rows']) for r in buckets]
    xlabels=['≤0.1%','0.1–1%','1–5%','5–10%','10–100%','100–1,000%','>1,000%']
    axes[1].bar(range(7),pct,color=['#215a88']*4+['#98623b']*3)
    axes[1].set_xticks(range(7),xlabels,rotation=40,ha='right');axes[1].set_ylim(0,44)
    axes[1].set_ylabel('Retained rows (%)');axes[1].set_xlabel('Relative difference from reported P')
    for i,p in enumerate(pct):axes[1].text(i,p+.8,('%.2f%%'%p if p>=.01 else '<0.01%'),ha='center',fontsize=8)
    axes[1].set_title('B  Reported P versus CI-derived Z',loc='left',fontweight='bold',pad=13)
    fig.suptitle('MVP GIA clinical-insomnia source preprocessing · no correlation tested',x=.01,ha='left',fontsize=13,fontweight='bold')
    fig.text(.01,.01,'Historical-rule CI-derived Z assumes 95% log-scale Wald intervals; discrepancy is descriptive and does not filter rows.\n19,703,815 source rows; 826,027 retained. All source/reference hashes unchanged. No heritability or power gate passed.',fontsize=8,color='#4f5961')
    fig.tight_layout(rect=(0,.12,1,.93));output=[]
    for suffix in ['png','pdf','svg']:
        p=figures/('MVP_source_preprocessing_QC_v3.'+suffix);fig.savefig(p,dpi=240,bbox_inches='tight',facecolor='white');output.append(str(p.relative_to(PACKAGE)))
    plt.close(fig)
    manifest={'source_receipt_sha256':sha(receipt_path),'script_sha256':sha(__file__),
              'tables':{str(p.relative_to(PACKAGE)):sha(p) for p in [count_path,p_path,c_path]},
              'outputs':{p:sha(PACKAGE/p) for p in output},'genetic_correlation_or_mechanism_claims':False,
              'matplotlib_version':matplotlib.__version__}
    with (PACKAGE/'manifests/mvp_QC_figure_receipt_v3.json').open('x') as f:json.dump(manifest,f,indent=2);f.write('\n')
    print(json.dumps({'status':'FIGURE_GENERATED_FOR_VISUAL_QA','outputs':output}),flush=True)

if __name__=='__main__':main()
