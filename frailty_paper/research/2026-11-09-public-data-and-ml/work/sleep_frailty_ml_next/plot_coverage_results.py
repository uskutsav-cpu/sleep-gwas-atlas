"""Compact scientific display of held-study errors, coverage and answer availability."""
import sys
from pathlib import Path
ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT/'work/sleep_frailty_ml_v9'))
from v9_common import *
DEST=ROOT/'outputs/sleep_frailty_ml_next';FINAL=DEST/'validation'

def main():
    audit=json.loads((FINAL/'reserved_results_independent_audit.json').read_text());assert audit['passed'];r=json.loads((FINAL/'reserved_outcome_evaluation.json').read_text());assert sha(FINAL/'reserved_outcome_evaluation.json')==audit['original_result_sha256']
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'axes.spines.top':False,'axes.spines.right':False})
    tests=['blood_K562','OL33_three_known_lines'];arms=['sequence','AVI'];colors=['#176B87','#A85525'];labels=['Independent publication\nBlood K562','Held library in same paper\nOL33 three cell lines'];fig,axes=plt.subplots(1,3,figsize=(12,4.3));xx=np.arange(2)
    for j,arm in enumerate(arms):
        values=[];lo=[];hi=[];coverage=[];accepted=[]
        for test in tests:
            a=r['primary_accuracy'][arm][test];c=r['confidence'][arm][test];v=100*a['relative_MAE_gain_vs_line_median'];ci=100*np.array(a['gain_vs_line_95CI']);values.append(v);lo.append(max(0,v-ci[0]));hi.append(max(0,ci[1]-v));coverage.append(100*c['observed_coverage']);accepted.append(100*c['accepted_fraction'])
        offset=(j-.5)*.25;axes[0].errorbar(xx+offset,values,yerr=[lo,hi],fmt='o',capsize=4,color=colors[j],label=arm);axes[1].bar(xx+offset,coverage,width=.24,color=colors[j],label=arm);axes[2].bar(xx+offset,accepted,width=.24,color=colors[j])
        for x,value in zip(xx+offset,accepted):axes[2].text(x,max(value,1)+1,f'{value:.1f}%',ha='center',fontsize=8)
    axes[0].axhline(0,color='#333333',linewidth=1);axes[0].axhline(5,color='#666666',linestyle=':',linewidth=1);axes[0].set_ylabel('Error improvement over typical effect (%)');axes[0].set_title('Accuracy with 95% window bootstrap CI');axes[0].legend(frameon=False,fontsize=9)
    axes[1].axhline(90,color='#333333',linestyle='--',linewidth=1);axes[1].set_ylim(0,105);axes[1].set_ylabel('Measurements inside interval (%)');axes[1].set_title('Coverage of nominal 90% intervals')
    axes[2].set_ylim(0,105);axes[2].set_ylabel('Examples accepted before target join (%)');axes[2].set_title('Answer availability before quality gates')
    for ax in axes:ax.set_xticks(xx,labels);ax.grid(axis='y',alpha=.15);ax.set_axisbelow(True)
    fig.suptitle('Reserved laboratory tests of DNA effect magnitude',fontsize=13);fig.text(.5,.015,'Reporter experiments only. Broad intervals or withholding every case do not establish useful prediction.',ha='center',fontsize=9);fig.tight_layout(rect=[0,.06,1,.94]);folder=FINAL/'figures';folder.mkdir(exist_ok=True);path=folder/'reserved_accuracy_coverage_availability.png';fig.savefig(path,dpi=180);plt.close(fig)
    write_json(FINAL/'figure_receipt.json',{'path':str(path),'sha256':sha(path),'result_sha256':sha(FINAL/'reserved_outcome_evaluation.json'),'all_numbers_derived_from_independently_audited_results':True,'visual_QA_pending':True});print(str(path),flush=True)

if __name__=='__main__':main()
