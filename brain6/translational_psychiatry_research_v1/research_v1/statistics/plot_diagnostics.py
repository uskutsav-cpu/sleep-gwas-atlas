"""Research diagnostic plots. All plotted estimates are truth-known simulations."""
import hashlib, json
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from calibration import HERE, sha

def main():
    table=HERE/'calibration_simulations.tsv'; x=pd.read_csv(table,sep='\t')
    fig,ax=plt.subplots(1,2,figsize=(10,4),layout='constrained')
    for i,h in enumerate(sorted(x.local_h2_both.unique())):
        for c in sorted(x.true_intercept.unique()):
            sub=x[(x.local_h2_both==h)&(x.true_intercept==c)].sort_values('ld_ar1')
            color=f'C{i}'; marker='o' if c==0 else 's'; style='-' if c==0 else ':'
            label=f'h²={h:g}, intercept={c:g}'
            ax[0].plot(sub.ld_ar1,100*sub['p_alpha0.05_rate'],style+marker,color=color,label=label)
            lo=100*sub['p_alpha0.05_rate']-100*sub['p_alpha0.05_wilson_lower']
            hi=100*sub['p_alpha0.05_wilson_upper']-100*sub['p_alpha0.05_rate']
            ax[0].errorbar(sub.ld_ar1,100*sub['p_alpha0.05_rate'],yerr=[lo,hi],fmt='none',color=color,capsize=2)
            ax[1].plot(sub.ld_ar1,sub.exact_over_reported_variance,style+marker,color=color,label=label)
    ax[0].axhline(5,color='black',linewidth=1,label='Nominal 5%')
    ax[0].set(ylabel='Null rejection at α=0.05 (%)',xlabel='AR(1) LD parameter',title='Archived implementation')
    ax[1].axhline(1,color='black',linewidth=1)
    ax[1].set(ylabel='Exact moment variance / mean reported variance',xlabel='AR(1) LD parameter',title='Point estimate / uncertainty mismatch')
    ax[0].legend(fontsize=7,loc='upper left')
    fig.suptitle('Truth-known Gaussian summary calibration — synthetic method diagnostics',fontsize=11)
    files=[]
    for ext in ['png','svg']:
        out=HERE/f'calibration_diagnostic.{ext}'; fig.savefig(out,dpi=160); files.append({'file':out.name,'sha256':sha(out),'bytes':out.stat().st_size})
    plt.close(fig)
    (HERE/'plot_provenance.json').write_text(json.dumps({'classification':'SYNTHETIC_RESEARCH_DIAGNOSTIC','input':table.name,'input_sha256':sha(table),'script_sha256':sha(__file__),'outputs':files},indent=2)+'\n')

if __name__=='__main__': main()
