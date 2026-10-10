#!/usr/bin/env python3
"""Replace crowded extension pages with measured variable-height rows."""
import datetime
import importlib.util
import json
import textwrap
from pathlib import Path

P=Path(__file__).resolve().parents[1]
BASE=P/'scripts/89_finish_research_figures_v4_2.py'
spec=importlib.util.spec_from_file_location('_verified_figure_export',BASE)
f=importlib.util.module_from_spec(spec);spec.loader.exec_module(f)
OUT=P/'figures/compact_extension_reflow_v4_3'
RECORD=P/'manifests/compact_research_figure_exports_v4_3.json'

def main():
    if RECORD.exists():raise RuntimeError('EXISTING_FIGURE_MANIFEST_PRESERVED')
    previous_path=P/'manifests/compact_research_figure_exports_v4_2.json'
    if f.sha(previous_path)!='a726fc56a2f2e92659649107d6da5beaa38be3cc73d53ee63dacd24159b36c30':
        raise RuntimeError('PREVIOUS_FIGURE_MANIFEST_CHANGED')
    previous=json.loads(previous_path.read_text())
    frozen=dict(previous['source_sha256'])
    frozen[str(previous_path)]=f.sha(previous_path)
    frozen[str(Path(__file__))]=f.sha(Path(__file__))
    for path,digest in frozen.items():
        if f.sha(path)!=digest:raise RuntimeError('VERIFIED_FIGURE_SOURCE_CHANGED')
    for item in previous['figures']:
        for path,digest in item['output_sha256'].items():
            if f.sha(path)!=digest:raise RuntimeError('PREVIOUS_FIGURE_CHANGED')
    source=P/'tables/extension_native_full_precision_rg_v4.tsv'
    data=f.rows(source);outcomes=list(dict.fromkeys(r['outcome_trait'] for r in data))
    lookup={(r['sleep_trait'],r['outcome_trait']):r for r in data}
    if len(data)!=1200 or len(lookup)!=1200 or len(outcomes)!=100:
        raise RuntimeError('COMPLETE_1200_MATRIX_REQUIRED')
    if sum(r['native_fdr_pass_0_05']=='True' for r in data)!=603:
        raise RuntimeError('FROZEN_603_BH_BOUNDARY_CHANGED')
    lim=f.np.ceil(max(1.,max(abs(float(r['rg'])) for r in data))*10)/10
    semax=max(float(r['se']) for r in data)
    f.plt.rcParams.update({'font.family':'DejaVu Sans','font.size':9,'pdf.fonttype':42,
                          'svg.fonttype':'none','axes.spines.top':False,'axes.spines.right':False})
    OUT.mkdir(exist_ok=False);f.OUT=OUT;records=[]
    for page,start in enumerate(range(0,100,13),1):
        selected=outcomes[start:start+13]
        labels=['\n'.join(textwrap.wrap(lookup[(f.SLEEPS[0],t)]['original_phenotype_name'],27)) for t in selected]
        heights=f.np.array([len(label.splitlines())*11.3+9 for label in labels])
        heights*=max(1.,430/heights.sum())
        edges=f.np.r_[0.,f.np.cumsum(heights)];centers=(edges[:-1]+edges[1:])/2
        height_inches=(float(edges[-1])+310)/72
        fig,axes=f.plt.subplots(1,2,figsize=(7.2,height_inches))
        total_points=height_inches*72
        fig.subplots_adjust(left=.385,right=.97,bottom=180/total_points,
                            top=(180+edges[-1])/total_points,wspace=.12)
        points=f.np.array([[float(lookup[(sleep,t)]['rg']) for sleep in f.SLEEPS] for t in selected])
        se=f.np.array([[float(lookup[(sleep,t)]['se']) for sleep in f.SLEEPS] for t in selected])
        if not f.np.isfinite(points).all() or not f.np.isfinite(se).all() or not (se>0).all():
            raise RuntimeError('NONFINITE_POINT_OR_SE')
        plots=[axes[0].pcolormesh(f.np.arange(13)-.5,edges,points,cmap='RdBu_r',vmin=-lim,vmax=lim,shading='flat'),
               axes[1].pcolormesh(f.np.arange(13)-.5,edges,se,cmap='Blues',vmin=0,vmax=semax,shading='flat')]
        for ax in axes:
            ax.set_ylim(edges[-1],0);ax.set_xlim(-.5,11.5)
            ax.set_xticks(range(12),f.SLEEP_LABELS,rotation=78,ha='right',fontsize=9)
            ax.tick_params(length=0)
        axes[0].set_yticks(centers,labels,fontsize=9)
        axes[1].set_yticks(centers,['']*len(selected),fontsize=9)
        for j,ax in enumerate(axes):
            ax.set_title('A  Genetic correlation' if j==0 else 'B  Standard error',fontsize=10,loc='left',pad=10)
            position=ax.get_position()
            cax=fig.add_axes([position.x0+.015,67/total_points,position.width-.03,9/total_points])
            bar=fig.colorbar(plots[j],cax=cax,orientation='horizontal')
            bar.set_ticks([-lim,0,lim] if j==0 else [0,semax/2,semax])
            bar.set_label('rg' if j==0 else 'SE of rg',fontsize=9)
            bar.ax.tick_params(labelsize=9)
        for i,t in enumerate(selected):
            for j,sleep in enumerate(f.SLEEPS):
                if lookup[(sleep,t)]['native_fdr_pass_0_05']=='True':
                    axes[0].plot(j,centers[i],'o',color='#111111',ms=2.5)
        fig.text(.045,1-22/total_points,'Phenome extension: points and uncertainty',fontsize=12,weight='bold',va='top')
        fig.text(.045,1-55/total_points,'Complete frozen family; page '+str(page)+' of 8\n'+str(12*len(selected))+' marginal estimates; all rows retained.',fontsize=9,va='top')
        note='Dots: original all-1,200 BH q < 0.05. Fixed scales across pages.\nExact values and marginal intervals are in the table. No calibrated\ncontrast is shown. Act. = actigraphy; source qualifications remain.'
        fig.text(.045,10/total_points,note,fontsize=9,va='bottom')
        fig.canvas.draw();renderer=fig.canvas.get_renderer()
        boxes=sorted((t.get_window_extent(renderer) for t in axes[0].get_yticklabels()),key=lambda b:b.y0)
        if any(a.y1+1>b.y0 for a,b in zip(boxes,boxes[1:])):
            raise RuntimeError('EXTENSION_LABELS_OVERLAP')
        alt='Extension page '+str(page)+': '+str(len(selected))+' outcomes and all 12 sleep traits, with marginal rg and SE; row height follows full label length. '+note.replace('\n',' ')
        f.export(fig,'03_extension_native_reflow_v4_page%02d'%page,alt,source,12*len(selected),records)
        records[-1]['label_bounding_boxes_nonoverlapping']=True
        records[-1]['row_height_points']=list(heights)
    selected_records=[];cursor=iter(records)
    for item in previous['figures']:
        if item['figure'].startswith('03_extension_'):item=next(cursor)
        else:
            item=dict(item);item['visual_QA_pending']=False
            item['visual_QA']='Root inspected Poppler-rendered contact sheets; labels, intervals and qualifications readable.'
        selected_records.append(item)
    if len(selected_records)!=21 or sum(r['supported_value_count'] for r in records)!=1200:
        raise RuntimeError('FINAL_FIGURE_SELECTION_CARDINALITY_CHANGED')
    for path,digest in frozen.items():
        if f.sha(path)!=digest:raise RuntimeError('FIGURE_SOURCE_CHANGED_AFTER_EXPORT')
    result=dict(schema='selected_scientific_figure_exports_v4_3',created_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),
                figure_page_count=21,export_count=63,source_sha256=frozen,figures=selected_records,
                prior_extension_pages_preserved=True,prior_extension_QA='Eight uniform-row pages replaced because long labels overlapped.',
                replacement_page_count=8,visual_QA_pending=True,new_estimator_fits=0,new_tests=0,
                calibrated_contrast_P_values=0,manuscript_legends_authored=False)
    with RECORD.open('x') as out:json.dump(result,out,indent=2,allow_nan=False);out.write('\n')
    print(json.dumps({'selected_pages':21,'replacement_pages':8,'manifest_sha256':f.sha(RECORD)}),flush=True)

if __name__=='__main__':main()
