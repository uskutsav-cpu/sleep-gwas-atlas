#!/usr/bin/env python3
"""Scientific figures from verified full-precision tables; no new inference."""
import csv
import datetime
import hashlib
import json
import os
from pathlib import Path
import platform
import resource
import shutil
import textwrap
import xml.etree.ElementTree as ET

P=Path(__file__).resolve().parents[1]
SSD=Path('/Volumes/Extreme SSD/sleep-unified-research-v1/research-completion-2026-10-09/sleep_unified_research_v4')
os.environ.setdefault('MPLCONFIGDIR',str(SSD/'cache/matplotlib'))
os.environ.setdefault('OPENBLAS_NUM_THREADS','1')
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle,FancyBboxPatch
import numpy as np
from PIL import Image

OUT=P/'figures/compact_research_v4_1'
RECORD=P/'manifests/compact_research_figure_exports_v4_1.json'
SLEEPS=['insomnia','sleepdur','shortsleep','longsleep','chronotype','sleepiness','napping','snoring','sleep_apnea','sleep_efficiency','accel_sleep_duration','sleep_timing']
SLEEP_LABELS=['Insomnia','Duration','Short <7 h','Long >=9 h','Chronotype','Sleepiness','Napping','Snoring','Apnea','Act. efficiency','Act. duration','Act. timing']


def sha(path):
    h=hashlib.sha256()
    with Path(path).open('rb') as f:
        for b in iter(lambda:f.read(65536),b''):h.update(b)
    return h.hexdigest()


def rows(path):
    with Path(path).open() as f:return list(csv.DictReader(f,delimiter='\t'))


def guard():
    if shutil.disk_usage('/System/Volumes/Data').free<3<<30 or shutil.disk_usage(SSD).free<5<<30:
        raise RuntimeError('UNCHANGED_FIGURE_STORAGE_FLOOR_FAILS')
    if resource.getrusage(resource.RUSAGE_SELF).ru_maxrss>2<<30:
        raise RuntimeError('FIGURE_MEMORY_LIMIT_FAILS')


def export(fig,name,alt,table,supported_values,records):
    guard();fig.canvas.draw();renderer=fig.canvas.get_renderer();outside=[]
    for text in fig.findobj(matplotlib.text.Text):
        if not text.get_visible() or not text.get_text():continue
        box=text.get_window_extent(renderer)
        if box.x0 < -1 or box.y0 < -1 or box.x1>fig.bbox.width+1 or box.y1>fig.bbox.height+1:outside.append(text.get_text())
    if outside:raise RuntimeError('TEXT_OUTSIDE_CANVAS: '+repr(outside))
    paths=[OUT/(name+'.'+ext) for ext in ['pdf','svg','png']]
    if any(path.exists() for path in paths):raise RuntimeError('PRIOR_FIGURE_PRESERVED_NO_OVERWRITE')
    for path in paths:fig.savefig(path,dpi=300,facecolor='white')
    ns='http://www.w3.org/2000/svg';ET.register_namespace('',ns)
    tree=ET.parse(paths[1]);root=tree.getroot()
    title=ET.Element('{'+ns+'}title',{'id':name+'-title'});title.text=name.replace('_',' ')
    desc=ET.Element('{'+ns+'}desc',{'id':name+'-desc'});desc.text=alt
    root.insert(0,title);root.insert(1,desc);root.set('role','img');root.set('aria-labelledby',name+'-title '+name+'-desc')
    tree.write(paths[1],encoding='utf-8',xml_declaration=True)
    with Image.open(paths[2]) as png:
        pixels=png.size;dpi=png.info['dpi']
    if any(abs(x-300)>.02 for x in dpi):raise RuntimeError('300DPI_EXPORT_FAILED')
    records.append(dict(figure=name,output_sha256={str(path):sha(path) for path in paths},
                        export_size_inches=list(fig.get_size_inches()),PNG_dimensions_pixels=pixels,PNG_dpi=dpi,
                        minimum_planned_text_points=9,all_text_inside_canvas=True,
                        factual_description_metadata=alt,table_alternative=str(table),supported_value_count=supported_values,
                        SVG_accessibility_title_description_role=True,PDF_tagged_accessibility_claimed=False,
                        scientific_inference_scope='Verified marginal estimates and frozen QC; no comparative P or causal claim',
                        visual_QA_pending=True))
    plt.close(fig);print(json.dumps({'figure':name,'values':supported_values,'text_bounds_pass':True}),flush=True)


def maps(data,label_lookup,stem,title,nper,records,source):
    outcomes=list(dict.fromkeys(r['outcome_trait'] for r in data));lookup={(r['sleep_trait'],r['outcome_trait']):r for r in data}
    if len(lookup)!=len(data) or len(data)!=12*len(outcomes):raise RuntimeError('COMPLETE_MATRIX_CARDINALITY_REQUIRED')
    lim=max(1.,max(abs(float(r['rg'])) for r in data));lim=np.ceil(lim*10)/10
    semax=max(float(r['se']) for r in data)
    for page,start in enumerate(range(0,len(outcomes),nper),1):
        selected=outcomes[start:start+nper]
        point=np.array([[float(lookup[(sleep,t)]['rg']) for sleep in SLEEPS] for t in selected])
        se=np.array([[float(lookup[(sleep,t)]['se']) for sleep in SLEEPS] for t in selected])
        if not np.isfinite(point).all() or not np.isfinite(se).all() or not (se>0).all():raise RuntimeError('FINITE_COMPLETE_POINTS_AND_UNCERTAINTY_REQUIRED')
        fig,axes=plt.subplots(1,2,figsize=(7.2,9.2),gridspec_kw={'width_ratios':[1,1]})
        fig.subplots_adjust(left=.385,right=.97,top=.845,bottom=.255,wspace=.12)
        plots=[axes[0].imshow(point,cmap='RdBu_r',vmin=-lim,vmax=lim,aspect='auto',interpolation='none'),
               axes[1].imshow(se,cmap='Blues',vmin=0,vmax=semax,aspect='auto',interpolation='none')]
        for i,t in enumerate(selected):
            historical_sensitivity=lookup[(SLEEPS[0],t)].get('original_analysis_tier','PRIMARY_PHASE1')!='PRIMARY_PHASE1'
            for ax in axes:
                if historical_sensitivity:ax.add_patch(Rectangle((-.5,i-.5),12,1,fill=False,hatch='////',edgecolor='#59636c',lw=0))
            for j,sleep in enumerate(SLEEPS):
                if lookup[(sleep,t)]['native_fdr_pass_0_05']=='True':axes[0].plot(j,i,'o',color='#111111',ms=2.5)
        for j,ax in enumerate(axes):
            ax.set_xticks(range(12),SLEEP_LABELS,rotation=78,ha='right',fontsize=9)
            ax.set_yticks(range(len(selected)),['\n'.join(textwrap.wrap(label_lookup[t],27)) for t in selected] if j==0 else ['']*len(selected),fontsize=9)
            ax.tick_params(length=0);ax.set_title('A  Genetic correlation' if j==0 else 'B  Standard error',fontsize=10,loc='left',pad=10)
        for ax,plot,label in zip(axes,plots,['rg','SE of rg']):
            bar=fig.colorbar(plot,ax=ax,orientation='horizontal',fraction=.035,pad=.245,shrink=.93)
            bar.set_ticks([-lim,0,lim] if label=='rg' else [0,semax/2,semax])
            bar.set_label(label,fontsize=9);bar.ax.tick_params(labelsize=9)
        fig.text(.045,.971,title,fontsize=12,weight='bold',va='top')
        fig.text(.045,.925,'Complete frozen family; page '+str(page)+' of '+str((len(outcomes)+nper-1)//nper)+'\n'+str(len(selected)*12)+' marginal estimates on this page; all rows retained.',fontsize=9,va='top')
        note='Dots: original complete-family BH q < 0.05. Both panels use\nfixed scales across pages. Exact values and marginal 95% intervals\nare in the table. No covariance-adjusted contrast is shown.'
        if stem.startswith('02'):note+='\nHatching: original sensitivity tier. Lipid source QC remains qualified.'
        else:note+='\nAct. = actigraphy. Source and ancestry qualifications remain.'
        fig.text(.045,.035,note,fontsize=9,va='bottom')
        alt=title+'. Page '+str(page)+' contains '+str(len(selected))+' outcome rows and all12 sleep columns in frozen order, showing rg and SE. '+note.replace('\n',' ')
        export(fig,stem+'_page%02d'%page,alt,source,len(selected)*12,records)


def forests(data,stem,title,labels,nper,records,source,sensitivity=False):
    for page,start in enumerate(range(0,len(data),nper),1):
        selected=data[start:start+nper];fig,ax=plt.subplots(figsize=(7.2,9.2))
        fig.subplots_adjust(left=.50,right=.965,top=.845,bottom=.185)
        for i,row in enumerate(selected):
            if sensitivity:
                for arm,offset,color,marker in [('baseline',-.12,'#607887','o'),('sensitivity',.12,'#b54e25','s')]:
                    point=float(row[arm+'_rg_ratio']);lo=float(row[arm+'_CI95_lower']);hi=float(row[arm+'_CI95_upper'])
                    ax.errorbar(point,i+offset,xerr=[[point-lo],[hi-point]],fmt=marker,color=color,ms=3.2,lw=.9,capsize=2)
            else:
                point=float(row['rg']);lo=float(row['ci_lower_95']);hi=float(row['ci_upper_95'])
                ax.errorbar(point,i,xerr=[[point-lo],[hi-point]],fmt='o',color='#176c92',ms=3.5,lw=1,capsize=2)
        ax.axvline(0,color='#505a63',lw=.75);ax.set_yticks(range(len(selected)),['\n'.join(textwrap.wrap(label,39)) for label in labels[start:start+nper]],fontsize=9)
        ax.set_ylim(len(selected)-.5,-.5);ax.set_xlabel('rg and marginal 95% interval',fontsize=9)
        ax.grid(axis='x',color='#dee3e6',lw=.5);ax.tick_params(axis='y',length=0)
        fig.text(.045,.971,title,fontsize=12,weight='bold',va='top')
        fig.text(.045,.925,'Page '+str(page)+' of '+str((len(data)+nper-1)//nper)+'; all '+str(len(data))+' predefined estimates retained across pages.',fontsize=9,va='top')
        if sensitivity:
            fig.text(.045,.88,'Blue circles: baseline. Orange squares: frozen sensitivity.',fontsize=9)
            note='Marginal intervals compare estimators descriptively. Shared-source\ncovariance is uncalibrated; no comparative P, new BH family, or\nsource clearance is inferred. Original q values/exclusions are retained.'
        else:
            note='These41 historical external-outcome estimates reuse discovery sleep\nGWAS. Zero fully independent two-trait replications are established.\nThe full217-candidate classification and source/QC limits are in the tables.'
        fig.text(.045,.055,note,fontsize=9,va='bottom')
        export(fig,stem+'_page%02d'%page,title+'. '+note.replace('\n',' '),source,len(selected),records)


def main():
    if RECORD.exists() or OUT.exists():raise RuntimeError('PRIOR_FIGURE_BATCH_PRESERVED_NO_OVERWRITE')
    inputs=[P/'tables'/name for name in ['core_native_full_precision_rg_v4.tsv','extension_native_full_precision_rg_v4.tsv',
        'validation_native_full_precision_rg_v4.tsv','independent_sensitivity_62_full_precision_with_CI_v4_6.tsv',
        'independent_sensitivity_liability_presentations_v4_6.tsv']]
    metadata=P.parent/'sleep_unified_research_v1/tables/phenotype_and_source_metadata.tsv';inputs += [metadata,Path(__file__)]
    expected_reviews={'independent_whole_extension_v4.sha256':'eb4051e7f26e76de71883c42ab0c4f4567f9c54749fc3ab9296e2ab23accaa22',
        'independent_whole_validation_v4.sha256':'88886794f090fa12cca1b1ffdb426fecb2086e53a3fbe9cca4c8c67b6b959101',
        'independent_sensitivity_whole62_review_seal_v4_6.json':'b85752b4a6ab74c2d8e9a69177b66738e3877194a0c065635afddd93ddc65619'}
    for name,digest in expected_reviews.items():
        path=P/'reviews'/name
        if sha(path)!=digest:raise RuntimeError('VERIFIED_RESEARCH_REVIEW_CHANGED')
        inputs.append(path)
    inputs += [P/'logs/native190_root_independent_completion_addendum_v4.json',P/'logs/sensitivity_root_independent_adjudication_v4_6.json']
    frozen={str(path):sha(path) for path in inputs}
    core,extension,validation,sens,liability=[rows(path) for path in inputs[:5]]
    if (len(core),len(extension),len(validation),len(sens),len(liability))!=(396,1200,41,62,10):raise RuntimeError('FROZEN_COMPLETE_TABLE_CARDINALITY_DIFFERS')
    if sum(x['native_fdr_pass_0_05']=='True' for x in core)!=161 or sum(x['native_fdr_pass_0_05']=='True' for x in extension)!=603:
        raise RuntimeError('ORIGINAL_COMPLETE_FAMILY_BOUNDARY_DIFFERS')
    labels={x['trait_id']:x['phenotype'] for x in rows(metadata)}
    extension_labels={x['outcome_trait']:x['original_phenotype_name'] for x in extension}
    source_design=json.loads((P/'statistical_validation/validation_raw_replay_design_v1.json').read_text())
    source_labels={m['source_id']:m['frozen_source_metadata']['replication_phenotype_definition'] for m in source_design['members']}
    frozen[str(P/'statistical_validation/validation_raw_replay_design_v1.json')]=sha(P/'statistical_validation/validation_raw_replay_design_v1.json')
    OUT.mkdir();records=[]
    plt.rcParams.update({'font.family':'DejaVu Sans','font.size':9,'axes.titlesize':10,'pdf.fonttype':42,'svg.fonttype':'none','axes.spines.top':False,'axes.spines.right':False})
    maps(core,labels,'02_core_native_compact_v4','Core atlas: points and uncertainty',11,records,inputs[0])
    maps(extension,extension_labels,'03_extension_native_compact_v4','Phenome extension: points and uncertainty',13,records,inputs[1])
    forests(validation,'04_qualified_external_validation_v4','Qualified external-outcome validation',
            [labels[x['sleep_trait']]+' x '+source_labels[x['outcome_trait']] for x in validation],14,records,inputs[2])
    for group,stem,title in [('lipid','05a_lipid_sensitivity_v4','Lipid estimator sensitivity'),('binary','05b_binary_N_sensitivity_v4','MS/melanoma sample-size sensitivity')]:
        selected=[x for x in sens if x['kind']=='rg' and (x['outcome_trait'] in ['ldl','hdl','triglycerides'])==(group=='lipid')]
        if len(selected)!=(36 if group=='lipid' else 24):raise RuntimeError('ALL_FROZEN_SENSITIVITY_PAIRS_REQUIRED')
        forests(selected,stem,title,[labels[x['sleep_trait']]+' x '+labels[x['outcome_trait']] for x in selected],12,records,inputs[3],True)
    fig,ax=plt.subplots(figsize=(7.2,8.0));fig.subplots_adjust(left=.50,right=.965,top=.84,bottom=.21)
    presentation_labels=[]
    for i,x in enumerate(liability):
        point=float(x['presented_h2']);lo=float(x['presented_CI95_lower']);hi=float(x['presented_CI95_upper'])
        ax.errorbar(point,i,xerr=[[point-lo],[hi-point]],fmt='o',color='#176c92' if x['point_in_0_1']=='True' else '#b54e25',ms=4,lw=1,capsize=2)
        presentation_labels.append('\n'.join(textwrap.wrap(labels[x['trait']]+': '+x['N_arm'].replace('_',' ')+'; '+x['P_convention'].replace('_',' '),38)))
    ax.set_yticks(range(10),presentation_labels,fontsize=9);ax.set_ylim(9.5,-.5);ax.axvspan(0,1,color='#edf5f4');ax.set_xlabel('Liability presentation and marginal 95% interval')
    ax.grid(axis='x',color='#dee3e6',lw=.5);ax.tick_params(axis='y',length=0)
    fig.text(.045,.965,'Binary liability presentations',fontsize=12,weight='bold',va='top')
    fig.text(.045,.915,'Original K values retained;10 descriptive transformations\nfrom original or total-N h2. Background marks zero to one.',fontsize=9,va='top')
    note='A physical value does not validate population prevalence or ascertainment.\nMS h2 Z remains about5.00; melanoma about3.57 and its DROP remains.\nNo new rg, population variance-fraction claim, or source admission follows.'
    fig.text(.045,.04,note,fontsize=9,va='bottom')
    export(fig,'05c_liability_presentations_v4','Ten liability-scale diagnostic presentations; '+note.replace('\n',' '),inputs[4],10,records)
    fig,ax=plt.subplots(figsize=(7.2,7.8));ax.set_axis_off()
    boxes=[(.055,.72,.40,.16,'Core\n45 traits /57 native jobs\n396 rg +45 h2\n161 historical BH positives'),
           (.545,.72,.40,.16,'Extension\n100 traits /112 native jobs\n1,200 rg +100 h2\n603 historical BH positives'),
           (.055,.46,.89,.16,'Validation\n217 candidates /21 native jobs /13 external sources\n41 rg +13 h2\n23 qualified historical positives;18 directional'),
           (.055,.20,.89,.16,'Validation limitations\n17 QC/power failures +159 unavailable candidates\n0 fully independent two-trait replications\nNo primary novelty or causal claim admitted')]
    for x,y,w,h,t in boxes:
        ax.add_patch(FancyBboxPatch((x,y),w,h,boxstyle='round,pad=.012',facecolor='#eef3f6',edgecolor='#547083',lw=1,transform=ax.transAxes))
        ax.text(x+w/2,y+h/2,t,transform=ax.transAxes,ha='center',va='center',fontsize=10,linespacing=1.55)
    ax.annotate('',xy=(.5,.65),xytext=(.5,.71),xycoords='axes fraction',arrowprops={'arrowstyle':'->','color':'#547083'})
    ax.annotate('',xy=(.5,.39),xytext=(.5,.45),xycoords='axes fraction',arrowprops={'arrowstyle':'->','color':'#547083'})
    fig.text(.065,.96,'Frozen native evidence flow',fontsize=13,weight='bold',va='top')
    fig.text(.065,.90,'190 completed native jobs;158 h2 and1,637 rg estimates.\nIndependent arithmetic verification retains scientific qualifications.',fontsize=9,va='top')
    fig.text(.065,.035,'The12 sleep constructs span self-report, clinical labels and actigraphy.\nThe test families remain396,1,200 and217. Numerical reproduction\nand historical significance do not independently clear source QC.',fontsize=9,va='bottom')
    export(fig,'01_native_evidence_flow_v4','Frozen native evidence flow: '+str(190)+' jobs and1795 estimates; historical classifications and zero fully independent two-trait replications remain.',inputs[2],1795,records)
    if len(records)!=21:raise RuntimeError('EXPECTED21_FIGURE_PAGES_REQUIRED')
    for path,digest in frozen.items():
        if sha(path)!=digest:raise RuntimeError('FIGURE_SOURCE_CHANGED')
    guard()
    record=dict(schema='complete_scientific_figure_exports_v4',created_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),
                figure_page_count=21,export_count=63,source_sha256=frozen,figures=records,
                environment=dict(python=platform.python_version(),matplotlib=matplotlib.__version__,numpy=np.__version__,Pillow=Image.__version__),
                peak_self_RSS_bytes=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,visual_QA_pending=True,
                new_estimator_fits=0,new_tests=0,calibrated_contrast_P_values=0,manuscript_legends_authored=False)
    with RECORD.open('x') as f:json.dump(record,f,indent=2,allow_nan=False);f.write('\n')
    print(json.dumps({'pages':21,'exports':63,'manifest_sha256':sha(RECORD),'visual_QA_pending':True}),flush=True)


if __name__=='__main__':main()
