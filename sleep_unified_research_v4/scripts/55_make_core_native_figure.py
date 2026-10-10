#!/usr/bin/env python3
"""Core point and uncertainty maps from independently reviewed native estimates."""
import csv
from datetime import datetime, timezone
import hashlib
import json
import os
from pathlib import Path
import platform
import textwrap
import xml.etree.ElementTree as ET

SSD = Path('/Volumes/Extreme SSD/sleep-unified-research-v1/research-completion-2026-10-09/sleep_unified_research_v4')
os.environ.setdefault('MPLCONFIGDIR', str(SSD/'cache/matplotlib'))
os.environ.setdefault('OPENBLAS_NUM_THREADS', '1')
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
import numpy as np
from PIL import Image

PACKAGE = Path(__file__).resolve().parents[1]
ROOT = PACKAGE.parent


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def rows(path):
    with path.open(newline='') as handle:
        return list(csv.DictReader(handle, delimiter='\t'))


def main():
    source = PACKAGE/'tables/core_native_full_precision_rg_v4.tsv'
    metadata = ROOT/'sleep_unified_research_v1/tables/phenotype_and_source_metadata.tsv'
    review = PACKAGE/'reviews/independent_core_numerical_adjudication_v4.json'
    reviewed = json.loads(review.read_text())
    assert reviewed['verdict'] == 'PROCESSED_INPUT_CORE_NUMERICAL_PASS_WITH_PRECISION_ADJUDICATION'
    assert reviewed['fdr_count'] == 161 and reviewed['primary_fdr_count'] == 153 and reviewed['sensitivity_fdr_count'] == 8
    inputs = {str(path):sha(path) for path in (source, metadata, review, Path(__file__))}
    data = rows(source)
    labels = {row['trait_id']:row['phenotype'] for row in rows(metadata)}
    sleeps = list(dict.fromkeys(row['sleep_trait'] for row in data))
    outcomes = list(dict.fromkeys(row['outcome_trait'] for row in data))
    lookup = {(row['sleep_trait'],row['outcome_trait']):row for row in data}
    assert (len(data),len(lookup),len(sleeps),len(outcomes)) == (396,396,12,33)
    assert all(row['reproduction_status']=='PRINTED_PRECISION_CONCORDANT' for row in data)
    matrix = np.array([[float(lookup[(sleep,trait)]['rg']) for sleep in sleeps] for trait in outcomes])
    se = np.array([[float(lookup[(sleep,trait)]['se']) for sleep in sleeps] for trait in outcomes])
    assert np.isfinite(matrix).all() and np.isfinite(se).all() and (se>0).all()
    assert sum(row['native_fdr_pass_0_05']=='True' for row in data) == 161
    plt.rcParams.update({'font.family':'DejaVu Sans','font.size':9,'axes.titlesize':12,
                         'pdf.fonttype':42,'svg.fonttype':'none','axes.spines.top':False,'axes.spines.right':False})
    fig,axes = plt.subplots(1,2,figsize=(15,12),gridspec_kw={'width_ratios':[1,1]})
    plt.subplots_adjust(left=.22,right=.94,top=.93,bottom=.22,wspace=.16)
    image_a = axes[0].imshow(matrix,cmap='RdBu_r',vmin=-1,vmax=1,aspect='auto',interpolation='none')
    image_b = axes[1].imshow(se,cmap='Blues',vmin=0,vmax=float(se.max()),aspect='auto',interpolation='none')
    for i,trait in enumerate(outcomes):
        sensitivity = lookup[(sleeps[0],trait)]['original_analysis_tier'] != 'PRIMARY_PHASE1'
        for ax in axes:
            if sensitivity:
                ax.add_patch(Rectangle((-.5,i-.5),12,1,fill=False,hatch='////',edgecolor='#727272',lw=0))
        for j,sleep in enumerate(sleeps):
            row = lookup[(sleep,trait)]
            if row['native_fdr_pass_0_05']=='True':
                axes[0].plot(j,i,'o',ms=3.2,color='#727272' if sensitivity else '#111111')
    ylabels = [textwrap.fill(labels[trait],38)+(' †' if trait in ('ldl','hdl','triglycerides') else '') for trait in outcomes]
    for index,ax in enumerate(axes):
        ax.set_xticks(range(12),[labels[sleep] for sleep in sleeps],rotation=58,ha='right',fontsize=9)
        ax.set_yticks(range(33),ylabels if index==0 else ['']*33,fontsize=9)
        ax.tick_params(length=0)
        ax.set_title('A   Global genetic correlation',loc='left',pad=12) if index==0 else ax.set_title('B   Standard error',loc='left',pad=12)
    bar_a = fig.colorbar(image_a,ax=axes[0],orientation='horizontal',fraction=.035,pad=.19,shrink=.82)
    bar_a.set_label('rg');bar_a.set_ticks([-1,-.5,0,.5,1])
    bar_b = fig.colorbar(image_b,ax=axes[1],orientation='horizontal',fraction=.035,pad=.19,shrink=.82)
    bar_b.set_label('SE of rg')
    fig.suptitle('Core atlas: native points and uncertainty · frozen 396-test family',x=.22,ha='left',fontsize=15,y=.978)
    fig.text(.22,.956,'45 traits · 153 historical primary positives + 8 sensitivity positives · printed-precision reproduction',ha='left',fontsize=10,color='#44515b')
    fig.text(.22,.065,'Dots: BH q < 0.05 within all 396 tests. Hatching: historical sensitivity rows.\n'
             '† All 36 sleep–lipid pairwise outcome intercepts exceed 1.2; scientific clearance remains unresolved.\n'
             'Marginal 95% intervals and source/QC details are in the complete numerical table; measurement contrasts are untested.',
             ha='left',va='bottom',fontsize=9,color='#283540')
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    outside = []
    for text in fig.findobj(matplotlib.text.Text):
        if not text.get_visible() or not text.get_text():continue
        bounds = text.get_window_extent(renderer)
        if bounds.x0 < -1 or bounds.y0 < -1 or bounds.x1 > fig.bbox.width+1 or bounds.y1 > fig.bbox.height+1:
            outside.append(text.get_text())
    assert not outside, 'TEXT_OUTSIDE_CANVAS: '+repr(outside)
    destination = PACKAGE/'figures'
    destination.mkdir(exist_ok=True)
    name = '02_core_native_points_and_uncertainty_v4'
    paths = [destination/(name+'.'+suffix) for suffix in ('pdf','svg','png')]
    assert not any(path.exists() for path in paths)
    for path in paths:
        fig.savefig(path,dpi=300,facecolor='white')
    plt.close(fig)
    alt = ('Two maps show all396 native sleep-by-outcome genetic correlations and their standard errors, with the same12 sleep columns and33 outcome rows. '
           'Dots identify161 positives under the frozen396-test BH family, including153 historical primary and8 sensitivity results. '
           'Type2 diabetes and melanoma retain historical sensitivity hatching. LDL, HDL and triglycerides carry an intercept warning. '
           'These are reproduced marginal estimates; no direct measurement contrast, causal inference or corrected shared-source heterogeneity is shown. '
           'The accompanying complete table supplies exact values, marginal95% intervals, source identities and QC.')
    namespace = 'http://www.w3.org/2000/svg'
    ET.register_namespace('',namespace)
    tree = ET.parse(paths[1]);svg = tree.getroot()
    title = ET.Element('{'+namespace+'}title',{'id':'core-native-title'});title.text='Core native genetic correlations and uncertainty'
    desc = ET.Element('{'+namespace+'}desc',{'id':'core-native-desc'});desc.text=alt
    svg.insert(0,title);svg.insert(1,desc);svg.set('role','img');svg.set('aria-labelledby','core-native-title core-native-desc')
    tree.write(paths[1],encoding='utf-8',xml_declaration=True)
    with Image.open(paths[2]) as png:
        dimensions,dpi=png.size,png.info['dpi']
    assert dimensions == (4500,3600) and all(abs(value-300)<.02 for value in dpi)
    assert all(sha(path)==expected for path,expected in inputs.items())
    record = {'created_utc':datetime.now(timezone.utc).isoformat(),'figure':name,
              'source_sha256':inputs,'output_sha256':{str(path):sha(path) for path in paths},
              'environment':{'python':platform.python_version(),'matplotlib':matplotlib.__version__,
                             'numpy':np.__version__,'Pillow':Image.__version__},
              'PNG_dimensions_pixels':dimensions,'PNG_dpi':dpi,'export_size_inches':[15,12],
              'all_text_inside_canvas':True,'minimum_planned_label_font_points':9,
              'factual_description_metadata':alt,'table_alternative':str(source),
              'SVG_title_description_role_present':True,'PDF_tagged_accessibility_claimed':False,
              'selection_reason':'Original figure inspected:220dpi and historical rounded values; new plot adds native precision and SE panel.',
              'scientific_scope':'Processed-input native printed-precision reproduction; source/QC qualifications retained.',
              'new_biological_inference':False,'visual_QA_pending':True}
    out = PACKAGE/'manifests/core_native_figure_export_v4.json'
    with out.open('x') as handle:json.dump(record,handle,indent=2);handle.write('\n')
    print(json.dumps({'figure':name,'exports':3,'dpi':dpi,'pixels':dimensions,'text_bounds_pass':True}))


if __name__=='__main__':main()
