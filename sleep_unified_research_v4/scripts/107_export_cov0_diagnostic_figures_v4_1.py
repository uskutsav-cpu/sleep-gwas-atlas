"""Publication-quality diagnostic figures; the Cov=0 intervals are unvalidated."""
from pathlib import Path
import csv
import hashlib
import json
import os
import platform
import resource
import shutil
import sys
import textwrap
import xml.etree.ElementTree as ET
from datetime import datetime, timezone

P = Path(__file__).resolve().parents[1]
CACHE = P.parent.parent.parent/'work/cov0_figure_cache_v4_1'
os.environ.setdefault('MPLCONFIGDIR', str(CACHE))
for key in ('OPENBLAS_NUM_THREADS', 'OMP_NUM_THREADS', 'MKL_NUM_THREADS'):
    os.environ[key] = '1'
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from PIL import Image

TABLE = P/'tables/shared_sleep_cov0_diagnostic41_v4_1.tsv'
TABLE_SHA = '0698b9493637176c03b975bc36e9fbd6e2c5baf6dfdc545486c004df693149bb'
OUT = P/'figures/shared_sleep_cov0_diagnostic_v4_1'
MANIFEST = P/'manifests/shared_sleep_cov0_diagnostic_figure_exports_v4_1.json'
SSD = Path('/Volumes/Extreme SSD/sleep-unified-research-v1/research-completion-2026-10-09/sleep_unified_research_v4')
LABELS = {'insomnia':'Insomnia', 'sleepdur':'Sleep duration', 'shortsleep':'Short sleep',
          'longsleep':'Long sleep', 'chronotype':'Chronotype', 'sleepiness':'Sleepiness',
          'napping':'Napping', 'snoring':'Snoring', 'sleep_apnea':'Sleep apnea',
          'sleep_efficiency':'Actigraphy efficiency', 'accel_sleep_duration':'Actigraphy duration',
          'sleep_timing':'Sleep timing'}


def sha(q):
    return hashlib.sha256(Path(q).read_bytes()).hexdigest()


def guard():
    assert shutil.disk_usage('/System/Volumes/Data').free >= 3 << 30
    assert shutil.disk_usage(SSD).free >= 5 << 30
    assert resource.getrusage(resource.RUSAGE_SELF).ru_maxrss <= 1 << 30


def main():
    assert sha(TABLE) == TABLE_SHA and not MANIFEST.exists() and not OUT.exists()
    guard()
    with TABLE.open(newline='') as stream:
        rows = list(csv.DictReader(stream, delimiter='\t'))
    assert len(rows) == len({r['pair_id'] for r in rows}) == 41
    for r in rows:
        assert r['interpretation'] == 'UNVERIFIED_ZERO_COVARIANCE_DIAGNOSTIC_ONLY'
        assert r['confirmed_effect_difference'] == r['sampling_covariance_established'] == 'False'
    OUT.mkdir()
    plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'pdf.fonttype':42,
                         'svg.fonttype':'none','axes.spines.top':False,'axes.spines.right':False})
    low = min(float(r['difference_CI95_low_cov0_unvalidated']) for r in rows)
    high = max(float(r['difference_CI95_high_cov0_unvalidated']) for r in rows)
    span = max(high-low, 0.5)
    limits = (min(low-0.06*span, -0.03), max(high+0.06*span, 0.03))
    exports = []
    for page, start in enumerate(range(0, len(rows), 21), 1):
        selected = rows[start:start+21]
        wrapped = [textwrap.fill(LABELS[r['sleep_trait']]+' | '+r['phenotype'], width=45,
                                 break_long_words=False, break_on_hyphens=False) for r in selected]
        heights = [0.25*len(s.splitlines())+0.16 for s in wrapped]
        panel_h = sum(heights)+0.2
        fig_h = panel_h+2.0
        fig = plt.figure(figsize=(14, fig_h), facecolor='white')
        fig.text(0.04, 1-0.36/fig_h, 'Discovery–validation effect differences',
                 fontsize=15, fontweight='bold', color='#172b45')
        fig.text(0.04, 1-0.66/fig_h, 'Unverified Cov=0 diagnostic; the same sleep GWAS is reused', fontsize=11)
        fig.text(0.96, 1-0.36/fig_h, f'Page {page}/2', ha='right', fontsize=10)
        bottom = 0.92/fig_h
        axes = [fig.add_axes([x, bottom, w, panel_h/fig_h])
                for x,w in ((0.04,0.36),(0.415,0.355),(0.79,0.17))]
        left, forest, numbers = axes
        for ax in axes:
            ax.set_ylim(panel_h, 0)
            ax.set_yticks([])
        left.set_axis_off(); numbers.set_axis_off()
        left.set_xlim(0,1); numbers.set_xlim(0,1)
        forest.set_xlim(*limits)
        forest.axvline(0, color='#9aa3ad', linewidth=1, zorder=0)
        forest.grid(axis='x', color='#e3e6ea', linewidth=.6)
        forest.set_xlabel('External rg − discovery rg; interval assumes Cov=0', fontsize=10, labelpad=9)
        forest.tick_params(axis='x',labelsize=9.5)
        left.text(0,-.09,'Sleep trait | outcome', fontsize=10, fontweight='bold', clip_on=False)
        numbers.text(0,-.09,'Difference [95% CI0]', fontsize=10, fontweight='bold', clip_on=False)
        cursor=.1
        for r,label,h in zip(selected,wrapped,heights):
            y=cursor+h/2
            d,lo,hi = (float(r[k]) for k in ('difference_external_minus_discovery',
                       'difference_CI95_low_cov0_unvalidated','difference_CI95_high_cov0_unvalidated'))
            left.text(0,y,label,fontsize=9.5,va='center',linespacing=1.3)
            forest.plot([lo,hi],[y,y],color='#263f58',linewidth=1.25,zorder=2)
            forest.plot(d,y,'o',color='#263f58',markersize=4,zorder=3)
            numbers.text(0,y,f'{d:.3f} [{lo:.3f}, {hi:.3f}]',fontsize=9.5,va='center')
            cursor+=h
        fig.text(.04,.28/fig_h,'All 41 estimated locked pairs retained. Shared sampling covariance is unestimated; these intervals cannot establish heterogeneity.',fontsize=9.5)
        fig.text(.04,.08/fig_h,'Table: shared_sleep_cov0_diagnostic41_v4_1.tsv | No new test family, primary claim or fully independent two-trait replication.',fontsize=9.5)
        guard();fig.canvas.draw();renderer=fig.canvas.get_renderer()
        for t in fig.findobj(matplotlib.text.Text):
            if t.get_visible() and t.get_text():
                box=t.get_window_extent(renderer)
                assert box.x0 >= -1 and box.y0 >= -1 and box.x1<=fig.bbox.width+1 and box.y1<=fig.bbox.height+1, t.get_text()
        stem=f'shared_sleep_cov0_diagnostic_page{page}'
        paths=[OUT/(stem+'.'+ext) for ext in ('pdf','svg','png')]
        for q in paths:
            assert not q.exists();fig.savefig(q,dpi=300,facecolor='white')
        ns='http://www.w3.org/2000/svg';ET.register_namespace('',ns)
        tree=ET.parse(paths[1]);svg=tree.getroot()
        title=ET.Element('{'+ns+'}title',id=stem+'-title');title.text='Qualified zero-covariance effect-difference diagnostic'
        desc=ET.Element('{'+ns+'}desc',id=stem+'-desc')
        desc.text=f'Page {page} shows {len(selected)} of all41 existing estimated discovery/external-outcome pairs. Dots are external minus discovery rg; bars use the unverified Cov=0 standard error. The sleep GWAS is reused. No covariance-adjusted inference or confirmed heterogeneity is supported. Exact pair values and qualifications are in the linked41-row TSV.'
        svg.insert(0,title);svg.insert(1,desc);svg.set('role','img');svg.set('aria-labelledby',stem+'-title '+stem+'-desc');tree.write(paths[1],encoding='utf-8',xml_declaration=True)
        with Image.open(paths[2]) as im:
            pixels=im.size;dpi=im.info['dpi']
        assert all(abs(x-300)<.02 for x in dpi)
        exports.append({'page':page,'pair_ids':[r['pair_id'] for r in selected],
            'supported_value_count':len(selected),'output_sha256':{str(q):sha(q) for q in paths},
            'export_size_inches':list(fig.get_size_inches()),'PNG_dimensions_pixels':pixels,'PNG_dpi':dpi,
            'minimum_text_points':9.5,'all_text_inside_canvas':True,
            'SVG_title_description_role':True,'PDF_tagged_accessibility_claimed':False,
            'visual_QA_pending':True,'scientific_scope':'UNVERIFIED_ZERO_COVARIANCE_DIAGNOSTIC_ONLY'})
        plt.close(fig)
        print(json.dumps({'figure_page_exported':page,'pairs':len(selected)}),flush=True)
    assert sha(TABLE)==TABLE_SHA and sum(r['supported_value_count'] for r in exports)==41
    manifest={'schema':'shared_sleep_cov0_diagnostic_figure_exports_v1','created_utc':datetime.now(timezone.utc).isoformat(),
        'script_sha256':sha(Path(__file__)),'table_sha256':TABLE_SHA,'figure_page_count':2,'export_count':6,
        'figures':exports,'environment':{'python':platform.python_version(),'executable':sys.executable,
        'matplotlib':matplotlib.__version__,'numpy':np.__version__,'Pillow':Image.__version__},
        'pycache_writes_disabled':sys.dont_write_bytecode,'peak_self_RSS_bytes':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
        'manuscript_legends_authored':False,'new_estimator_fits':0,'new_BH_or_other_correction':False,
        'sampling_covariance_established':False,'confirmed_heterogeneity':False,'visual_QA_pending':True}
    with MANIFEST.open('x') as stream:stream.write(json.dumps(manifest,indent=2,allow_nan=False)+'\n')
    print(json.dumps({'manifest_sha256':sha(MANIFEST),'pages':2,'exports':6,'visual_QA_pending':True}))


if __name__=='__main__':
    main()
