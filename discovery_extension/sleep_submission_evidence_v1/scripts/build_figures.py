#!/usr/bin/env python3
"""Deterministic figures from historical observations; never writes original outputs."""
from pathlib import Path
import hashlib
import json
import textwrap
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

PACKAGE = Path(__file__).resolve().parents[1]
ROOT = PACKAGE.parents[1]
OUT = PACKAGE / 'figures'
SRC = OUT / 'source_data'
plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 9,
                     'pdf.fonttype': 42, 'svg.fonttype': 'none', 'svg.hashsalt': 'sleep-phenome-evidence-v1',
                     'axes.spines.top': False, 'axes.spines.right': False})
BLUE, ORANGE = '#0072B2', '#D55E00'
META = {'CreationDate': None, 'ModDate': None, 'Creator': 'sleep-phenome-evidence-v1'}

def read(path):
    return pd.read_csv(ROOT / path, sep='\t')

def save(fig, name):
    fig.savefig(OUT / (name + '.pdf'), metadata=META)
    fig.savefig(OUT / (name + '.svg'), metadata={'Date': None})
    fig.savefig(OUT / (name + '.png'), dpi=300)
    plt.close(fig)

def table(frame, name):
    frame.to_csv(SRC / (name + '.tsv'), sep='\t', index=False, na_rep='NA', float_format='%.17g')

def label(row):
    return row['sleep_trait'].replace('_', ' ') + ' / ' + row['external_phenotype_name']

def forest(rows, name, title):
    rows = rows.sort_values(['external_phenotype_name', 'sleep_trait']).reset_index(drop=True)
    fig, ax = plt.subplots(figsize=(8.4, max(4.6, 2.0 + .40 * len(rows))))
    fig.subplots_adjust(left=.52, right=.97, top=.90, bottom=.23)
    y = np.arange(len(rows))
    ax.errorbar(rows.discovery_rg, y-.10, xerr=1.959963984540054*rows.discovery_se,
                fmt='o', markersize=3.8, color=BLUE, capsize=2, label='Discovery')
    ax.errorbar(rows.replication_rg, y+.10, xerr=1.959963984540054*rows.replication_se,
                fmt='s', markersize=3.8, color=ORANGE, capsize=2, label='Outcome-side validation')
    labels = [textwrap.fill(label(r), 44) + (' [C]' if r['replication_phenotype_definition'].startswith('COMPARABLE') else '')
              for _, r in rows.iterrows()]
    # Definition labels in the original manifest are joined below, not inferred from significance.
    if 'phenotype_match_status' in rows:
        labels = [textwrap.fill(label(r), 44) + (' [C]' if 'COMPARABLE' in str(r['phenotype_match_status']) else '')
                  for _, r in rows.iterrows()]
    ax.set_yticks(y, labels, fontsize=8)
    ax.invert_yaxis()
    ax.axvline(0, color='#777777', linewidth=.7)
    ax.set_xlabel('Genetic correlation (rg), 95% normal interval')
    ax.set_xlim(min(0, (rows.discovery_rg-1.96*rows.discovery_se).min(), (rows.replication_rg-1.96*rows.replication_se).min())-.06,
                max((rows.discovery_rg+1.96*rows.discovery_se).max(), (rows.replication_rg+1.96*rows.replication_se).max())+.06)
    ax.grid(axis='x', color='#dddddd', linewidth=.5)
    fig.suptitle(title, fontsize=12, x=.04, ha='left')
    fig.legend(loc='center', bbox_to_anchor=(.73, .11), ncol=2, frameon=False, fontsize=8)
    footer = 'Sleep GWAS reused. [C] = comparable definition. Estimates have archived display precision.'
    if name == 'fig4_heterogeneity':
        footer += '\nSelection: nominal heterogeneity P < 0.05, with zero covariance assumed.'
    fig.text(.04, .025, footer, fontsize=8)
    rows['discovery_ci_low'] = rows.discovery_rg-1.959963984540054*rows.discovery_se
    rows['discovery_ci_high'] = rows.discovery_rg+1.959963984540054*rows.discovery_se
    rows['replication_ci_low'] = rows.replication_rg-1.959963984540054*rows.replication_se
    rows['replication_ci_high'] = rows.replication_rg+1.959963984540054*rows.replication_se
    rows['display_order'] = y
    table(rows, name)
    save(fig, name)

def main():
    OUT.mkdir(parents=True, exist_ok=True); SRC.mkdir(exist_ok=True)
    d = read('discovery_extension/results/ldsc/extension_rg_matrix.tsv')
    r = read('discovery_extension/results/replication/replication_results.tsv')
    panel = read('discovery_extension/config/candidate_traits.tsv')
    manifest = read('discovery_extension/config/replication_manifest.tsv')
    matching = [c for c in manifest if 'match' in c or 'compatib' in c]
    if matching:
        chosen = 'phenotype_match_status' if 'phenotype_match_status' in matching else matching[0]
        r = r.merge(manifest[['pair_id', chosen]].rename(columns={chosen:'phenotype_match_status'}), on='pair_id', how='left', validate='one_to_one')
        assert r.phenotype_match_status.notna().all(), 'Missing phenotype comparison'
    positives = r[r.replication_class=='REPLICATED'].copy()
    sleeps = read('config/analysis_panel.tsv').query("domain == 'sleep'").trait_id.tolist()
    order = panel.sort_values(['phenotype_domain','phenotype_name']).extension_trait_id.tolist()
    assert len(d)==len(sleeps)*len(panel)
    assert not d.duplicated(['sleep_trait','extension_trait_id']).any()
    assert set(zip(d.sleep_trait,d.extension_trait_id)) == {(s,e) for s in sleeps for e in panel.extension_trait_id}
    assert np.isfinite(d.rg.to_numpy()).all() and d.rg.between(-.6,.6).all()
    mat = d.pivot(index='extension_trait_id', columns='sleep_trait', values='rg').loc[order,sleeps]
    sig = d.pivot(index='extension_trait_id', columns='sleep_trait', values='extension_fdr').loc[order,sleeps] < .05
    table(d.merge(panel[['extension_trait_id','phenotype_definition']],on='extension_trait_id',validate='many_to_one'), 'fig1_atlas')
    descriptions = []
    with PdfPages(OUT/'fig1_atlas.pdf', metadata=META) as pdf:
        for page, start in enumerate(range(0, len(panel), 50), 1):
            ids = order[start:start+50]
            fig, ax = plt.subplots(figsize=(8.4, 11.7))
            fig.subplots_adjust(left=.53, right=.96, top=.94, bottom=.18)
            im = ax.imshow(mat.loc[ids], cmap='RdBu_r', vmin=-.6, vmax=.6, aspect='auto')
            named = panel.set_index('extension_trait_id').loc[ids]
            duplicates = set(panel.loc[panel.phenotype_name.duplicated(False), 'phenotype_name'])
            names = [str(row.phenotype_name) + (' / ' + str(row.phenotype_definition).split(' | ')[1]
                     if row.phenotype_name in duplicates else '') for _, row in named.iterrows()]
            ax.set_yticks(range(len(ids)), [textwrap.shorten(str(x), width=57, placeholder='...') for x in names], fontsize=8)
            ax.set_xticks(range(len(sleeps)),[x.replace('_',' ') for x in sleeps],rotation=65,ha='right',fontsize=8)
            yy, xx = np.where(sig.loc[ids].to_numpy()); ax.scatter(xx, yy, s=6, c='black', marker='o')
            ax.set_title(f'Phenome-wide genetic correlation atlas ({page}/2)', fontsize=12, loc='left', pad=16)
            fig.colorbar(im, ax=ax, orientation='horizontal',fraction=.03,pad=.14,label='Genetic correlation (rg)')
            fig.text(.04,.025,'Black dots: BH FDR < 0.05 across the fixed 1,200-test family.\n'
                     'Rows grouped by phenotype domain; correlated phenotypes are retained. Full labels in source table.',fontsize=8)
            pdf.savefig(fig)
            fig.savefig(OUT/f'fig1_atlas_page{page}.png',dpi=300)
            fig.savefig(OUT/f'fig1_atlas_page{page}.svg',metadata={'Date':None})
            plt.close(fig)
    forest(positives,'fig2_outcome_replication',f'{len(positives)} sleep–phenotype pairs with outcome-side validation')
    forest(positives[positives.heterogeneity_p<.05], 'fig4_heterogeneity', 'Nominally heterogeneous positive pairs')
    counts = r.replication_class.value_counts().to_dict()
    n_candidate = len(read('discovery_extension/results/candidate_pool.tsv'))
    n_sig = int((d.extension_fdr<.05).sum())
    flow = pd.DataFrame([('Candidate phenotypes',n_candidate,'Prospective source/eligibility pool'),
         ('Selected phenotypes',len(panel),'Ordered frozen panel'),('Discovery tests',len(d),'12 sleep traits × 100 phenotypes'),
         ('FDR-significant pairs',n_sig,'BH within 1,200'),('Replication candidates',len(r),'Selected 217-pair family'),
         ('Outcome-side positives',len(positives),'Bonferroni 0.05 / 217; same sleep input')],columns=['stage','count','definition'])
    table(flow,'fig3_flow')
    fig,ax=plt.subplots(figsize=(8.4,4.3)); ax.axis('off')
    for i,row in flow.iterrows():
        x=.02+(i%3)*.335; y=.75-(i//3)*.40
        ax.text(x,y,f"{row['count']:,}\n{row['stage']}", transform=ax.transAxes,fontsize=13,
                ha='left',va='top',bbox={'boxstyle':'round,pad=.55','facecolor':'#eef3f7','edgecolor':BLUE})
    ax.set_title('Discovery and validation accounting',loc='left',fontsize=13)
    outcome_labels={'REPLICATED':'outcome-side positives','DIRECTIONALLY_CONCORDANT':'directionally concordant',
                    'UNDERPOWERED':'underpowered / QC-ineligible','NO_INDEPENDENT_DATASET':'without an eligible source'}
    ax.text(.02,-.04,'Replication outcomes: '+', '.join(f'{v} {outcome_labels.get(k,k)}' for k,v in sorted(counts.items())),
            transform=ax.transAxes,fontsize=8,wrap=True)
    fig.subplots_adjust(left=.05,right=.96,top=.86,bottom=.19); save(fig,'fig3_flow')
    domain = d.assign(significant=d.extension_fdr<.05).groupby('phenotype_domain').agg(
        tests=('rg','size'),significant_pairs=('significant','sum'),phenotypes=('extension_trait_id','nunique')).reset_index()
    domain['fraction_significant']=domain.significant_pairs/domain.tests
    table(domain,'fig5_domains')
    fig,ax=plt.subplots(figsize=(8.4,4.8)); y=np.arange(len(domain))
    ax.barh(y,domain.tests,color='#dfe6ed',label='All planned comparisons')
    ax.barh(y,domain.significant_pairs,color=BLUE,label='BH FDR < 0.05')
    ax.set_yticks(y,[x.replace('_',' ') for x in domain.phenotype_domain]);ax.invert_yaxis()
    ax.set_xlabel('Number of sleep–phenotype comparisons');ax.set_title('Discovery by phenotype domain',loc='left')
    ax.legend(frameon=False,fontsize=8);fig.subplots_adjust(left=.30,bottom=.19,right=.96,top=.9)
    fig.text(.04,.035,'Counts describe comparisons. Correlated phenotypes and reused sleep traits limit biological independence.',fontsize=8)
    save(fig,'fig5_domains')
    table(flow,'graphical_abstract')
    fig,ax=plt.subplots(figsize=(8.4,4.0));ax.axis('off')
    ax.text(.03,.95,'Sleep GWAS phenome-wide extension',fontsize=18,va='top')
    cells=[('Discovery',f'{len(sleeps)} sleep traits\n{len(panel)} phenotypes\n{len(d):,} genetic correlations'),
           ('Global associations',f'{n_sig} BH FDR-positive pairs\nSeparate exploratory family'),
           ('Outcome-side validation',f'{len(positives)} positive pairs\n{positives.extension_trait_id.nunique()} external phenotypes\nSleep GWAS reused')]
    for i,(title,body) in enumerate(cells):
        ax.text(.03+i*.33,.65,title+'\n\n'+body,fontsize=12,va='top',linespacing=1.45,
                bbox={'boxstyle':'round,pad=.6','facecolor':'#eef3f7','edgecolor':BLUE})
    ax.text(.03,.04,'Global genetic sharing. No causal, locus, gene, mechanism or first-ever novelty claim.',fontsize=10)
    fig.subplots_adjust(left=.03,right=.97,top=.96,bottom=.08);save(fig,'graphical_abstract')
    descriptions=[
      ('fig1_atlas','Two-page matrix of all 1,200 archived genetic correlations. Black dots indicate BH FDR < 0.05. The color scale spans -0.6 to 0.6 and contains every estimate.'),
      ('fig2_outcome_replication','Discovery and outcome-side validation point estimates with 95% normal intervals for the 23 historically replicated pairs. Sleep GWAS inputs are reused. Comparable definitions are marked C.'),
      ('fig3_flow','Accounting of the prospective candidate pool, frozen phenotype panel, discovery family, FDR-positive comparisons, selected replication family and outcome-side positives. Counts have different units and denominators.'),
      ('fig4_heterogeneity','Discovery and outcome-side estimates for seven positive pairs with nominal heterogeneity P < 0.05 under the historical zero-covariance approximation. This is not a family-corrected heterogeneity claim.'),
      ('fig5_domains','Planned comparison counts and BH FDR-positive comparison counts by phenotype domain. Counts are not independent diseases, loci or mechanisms.'),
      ('graphical_abstract','Separate preliminary graphical abstract showing study size and outcome-side validation counts, explicitly retaining reused sleep GWAS and global association limitations.')]
    pd.DataFrame(descriptions,columns=['figure_id','technical_description_and_alt_text']).to_csv(OUT/'FIGURE_METADATA.tsv',sep='\t',index=False)
    receipt={'figure_types':len(descriptions),'atlas_pages':2,'source_rows':len(d),'replicated_source_rows':len(positives),
             'raster_dpi':300,'native_gwas_rerun':False,'source_hashes':{str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in SRC.glob('*.tsv')}}
    receipt['original_input_hashes'] = {name:hashlib.sha256((ROOT/name).read_bytes()).hexdigest() for name in
        ['discovery_extension/results/ldsc/extension_rg_matrix.tsv','discovery_extension/results/replication/replication_results.tsv',
         'discovery_extension/config/candidate_traits.tsv','discovery_extension/config/replication_manifest.tsv','config/analysis_panel.tsv']}
    (OUT/'FIGURE_EXECUTION_RECEIPT.json').write_text(json.dumps(receipt,indent=2)+'\n')

if __name__=='__main__': main()
