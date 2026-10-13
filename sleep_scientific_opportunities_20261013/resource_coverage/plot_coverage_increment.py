"""Plot authored snapshot classifications, never third-party rg estimates."""
import json
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parent
summary = json.loads((ROOT / 'results_corrected_v3/coverage_summary.json').read_text())
OUT = ROOT / 'figures'
OUT.mkdir(exist_ok=True)
plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 11,
                     'svg.fonttype': 'none', 'axes.spines.top': False,
                     'axes.spines.right': False})
fig, axes = plt.subplots(1, 2, figsize=(11.5, 5.9), gridspec_kw={'width_ratios': [1.4, 1]})
classes = [('B', 'Related prior coverage', '#39789c'),
           ('D', 'Absent within mapped scope', '#d5e4ea'),
           ('E', 'Unresolved', '#b7b4ad')]
for y, field in enumerate(['before_four_class_counts', 'after_five_class_counts']):
    left = 0
    for code, label, color in classes:
        n = summary[field].get(code, 0)
        axes[0].barh(y, n, left=left, color=color, height=.5,
                     label=label if y == 0 else None)
        axes[0].text(left+n/2, y, str(n), ha='center', va='center',
                     color='white' if code == 'B' else '#243742', fontsize=10)
        left += n
axes[0].set(yticks=[0, 1], yticklabels=['Four snapshots', 'Five snapshots'],
            xlim=(0, 1637), xlabel='Atlas records (fixed population: 1,637)',
            title='A  Dated comparator coverage')
axes[0].invert_yaxis()
fig.legend(*axes[0].get_legend_handles_labels(), frameon=False,
           loc='lower left', bbox_to_anchor=(.12, .12), ncol=3, fontsize=9)
vals = [summary['qualified_absent_before']-3,
        summary['qualified_absent_after_outside_hold']]
axes[1].bar([0, 1], vals, width=.58, color=['#bbcdd5', '#39789c'])
for x, n in enumerate(vals):
    axes[1].text(x, n+6, str(n), ha='center', fontsize=12)
axes[1].set(xticks=[0, 1], xticklabels=['Four snapshots', 'Five snapshots'],
            ylim=(0, 410), ylabel='Historically qualified core/extension records',
            title='B  Mapped absence\noutside direction holds')
axes[1].title.set_fontsize(11)
fig.suptitle('Adding SleepChart changes ten coverage classifications',
             x=.12, ha='left', fontsize=15, fontweight='bold')
fig.text(.12, .885, '396 core + 1,200 extension + 41 estimated validation records; no new genetic fits', fontsize=10)
fig.text(.12, .045, 'B is related context, not exact equivalence. D is bounded snapshot/alias absence, not global novelty.\n'
                   'Four comparator annotations are inherited. Counts are not discoveries or significant-result counts.', fontsize=10)
fig.subplots_adjust(left=.12, right=.98, bottom=.30, top=.74, wspace=.46)
fig.savefig(OUT / 'five_snapshot_coverage_increment.svg')
svg = OUT / 'five_snapshot_coverage_increment.svg'
svg.write_text('\n'.join(line.rstrip() for line in svg.read_text().splitlines())+'\n')
fig.savefig(OUT / 'five_snapshot_coverage_increment.png', dpi=180)
plt.close(fig)
(OUT / 'figure_data.json').write_text(json.dumps({
    'source': '../results_corrected_v3/coverage_summary.json',
    'fixed_population': 1637,
    'four_snapshot_union': summary['before_four_class_counts'],
    'five_snapshot_union': summary['after_five_class_counts'],
    'qualified_mapped_absence_outside_holds': dict(zip(['four', 'five'], vals)),
    'four_count_basis': '361 qualified absent minus three preserved direction-held records',
    'changed_records': {'core_colorectal_cancer_duration': 4,
                        'estimated_validation_related_outcomes': 6},
    'new_genetic_estimates': 0,
    'interpretation': 'Bounded retrospective metadata extension; no novelty/discovery certification'
}, indent=2)+'\n')
