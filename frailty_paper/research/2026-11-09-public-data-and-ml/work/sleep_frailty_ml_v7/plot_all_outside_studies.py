"""Display all three outside cell contexts without changing inference."""
from common import *
os.environ.setdefault('MPLCONFIGDIR', str(DATA/'runtime_tmp/matplotlib'))
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt


def main():
    primary = json.loads((OUT/'reserved_study_results.json').read_text())
    extra = json.loads((OUT/'additional_neural_results.json').read_text())
    plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 9, 'axes.spines.top': False, 'axes.spines.right': False})
    labels = ['Neural progenitors\nMcAfee\nn = 4,817', 'HEK293\nMcAfee\nn = 4,785', 'Baseline neurons\nRetallick-Townsley\nn = 547']
    xx = np.arange(3)
    fig, axes = plt.subplots(1, 3, figsize=(13, 4.4), layout='constrained')
    source = []
    for i, (arm, label, color) in enumerate([('sequence_only', 'DNA and cell inputs', '#2166ac'), ('AVI_added', 'AVI added', '#d95f02')]):
        metrics = [primary['arms'][arm]['cells'][c] for c in ['HNP', 'HEK293']] + [extra['conditions']['baseline'][arm]]
        errors = [m['MAE'] for m in metrics]
        bounds = np.maximum([[m['MAE']-m['MAE_genomic_window_CI95'][0] for m in metrics], [m['MAE_genomic_window_CI95'][1]-m['MAE'] for m in metrics]], 0)
        x = xx+(i-.5)*.34
        axes[0].bar(x, errors, .34, color=color, label=label, yerr=bounds, capsize=3)
        axes[1].bar(x, [m['interval_coverage']*100 for m in metrics], .34, color=color)
        axes[2].bar(x, [m['ranked_fraction']*100 for m in metrics], .34, color=color)
        for j, m in enumerate(metrics):
            source.append({'context': labels[j].replace('\n', ' '), 'arm': arm, 'n': m['n'], 'MAE': m['MAE'], 'baseline_MAE': m['baseline_MAE'], 'interval_coverage': m['interval_coverage'], 'median_interval_width': m['median_interval_width'], 'frozen_input_policy_accepted_fraction': m['ranked_fraction'], 'final_delivery_accepted_fraction': 0})
    axes[0].scatter(xx, [r['baseline_MAE'] for r in source[:3]], color='black', marker='_', s=260, zorder=3, label='Training-cell median')
    axes[0].set(title='Effect-size error (lower is better)', ylabel='Mean absolute error, log2 units', ylim=(0, .12))
    axes[0].legend(fontsize=8)
    axes[1].axhline(90, color='black', ls='--', lw=1)
    axes[1].set(title='Coverage of broad intervals', ylabel='Observed effects inside interval, %', ylim=(0, 106))
    axes[2].scatter(xx, [0, 0, 0], marker='x', s=65, color='black', zorder=5, label='Final answers released: 0%')
    axes[2].set(title='Original input-policy acceptance', ylabel='Accepted measurements, %', ylim=(-5, 106))
    axes[2].legend(fontsize=8, loc='upper left')
    for ax in axes:
        ax.set_xticks(xx, labels, fontsize=8)
    fig.suptitle('Outside tests: useful prediction was not established', fontsize=13)
    fig.supxlabel('Genomic-window bootstrap error bars; broad coverage and abstention do not demonstrate useful precision.', fontsize=9)
    figures = OUT/'figures'
    image_path = figures/'all_studies_error_and_answer_availability.png'
    fig.savefig(image_path, dpi=200)
    plt.close(fig)
    data_path = figures/'all_studies_error_and_answer_availability.csv'
    with data_path.open('w') as f:
        w = csv.DictWriter(f, fieldnames=list(source[0])); w.writeheader(); w.writerows(source)
    notes = figures/'figure_notes.txt'
    notes.write_text('The original outside-study panels show only the McAfee primary test. The all_studies panel adds the fixed additional neuronal test. All use original fitted magnitudes, not disease-risk outcomes. Error bars are 95% genomic 1-Mb-window bootstrap intervals, not ancestry-specific LD blocks. Measurements within a study or region are correlated. The 90% line is a nominal calibration target; broad ranges can yield high coverage without useful precision. Original input-policy acceptance differs from final delivery: all 547 baseline neuronal measurements passed the former, but no cell group passed the outside-performance requirement for released answers. Both primary cell groups lacked support from two development studies. Learning-curve shading in the separate development figure is the min–max across model seeds; seeds are not biological replicates. No model fitting or retention changes were made to create these figures.\n')
    receipt = json.loads((OUT/'figure_receipt.json').read_text())
    for path in [image_path, data_path, notes]:
        receipt['figures_and_sources'][str(path)] = sha(path)
    receipt['additional_locked_result_sha256'] = sha(OUT/'additional_neural_results.json')
    write_json(OUT/'figure_receipt.json', receipt)
    print('ALL OUTSIDE STUDIES PLOTTED', flush=True)


if __name__ == '__main__':
    main()
