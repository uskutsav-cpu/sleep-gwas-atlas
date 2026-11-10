"""Write only results that are present; never infer successful training."""
from common import *
from uncertainty import estimate
estimate()
sel=json.loads((OUT/'model_selection.json').read_text());records=json.loads((OUT/'search_results.json').read_text());bench=json.loads((OUT/'reused_benchmark_results.json').read_text());ext=json.loads((OUT/'external_assay_results.json').read_text()) if (OUT/'external_assay_results.json').exists() else None
student=json.loads((OUT/'student_status.json').read_text());avi=json.loads((OUT/'avi_retrieval_audit.json').read_text());selected=sel['selected_model']
student_config=json.loads((OUT/'student_config.json').read_text()) if (OUT/'student_config.json').exists() else {}
extension_notice='The auxiliary native-weight experiment was extended after observing round1 benchmark and separate-study results, despite its original validation trigger not firing. The architecture/loss settings were already specified, but this execution decision is exploratory. Round1 is preserved under round1/. Both existing evaluation sets are reused development comparisons for this extension; a further separate source is required for a fresh evaluation.' if student_config.get('post_round1_extension') else 'The auxiliary native-weight experiment followed its original validation-only trigger.'
table=['| Model | Validation AP | Reused benchmark AUC | Reused benchmark AP | Separate assay AUC | Separate assay AP |','|---|---:|---:|---:|---:|---:|']
for r in sorted(records,key=lambda r:r['validation']['average_precision'],reverse=True):
    n=r['name'];b=bench['models'][n];e=ext['models'].get(n) if ext else None
    table.append(f"| {n}{' **selected**' if n==selected else ''} | {r['validation']['average_precision']:.4f} | {b['auroc']:.3f} | {b['average_precision']:.4f} | {format(e['auroc'],'.3f') if e else 'pending'} | {format(e['average_precision'],'.4f') if e else 'pending'} |")
mapping=json.loads((OUT/'external_mapping_audit.json').read_text())
text=f'''# AlphaGenome-assisted Sleep–Frailty research models, v2

The selected model is **{selected}**, chosen by the largest average precision (AP) on Salomon chromosomes 17–18. All {len(records)} tested model configurations are retained. A larger model or an AlphaGenome input is not assumed to be better. The separate assay and reused benchmark never enter the selection rule. The finite multi-model search used only 78 positive validation calls; validation optimism is an important limitation.

{extension_notice}

## What AlphaGenome contributes

The training input is the downloadable AVI SNV raw/PHRED score, not the detailed API predictions. The official [downloads page](https://deepmind.google.com/science/alphagenome/downloads) lists AVI SNV scores as a permissive-use artifact. The [service terms](https://deepmind.google.com/science/alphagenome/terms), inspected on 6 November 2026, explicitly except those artifacts from the similar-model-training restriction. This does not authorize training from the other API outputs or AVI feature breakdowns. Exact-allele indexed retrieval from the [Genomic Resource Repository](https://grr.iossifovlab.com/hg38/scores/AVI/index.html) found {avi['exact_matches']:,} scores. Missing scores remain unavailable, not zero biological effects. The mirror index checksum was verified; the whole 88 GB table was not independently hashed.

AVI combines predicted functional information and other annotations; it is neither a measured neuronal effect nor a sleep/frailty label. Its model-training and benchmark exposure to the evaluated assays has not been established. Adding AVI can transfer information from a larger model, but does not create new experimental validation.

The paired search compares identical fitting settings with and without AVI: regularized logistic models, tree boosting and small neural heads on frozen pretrained DNA representations. Neural heads were optimized by backpropagation across three seeds. A separate native-weight student was conditionally tested if the AVI-assisted candidates improved validation. Student status: **{student['status']}**. Its primary target is the measured reporter call, and AVI is an auxiliary target without being an input to its assay prediction head. Learning an AVI target is imitation; only measured-assay evaluation can assess experimental usefulness. See `student_config.json` and `student_status.json` when trained.

## Data and evaluation

Primary fitting uses the v1 audited 37,349 Salomon human excitatory-neuron 270-base DNA pairs: 31,333 on training chromosomes 1–16, 2,474 validation pairs on 17–18, and 3,542 pairs on chromosomes 19–22. The last group was already inspected in v1. It is a **reused comparison benchmark**, not a fresh independent test. V1 results and weights are preserved. Only paired DNA and permitted AVI enter new model inputs; gene labels, disease GWAS statistics, detailed AlphaGenome outputs and outcome-derived source annotations do not.

For a separate-study check we reconstructed Retallick-Townsley et al.'s baseline neuronal reporter assay ([paper](https://doi.org/10.1038/s41467-026-77114-x)). All 962 source DNA constructs matched hg38 exactly, including one observed allele. Removing variants within 1 kb of any Salomon variant retained {mapping['retained_after_1kb_region_exclusion']} pairs with {mapping['detected_calls_retained']} detected calls. Public native 270-base contexts were used for consistency with model input length, after confirming the measured 199-base constructs inside them. The extra native flanks were not tested by that experiment, and biological state/assay selection differ. Region exclusion does not prove no sequence homology or foundation-model exposure. This study was previously inspected elsewhere in the project; it is external to current fitting, not an entirely unseen publication. No external labels train or select these models.

## Results

{chr(10).join(table)}

AUC 0.5 denotes chance ranking, not 50% clinical accuracy. AP depends on positive prevalence; compare models within the same assay rather than comparing AP values across different datasets. Binary calls reflect statistical detection versus nonsignificance, not proven biological active/inactive truth. Class weighting and source selection mean raw neural outputs are uncalibrated. Direct AVI is ranked through a monotonic transform solely to keep score ranges finite; its Brier statistic is not a probability-calibration result.

`conditional_uncertainty.json` reports paired chromosome-resampling percentile ranges for the selected model against the DNA-only model chosen by validation and against direct AVI. These are conditional descriptive sensitivities, not complete uncertainty guarantees. Only four reused benchmark chromosomes exist. Training, search/selection, study sampling and foundation-model exposure uncertainty are not included. Point improvements alone are insufficient to establish that adding AlphaGenome reliably helps.

The selected model's reused benchmark AUC is {bench['models'][selected]['auroc']:.3f}, AP {bench['models'][selected]['average_precision']:.4f}. {'Its separate-study AUC is '+format(ext['models'][selected]['auroc'],'.3f')+' and AP '+format(ext['models'][selected]['average_precision'],'.4f')+'.' if ext else 'Separate-study evaluation remains pending.'} See all model results, including failures to improve. No recommendation was changed to whichever model won the separate-study comparison.

## Verification, scope and remaining work

Frozen encoder caches were replayed on sampled rows against the pinned original weights with exact agreement; source/checkpoint/feature hashes are retained. All selected heads and comparators are reloaded and checked against recorded validation predictions in `checkpoint_replay.json`. Independent R ranking calculations are compared with the Python metrics in `independent_metric_agreement.json`. This research code uses the previously checked custom NumPy/MLX HyenaDNA port; official PyTorch-executable parity remains untested.

The 33 project candidate predictions are exploratory model outputs without matching experimental candidate labels. They do not establish disease association, causality, a causal gene, cell specificity across tissues, genetic correlation, or validation of 15-/22-variant scores. None of these models predicts an individual's insomnia, frailty or circadian disorder. Stronger claims need additional independent assay data, appropriate uncertainty and cell-state testing, and separate genotype/outcome evidence for clinical prediction.

Source sequence-level extracts and external predictions remain local. The Retallick-Townsley article is CC BY-NC-ND 4.0; do not redistribute adapted sequence-level material. Salomon source/data and HyenaDNA licenses remain as documented in the v1 report. Preserve upstream notices. No model package has been published externally.
'''
(OUT/'README.md').write_text(text)
with (OUT/'model_comparison.csv').open('w',newline='') as f:
    fields=['model','selected','validation_ap','reused_benchmark_auc','reused_benchmark_ap','external_auc','external_ap'];writer=csv.DictWriter(f,fieldnames=fields);writer.writeheader()
    for r in records:
        n=r['name'];b=bench['models'][n];e=ext['models'][n] if ext else {}
        writer.writerow(dict(zip(fields,[n,n==selected,r['validation']['average_precision'],b['auroc'],b['average_precision'],e.get('auroc',''),e.get('average_precision','')])))
print('Wrote v2 model report',flush=True)
