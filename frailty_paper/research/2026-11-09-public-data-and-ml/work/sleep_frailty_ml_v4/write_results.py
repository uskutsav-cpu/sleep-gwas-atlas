"""Write the bounded additional-data experiment's evidence and limitations."""
from common import *

def main():
    status=json.loads((OUT/'teacher_training_status.json').read_text())
    probe=json.loads((OUT/'assay_transfer_results.json').read_text())
    audit=json.loads((OUT/'assay_probe_replay_audit.json').read_text())
    replay=json.loads((OUT/'teacher_checkpoint_replay.json').read_text())
    bounds=json.loads((OUT/'assay_probe_block_uncertainty.json').read_text())
    assert audit['passed'] and replay['passed'] and audit['old_V3_registry_unchanged']
    ensembles={r['variant']:r['metrics'] for r in probe['results'] if r['seed']=='ensemble'}
    rows=[]
    for s,study in enumerate(['Salomon','Siraj']):
        for subset,key in [('all measured pairs','studies'),('pairs in active sequences','within_active')]:
            a=ensembles['original_control'][key][str(s)];b=ensembles['AVI_taught'][key][str(s)]
            ci=bounds[study.lower()+('_all' if key=='studies' else '_active')]['paired_differences']['AVI_taught']['original_control']
            rows.append({'study':study,'subset':subset,'n':a['n'],'calls':a['calls'],'original_AUROC':a['AUROC'],'taught_AUROC':b['AUROC'],'original_AP':a['AP'],'taught_AP':b['AP'],'AP_difference':b['AP']-a['AP'],'paired_AP_difference_CI95':ci['paired_AP_difference_CI95']})
    positive=sum(r['AP_difference']>0 for r in rows)
    decision='Retain the previously evaluated V3 model. This V4 experiment is a development candidate; no independent biological validation or clinical prediction is established.'
    if positive<4:decision+=' The measured development comparisons are mixed; more AlphaGenome imitation has not consistently strengthened the model.'
    else:decision+=' Development comparisons favor the added training, but a new independent measured study remains necessary before promotion.'
    result={'status':'Completed bounded additional-AVI development comparison; native run interrupted after one completed epoch',
        'training_rows':status['training_rows'],'new_training_pairs':25203,'new_teacher_validation_pairs':5580,
        'new_teacher_validation_regions':60,'teacher_checkpoint_sha256':status['checkpoint_sha256'],
        'teacher_metrics':status['selected_validation'],'measured_reused_development_comparisons':rows,
        'decision':decision,'V3_registry_sha256':sha(V3/'final_model_registry.json'),
        'not_established':['Independent biological replication','Patient insomnia or frailty prediction','Genetic-score validation','Causation'],
        'all_available_outputs_used':False,
        'source_scope':'Only designated permissive AVI SNV raw/PHRED predictions; checking examples withheld, detailed API/Atlas outputs excluded.',
        'comparison_scope':'Whole AVI teaching procedure versus an untaught matched control. Does not isolate the added 25,203 examples from other teaching changes, and is not a head-to-head fresh test against frozen V3.',
        'next_step':'If expanding, make a separately versioned larger region sample and compare on new measured laboratory data. Stop additional epochs when checking performance worsens; do not treat prior V3 tests as fresh.'}
    write_json(OUT/'additional_training_assessment.json',result)
    lines=[
        'ADDITIONAL ALPHAGENOME TRAINING — RESULTS',
        '',
        'We trained a separate small DNA model on 110,070 examples. This includes 25,203 newly collected training examples. We kept 5,580 additional examples, from 60 other DNA regions, for checking and choosing the saved model.',
        'These new labels are AlphaGenome predictions, not newly performed experiments. Nearby examples are related; they are not 30,783 independent biological observations.',
        '',
        'WHAT ACTUALLY RAN',
        'The prediction head completed 14 passes, with its best saved result at pass 10. A later pass updated 217,152 values inside the DNA encoder and was saved. The run was interrupted during a second encoder pass; that incomplete pass was not used. No claim is made that the planned two passes finished.',
        'The saved encoder passed a separate NumPy-versus-MLX arithmetic check. Six matched laboratory readout models, across three fixed random seeds per comparison, were then fitted and checked.',
        '',
        'HOW WELL IT COPIES ALPHAGENOME',
        f"For separating high AlphaGenome scores from other scores, the ranking measure was {status['selected_validation']['high_AVI_AUROC']:.3f} (0.5 is chance; 1 is perfect). This is a teacher-copying result, not a disease prediction accuracy.",
        'Distinguishing neighboring DNA changes remained weak: average local score correlation was 0.103. Numerical score errors were worse than a constant based on the checking-set mean (negative R-squared). Updating the encoder slightly improved ranking but worsened numerical score matching.',
        '',
        'CHECK AGAINST REAL LABORATORY RESULTS',
        'Both comparison models received the same permitted AlphaGenome inputs, 579 input values, training examples, three random seeds, and fitting settings. The adapted version also used the newly taught encoder and its learned representation. Real laboratory outcomes were targets only, never prediction inputs.',
        'This checks the whole added teaching procedure. It does not isolate the value of the 25,203 new examples by themselves, and it is not a fresh test against the frozen V3 model.',
        'These laboratory examples have been used for development before. They were also used to choose when fitting stopped. The results below are development comparisons, not new independent tests.',
    ]
    for r in rows:
        lo,hi=r['paired_AP_difference_CI95']
        lines.append(f"{r['study']}, {r['subset']}: {r['n']:,} pairs, {r['calls']:,} detected effects; ranking area {r['original_AUROC']:.4f} -> {r['taught_AUROC']:.4f}; precision-recall summary {r['original_AP']:.4f} -> {r['taught_AP']:.4f}; paired precision-recall change {r['AP_difference']:+.4f} (development block interval {lo:+.4f} to {hi:+.4f}).")
    lines.extend([
        'Higher ranking and precision-recall summaries are better. These numbers are not the percentage of patients correctly diagnosed. Intervals resample 1 Mb DNA blocks; they do not remove bias from using development examples.',
        'Saved laboratory ensembles were independently replayed, and 16 ranking metrics were checked in R.',
        '',
        'DECISION',decision,
        '',
        'CAN WE USE EVERYTHING AND KEEP TRAINING?',
        'We can add eligible public scores in batches. We must keep some DNA regions out of fitting to detect memorization and check against measured experiments. Repeating training after checking results get worse can reduce usefulness.',
        'Only the designated permissive downloadable AVI SNV scores were used for this training. Detailed AlphaGenome API/Atlas predictions, embeddings, splicing downloads and feature-importance outputs were excluded under the current training-use terms.',
        'The full chromosome 1–3 reserved teacher-region sample was not collected or evaluated. Those would check copying AlphaGenome, not real biological accuracy.',
        '',
        'WHAT THIS DOES NOT SHOW',
        'This does not predict who develops insomnia or frailty, validate a patient genetic score, establish genetic correlation, identify proven causal genes, or establish causation. The original scientific analyses and V3 model remain unchanged.',
        '',
        'SOURCES',
        'https://deepmind.google.com/science/alphagenome/downloads',
        'https://deepmind.google.com/science/alphagenome/terms',
        'https://github.com/google-deepmind/alphagenome',
        '',
        'AUDIT FILES',
        'teacher_data_receipt.json; teacher_training_history.json; teacher_training_status.json; teacher_checkpoint_replay.json; assay_transfer_results.json; assay_probe_replay_audit.json; assay_probe_block_uncertainty.json; additional_training_assessment.json.',
    ])
    (OUT/'PLAIN_LANGUAGE_RESULTS.txt').write_text('\n'.join(lines)+'\n')
    sources={p.name:sha(p) for p in WORK.glob('*.py')}
    write_json(OUT/'completion_receipt.json',{'assessment_sha256':sha(OUT/'additional_training_assessment.json'),'plain_report_sha256':sha(OUT/'PLAIN_LANGUAGE_RESULTS.txt'),'source_hashes':sources,'teacher_checkpoint_sha256':sha(status['checkpoint']),'teacher_replay_passed':replay['passed'],'lab_replay_and_R_checks_passed':audit['passed'],'V3_registry_unchanged':audit['old_V3_registry_unchanged'],'no_fresh_test_or_disease_claim':True})
    print(json.dumps({'decision':decision,'comparisons':rows},indent=2))
if __name__=='__main__':main()
