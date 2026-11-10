"""Lock final roles and refit recipes from completed development folds only."""
import json,math
from common import *
def main():
    assert (OUT/'development_robustness_finished.json').is_file()
    comparison=json.loads((OUT/'development_comparison.json').read_text());assert len(comparison['aggregate'])==12 and all(r['completed_folds']==3 for r in comparison['aggregate'])
    protocol=json.loads((OUT/'final_selection_protocol.json').read_text());assert not comparison['test_evaluated']
    choose=lambda key:sorted(comparison['aggregate'],key=lambda r:(-r[key],r['family']))[0]['family']
    roles={'overall_assay_call_detector':choose('all_fold_mean_all_call_AP_skill'),'within_active_variant_ranker':choose('all_fold_mean_conditional_AP_skill')}
    # Controls use exactly the existing fixed screening-selected baselines.
    controls={r['name']:'baseline_'+r['name'] for r in json.loads((OUT/'development_robustness_locked_config.json').read_text())['classical_records']};ablations={'ablation_'+f:f for f in protocol['matched_neural_ablation_refits']['families']};stacking=json.loads((OUT/'conditional_stacking_final_artifacts.json').read_text());assert stacking['development_comparison_sha256']==sha(OUT/'development_comparison.json');stacking_controls={'combination_'+f:f for f in stacking['families']};families=sorted(set(roles.values())|set(controls.values())|set(ablations.values())|set(stacking_controls.values()));recipes=[]
    for family in families:
        folds=[r for r in comparison['fold_results'] if r['family']==family];assert len(folds)==3
        if family in stacking['families']:
            info=stacking['families'][family];recipe={'family_id':family,'kind':'conditional_stack','source':info,'base_family':info['base_family'],'include_DNA_expert':info['include_DNA_expert']}
        elif family.startswith('baseline_'):
            source=next(r for r in json.loads((OUT/'development_robustness_locked_config.json').read_text())['classical_records'] if 'baseline_'+r['name']==family);recipe={'family_id':family,'source_record':source,'kind':'classical'}
        else:
            members=[]
            for fold in folds:
                path=next(p for p in OUT.glob('*_search_results.json') if any(r['name']==fold['record_name'] for r in json.loads(p.read_text())))
                records={r['name']:r for r in json.loads(path.read_text())}
                for name in fold['members']:members.append(records[name])
            seeds=[42,73,101];epochs={str(s):max(1,int(math.floor(float(np.median([r['config']['best_epoch'] for r in members if r['config']['seed']==s]))+.5))) for s in seeds};assert len(members)==9
            sample=members[0];config=json.loads((OUT/sample['config']['locked_config']).read_text());recipe={'family_id':family,'kind':'neural','network_family':sample['family'],'features':sample['features'],'task':sample['config']['task'],'AVI_as_input':config.get('AVI_as_input',False),'teacher':'full','seeds':seeds,'fixed_epochs_by_seed':epochs,'source_member_names':[r['name'] for r in members],'original_locked_config':sample['config']['locked_config']}
        recipes.append(recipe)
    plan={'roles':roles,'controls':controls,'matched_neural_ablations':ablations,'conditional_stacking_controls':stacking_controls,'recipes':recipes,'training_chromosomes':list(range(4,23)),'reserved_test_chromosomes':[1,2,3],'protocol_sha256':sha(OUT/'final_selection_protocol.json'),'development_comparison_sha256':sha(OUT/'development_comparison.json'),'conditional_stacking_specification_sha256':sha(OUT/'conditional_stacking_prespecification.json'),'data_audit_sha256':sha(OUT/'data_audit.json'),'AVI_audit_sha256':sha(OUT/'avi_retrieval_audit.json'),'source_sequence_sha256':{s:sha(WORK/(s+'_sequences_local.npz')) for s in ['salomon','siraj']},'scope':'Final research ranking models, classical controls, matched neural AVI/no-AVI ablations and conditional combination models. No reserved-test predictions made or used for selection.','blind_test_evaluated':False}
    path=OUT/'final_training_locked_plan.json'
    if path.exists():assert json.loads(path.read_text())==plan,'Do not replace a different final plan'
    else:path.write_text(json.dumps(plan,indent=2)+'\n')
    print('Locked final model roles',roles,'and',len(recipes),'refit recipes',flush=True)
if __name__=='__main__':main()
