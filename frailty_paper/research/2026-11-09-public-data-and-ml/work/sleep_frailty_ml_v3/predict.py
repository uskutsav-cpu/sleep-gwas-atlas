"""Score paired experimental DNA with a completed exploratory model.

Input JSON list: ref_sequence,alt_sequence,chromosome,position_hg38,assay_context
(Salomon_NGN2 or Siraj_OL42_SKNSH), optional id/genomic_strand. GRCh38 coordinates
must describe the substituted base; the caller supplies and verifies coordinates.
No lab outcome, sleep/frailty measurement or actual reporter activity is accepted.
"""
import argparse,json,numpy as np
from common import *
from raw_features import RawDNAInference
from validation_registry import completed_candidates
def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--input',type=Path,required=True);p.add_argument('--output',type=Path,required=True);p.add_argument('--model',required=True);a=p.parse_args();registry={};status='Exploratory trained checkpoint; robust family selection and blind validation pending'
    if (OUT/'final_model_registry.json').is_file():
        final=json.loads((OUT/'final_model_registry.json').read_text());assert final['plan_sha256']==sha(OUT/'final_training_locked_plan.json')
        for name,r in final['records'].items():
            if r['artifact']:assert sha(OUT/r['artifact'])==r['checkpoint_sha256']
            assert name not in registry;registry[name]=r
        if a.model in final['records']:status='Final development-selected/refitted research model; consult held-out evaluation report for performance'
    if a.model not in registry:
        development,_=completed_candidates();registry.update(development)
    assert a.model in registry,'Specify a completed model name from completed_validation_registry.json or final_model_registry.json';records=json.loads(a.input.read_text());engine=RawDNAInference(records);pred=engine.predict(registry[a.model],np.arange(len(records)),registry);scores=[]
    for i,r in enumerate(records):
        row={'id':r.get('id',str(i)),'variant_GRCh38':engine.manifest[i]['variant_GRCh38'],'assay_context':r['assay_context'],'assay_call_score_uncalibrated':float(pred['overall'][i])}
        if 'conditional_scores' in pred:row['active_region_expert_score_uncalibrated']=float(pred['conditional_scores'][i,1])
        if 'predicted_active' in pred:row['predicted_activity_gate']=float(pred['predicted_active'][i])
        if 'absolute_effect' in pred:row['effect_head_scaled_output']=float(pred['absolute_effect'][i])
        scores.append(row)
    output={'model':a.model,'model_status':status,'input_sha256':sha(a.input),'coordinate_verification':'Caller-supplied GRCh38 coordinates; allele identity checked against paired sequences and strand. Native flanks are not independently checked by this CLI.','scope':'Laboratory reporter detection ranking in two study/cell contexts. Not patient insomnia/frailty prediction, causal evidence, gene assignment or genetic-score validation.','AVI_usage':'Only exact matching rows of the audited designated permissive downloadable AVI table; unavailable input scores are imputed using the model training transform. No API or detailed-output derivatives.','scores':scores};a.output.write_text(json.dumps(output,indent=2,allow_nan=False)+'\n');print('Wrote research predictions',a.output)
if __name__=='__main__':main()
