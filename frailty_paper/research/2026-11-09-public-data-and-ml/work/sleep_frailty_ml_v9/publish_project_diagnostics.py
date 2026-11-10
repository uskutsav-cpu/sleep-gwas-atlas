"""Preserve project evidence; do not turn raw projections into authorized ranks."""
from final_validation_common import *


def main():
    result=json.loads((FINAL/'fresh_outcome_evaluation.json').read_text());rows=list(csv_rows(V8/'blind_inference_row_metadata.csv'));raw={}
    for arm in ['sequence','AVI']:
        with np.load(FINAL/f'{arm}_blinded_predictions.npz') as z:p={k:z[k] for k in z.files}
        raw[arm]={int(index):{k:v[j] for k,v in p.items() if k!='row_index'} for j,index in enumerate(p['row_index'])}
    decisions=[];diagnostics=[]
    for index,r in enumerate(rows):
        if r['role']!='project':continue
        reason=r['unscorable_reason'] or 'Unmeasured native genomic context: applicability of assay-active-only external validation not established.'
        decision={**r,'authorized':False,'authorized_prediction':'','authorized_rank':'','reason':reason}
        decisions.append(decision)
        if index in raw['sequence']:
            diagnostic={**r,'interpretation':'UNVALIDATED molecular projections; no authorized project ranking, disease prediction or target-gene claim.'}
            for arm in ['sequence','AVI']:
                p=raw[arm][index]
                for key in ['magnitude','signed_effect','interval_halfwidth','estimated_error','estimated_baseline_error','maximum_standardized_error_feature']:
                    diagnostic[f'{arm}_{key}']=float(p[key])
            diagnostic['AVI_retained']=bool(result['AVI_comparison']['retain_AVI']);diagnostics.append(diagnostic)
    assert len(decisions)==2490 and len(diagnostics)==2350
    assert all(not r['authorized'] and r['authorized_prediction']==r['authorized_rank']=='' for r in decisions)
    write_csv(FINAL/'project_answer_availability.csv',decisions);write_csv(FINAL/'project_raw_diagnostics_NOT_VALIDATED.csv',diagnostics)
    write_json(FINAL/'project_output_audit.json',{'total_contexts':len(decisions),'supported_SNV_contexts':len(diagnostics),'unsupported_indel_contexts':140,
               'authorized_predictions':0,'authorized_ranks':0,'original_gene_evidence_modified':False,'reason':'Project DNA has no measured assay-activity proof. One-context active-element validation cannot establish its prediction accuracy. Raw diagnostics are separate.'})
    print('PROJECT ANSWERS WITHHELD; RAW DIAGNOSTICS SEPARATE',flush=True)


if __name__=='__main__':main()
