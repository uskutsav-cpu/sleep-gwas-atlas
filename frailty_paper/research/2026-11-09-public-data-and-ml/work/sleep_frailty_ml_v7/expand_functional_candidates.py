"""Score the complete existing molecular candidate set with honest nonanswers."""
import gzip,subprocess,sys,bisect,collections
from common import *
def write_csv(path,rows):
    with Path(path).open('w') as f:w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
def main():
    audit=json.loads((OUT/'all_project_variant_input_audit.json').read_text());assert sha(DATA/'all_project_variant_inputs.csv')==audit['input_sha256'];metadata=json.loads((DATA/'all_project_variant_metadata.json').read_text());byid={r['candidate_id']:r for r in metadata};registry=json.loads((OUT/'inference_registry.json').read_text());inputs=list(csv.DictReader((DATA/'all_project_variant_inputs.csv').open()))
    if registry['recommended_arm']=='AVI_added':
        receipt=json.loads((OUT/'extra_AVI_inference_receipt.json').read_text());assert sha(receipt['output'])==receipt['sha256'];avi={r['variant_GRCh38']:r for r in csv.DictReader(open(receipt['output']))};path=DATA/'AVI_inference_lookup.csv.gz'
        with gzip.open(path,'rt') as f:lookups=list(csv.DictReader(f))
        seen={(r['variant_GRCh38'],r['reference51'],int(r['changed_position51'])) for r in lookups}
        for r in inputs:
            m=byid[r['id'].split('__')[0]]
            if m['unscorable_reason']:continue
            key=r['variant_GRCh38'];value=avi[key];identity=(key,r['reference'][50:101],25)
            if identity in seen or value['available']!='1':continue
            seen.add(identity);lookups.append({'variant_GRCh38':key,'reference51':identity[1],'changed_position51':25,'raw':value['raw'],'phred':value['phred']})
        with gzip.open(path,'wt') as f:w=csv.DictWriter(f,fieldnames=list(lookups[0]));w.writeheader();w.writerows(lookups)
    out=OUT/'all_project_variant_predictions.csv';subprocess.run([sys.executable,str(WORK/'score.py'),'--input',str(DATA/'all_project_variant_inputs.csv'),'--output',str(out)],env=environment(),check=True);predictions=list(csv.DictReader(out.open()));assert len(predictions)==audit['rows'];positions=collections.defaultdict(set)
    with np.load(DATA/'training_arrays.npz') as z:fitbank=set(z['bank_index'][z['split']<2].tolist())
    with gzip.open(DATA/'sequence_bank.csv.gz','rt') as f:
        for r in csv.DictReader(f):
            if int(r['bank_id']) in fitbank:positions[r['chromosome']].add(int(r['position_hg38']))
    positions={c:sorted(v) for c,v in positions.items()};byvariant={}
    for row in predictions:
        m=byid[row['id'].split('__')[0]];c,p,_,_=m['variant_GRCh38'].split(':');p=int(p);ps=positions.get(c,[]);j=bisect.bisect_left(ps,p-2000);overlap=j<len(ps) and ps[j]<=p+2000
        if m['unscorable_reason']:row['unranked_reason']=m['unscorable_reason'];assert row['raw_model_projection']=='' and row['predicted_absolute_log2_effect']==''
        row.update({'variant_GRCh38':m['variant_GRCh38'],'input_labels':';'.join(m['input_labels']),'variant_class':m['variant_class'],'in_existing_488_molecular_candidates':int(m['in_existing_molecular_candidate_set']),'in_original_33_score_markers':int(m['in_original_33_score_markers']),'within_2kb_of_final_training_DNA':int(overlap),'gene_context_associations':';'.join(m['gene_context_associations']),'interpretation':'Genomic reporter-magnitude hypothesis only. Raw projections are unvalidated when authorized fields are blank; no signed regulation, causal gene or clinical genetic-score validation.'});byvariant.setdefault(m['variant_GRCh38'],{})[row['cell_group']]=row
    write_csv(out,predictions);functional=[r for r in predictions if r['in_existing_488_molecular_candidates']==1];assert len(functional)==488*5;write_csv(OUT/'functional_candidate_predictions.csv',functional)
    # Confirm that the original33marker scores reproduce in the independently
    # verified151bpgenomic inputs; the original270bpinputs remain unchanged.
    old=list(csv.DictReader((OUT/'project_molecular_predictions.csv').open()));maximum=0.
    for r in old:
        other=byvariant[r['variant_GRCh38']][r['cell_group']];assert r['status']==other['status'];maximum=max(maximum,abs(float(r['raw_model_projection'])-float(other['raw_model_projection'])))
    assert maximum<1e-5
    source=ROOT/'outputs/stage5_variant_gene_tissue_phenotype_integrated.csv';atlas=list(csv.DictReader(source.open()));credible=collections.defaultdict(set)
    for row in atlas:
        key=row['variant_GRCh38_ref_alt'].replace('>',':');rows=byvariant[key];row['V7_raw_model_projections_unvalidated']=';'.join(cell+'='+p['raw_model_projection'] for cell,p in rows.items() if p['raw_model_projection']);row['V7_authorized_research_cells']=';'.join(cell for cell,p in rows.items() if p['status']=='research_candidate');row['V7_unanswered_reasons']=';'.join(cell+':'+p['unranked_reason'] for cell,p in rows.items() if p['status']=='unranked');row['V7_interpretation']='Prediction annotation only. QTL/nearby/expressed genes and raw projections do not establish a disease-causal gene.'
        if row['credible_set_id']:credible[key].add(tuple(row.get(k,'') for k in ['qtl_source_collection','qtl_dataset','gene_id','tissue_or_cell','molecular_measure','credible_set_id']))
    assert len(atlas)==1351;write_csv(OUT/'complete_gene_evidence_with_model_annotations.csv',atlas)
    # Same credible-set source, same modeled cell, <=25bp. Non-overlapping
    # nominal marginal ranges are conservative research diagnostics only.
    keys=sorted(credible);comparisons=[]
    for i,a in enumerate(keys):
        ca,pa,_,_=a.split(':');pa=int(pa)
        for b in keys[i+1:]:
            cb,pb,_,_=b.split(':');distance=abs(pa-int(pb));shared=credible[a]&credible[b]
            if ca!=cb or distance>25 or not shared:continue
            for cell in audit['cell_contexts']:
                aa,bb=byvariant[a][cell],byvariant[b][cell];winner='';reason='individual_prediction_withheld';difference=''
                if aa['raw_model_projection'] and bb['raw_model_projection']:difference=float(aa['raw_model_projection'])-float(bb['raw_model_projection'])
                if aa['status']==bb['status']=='research_candidate':
                    if float(aa['lower90'])>float(bb['upper90']):winner=a;reason='disjoint_nominal_ranges'
                    elif float(bb['lower90'])>float(aa['upper90']):winner=b;reason='disjoint_nominal_ranges'
                    else:reason='uncertainty_ranges_overlap'
                comparisons.append({'variant_A':a,'variant_B':b,'distance_bp':distance,'shared_molecular_credible_set_sources':len(shared),'cell_group':cell,'raw_projection_difference_unvalidated':difference,'research_order_winner':winner,'status':'exploratory_order_supported_by_disjoint_nominal_ranges' if winner else 'unanswered','reason':reason,'scope':'No joint90%coverage guarantee, signed regulation, disease causation or causal target-gene conclusion.'})
    if comparisons:write_csv(OUT/'nearby_credible_variant_comparisons.csv',comparisons)
    report={'union_variants':len(metadata),'existing_molecular_candidates':488,'score_markers':33,'variant_cell_rows':len(predictions),'functional_candidate_cell_rows':len(functional),'supported_SNV_variants':sum(not r['unscorable_reason'] for r in metadata),'unsupported_variants_retained':sum(bool(r['unscorable_reason']) for r in metadata),'authorized_research_answers':sum(r['status']=='research_candidate' for r in predictions),'gene_evidence_rows_preserved':len(atlas),'nearby_same_credible_set_comparisons':len(comparisons),'nearby_research_orders_with_nonoverlapping_ranges':sum(bool(r['research_order_winner']) for r in comparisons),'original33marker_raw_projection_replay_max_error':maximum,'original_gene_evidence_unchanged_sha256':sha(source),'inputs_audit_sha256':sha(OUT/'all_project_variant_input_audit.json'),'outputs_sha256':{n:sha(OUT/n) for n in ['all_project_variant_predictions.csv','functional_candidate_predictions.csv','complete_gene_evidence_with_model_annotations.csv']},'clinical_or_genetic_score_validation':False};write_json(OUT/'complete_project_model_integration_audit.json',report);print('COMPLETE MOLECULAR CANDIDATE SET INTEGRATED',json.dumps(report),flush=True)
if __name__=='__main__':main()
