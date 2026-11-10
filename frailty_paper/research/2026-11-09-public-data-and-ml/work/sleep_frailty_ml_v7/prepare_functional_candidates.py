"""Prepare every existing molecular candidate plus the33score markers.

These are computational genomic windows, not invented experimental oligos.
Indels remain explicitly unsupported by the single-substitution model.
"""
import collections
from common import *
from reference import get_reference
CONTEXTS={'induced_excitatory_neuron':'WTC11_NGN2','neuronal_tumor':'SK-N-SH','neural_progenitor':'H1_N2_NPC','kidney_derived':'HEK293T','hepatocyte_tumor':'HepG2'}
def main():
    source=ROOT/'outputs/stage5_variant_gene_tissue_phenotype_integrated.csv';atlas=list(csv.DictReader(source.open()));by=collections.defaultdict(list)
    for r in atlas:by[r['variant_GRCh38_ref_alt'].replace('>',':')].append(r)
    original=json.loads((ROOT/'outputs/sleep_frailty_ai_v3/project_raw_inputs.json').read_text());markers=collections.defaultdict(set)
    for r in original:markers[r['variant_GRCh38']].add(r['rsid'])
    keys=sorted(set(by)|set(markers));genome=get_reference();metadata=[];inputs=[];classes=collections.Counter();statuses=collections.Counter()
    for number,key in enumerate(keys):
        c,position,ref,alt=key.split(':');position=int(position);kind='SNV' if len(ref)==len(alt)==1 else 'indel' if len(ref)!=len(alt) else 'multiple_substitution';classes[kind]+=1;reason='';verified=False;reference=alternate='';observed=''
        try:observed=genome.sequence(c,position-1,position-1+len(ref));verified=observed==ref
        except (KeyError,AssertionError):reason='unavailable_reference_coordinate'
        if not reason and not verified:reason='source_reference_allele_does_not_match_official_hg38'
        if not reason and kind!='SNV':reason='unsupported_'+kind+'_by_single_substitution_model'
        if not reason:
            try:reference=genome.sequence(c,position-76,position+75)
            except (KeyError,AssertionError):reason='insufficient_supported_flanking_sequence'
            if reference and set(reference)-set('ACGT'):reason='ambiguous_flanking_reference_bases'
            if not reason:assert reference[75]==ref;alternate=reference[:75]+alt+reference[76:];assert sum(a!=b for a,b in zip(reference,alternate))==1
        ident='candidate_'+hashlib.sha256(key.encode()).hexdigest()[:12];rows=by.get(key,[]);labels=sorted(set(markers.get(key,set()))|{v for r in rows for v in r.get('input_labels','').split('|') if v});genes=sorted({r['gene_symbol'] or r['gene_id'] for r in rows if r['gene_symbol'] or r['gene_id']});sets=sorted({r['credible_set_id'] for r in rows if r['credible_set_id']});m={'candidate_id':ident,'variant_GRCh38':key,'input_labels':labels,'gene_evidence_rows':len(rows),'gene_context_associations':genes,'credible_set_ids':sets,'in_existing_molecular_candidate_set':key in by,'in_original_33_score_markers':key in markers,'variant_class':kind,'source_reference_verified':verified,'observed_hg38_reference':observed,'unscorable_reason':reason,'changed_position':75 if not reason else None};metadata.append(m);statuses[reason or 'supported_genomic_SNV_input']+=1
        for cell,line in CONTEXTS.items():inputs.append({'id':ident+'__'+cell,'reference':reference if not reason else '','alternate':alternate if not reason else '','cell_group':cell,'cell_line':line,'variant_GRCh38':key})
    path=DATA/'all_project_variant_inputs.csv'
    with path.open('w') as f:w=csv.DictWriter(f,fieldnames=list(inputs[0]));w.writeheader();w.writerows(inputs)
    write_json(DATA/'all_project_variant_metadata.json',metadata);report={'existing_molecular_candidates':len(by),'score_markers':len(markers),'overlap_between_sets':len(set(by)&set(markers)),'union_variants':len(keys),'cell_contexts':CONTEXTS,'rows':len(inputs),'variant_classes':dict(classes),'input_status':dict(statuses),'source_gene_evidence_sha256':sha(source),'source_marker_raw_inputs_sha256':sha(ROOT/'outputs/sleep_frailty_ai_v3/project_raw_inputs.json'),'input_sha256':sha(path),'metadata_sha256':sha(DATA/'all_project_variant_metadata.json'),'reference_receipt_sha256':sha(OUT/'reference_genome_receipt.json'),'interpretation':'151bp genomic sequence hypotheses verified against official hg38. Not experimentally measured constructs, a causal gene map, GWAS allele-convention validation or independent genetic-score validation. Indels/missing inputs are retained as unanswered cases. No model fitting or target-label input.'};write_json(OUT/'all_project_variant_input_audit.json',report);print('COMPLETE PROJECT CANDIDATE INPUTS PREPARED',json.dumps(report),flush=True)
if __name__=='__main__':main()
