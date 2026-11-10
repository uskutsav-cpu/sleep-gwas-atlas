"""Reconcile real allele identities and source assay QC without model predictions."""
import csv,json,hashlib,collections,math
from pathlib import Path

ROOT=Path(__file__).resolve().parents[2]
OUT=ROOT/'outputs/sleep_frailty_ai_v9_validation/next_source_screen'


def rows(path):
    with Path(path).open() as f:return list(csv.DictReader(f,delimiter='\t' if str(path).endswith('.tsv') else ','))


def main():
    for receipt_name in ['targeted_archive_retrieval_receipt.json','author_counts_and_code_retrieval_receipt.json','author_comparator_and_filter_code_receipt.json']:
        for f in json.loads((OUT/receipt_name).read_text())['files']:
            assert hashlib.sha256(Path(f['path']).read_bytes()).hexdigest()==f['sha256']
    accepted=rows(OUT/'Tardaguila_reserved_no_historical_overlap.csv');design={ (r['SNP'],r['comparison'],r['window']):r for r in rows(OUT/'Tardaguila2026_sequence_information.tsv')}
    source=rows(OUT/'Tardaguila2026_TableS3.tsv');lookup={ (r['SNP'],r['comparison'],r['window'],r['Cell_Type']):r for r in source};assert len(lookup)==len(source)
    attrs={}
    for comparison in ['ref_vs_alt','ref_ref_vs_alt_ref','ref_ref_vs_ref_alt']:
        stem='single_variants' if comparison=='ref_vs_alt' else 'diplotypes_1_'+comparison
        rr=rows(OUT/f'author_attributes_{stem}.tsv');attrs[comparison]={r['ID']:r for r in rr};assert len(attrs[comparison])==len(rr)
    model=(OUT/'324_MPRA_bc_synthesis_MPRA_model_adapted_v3.R').read_text()
    for text in ['levels = c("ref","alt")','~material + allele + material:allele','name="materialRNA.allelealt"','paired=F']:
        assert text in model,text
    assert 'out2$Log2Skew <- out2$B_log2FC - out2$A_log2FC' in model
    per_cell={};qc=[];identity_count=0
    for cell in ['K562','HL60','CHRF','THP1']:
        counts={r['Oligo']:r for r in rows(OUT/f'{cell}_independent_oligo_count_sums.tsv')}
        with (OUT/f'{cell}_author_sample_metadata.tsv').open() as f:metadata=list(csv.reader(f,delimiter='\t'))
        names={r[0] for r in metadata};assert len(names)==len(metadata)==14
        rna=[r[0] for r in metadata if r[1]=='RNA'];dna=[r[0] for r in metadata if r[1]=='DNA'];assert len(rna)==len(dna)==7
        rep_names={x.removesuffix('_cDNA') for x in rna};assert rep_names=={x.removesuffix('_gDNA') for x in dna}
        min_barcodes=[];min_good_dna=[];min_good_rna=[];nonpositive_se=0
        for r in accepted:
            key=(r['source_SNP'],r['comparison'],r['window']);physical=design[key];table=lookup[(*key,cell)];assert table['ID']==physical['ID']==physical['Oligo_ID_REF']
            a,b=physical['Oligo_ID_REF'],physical['Oligo_ID_ALT'];aa,bb=attrs[r['comparison']][a],attrs[r['comparison']][b]
            assert aa['Allele']=='ref' and bb['Allele']=='alt';assert aa['SNP']==bb['SNP']==r['source_SNP'];assert aa['Window']==bb['Window']==r['window'];assert aa['strand']==bb['strand']=='+'
            assert a in counts and b in counts;identity_count+=1
            # Source Table S3's original GLM allelic coefficient is the target;
            # no activity-difference or simple RNA/DNA ratio replaces it.
            effect=float(table['log2Skew']);se=float(table['Skew_SE']);assert math.isfinite(effect) and math.isfinite(se) and se>=0
            nonpositive_se+=int(se<=0)
            nb=[int(counts[z]['distinct_barcodes']) for z in [a,b]];assert min(nb)>0
            raw={}
            for allele,ident in enumerate([a,b]):
                assert names.issubset(counts[ident])
                for name in names:
                    value=float(counts[ident][name])-nb[allele]
                    assert value>=0 and abs(value-round(value))<1e-8
                    raw[(allele,name)]=value
            dna_reps=sum(all(raw[(allele,rep+'_gDNA')]>=5 for allele in [0,1]) for rep in rep_names)
            rna_reps=sum(all(raw[(allele,rep+'_cDNA')]>=5 for allele in [0,1]) for rep in rep_names)
            qc.append({'source_SNP':r['source_SNP'],'comparison':r['comparison'],'window':r['window'],'variant_GRCh38':r['variant_GRCh38'],'Cell_Type':cell,
                       'source_replicates':7,'minimum_distinct_barcodes_across_alleles':min(nb),'replicates_both_alleles_DNA_total_ge5':dna_reps,'replicates_both_alleles_RNA_total_ge5':rna_reps,
                       'source_GLM_standard_error_positive':se>0,'excluded_by_this_QC_profile':False})
            min_barcodes.append(min(nb));min_good_dna.append(dna_reps);min_good_rna.append(rna_reps)
        per_cell[cell]={'accepted_tile_comparisons':len(accepted),'matched_RNA_DNA_replicates':7,'minimum_barcodes_over_accepted_allele_pairs':min(min_barcodes),
                        'pairs_with_at_least_15_barcodes_per_allele':sum(n>=15 for n in min_barcodes),'pairs_with_all_7_DNA_replicates_ge5_raw_counts_per_allele':sum(n==7 for n in min_good_dna),
                        'pairs_with_all_7_RNA_replicates_ge5_raw_counts_per_allele':sum(n==7 for n in min_good_rna),'nonpositive_original_GLM_SE':nonpositive_se}
    with (OUT/'candidate_source_assay_QC_profile.csv').open('w') as f:
        w=csv.DictWriter(f,fieldnames=list(qc[0]));w.writeheader();w.writerows(qc)
    audit={'passed':True,'reference_and_alternate_source_attributes_reconciled':identity_count,'full_148bp_reference_verified_designs':len(accepted),'eligible_original_GLM_tile_cell_measurements':len(qc),
           'signed_target_convention':'Original GLM RNA-by-alt interaction with ref and DNA reference levels: positive indicates greater alternate-versus-reference RNA/DNA activity. Exact original comparator allele IDs and physical single-change pairs verified.',
           'primary_magnitude_target':'Absolute original Table S3 log2Skew; not log2FC enhancer activity and not a multi-tile meta-analysis assigned to one construct.',
           'by_cell':per_cell,'source_filtering':'Author code selects oligos that pass both RNA and DNA total-count >=5 in at least one selected replicate, then retains their rows across replicates; adds pseudocount 1 and fills missing barcode/sample entries with 1. This does not establish adequate counts in every replicate.',
           'source_GLM':'DESeq2, unpaired material+allele+material:allele; actual pairing exists in sample metadata but original coefficient/SE were not replaced by a paired reanalysis.',
           'original_barcode_to_counts_pipeline_reconstructed':False,'original_GLM_coefficients_and_SE_independently_refitted':False,
           'independent_R_filtered_count_integrity_and_oligo_aggregation':True,'new_significance_or_activity_filter_applied':False,'additional_QC_profile_used_to_select_favorable_effects':False,
           'model_predictions_generated_or_compared':False,'future_pretraining_overlap_audit_and_frozen_test_protocol_required':True}
    (OUT/'original_comparator_and_assay_QC_audit.json').write_text(json.dumps(audit,indent=2)+'\n');print(json.dumps(audit),flush=True)


if __name__=='__main__':main()
