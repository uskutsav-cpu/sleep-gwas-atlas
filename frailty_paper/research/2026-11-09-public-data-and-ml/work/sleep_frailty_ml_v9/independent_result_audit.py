"""Independent source XML extraction, metrics, weighted-block bootstrap and gates."""
from v9_common import *
import zipfile,xml.etree.ElementTree as ET,posixpath,re
from scipy.stats import rankdata

FINAL=ROOT/'outputs/sleep_frailty_ai_v9_validation'


def original_effects(path):
    ns={'s':'http://schemas.openxmlformats.org/spreadsheetml/2006/main'}
    with zipfile.ZipFile(path) as z:
        shared=[]
        if 'xl/sharedStrings.xml' in z.namelist():
            root=ET.fromstring(z.read('xl/sharedStrings.xml'));shared=[''.join(item.itertext()) for item in root]
        wb=ET.fromstring(z.read('xl/workbook.xml'));rel={r.attrib['Id']:r.attrib['Target'] for r in ET.fromstring(z.read('xl/_rels/workbook.xml.rels'))}
        sheet=next(s for s in wb.find('s:sheets',ns) if s.attrib['name']=='b. Allelic effect (Veh)');target=rel[sheet.attrib['{http://schemas.openxmlformats.org/officeDocument/2006/relationships}id']]
        target=target.lstrip('/') if target.startswith('/') else posixpath.normpath('xl/'+target);root=ET.fromstring(z.read(target));table=[]
        for row in root.find('s:sheetData',ns):
            values={}
            for cell in row:
                col=re.sub(r'[0-9]','',cell.attrib['r']);kind=cell.attrib.get('t','n');value=cell.find('s:v',ns)
                if kind=='inlineStr':text=''.join(cell.find('s:is',ns).itertext())
                elif value is None:text=''
                elif kind=='s':text=shared[int(value.text)]
                else:text=value.text
                values[col]=text
            table.append(values)
        header=table[0];idcol=next(c for c,v in header.items() if v=='SNP');effectcol=next(c for c,v in header.items() if v=='logFC');effects={}
        for row in table[1:]:
            if not row.get(idcol):continue
            identity=row[idcol];assert identity not in effects;effects[identity]=abs(float(row[effectcol]))
        assert len(effects)==661;return effects


def main():
    lock=json.loads((FINAL/'before_fresh_outcomes_lock.json').read_text())
    for name,digest in lock['files'].items():assert sha(name)==digest,name
    expected=json.loads((FINAL/'fresh_outcome_evaluation.json').read_text());rows=list(csv_rows(FINAL/'fresh_study_predictions_and_observations.csv'));source=original_effects(V8/'sources/Matoba_journal_Supplementary_Table3.xlsx')
    y=np.array([source[r['rsid']] for r in rows]);stored=np.array([float(r['measured_absolute_original_log2_effect']) for r in rows]);assert np.array_equal(y,stored)
    assert len(rows)==610 and len(set(r['rsid'] for r in rows))==610
    groups=np.array([r['bootstrap_1Mb_window'] for r in rows]);u,inv=np.unique(groups,return_inverse=True);n=np.bincount(inv);rng=np.random.default_rng(20261108)
    draws=rng.integers(len(u),size=(2000,len(u)));multiplicities=np.stack([np.bincount(row,minlength=len(u)) for row in draws]);denominator=multiplicities@n;checks={}
    for arm in ['sequence','AVI']:
        p=np.array([float(r[f'{arm}_magnitude']) for r in rows]);b=np.array([float(r[f'{arm}_baseline']) for r in rows]);cb=np.array([float(r[f'{arm}_cell_baseline']) for r in rows]);hw=np.array([float(r[f'{arm}_interval_halfwidth']) for r in rows]);e=abs(y-p);be=abs(y-b)
        m=expected['arms'][arm]['all_cases'];assert abs(e.mean()-m['MAE'])<1e-12;assert abs(be.mean()-m['line_or_cell_baseline_MAE'])<1e-12
        ranked=float(np.corrcoef(rankdata(y),rankdata(p))[0,1]);assert abs(ranked-m['magnitude_Spearman'])<1e-12
        group_gain=np.bincount(inv,weights=be-e);boot=multiplicities@group_gain/denominator;ci=np.quantile(boot,[.025,.975]);assert np.max(abs(ci-np.array(m['absolute_error_improvement_vs_line_baseline_CI95'])))<1e-10
        coverage=float(np.mean((y>=np.maximum(0,p-hw))&(y<=p+hw)));assert coverage==expected['arms'][arm]['all_case_interval_coverage']
        assert not expected['arms'][arm]['full_model_quality_gate_passed'] and not expected['arms'][arm]['uncertainty_and_delivery_gate_passed']
        checks[arm]={'source_values_identical':True,'MAE_recomputed':float(e.mean()),'rank_correlation_independently_recomputed':ranked,'alternative_aggregate_bootstrap_CI95':ci.tolist(),'coverage_from_explicit_interval_bounds':coverage,'failure_gates_consistent':True}
    projects=list(csv_rows(FINAL/'project_answer_availability.csv'));assert len(projects)==2490
    assert all(r['authorized']=='False' and r['authorized_prediction']==r['authorized_rank']=='' for r in projects)
    assert not expected['AVI_comparison']['retain_AVI'];assert expected['AVI_comparison']['relative_MAE_improvement_vs_sequence']<0
    write_json(FINAL/'independent_result_audit.json',{'passed':True,'source_reader':'Separate direct XML tree parser; did not import the primary worksheet extraction helper.','checks':checks,'all_project_authorization_and_rank_fields_withheld':True,'all_pre_outcome_locked_files_unchanged':True,'does_not_establish_scientific_success':True})
    assert_previous_sources();print('INDEPENDENT SOURCE, METRIC, BOOTSTRAP AND WITHHOLDING AUDIT PASS',flush=True)


if __name__=='__main__':main()
