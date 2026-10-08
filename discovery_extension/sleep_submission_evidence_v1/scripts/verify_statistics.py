#!/usr/bin/env python3
"""Versioned source-free arithmetic audit; never changes historical results.

Uses standard-library erfc and reverse-rank BH. Output P values are calculations
from printed estimates, not recovered unrounded jackknife estimates or native
LDSC reruns. All sensitivity analyses are retrospective, descriptive audits.
"""
import argparse
import csv
import hashlib
import json
import math
import re
from collections import Counter
from pathlib import Path

Z95 = 1.959963984540054
NA = 'NA'

def read(path):
    with Path(path).open(newline='') as f:
        return list(csv.DictReader(f, delimiter='\t'))

def write(path, rows):
    if not rows:
        raise ValueError('empty table requires explicit schema')
    fields = list(dict.fromkeys(k for row in rows for k in row))
    with Path(path).open('w', newline='') as f:
        w = csv.DictWriter(f, fields, delimiter='\t', lineterminator='\n', restval=NA)
        w.writeheader()
        w.writerows(rows)

def numeric(value):
    if value in ('', 'NA', 'NOT_AVAILABLE', None):
        return None
    x = float(value)
    if not math.isfinite(x):
        raise ValueError('nonfinite numerical input')
    return x

def normal_p(z):
    if not math.isfinite(z):
        raise ValueError('nonfinite Z')
    return math.erfc(abs(z) / math.sqrt(2))

def bh(p):
    if not p or any(not math.isfinite(x) or not 0 <= x <= 1 for x in p):
        raise ValueError('BH requires a complete finite P family')
    order = sorted(range(len(p)), key=lambda i: p[i])
    out = [None] * len(p)
    q = 1.0
    for rank in reversed(range(len(p))):
        i = order[rank]
        q = min(q, p[i] * len(p) / (rank + 1))
        out[i] = q
    return out

def estimates(effect, se):
    if not all(math.isfinite(x) for x in (effect, se)) or se <= 0:
        raise ValueError('estimate needs finite effect and positive SE')
    return effect / se, normal_p(effect / se), effect - Z95 * se, effect + Z95 * se

def difference(a, sa, b, sb, rho=0):
    estimates(a, sa); estimates(b, sb)
    if not -1 <= rho <= 1:
        raise ValueError('invalid assumed estimator correlation')
    variance = sa**2 + sb**2 - 2 * rho * sa * sb
    if variance <= 0:
        raise ValueError('nonpositive difference variance')
    delta = b - a
    se = math.sqrt(variance)
    z = delta / se
    return delta, se, z, z*z, normal_p(z), delta-Z95*se, delta+Z95*se

def unique(rows, field):
    vals = [r[field] for r in rows]
    if len(set(vals)) != len(vals):
        raise ValueError('duplicate identity: ' + field)
    return set(vals)

def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()

def printed_logs(root):
    """Independent lightweight parser of the retained replication summary tables."""
    out = {}
    for path in sorted((root/'discovery_extension/logs/replication/rg').glob('*.log')):
        lines = path.read_text().splitlines()
        if any('SYNTHETIC' in line for line in lines):
            raise ValueError('synthetic log in real result namespace')
        start = next(i for i,l in enumerate(lines) if l == 'Summary of Genetic Correlation Results')
        header = lines[start+1].split()
        scalar_p = [line.split(':',1)[1].strip() for line in lines[:start] if line.startswith('P:')]
        sleep_block=path.read_text().split('Heritability of phenotype 1\n')[1].split('Heritability of phenotype ')[0]
        hm=re.search(r'Total (\w+) scale h2: ([\d.eE+-]+) \(([\d.eE+-]+)\)',sleep_block)
        im=re.search(r'Intercept: ([\d.eE+-]+)',sleep_block)
        sleep_h2=float(hm[2]);sleep_se=float(hm[3]);sleep_intercept=float(im[1])
        if sleep_se<=0:raise ValueError('invalid sleep h2 SE in retained log')
        body_index=0
        for line in lines[start+2:]:
            values = line.split()
            if len(values) != len(header):
                break
            r = dict(zip(header, values))
            sleep = Path(r['p1']).name.removesuffix('.sumstats.gz')
            source = Path(r['p2']).name.removesuffix('.sumstats.gz')
            key = (sleep, source)
            if key in out:
                raise ValueError('duplicate printed log result')
            out[key] = r | {'log_path': str(path.relative_to(root)), 'scalar_p_display':scalar_p[body_index],
                            'sleep_h2':sleep_h2,'sleep_h2_se':sleep_se,'sleep_h2_z':sleep_h2/sleep_se,
                            'sleep_h2_intercept':sleep_intercept}
            body_index+=1
        if body_index!=len(scalar_p):raise ValueError('printed scalar P/result table count mismatch')
    return out

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--root', type=Path, default=Path(__file__).resolve().parents[3])
    ap.add_argument('--out', type=Path)
    args = ap.parse_args(); root = args.root.resolve()
    out = args.out or root/'discovery_extension/sleep_submission_evidence_v1'
    (out/'tables').mkdir(parents=True, exist_ok=True); (out/'logs').mkdir(parents=True, exist_ok=True)
    paths = {
        'global': 'discovery_extension/results/ldsc/extension_rg_matrix.tsv',
        'universe': 'discovery_extension/results/ldsc/extension_pair_universe.tsv',
        'replication': 'discovery_extension/results/replication/replication_results.tsv',
        'manifest': 'discovery_extension/config/replication_manifest.tsv',
        'lock': 'discovery_extension/config/replication_manifest.lock.json',
        'candidate_lock': 'discovery_extension/config/replication_candidate_family.lock.json',
        'sleep_panel': 'config/analysis_panel.tsv',
        'external_panel': 'discovery_extension/config/candidate_traits.tsv',
        'h2': 'discovery_extension/results/ldsc/extension_trait_readiness.tsv',
        'local_queue':'discovery_extension/results/local/local_analysis_queue.tsv',
        'replication_h2':'discovery_extension/results/replication/replication_source_h2.tsv',
    }
    g = read(root/paths['global']); rep = read(root/paths['replication'])
    manifest = {r['pair_id']: r for r in read(root/paths['manifest'])}
    panel = {r['extension_trait_id']: r for r in read(root/paths['external_panel'])}
    sleep = {r['trait_id']: r for r in read(root/paths['sleep_panel']) if r['domain']=='sleep'}
    expected = {(s,p) for s in sleep for p in panel}
    actual = {(r['sleep_trait'],r['extension_trait_id']) for r in g}
    if len(g)!=len(actual) or actual != expected:
        raise ValueError('incomplete or duplicate discovery universe')
    universe = {(r['sleep_trait'],r['extension_trait_id']) for r in read(root/paths['universe'])}
    if universe != expected:
        raise ValueError('historical pair universe disagreement')
    lock = json.loads((root/paths['lock']).read_text())
    family = lock['pair_ids_in_locked_order']
    if unique(rep,'pair_id') != set(family) or set(manifest) != set(family):
        raise ValueError('incomplete replication universe')
    alpha = .05/len(family)
    if not math.isclose(alpha,float(lock['bonferroni_alpha']),rel_tol=1e-12):
        raise ValueError('wrong locked correction denominator')
    q = bh([float(r['p']) for r in g])
    # Four-decimal display uncertainty: conservative absolute half-unit bounds.
    rounding_low_p = [normal_p((abs(float(r['rg']))+.00005)/(float(r['se'])-.00005)) for r in g]
    rounding_high_p = [normal_p(max(0,abs(float(r['rg']))-.00005)/(float(r['se'])+.00005)) for r in g]
    rounding_low_q, rounding_high_q = bh(rounding_low_p), bh(rounding_high_p)
    qz = bh([normal_p(float(r['z'])) for r in g])
    qr = bh([estimates(float(r['rg']),float(r['se']))[1] for r in g])
    sensitivities = []
    def sens(kind, scope, n, hits, value=NA, details=''):
        sensitivities.append({'analysis_id':kind,'scope':scope,'denominator':n,'n_significant':hits,'value':value,'details':details,'design':'RETROSPECTIVE_SOURCE_FREE','native_GWAS_rerun':'False'})
    for i,r in enumerate(g):
        est,se=float(r['rg']),float(r['se']); z,p,lo,hi=estimates(est,se)
        if not math.isclose(q[i],float(r['extension_fdr']),rel_tol=2e-13,abs_tol=1e-300):
            raise ValueError('historical BH mismatch')
        r.update({'pair_id':r['sleep_trait']+'__'+r['extension_trait_id'],
                  'p_original_printed':r['p'],'z_from_printed_rg_se':z,'p_from_printed_rg_se':p,
                  'p_from_printed_z':normal_p(float(r['z'])),'fdr_recomputed_original_p':q[i],
                  'fdr_from_printed_z':qz[i],'fdr_from_printed_rg_se':qr[i],
                  'ci95_low':lo,'ci95_high':hi,'historical_significant':q[i]<.05,
                  'p_rounding_lower_bound':rounding_low_p[i],'p_rounding_upper_bound':rounding_high_p[i],
                  'fdr_rounding_lower_bound':rounding_low_q[i],'fdr_rounding_upper_bound':rounding_high_q[i],
                  'precision_sensitivity_changed':(q[i]<.05)!=(qz[i]<.05) or (q[i]<.05)!=(qr[i]<.05),
                  'sleep_source_id':sleep[r['sleep_trait']]['source_id'],
                  'external_source_accession':panel[r['extension_trait_id']]['study_accession'],
                  'sleep_measurement':'ACCELEROMETER' if r['sleep_trait'] in ('sleep_efficiency','accel_sleep_duration','sleep_timing') else 'REGISTER_DIAGNOSIS' if r['sleep_trait']=='sleep_apnea' else 'SELF_REPORT',
                  'source_result_path':paths['global'], 'numerical_verification':'RECOMPUTED_FROM_PRINTED_ESTIMATES',
                  'native_reproduction':'BLOCKED_MISSING_INPUTS','full_precision_estimates_available':False})
    global_lookup={(r['sleep_trait'],r['extension_trait_id']):r for r in g}
    h2=read(root/paths['h2'])
    if unique(h2,'extension_trait_id')!=set(panel):raise ValueError('incomplete h2 source family')
    for r in h2:
        if not math.isclose(float(r['h2_z']),estimates(float(r['h2']),float(r['h2_se']))[0],rel_tol=1e-12):raise ValueError('external h2 Z mismatch')
    if any(float(r['h2_z'])<4 or float(r['LDSC_intercept'])>1.2 or r['primary_rg_eligibility']!='PRIMARY_PASS' for r in h2):raise ValueError('primary discovery h2 eligibility disagreement')
    sig=[r for r in g if r['historical_significant']]
    write(out/'tables/global_1200.tsv',g); write(out/'tables/significant_603.tsv',sig)
    local=read(root/paths['local_queue'])
    null_ids={r['pair_id'] for r in g if not r['historical_significant']}
    if unique(local,'pair_id')!={r['pair_id'] for r in local}:raise ValueError('duplicate local queue')
    if {r['pair_id'] for r in local if r['selection_status']=='PENDING_GLOBAL_NULL_CURATION'}!=null_ids:raise ValueError('local global-null candidate set differs from non-BH universe')
    if {r['pair_id'] for r in local if r['selection_status']=='SELECTED_PRIORITY_DISCOVERY'}!=set(family):raise ValueError('local priority set differs from217 candidate universe')
    rh2=read(root/paths['replication_h2']);unique(rh2,'replication_source_id')
    for r in rh2:
        htext=(root/r['input_log']).read_text()
        hm=re.search(r'Total (\w+) scale h2: ([\d.eE+-]+) \(([\d.eE+-]+)\)',htext)
        im=re.search(r'^Intercept: ([\d.eE+-]+)',htext,re.M)
        if float(r['h2'])!=float(hm[2]) or float(r['h2_se'])!=float(hm[3]) or float(r['LDSC_intercept'])!=float(im[1]):raise ValueError('replication standalone h2 log/table mismatch')
        hz=estimates(float(r['h2']),float(r['h2_se']))[0]
        if not math.isclose(hz,float(r['h2_z']),rel_tol=1e-12):raise ValueError('replication h2 Z mismatch')
        if (hz>=4 and float(r['LDSC_intercept'])<=1.2)!=(r['primary_status']=='PASS'):raise ValueError('replication standalone h2 gate mismatch')
    logs = printed_logs(root)
    tested = []
    for r in rep:
        m=manifest[r['pair_id']]
        if (r['sleep_trait'],r['extension_trait_id']) not in actual:
            raise ValueError('replication pair outside discovery')
        discovery=global_lookup[(r['sleep_trait'],r['extension_trait_id'])]
        for rk,gk in [('discovery_rg','rg'),('discovery_se','se'),('discovery_fdr','extension_fdr')]:
            if not math.isclose(float(r[rk]),float(discovery[gk]),rel_tol=1e-12):raise ValueError('replication discovery crosswalk disagreement')
        for key in ('discovery_rg','discovery_se','discovery_fdr'):
            if float(r[key])!=float(m[key]):
                raise ValueError('replication discovery estimate differs from locked manifest')
        a,sa=float(r['discovery_rg']),float(r['discovery_se'])
        az,ap,al,ah=estimates(a,sa)
        r.update({'replication_alpha_exact':alpha,'phenotype_match_status':m['phenotype_match_status'],
                  'discovery_ci95_low':al,'discovery_ci95_high':ah,
                  'sleep_source_id':m['discovery_sleep_source_id'],'sleep_GWAS_reused':True,
                  'two_trait_independence':'NOT_FULLY_INDEPENDENT_SLEEP_REUSED',
                  'source_result_path':paths['replication'],'source_manifest_path':paths['manifest'],
                  'full_precision_estimates_available':False,'native_reproduction':'BLOCKED_MISSING_INPUTS'})
        if r['replication_rg']=='NA':
            if r['replication_class'] not in ('UNDERPOWERED','NO_INDEPENDENT_DATASET'):
                raise ValueError('unestimated row incorrectly classified')
            r.update({'audit_replication_class':r['replication_class'],'numerical_verification':'NOT_ESTIMABLE','heterogeneity_assumption':'NOT_ESTIMABLE'})
            continue
        b,sb=float(r['replication_rg']),float(r['replication_se']); z,p,lo,hi=estimates(b,sb)
        if float(r['replication_h2_z'])<4 or float(r['replication_LDSC_intercept'])>1.2 or r['replication_h2_pass']!='True':raise ValueError('tested replication source fails historical h2 gate')
        r['replication_p_rounding_lower_bound']=normal_p((abs(b)+.00005)/(sb-.00005))
        r['replication_p_rounding_upper_bound']=normal_p(max(0,abs(b)-.00005)/(sb+.00005))
        pz=normal_p(float(r['replication_z'])); delta,dse,dz,dq,dp,dl,dh=difference(a,sa,b,sb)
        log=logs.pop((r['sleep_trait'],r['replication_source_id']))
        for field, key in [('rg','replication_rg'),('se','replication_se'),('z','replication_z')]:
            if float(log[field])!=float(r[key]):
                raise ValueError('retained replication log/table disagreement')
        original_class='REPLICATED' if a*b>0 and pz<alpha else 'DIRECTIONALLY_CONCORDANT' if a*b>0 else 'FAILED_REPLICATION'
        if original_class!=r['replication_class'] or not math.isclose(pz,float(r['replication_p']),rel_tol=2e-13):
            raise ValueError('historical replication classification/P disagreement')
        if not math.isclose(dp,float(r['heterogeneity_p']),rel_tol=2e-13):
            raise ValueError('historical zero-covariance difference calculation differs')
        r.update({'z_from_printed_rg_se':z,'p_from_printed_rg_se':p,'p_from_printed_z':pz,
                  'replication_p_ldsc_display':log['p'],'replication_p_ldsc_scalar_display':log['scalar_p_display'],
                  'sleep_h2_validation_log':log['sleep_h2'],'sleep_h2_se_validation_log':log['sleep_h2_se'],
                  'sleep_h2_z_validation_log':log['sleep_h2_z'],'sleep_h2_intercept_validation_log':log['sleep_h2_intercept'],
                  'sleep_h2_qc_evidence':'RETAINED_REPLICATION_FIRST_PAIR_PRINTED_LOG_NOT_ORIGINAL_DISCOVERY_QC','replication_ci95_low':lo,'replication_ci95_high':hi,
                  'difference_replication_minus_discovery':delta,'difference_se_cov0':dse,
                  'difference_ci95_low_cov0':dl,'difference_ci95_high_cov0':dh,
                  'heterogeneity_p_recomputed_cov0':dp,'heterogeneity_assumption':'ZERO_COVARIANCE_ASSUMPTION_UNVERIFIED',
                  'direction_concordant_recomputed':a*b>0,'p_bonferroni_217':min(1,pz*len(family)),
                  'precision_sensitivity_changed':(pz<alpha)!=(p<alpha),
                  'audit_replication_class':'EXTERNAL_OUTCOME_SIDE_REPLICATION' if original_class=='REPLICATED' else original_class,
                  'numerical_verification':'TABLE_AND_RETAINED_PRINTED_LOG_VERIFIED','replication_printed_log':log['log_path'],
                  'heterogeneity_decision':'DIRECTION_REPLICATED_EFFECT_MAGNITUDE_QUALIFIED' if dp<.05 else 'DIRECTION_REPLICATED_EFFECT_CONSISTENCY_UNRESOLVED',
                  'source_comparison':f"Pan-UKB {panel[r['extension_trait_id']]['phenotype_definition']} versus {m['replication_phenotype_definition']}; sleep source reused",
                  'explanation_supported':'NO_SPECIFIC_CAUSE_ESTABLISHED',
                  'ascertainment_limitation':'MVP veteran cohort' if m['replication_source_id'].startswith('gwas_catalog') else 'FinnGen register endpoint versus UKB phenotype'})
        for rho in (-.9,-.5,0,.5,.9):
            rp=difference(a,sa,b,sb,rho)[4]
            r['heterogeneity_p_rho_'+str(rho).replace('-','minus').replace('.','p')]=rp
            sens('ESTIMATOR_COVARIANCE_GRID',r['pair_id'],1,int(rp<.05),rp,f'assumed rho={rho}; not an estimated sample-overlap correction')
        tested.append(r)
    if logs:
        raise ValueError('unaccounted printed log pairs')
    hetq=bh([float(r['heterogeneity_p_recomputed_cov0']) for r in tested])
    for r,qh in zip(tested,hetq):
        r['heterogeneity_bh_41_cov0']=qh
        r['heterogeneity_bonferroni_41_cov0']=min(1,float(r['heterogeneity_p_recomputed_cov0'])*len(tested))
    replicated=[r for r in rep if r['replication_class']=='REPLICATED']
    heterogeneous=[r for r in replicated if float(r['heterogeneity_p_recomputed_cov0'])<.05]
    write(out/'tables/replication_family_217.tsv',rep); write(out/'tables/replicated_23.tsv',replicated)
    write(out/'tables/heterogeneity_7.tsv',heterogeneous); write(out/'tables/HETEROGENEITY_MASTER.tsv',tested)
    for label,qs in [('ORIGINAL_PRINTED_P',q),('PRINTED_Z_NORMAL_TAIL',qz),('PRINTED_RG_SE_NORMAL_TAIL',qr)]:
        sens('BH_PRECISION',label,len(g),sum(x<.05 for x in qs),details='Family fixed at1200; historical decisions preserved')
    sens('BH_ROUNDING_BOUNDS','MOST_SIGNIFICANT_COMPATIBLE_RG_SE',len(g),sum(x<.05 for x in rounding_low_q),details='Assumes nearest four-decimal rg/se print; all1200 family fixed')
    sens('BH_ROUNDING_BOUNDS','LEAST_SIGNIFICANT_COMPATIBLE_RG_SE',len(g),sum(x<.05 for x in rounding_high_q),details='Assumes nearest four-decimal rg/se print; all1200 family fixed')
    sens('DISCOVERY_BONFERRONI','ALL1200',len(g),sum(float(r['p'])<.05/len(g) for r in g),.05/len(g),'Retrospective stricter threshold, not replacement primary analysis')
    for group in ('SELF_REPORT','ACCELEROMETER','REGISTER_DIAGNOSIS'):
        rows=[r for r in g if r['sleep_measurement']==group]
        sens('MEASUREMENT_STRATIFICATION',group,len(rows),sum(r['historical_significant'] for r in rows),details='Uses full1200 BH labels; measurement categories from source panel')
    for domain in sorted({r['phenotype_domain'] for r in g}):
        rows=[r for r in g if r['phenotype_domain']==domain]
        sens('DOMAIN_STRATIFICATION',domain,len(rows),sum(r['historical_significant'] for r in rows),len({r['extension_trait_id'] for r in rows}),'Counts correlated pairs, not independent biological discoveries')
    for match in ('EXACT','COMPARABLE_WITH_DOCUMENTED_DIFFERENCES'):
        rows=[r for r in rep if r['phenotype_match_status']==match]
        sens('PHENOTYPE_DEFINITION',match,len(rows),sum(r['replication_class']=='REPLICATED' for r in rows),details='Same frozen217 alpha; descriptive subgroup only')
    for s in sleep:
        rows=[r for r in g if r['sleep_trait']==s]
        iz=[float(r['cross_trait_LDSC_intercept'])/float(r['cross_trait_LDSC_intercept_se']) for r in rows]
        sens('CROSS_TRAIT_INTERCEPT_DIAGNOSTIC',s,len(rows),sum(normal_p(z)<.05 for z in iz),max(abs(float(r['cross_trait_LDSC_intercept'])) for r in rows),'Nominal intercept diagnostics; do not establish absence of overlap or estimator covariance')
    for rho in (-.9,-.5,0,.5,.9):
        ps=[difference(float(r['discovery_rg']),float(r['discovery_se']),float(r['replication_rg']),float(r['replication_se']),rho)[4] for r in tested]
        sens('HETEROGENEITY_GRID_SUMMARY',f'all_tested_rho={rho}',len(tested),sum(p<.05 for p in ps),sum(p<.05/len(tested) for p in ps),'value gives Bonferroni count across41; nominal n_significant; hypothetical covariance')
    extids=sorted({r['extension_trait_id'] for r in replicated});lookup={(r['sleep_trait'],r['extension_trait_id']):float(r['rg']) for r in g}
    profiles=[]
    for i,a in enumerate(extids):
        for b in extids[i+1:]:
            x=[lookup[(s,a)] for s in sleep];y=[lookup[(s,b)] for s in sleep]
            mx=sum(x)/len(x);my=sum(y)/len(y)
            corr=sum((u-mx)*(v-my) for u,v in zip(x,y))/math.sqrt(sum((u-mx)**2 for u in x)*sum((v-my)**2 for v in y))
            profiles.append(corr)
            sens('SLEEP_PROFILE_SIMILARITY',a+' versus '+b,len(sleep),NA,corr,'Pearson similarity across12 discovery sleep rg point estimates; descriptive, NOT phenotype genetic correlation; no P value or independence claim')
    sens('NEGATIVE_DISCOVERY','FDR_NOT_SIGNIFICANT',len(g),len(g)-len(sig),details='Failure to reject; neither proof of biological absence nor equivalence')
    sens('NEGATIVE_VALIDATION','TESTED_NOT_BONFERRONI_SIGNIFICANT',len(tested),sum(r['replication_class']!='REPLICATED' for r in tested),details='Concordant nonconfirmation; no direction failures observed')
    write(out/'tables/sensitivity_results.tsv',sensitivities)
    summary={
        'schema_version':'1.0','analysis_design':'RETROSPECTIVE_SOURCE_FREE_NUMERICAL_AUDIT',
        'discovery_tests':len(g),'sleep_traits':len(sleep),'external_traits':len(panel),'bh_significant':len(sig),
        'discovery_rounding_bound_counts':[sum(x<.05 for x in rounding_low_q),sum(x<.05 for x in rounding_high_q)],
        'discovery_precision_decision_changes':sum(r['precision_sensitivity_changed'] for r in g),
        'replication_family':len(rep),'replication_tests':len(tested),'replication_alpha_exact':alpha,
        'historical_replication_counts':dict(sorted(Counter(r['replication_class'] for r in rep).items())),
        'audit_replication_counts':dict(sorted(Counter(r['audit_replication_class'] for r in rep).items())),
        'replicated_external_phenotypes':len(extids),'replicated_sleep_traits':len({r['sleep_trait'] for r in replicated}),
        'replicated_exact':sum(r['phenotype_match_status']=='EXACT' for r in replicated),
        'replicated_comparable':sum(r['phenotype_match_status']=='COMPARABLE_WITH_DOCUMENTED_DIFFERENCES' for r in replicated),
        'replication_rounding_bound_counts':[sum(r['direction_concordant_recomputed'] and r['replication_p_rounding_lower_bound']<alpha for r in tested),sum(r['direction_concordant_recomputed'] and r['replication_p_rounding_upper_bound']<alpha for r in tested)],
        'replication_precision_decision_changes':sum(r['precision_sensitivity_changed'] for r in tested),
        'heterogeneity_replicated_nominal_cov0':len(heterogeneous),
        'heterogeneity_all_tested_nominal_cov0':sum(float(r['heterogeneity_p_recomputed_cov0'])<.05 for r in tested),
        'heterogeneity_all_tested_bonferroni_cov0':sum(float(r['heterogeneity_p_recomputed_cov0'])<.05/len(tested) for r in tested),
        'heterogeneity_all_tested_bh_cov0':sum(qh<.05 for qh in hetq),
        'heterogeneity_covariance':'UNKNOWN_SHARED_SLEEP_INPUT',
        'sensitivity_rows':len(sensitivities),'profile_similarity_range':[min(profiles),max(profiles)],
        'retained_replication_log_rows_verified':len(tested),
        'external_discovery_h2_gate_records_verified':len(h2),
        'replication_standalone_h2_logs_verified':len(rh2),'replication_standalone_h2_pass_sources':sum(r['primary_status']=='PASS' for r in rh2),
        'replication_log_sha256':{str(p.relative_to(root)):sha(p) for p in sorted((root/'discovery_extension/logs/replication').glob('*/*.log'))},
        'sleep_validation_log_h2_traits_checked':len({r['sleep_trait'] for r in tested}),
        'sleep_validation_log_h2_gate_passing_traits':len({r['sleep_trait'] for r in tested if float(r['sleep_h2_z_validation_log'])>=4 and float(r['sleep_h2_intercept_validation_log'])<=1.2}),
        'sleep_original_discovery_h2_qc':'BLOCKED_MISSING_ORIGINAL_CORE_H2',
        'local_queue_total':len(local),'local_queue_global_nulls_matched':len(null_ids),'local_queue_priority_matched':len(family),
        'source_level_native_ldsc_reruns':0,'full_precision_recovered_estimates':0,
        'source_hashes':{key:{'path':p,'sha256':sha(root/p)} for key,p in paths.items()},
        'table_sha256':{p.name:sha(p) for p in sorted((out/'tables').glob('*.tsv')) if p.name in ['global_1200.tsv','significant_603.tsv','replication_family_217.tsv','replicated_23.tsv','heterogeneity_7.tsv','HETEROGENEITY_MASTER.tsv','sensitivity_results.tsv']},
    }
    (out/'tables/numerical_summary.json').write_text(json.dumps(summary,indent=2,sort_keys=True)+'\n')
    (out/'logs/numerical_execution.json').write_text(json.dumps(summary,indent=2,sort_keys=True)+'\n')
    print(json.dumps({k:v for k,v in summary.items() if k not in ('source_hashes','table_sha256','replication_log_sha256')},sort_keys=True))

if __name__=='__main__':
    main()
