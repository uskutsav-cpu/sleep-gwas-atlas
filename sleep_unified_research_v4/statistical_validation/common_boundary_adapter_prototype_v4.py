#!/usr/bin/env python3
"""Result-free prototype for a separate pinned LDSC covariance experiment.

No GWAS runner/CLI is provided. Running this file executes small mock-interface
fixtures only. Real interval/source admission and calibration remain pending.
The adapter never edits pinned source; instrumentation is scoped and restored.
"""
from contextlib import contextmanager
from copy import deepcopy
from dataclasses import dataclass
from pathlib import Path
import hashlib
import json
from types import SimpleNamespace
import ast


class BoundaryAdmissionError(ValueError):
    pass


@dataclass(frozen=True)
class Interval:
    block_id: int
    chromosome: int
    start_inclusive: int
    end_exclusive: int


class CanonicalAdapter:
    """Thread-unsafe by design: one worker, one active pair, ordered callbacks."""

    def __init__(self, intervals, snp_coordinates, *, build, reference_sha256, fixture=False):
        self.intervals = tuple(intervals)
        self.coordinates = dict(snp_coordinates)
        self.build = build
        self.reference_sha256 = dict(reference_sha256)
        self.fixture = fixture
        self.pair = None
        self.frame = None
        self.audit = []
        if build != 'GRCh37/hg19':
            raise BoundaryAdmissionError('BUILD_NOT_ADMITTED')
        if not self.intervals or [x.block_id for x in self.intervals] != list(range(len(self.intervals))):
            raise BoundaryAdmissionError('NONCANONICAL_BLOCK_IDS')
        if len(self.intervals) != 200:
            raise BoundaryAdmissionError('PRIMARY_PROTOTYPE_REQUIRES_200_INTERVALS')
        for i, x in enumerate(self.intervals):
            if not 1 <= x.chromosome <= 22 or x.start_inclusive < 1 or x.end_exclusive <= x.start_inclusive:
                raise BoundaryAdmissionError('INVALID_GENOMIC_INTERVAL')
            if i and (self.intervals[i-1].chromosome, self.intervals[i-1].end_exclusive) > (x.chromosome, x.start_inclusive):
                raise BoundaryAdmissionError('UNORDERED_OR_OVERLAPPING_INTERVALS')
        if not fixture and (len(reference_sha256) != 22 or any(len(x) != 64 for x in reference_sha256.values())):
            raise BoundaryAdmissionError('EXACT_22_REFERENCE_HASHES_REQUIRED')

    def _assign(self, snp_ids):
        ids = tuple(snp_ids)
        if len(set(ids)) != len(ids):
            raise BoundaryAdmissionError('DUPLICATE_FINAL_SNP_IDS')
        try:
            coords = tuple(self.coordinates[x] for x in ids)
        except KeyError as e:
            raise BoundaryAdmissionError('FINAL_SNP_COORDINATE_UNAVAILABLE') from e
        if any(not (1 <= c <= 22 and p >= 1) for c, p in coords):
            raise BoundaryAdmissionError('INVALID_COORDINATE')
        if any(a > b for a, b in zip(coords, coords[1:])):
            raise BoundaryAdmissionError('FINAL_SNP_ORDER_NOT_GENOMIC')
        labels = []
        interval_index = 0
        for c, bp in coords:
            while interval_index < len(self.intervals):
                x = self.intervals[interval_index]
                if (c, bp) >= (x.chromosome, x.end_exclusive):
                    interval_index += 1
                else:
                    break
            if interval_index >= len(self.intervals):
                raise BoundaryAdmissionError('FINAL_SNP_OUTSIDE_MANIFEST')
            x = self.intervals[interval_index]
            if c != x.chromosome or not x.start_inclusive <= bp < x.end_exclusive:
                raise BoundaryAdmissionError('FINAL_SNP_OUTSIDE_MANIFEST')
            labels.append(interval_index)
        return ids, coords, tuple(labels)

    def _separators(self, labels):
        counts = [0]*len(self.intervals)
        for b in labels: counts[b] += 1
        # Strict prototype policy requested for both stages: fail on empty blocks.
        if any(c == 0 for c in counts):
            raise BoundaryAdmissionError('EMPTY_CANONICAL_BLOCK_IN_FIT_OR_STAGE')
        separators = [0]
        for count in counts: separators.append(separators[-1]+count)
        if len(labels) != separators[-1]:
            raise BoundaryAdmissionError('BOUNDARY_ROW_COUNT_MISMATCH')
        return separators, counts

    @contextmanager
    def pair_scope(self, snp_ids, n_blocks):
        if self.pair is not None:
            raise BoundaryAdmissionError('NESTED_OR_CONCURRENT_PAIR_UNSUPPORTED')
        if n_blocks != len(self.intervals):
            raise BoundaryAdmissionError('BLOCK_COUNT_NOT_FROZEN_200')
        ids, coords, labels = self._assign(snp_ids)
        s, counts = self._separators(labels)
        self.pair = dict(ids=ids, coords=coords, labels=labels, separators=s, counts=counts)
        try:
            yield
        finally:
            self.frame = None
            self.pair = None

    @contextmanager
    def regression_scope(self, n_rows, step1_ii, old_weights):
        if self.pair is None or self.frame is not None or n_rows != len(self.pair['ids']):
            raise BoundaryAdmissionError('BASE_REGRESSION_CONTEXT_MISMATCH')
        if old_weights:
            raise BoundaryAdmissionError('OLD_WEIGHTS_PATH_UNSUPPORTED')
        if step1_ii is None:
            stages = [('single_step', self.pair['labels'])]
        else:
            mask = step1_ii.reshape(-1).tolist() if hasattr(step1_ii, 'reshape') else list(step1_ii)
            if len(mask) != n_rows:
                raise BoundaryAdmissionError('TWO_STEP_MASK_LENGTH_MISMATCH')
            stages = [('two_step_filtered', tuple(b for b, keep in zip(self.pair['labels'], mask) if keep)),
                      ('two_step_full', self.pair['labels'])]
        # Validate both masks before any IRWLS call.
        prepared = [(name, labels, *self._separators(labels)) for name, labels in stages]
        self.frame = {'stages':prepared,'cursor':0}
        try:
            yield
            if self.frame['cursor'] != len(stages):
                raise BoundaryAdmissionError('UNEXPECTED_IRWLS_CALL_COUNT')
        finally:
            self.frame = None

    def irwls_separators(self, n_rows, incoming_separators):
        if self.frame is None or self.frame['cursor'] >= len(self.frame['stages']):
            raise BoundaryAdmissionError('IRWLS_OUTSIDE_EXPECTED_BASE_CONTEXT')
        name, labels, separators, counts = self.frame['stages'][self.frame['cursor']]
        if n_rows != len(labels):
            raise BoundaryAdmissionError('IRWLS_STAGE_ROW_COUNT_MISMATCH')
        self.frame['cursor'] += 1
        self.audit.append({'stage':name,'snp_count':n_rows,'block_counts':counts,
                           'canonical_separators':separators,
                           'incoming_separators_used':False,
                           'incoming_separators_present':incoming_separators is not None,
                           'ordered_final_snp_sha256':hashlib.sha256(('\n'.join(self.pair['ids'])+'\n').encode()).hexdigest()})
        return separators

    @contextmanager
    def instrument_modules(self, sumstats, reg, *, scope):
        if scope != 'FIXTURE' or not self.fixture:
            # Deliberate result-free guard. Remove only in a separately frozen,
            # independently checked method-validation executor, never here.
            raise BoundaryAdmissionError('PROTOTYPE_NATIVE_EXECUTION_NOT_ADMITTED')
        original_rg = sumstats._rg
        original_base = reg.LD_Score_Regression.__init__
        original_irwls = reg.IRWLS
        adapter = self

        def canonical_rg(rows, args, log, M_annot, ref_ld_cnames, w_ld_cname, i):
            # Pinned sumstats._rg:522–525 applies precisely this optional filter.
            # This adapter enters pair context AFTER the last filter, and passes
            # chisq_max=None to avoid performing it twice in the original call.
            final_rows = rows
            inner_args = deepcopy(args)
            if args.chisq_max is not None:
                keep = rows.Z1**2 * rows.Z2**2 < args.chisq_max**2
                final_rows = rows[keep]
                inner_args.chisq_max = None
            ids = final_rows.SNP.tolist()
            with adapter.pair_scope(ids, min(args.n_blocks, len(ids))):
                return original_rg(final_rows, inner_args, log, M_annot, ref_ld_cnames, w_ld_cname, i)

        def canonical_base(self, y, x, w, N, M, n_blocks, intercept=None, slow=False, step1_ii=None, old_weights=False):
            if n_blocks != len(adapter.intervals):
                raise BoundaryAdmissionError('BASE_BLOCK_COUNT_MISMATCH')
            with adapter.regression_scope(x.shape[0], step1_ii, old_weights):
                return original_base(self,y,x,w,N,M,n_blocks,intercept=intercept,slow=slow,step1_ii=step1_ii,old_weights=old_weights)

        class CanonicalIRWLS(original_irwls):
            def __init__(self, x, y, update_func, n_blocks, w=None, slow=False, separators=None):
                canonical = adapter.irwls_separators(x.shape[0], separators)
                super().__init__(x,y,update_func,n_blocks,w=w,slow=slow,separators=canonical)

        sumstats._rg = canonical_rg
        reg.LD_Score_Regression.__init__ = canonical_base
        reg.IRWLS = CanonicalIRWLS
        try:
            yield
        finally:
            sumstats._rg = original_rg
            reg.LD_Score_Regression.__init__ = original_base
            reg.IRWLS = original_irwls


def verify_pinned_interfaces():
    root=Path(__file__).resolve().parents[2]
    plan=json.loads((root/'sleep_unified_research_v4/manifests/ssd_native_execution_plan_v4_3.json').read_text())
    hashes={}
    trees={}
    for name in ['sumstats.py','regressions.py','irwls.py','jackknife.py']:
        path=root.parent/'ldsc-code/ldscore'/name
        source=path.read_bytes();actual=hashlib.sha256(source).hexdigest()
        assert actual==plan['dependencies_sha256'][str(path)]
        hashes[name]=actual;trees[name]=ast.parse(source)
    base=next(n for n in trees['regressions.py'].body if isinstance(n,ast.ClassDef) and n.name=='LD_Score_Regression')
    init=next(n for n in base.body if isinstance(n,ast.FunctionDef) and n.name=='__init__')
    assert [a.arg for a in init.args.args]==['self','y','x','w','N','M','n_blocks','intercept','slow','step1_ii','old_weights']
    irwls=next(n for n in trees['irwls.py'].body if isinstance(n,ast.ClassDef) and n.name=='IRWLS')
    init=next(n for n in irwls.body if isinstance(n,ast.FunctionDef) and n.name=='__init__')
    assert [a.arg for a in init.args.args]==['self','x','y','update_func','n_blocks','w','slow','separators']
    rg=next(n for n in trees['sumstats.py'].body if isinstance(n,ast.FunctionDef) and n.name=='_rg')
    assert [a.arg for a in rg.args.args]==['sumstats','args','log','M_annot','ref_ld_cnames','w_ld_cname','i']
    # Exact final filtering branch is audited, rather than inferred from logs.
    final_filter=next(n for n in rg.body if isinstance(n,ast.If))
    expected=ast.parse('if args.chisq_max is not None:\n ii = sumstats.Z1**2*sumstats.Z2**2 < args.chisq_max**2\n n_snp = np.sum(ii)\n sumstats = sumstats[ii]').body[0]
    assert ast.dump(final_filter,include_attributes=False)==ast.dump(expected,include_attributes=False)
    return {'matched_code_hashes':hashes,'base_and_irwls_signatures_verified':True,
            'exact_final_rg_filter_ast_verified':True,'native_import_or_regression_run':False}


def fixture_checks():
    intervals=[Interval(b,1,10*b+1,10*b+11) for b in range(200)]
    coords={f'rs{b}_{j}':(1,10*b+j) for b in range(200) for j in (1,2,3)}
    ids=list(coords)
    adapter=CanonicalAdapter(intervals,coords,build='GRCh37/hg19',reference_sha256={},fixture=True)
    # Mock interfaces match pinned Base/IRWLS signatures. No regression is fit.
    class Matrix:
        def __init__(self,n): self.shape=(n,1)
    class FakeIRWLS:
        def __init__(self,x,y,update_func,n_blocks,w=None,slow=False,separators=None):
            self.separators=separators
    class FakeBase:
        def __init__(self,y,x,w,N,M,n_blocks,intercept=None,slow=False,step1_ii=None,old_weights=False):
            if step1_ii is not None:
                filtered=Matrix(sum(step1_ii))
                first=reg.IRWLS(filtered,filtered,None,n_blocks)
                # Incoming mapped/default cuts intentionally wrong: wrapper
                # must replace them with full-row canonical coordinates.
                reg.IRWLS(x,y,None,n_blocks,separators=[0, x.shape[0]])
            else: reg.IRWLS(x,y,None,n_blocks)
    class SNPSeries:
        def __init__(self,data):self.data=data
        def tolist(self):return self.data
    rows=SimpleNamespace(SNP=SNPSeries(ids))
    reg=SimpleNamespace(LD_Score_Regression=FakeBase,IRWLS=FakeIRWLS)
    def fake_rg(rows,args,log,M_annot,ref_ld_cnames,w_ld_cname,i):
        m=Matrix(len(rows.SNP.tolist()))
        masks=[None,[j%3!=0 for j in range(m.shape[0])],[j%3!=1 for j in range(m.shape[0])]]
        for mask in masks: reg.LD_Score_Regression(m,m,None,None,None,200,step1_ii=mask)
        return 'FIXTURE_NO_ESTIMATE'
    sums=SimpleNamespace(_rg=fake_rg)
    old_base=FakeBase.__init__
    with adapter.instrument_modules(sums,reg,scope='FIXTURE'):
        result=sums._rg(rows,SimpleNamespace(chisq_max=None,n_blocks=200),None,None,None,None,0)
    assert result=='FIXTURE_NO_ESTIMATE'
    assert sums._rg is fake_rg and reg.IRWLS is FakeIRWLS and FakeBase.__init__ is old_base
    assert [x['stage'] for x in adapter.audit]==['single_step','two_step_filtered','two_step_full','two_step_filtered','two_step_full']
    assert [set(x['block_counts']) for x in adapter.audit]==[{3},{2},{3},{2},{3}]
    assert all(x['canonical_separators'][-1]==x['snp_count'] and len(x['canonical_separators'])==201 for x in adapter.audit)
    assert adapter.audit[2]['incoming_separators_present'] and adapter.audit[2]['canonical_separators']!=[0,600]
    # Exercise the optional final Z-product filter without importing NumPy.
    class Vector:
        def __init__(self,data):self.data=data
        def __pow__(self,p):return Vector([x**p for x in self.data])
        def __mul__(self,other):return Vector([x*y for x,y in zip(self.data,other.data)])
        def __lt__(self,value):return [x<value for x in self.data]
    class FinalRows:
        def __init__(self,ids,z1,z2):self.SNP=SNPSeries(ids);self.Z1=Vector(z1);self.Z2=Vector(z2)
        def __getitem__(self,keep):
            return FinalRows([x for x,k in zip(self.SNP.data,keep) if k],
                             [x for x,k in zip(self.Z1.data,keep) if k],
                             [x for x,k in zip(self.Z2.data,keep) if k])
    filtered_adapter=CanonicalAdapter(intervals,coords,build='GRCh37/hg19',reference_sha256={},fixture=True)
    captured={}
    def filtered_original(rows,args,log,M_annot,ref_ld_cnames,w_ld_cname,i):
        captured['ids']=rows.SNP.tolist();captured['chisq_max']=args.chisq_max
        m=Matrix(len(captured['ids']))
        reg.LD_Score_Regression(m,m,None,None,None,200)
        return 'FINAL_FILTER_FIXTURE_NO_ESTIMATE'
    sums2=SimpleNamespace(_rg=filtered_original)
    filter_input=FinalRows(ids,[10 if j%3==0 else 1 for j in range(600)],[1]*600)
    filter_args=SimpleNamespace(chisq_max=5,n_blocks=200)
    with filtered_adapter.instrument_modules(sums2,reg,scope='FIXTURE'):
        filter_result=sums2._rg(filter_input,filter_args,None,None,None,None,0)
    assert filter_result=='FINAL_FILTER_FIXTURE_NO_ESTIMATE' and len(captured['ids'])==400
    assert captured['chisq_max'] is None and filter_args.chisq_max==5 and len(filter_input.SNP.tolist())==600
    assert set(filtered_adapter.audit[0]['block_counts'])=={2}
    failures={}
    for name, altered_ids, mask in [('reversed_source_order',list(reversed(ids)),None),
                                    ('missing_source_coordinate',ids+['unknown'],None),
                                    ('duplicate_source_id',ids+ids[-1:],None),
                                    ('empty_final_block',ids[3:],None),
                                    ('empty_two_step_block',ids,[j>=3 for j in range(600)])]:
        try:
            with adapter.pair_scope(altered_ids,200):
                if mask is not None:
                    with adapter.regression_scope(len(altered_ids),mask,False): pass
        except BoundaryAdmissionError as e: failures[name]=str(e)
    assert len(failures)==5
    for name, operation in [
        ('wrong_genome_build',lambda:CanonicalAdapter(intervals,coords,build='GRCh38',reference_sha256={},fixture=True)),
        ('wrong_block_count',lambda:adapter.pair_scope(ids,100).__enter__())]:
        try:operation()
        except BoundaryAdmissionError as e:failures[name]=str(e)
    with adapter.pair_scope(ids,200):
        try:
            with adapter.regression_scope(600,None,True):pass
        except BoundaryAdmissionError as e:failures['unsupported_old_weights_path']=str(e)
        try:
            with adapter.regression_scope(600,None,False):pass
        except BoundaryAdmissionError as e:failures['unexpected_IRWLS_call_count']=str(e)
    # Module restoration must also occur after the final-pair context rejects.
    try:
        with adapter.instrument_modules(sums,reg,scope='FIXTURE'):
            sums._rg(SimpleNamespace(SNP=SNPSeries(list(reversed(ids)))),SimpleNamespace(chisq_max=None,n_blocks=200),None,None,None,None,0)
    except BoundaryAdmissionError:pass
    assert sums._rg is fake_rg and reg.IRWLS is FakeIRWLS and FakeBase.__init__ is old_base
    try:
        with adapter.instrument_modules(sums,reg,scope='METHOD_VALIDATION'): pass
    except BoundaryAdmissionError as e: failures['native_execution_guard']=str(e)
    assert 'native_execution_guard' in failures
    return {'status':'FIXTURE_INTERFACE_PASS_ONLY_NOT_NATIVE_ESTIMATOR_VALIDATION',
            'canonical_intervals':200,'fixture_final_snps':600,'stages_checked':len(adapter.audit),
            'nonempty_blocks_enforced':True,'restored_all_original_module_bindings':True,
            'restored_bindings_after_exception':True,
            'first_stage_mask_specific_canonical_counts':True,'second_stage_mapped_index_cuts_replaced':True,
            'optional_final_rg_filter_context_verified':True,'filtered_input_and_args_preserved':True,
            'pinned_source_interface_checks':verify_pinned_interfaces(),
            'fail_closed_cases':failures,'audit':adapter.audit,
            'native_estimator_executed':False,'real_boundary_manifest_created':False,
            'calibrated_covariance_established':False}


if __name__=='__main__':
    result=fixture_checks()
    out=Path(__file__).with_name('common_boundary_adapter_fixture_receipt_v4.json')
    result['script_sha256']=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    out.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({k:v for k,v in result.items() if k!='audit'},indent=2))
