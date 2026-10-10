"""Read-only identities at the pinned LDSC intersection; no alternative merging."""
import copy
import hashlib
import json
import math
import struct


def clean(value,np):
    if isinstance(value,np.generic):return clean(value.item(),np)
    if isinstance(value,np.ndarray):return clean(value.tolist(),np)
    if isinstance(value,(list,tuple)):return [clean(x,np) for x in value]
    if isinstance(value,dict):return {k:clean(v,np) for k,v in value.items()}
    if isinstance(value,float) and not math.isfinite(value):return None
    return value


def scalar_attrs(obj,np):
    fields=['tot','tot_se','tot_cov','intercept','intercept_se','mean_chisq','mean_z1z2','lambda_gc','ratio','ratio_se','z','p',
            'coef','coef_se','coef_cov','cat','cat_se','cat_cov','prop','prop_se','prop_cov','M','M_prop','enrichment',
            'constrain_intercept','n_snp','n_blocks','n_annot','twostep_filtered',
            'tot_delete_values','part_delete_values','intercept_delete_values']
    result={k:clean(getattr(obj,k,None),np) for k in fields}
    j=getattr(obj,'jknife',None)
    result['jknife']={k:clean(getattr(j,k,None),np) for k in
                       ['est','jknife_est','jknife_se','jknife_var','jknife_cov','delete_values','separators']}
    result['block_boundary_scope']='Recorded only when returned by stock object; absent separators are not inferred from array index.'
    return result


class IntersectionCapture:
    """Wrap stock read/merge/filter return paths and inspect their unmodified values."""
    def __init__(self,ss,np,pd):
        self.ss,self.np,self.pd=ss,np,pd
        self.rows=[]
        self.merged=None
        self.accepted=None
        self.compatibility=None
        self.current_p2=None
        self.h2_identity=None
        self.originals={n:getattr(ss,n) for n in ['_read_ld_sumstats','_read_other_sumstats','_merge_sumstats_sumstats','_filter_alleles']}

    @staticmethod
    def snp_digest(snps):
        h=hashlib.sha256()
        for snp in snps:h.update((str(snp)+'\n').encode('utf-8'))
        return h.hexdigest()

    @staticmethod
    def float_digest(columns):
        h=hashlib.sha256()
        for row in zip(*columns):
            for number in row:h.update(struct.pack('>d',float(number)))
        return h.hexdigest()

    def mask_info(self,mask):
        a=self.np.asarray(mask,dtype=self.np.uint8).reshape(-1)
        packed=self.np.packbits(a,bitorder='big').tobytes()
        return dict(count=int(a.sum()),total=len(a),packed_big_endian_sha256=hashlib.sha256(packed).hexdigest(),
                    serialization='np.packbits uint8 boolean vector, bitorder=big; total disambiguates terminal padding')

    def install(self):
        ss=self.ss
        originals=self.originals
        def read_ld(args,log,fh,alleles=False,dropna=True):
            result=originals['_read_ld_sumstats'](args,log,fh,alleles=alleles,dropna=dropna)
            if args.h2:
                M,w,ref,frame,novar=result
                self.h2_identity=self.h2_frame_identity(frame,args,len(ref),fh)
            return result
        def other(args,log,p2,sumstats,ref):
            self.current_p2=p2
            result=originals['_read_other_sumstats'](args,log,p2,sumstats,ref)
            self.merged=None
            return result
        def merge(args,s1,s2,log):
            result=originals['_merge_sumstats_sumstats'](args,s1,s2,log)
            self.merged=result
            return result
        def filt(alleles):
            mask=originals['_filter_alleles'](alleles)
            if self.merged is None:raise RuntimeError('no stock merge context')
            snps=self.merged.SNP.reindex(alleles.index)
            all_hash=hashlib.sha256()
            for snp,a,keep in zip(snps,alleles,mask):
                flip=str(int(ss.FLIP_ALLELES[a])) if keep else 'NA'
                all_hash.update((str(snp)+'\t'+str(a)+'\t'+str(int(keep))+'\t'+flip+'\n').encode('utf-8'))
            self.compatibility=dict(candidates_after_stock_dropna=len(mask),stock_compatible_count=int(mask.sum()),
                                    stock_rejected_count=int((~mask).sum()),all_candidates_ordered_allele_compatibility_sha256=all_hash.hexdigest())
            self.accepted=self.pd.DataFrame({'SNP':snps[mask],'alleles':alleles[mask]})
            self.merged=None
            return mask
        ss._read_ld_sumstats=read_ld
        ss._read_other_sumstats=other
        ss._merge_sumstats_sumstats=merge
        ss._filter_alleles=filt

    def pair_frame_identity(self,frame,args):
        if args.no_check_alleles:raise RuntimeError('frozen protocol requires stock allele checks')
        if self.accepted is None:raise RuntimeError('missing stock allele context')
        if not self.accepted.index.equals(frame.index) or not (self.accepted.SNP.values==frame.SNP.values).all():
            raise RuntimeError('accepted allele context does not match stock final loop')
        prefilter_count=len(frame)
        keep=self.np.ones(prefilter_count,dtype=bool)
        if args.chisq_max is not None:
            keep=self.np.asarray(frame.Z1**2*frame.Z2**2 < args.chisq_max**2)
        final=frame if keep.all() else frame[keep]
        meta=self.accepted if keep.all() else self.accepted[keep]
        h=hashlib.sha256();flips=0
        for snp,a in zip(meta.SNP,meta.alleles):
            flip=int(self.ss.FLIP_ALLELES[a]);flips+=flip
            h.update((str(snp)+'\t'+str(a)+'\t1\t'+str(flip)+'\n').encode('utf-8'))
        masks=None
        if args.two_step is not None:
            m1=self.np.asarray(final.Z1**2 < args.two_step)
            m2=self.np.asarray(final.Z2**2 < args.two_step)
            masks=dict(hsq1=self.mask_info(m1),hsq2=self.mask_info(m2),gencov=self.mask_info(self.np.logical_and(m1,m2)))
        row=dict(p1=args.rg.split(',')[0],p2=self.current_p2,
                 final_ordered_SNP_count=len(final),final_ordered_SNP_sha256=self.snp_digest(final.SNP),
                 final_ordered_SNP_allele_compatibility_sha256=h.hexdigest(),final_allele_flip_count=flips,
                 final_ordered_aligned_Z1_Z2_float64_big_endian_sha256=self.float_digest((final.Z1,final.Z2)),
                 final_ordered_N1_N2_float64_big_endian_sha256=self.float_digest((final.N1,final.N2)),
                 before_optional_final_Z_filter_count=prefilter_count,final_chisq_max=args.chisq_max,
                 final_Z_filter_removed=prefilter_count-len(final),effective_two_step=args.two_step,
                 n_blocks=min(args.n_blocks,len(final)),two_step_masks=masks,
                 stock_allele_compatibility=self.compatibility,
                 identity_serialization='UTF-8 SNP\\tA1+A2+A1x+A2x\\tcompatible_bit\\tstock_flip_bit\\n in exact final order',
                 source='instrumented stock read/merge/dropna/filter/align outputs; final _rg Z filter replayed read-only')
        self.accepted=None;self.compatibility=None
        self.rows.append(row)
        return row

    def h2_frame_identity(self,frame,args,n_annot,path):
        chisq_max=args.chisq_max
        if n_annot!=1 and chisq_max is None:chisq_max=max(.001*frame.N.max(),80)
        keep=self.np.ones(len(frame),dtype=bool) if chisq_max is None else self.np.asarray(frame.Z**2 < chisq_max)
        final=frame if keep.all() else frame[keep]
        effective=args.two_step
        if n_annot==1 and effective is None and args.intercept_h2 is None:effective=30
        return dict(input=path,final_ordered_SNP_count=len(final),final_ordered_SNP_sha256=self.snp_digest(final.SNP),
                    final_ordered_Z_float64_big_endian_sha256=self.float_digest((final.Z,)),
                    final_ordered_N_float64_big_endian_sha256=self.float_digest((final.N,)),
                    final_chisq_max=chisq_max,final_Z_filter_removed=len(frame)-len(final),effective_two_step=effective,
                    n_blocks=min(args.n_blocks,len(frame)),two_step_hsq_mask=self.mask_info(final.Z**2<effective) if effective is not None else None,
                    allele_scope='h2 does not use allele checking; source full-stream allele preservation proof is separate',
                    source='stock _read_ld_sumstats output; exact h2 final Z filter and automatic two-step branch replayed read-only')


def audit_without_fits(ss,np,args,log,capture):
    """Replay only the same pinned read/merge/allele functions, never construct regressions."""
    args=copy.deepcopy(args)
    if args.rg:
        paths,_=ss._parse_rg(args.rg)
        f=lambda x:ss._split_or_none(x,len(paths))
        args.intercept_h2,args.intercept_gencov,args.samp_prev,args.pop_prev=map(f,(args.intercept_h2,args.intercept_gencov,args.samp_prev,args.pop_prev))
        for x in (args.intercept_h2,args.intercept_gencov,args.samp_prev,args.pop_prev):
            if len(x)!=len(paths):raise ValueError('wrong stock argument vector length')
        if args.no_intercept:raise ValueError('not frozen: constrained intercept')
        M,w,ref,frame,_=ss._read_ld_sumstats(args,log,paths[0],alleles=True,dropna=True)
        if M.shape[1]==1 and args.two_step is None and args.intercept_h2 is None:args.two_step=30
        for p2 in paths[1:]:
            loop=ss._read_other_sumstats(args,log,p2,frame,ref)
            capture.pair_frame_identity(loop,args)
        return capture.rows
    if args.h2:
        if args.intercept_h2 is not None:args.intercept_h2=float(args.intercept_h2)
        if args.no_intercept:args.intercept_h2=1
        M,w,ref,frame,_=ss._read_ld_sumstats(args,log,args.h2)
        ss._check_ld_condnum(args,log,ref)
        ss._warn_length(log,frame)
        return [capture.h2_identity]
    raise ValueError('only frozen h2 or rg inputs are supported')
