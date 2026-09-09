"""Production-only effect-scale guard; binary LMM coefficients are not log odds."""
from .io import require


def validate_scale(method, cfg):
    if method=='susie':
        require(cfg.get('scale_reviewed') is True,'Fine-mapping phenotype/effect/N scales must be explicitly reviewed')
        require(cfg.get('n_semantics') in {'total','effective'} and cfg.get('sample_size_justification'),
                'Declare whether the reviewed RSS sample size is total or effective')
        kind=cfg.get('trait_type');scale=cfg.get('effect_scale')
        if kind=='cc':require(scale=='log_odds','Case-control beta must be log odds, not BOLT-LMM binary beta')
        elif kind=='quant':
            require(scale in {'quantitative_beta','standardized_beta','binary_lmm_beta'},'Unknown quantitative effect scale')
            require(cfg.get('sdY_source'),'Quantitative sdY needs a source, not an invented unit variance')
            if scale=='binary_lmm_beta':
                require(cfg.get('binary_lmm_approximation_reviewed') is True and cfg.get('binary_lmm_approximation_justification'),
                        'A binary-LMM-as-quantitative approximation requires an explicit scientific justification')
        else:require(False,'Unknown fine-mapping trait type')
    elif method=='mr':
        require(cfg.get('scale_reviewed') is True,'MR exposure/outcome effect scales and units require review')
        allowed={'log_odds','quantitative_beta','standardized_beta','binary_lmm_beta'}
        for side in ['exposure','outcome']:
            require(cfg.get(side+'_effect_scale') in allowed and cfg.get(side+'_effect_unit'),f'Declare MR {side} effect scale and unit')
        if 'binary_lmm_beta' in {cfg['exposure_effect_scale'],cfg['outcome_effect_scale']}:
            require(cfg.get('binary_lmm_interpretation_reviewed') is True,
                    'Binary-LMM MR cannot be automatically interpreted as a log-odds causal effect')
        if cfg.get('report_odds_ratios'):
            require(cfg['outcome_effect_scale']=='log_odds','Exponentiation is valid only for a log-odds outcome coefficient')
