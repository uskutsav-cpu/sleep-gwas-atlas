"""Small statistical utilities, never substitutes for named native methods."""
from __future__ import annotations
import math
from typing import Iterable
import numpy as np
from .io import require


def bh(pvalues: Iterable[float | None], family_size: int | None = None) -> list[float | None]:
    """BH with explicit denominator; invalid tests remain NA, not null results.

    Missing tests conservatively occupy family slots at p=1 for accounting.
    q-values alone do not make an analysis with excessive failures acceptable.
    """
    p = list(pvalues)
    m = len(p) if family_size is None else family_size
    require(m >= len(p) and m >= 1, "BH denominator cannot be smaller than supplied family")
    valid = []
    for i, x in enumerate(p):
        if x is None or (isinstance(x, (float, np.floating)) and math.isnan(x)):
            continue
        require(math.isfinite(x) and 0 <= x <= 1, f"Invalid p-value {x}")
        valid.append((float(x), i))
    valid.sort()
    out = [None] * len(p)
    q = 1.0
    for rank in range(len(valid), 0, -1):
        value, i = valid[rank - 1]
        q = min(q, value * m / rank)
        out[i] = q
    return out


def two_sided_normal(z: float) -> float:
    require(math.isfinite(z), "Nonfinite z-score")
    return math.erfc(abs(z) / math.sqrt(2))


def effective_n(cases: float, controls: float) -> float:
    require(cases > 0 and controls > 0 and math.isfinite(cases + controls), "Invalid case counts")
    return 4.0 / (1.0 / cases + 1.0 / controls)


def compare_replication(discovery_rg: float, replication_rg: float, replication_se: float,
                        replication_p: float, *, alpha: float, h2_z: float,
                        independent: bool, phenotype_match: str) -> str:
    if not independent or phenotype_match not in {"exact", "comparable"}:
        return "NO_VALID_REPLICATION"
    require(math.isfinite(discovery_rg) and abs(discovery_rg)<=1 and
            math.isfinite(h2_z) and math.isfinite(alpha) and 0<alpha<1,
            "Invalid discovery correlation, power statistic, or replication alpha")
    require(math.isfinite(replication_rg) and abs(replication_rg) <= 1,
            "Nonphysical replication correlation")
    require(replication_se > 0 and math.isfinite(replication_se), "Invalid replication SE")
    require(0 <= replication_p <= 1 and math.isfinite(replication_p), "Invalid replication P")
    if h2_z < 4:
        return "UNDERPOWERED"
    if discovery_rg * replication_rg <= 0:
        return "DISCORDANT" if replication_p <= alpha else "INCONCLUSIVE_OPPOSITE_SIGN"
    return "REPLICATED" if replication_p <= alpha else "DIRECTIONALLY_SUPPORTED"


def ivw(beta_x, se_x, beta_y, se_y) -> dict:
    """First-order, no-intercept, multiplicative random-effects IVW.

    Diagnostic utility for independent, aligned instruments only. Does not
    validate IV assumptions, account for overlap, or establish causality.
    Native TwoSampleMR remains the production MR adapter.
    """
    from scipy.stats import chi2
    bx, sx, by, sy = [np.asarray(a, dtype=float) for a in (beta_x, se_x, beta_y, se_y)]
    require(bx.ndim == 1 and bx.size >= 2 and all(a.shape == bx.shape for a in [sx, by, sy]),
            "Need at least two aligned instruments")
    require(all(np.isfinite(a).all() for a in [bx, sx, by, sy]), "Nonfinite MR inputs")
    require((sx > 0).all() and (sy > 0).all() and (bx != 0).all(), "Invalid MR effects/errors")
    w = 1 / sy**2
    denom = np.sum(w * bx**2)
    estimate = float(np.sum(w * bx * by) / denom)
    Q = float(np.sum(w * (by - estimate * bx)**2))
    df = bx.size - 1
    se = float(math.sqrt(max(1.0, Q / df) / denom))
    return {"method": "IVW_DIAGNOSTIC_NOT_CAUSAL_PROOF", "n_instruments": int(bx.size),
            "beta": estimate, "se": se, "p": two_sided_normal(estimate / se),
            "Q": Q, "Q_df": int(df), "Q_p": float(chi2.sf(Q, df)),
            "min_F": float(np.min((bx / sx)**2))}


def stratified_enrichment(scores, selected, strata, *, permutations=10000, seed=20260908) -> dict:
    """One-sided competitive permutation on independent loci / matched units.

    scores: one score per independent unit; selected: fixed membership.
    Labels permuted within predefined strata, preserving the selected count.
    Calling SNPs independent is NOT allowed: callers must collapse LD first.
    """
    scores = np.asarray(scores, dtype=float)
    selected = np.asarray(selected, dtype=bool)
    strata = np.asarray(strata, dtype=str)
    require(scores.ndim == 1 and scores.shape == selected.shape == strata.shape, "Shape mismatch")
    require(np.isfinite(scores).all() and selected.any() and (~selected).any(), "Invalid enrichment data")
    require(permutations >= 99, "At least 99 permutations required")
    groups = [np.where(strata == s)[0] for s in np.unique(strata)]
    require(all(0 < selected[g].sum() < len(g) for g in groups),
            "Each matched stratum needs selected and control units")
    obs = float(scores[selected].mean() - scores[~selected].mean())
    rng = np.random.default_rng(seed)
    exceed = 0
    for _ in range(permutations):
        labels = selected.copy()
        for g in groups:
            labels[g] = rng.permutation(labels[g])
        value = float(scores[labels].mean() - scores[~labels].mean())
        exceed += value >= obs - 1e-14
    return {"effect": obs, "p": (1 + exceed) / (permutations + 1),
            "permutations": permutations, "seed": seed, "n_units": len(scores)}
