"""DE-SWAN (Differential Expression - Sliding Window ANalysis) in Python.

Two variants:
  * `deswan_wilcoxon` - the "modified DE-SWAN" of Shen et al. 2024 (Nature Aging):
    Mann-Whitney U between the younger and older half of a window, BH-adjusted
    within each window centre.
  * `deswan_linear`  - the original Lehallier et al. 2019 formulation: linear
    model Feature ~ below/above indicator (+ covariates), type-II ANOVA p-value.

Window semantics follow Shen et al.: a window of total width `bucket` centred on
`midpoint`; young = [mid - bucket/2, mid), old = [mid, mid + bucket/2).
Note Shen et al. describe a "20-year window" whose halves are 15 years wide in the
first position; the replication code of Carbonneau et al. uses
`window=10` as the half-width, i.e. bucket=20. Use `half_width` here to be explicit.
"""
from __future__ import annotations
import numpy as np
from scipy.stats import mannwhitneyu
from statsmodels.stats.multitest import multipletests


def _bh(p):
    p = np.asarray(p, float)
    ok = np.isfinite(p)
    q = np.full(p.shape, np.nan)
    if ok.sum():
        q[ok] = multipletests(p[ok], method='fdr_bh')[1]
    return q


def deswan_wilcoxon(expression, ages, midpoints, half_width=10.0, min_n=3):
    """expression: (variables x samples). Returns (p, q) arrays of shape
    (n_variables, n_midpoints) plus per-midpoint group sizes."""
    expression = np.asarray(expression, float)
    ages = np.asarray(ages, float)
    midpoints = np.asarray(midpoints, float)
    P = np.full((expression.shape[0], midpoints.size), np.nan)
    sizes = np.zeros((midpoints.size, 2), int)

    for j, mid in enumerate(midpoints):
        y = (ages >= mid - half_width) & (ages < mid)
        o = (ages >= mid) & (ages < mid + half_width)
        sizes[j] = (y.sum(), o.sum())
        if y.sum() < min_n or o.sum() < min_n:
            continue
        for i in range(expression.shape[0]):
            a = expression[i, y]
            b = expression[i, o]
            a, b = a[np.isfinite(a)], b[np.isfinite(b)]
            if a.size < min_n or b.size < min_n:
                continue
            try:
                P[i, j] = mannwhitneyu(a, b, alternative='two-sided').pvalue
            except ValueError:      # all values identical
                P[i, j] = 1.0
    Q = np.column_stack([_bh(P[:, j]) for j in range(midpoints.size)])
    return P, Q, sizes


def deswan_linear(expression, ages, midpoints, half_width=10.0, covariates=None, min_n=3):
    """Original DE-SWAN test: OLS of feature on a below/above indicator."""
    import statsmodels.api as sm
    expression = np.asarray(expression, float)
    ages = np.asarray(ages, float)
    midpoints = np.asarray(midpoints, float)
    P = np.full((expression.shape[0], midpoints.size), np.nan)
    sizes = np.zeros((midpoints.size, 2), int)

    for j, mid in enumerate(midpoints):
        sel = (ages >= mid - half_width) & (ages < mid + half_width)
        grp = (ages[sel] >= mid).astype(float)
        sizes[j] = ((grp == 0).sum(), (grp == 1).sum())
        if sizes[j].min() < min_n:
            continue
        X = np.column_stack([np.ones(sel.sum()), grp])
        if covariates is not None:
            X = np.column_stack([X, np.asarray(covariates, float)[sel]])
        for i in range(expression.shape[0]):
            yv = expression[i, sel]
            m = np.isfinite(yv)
            if m.sum() < X.shape[1] + 1:
                continue
            res = sm.OLS(yv[m], X[m]).fit()
            P[i, j] = res.pvalues[1]
    Q = np.column_stack([_bh(P[:, j]) for j in range(midpoints.size)])
    return P, Q, sizes


def count_significant(Q, threshold=0.05):
    return np.nansum(Q < threshold, axis=0)
