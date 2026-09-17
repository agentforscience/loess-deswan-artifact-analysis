"""Pure-Python reimplementation of R's stats::loess for the 1-D case used by the
iPOP aging pipeline (family="gaussian", degree=2, tricube weights, surface="direct"),
plus the leave-one-out span selection from jaspershen-lab/ipop_aging.

R is not required. Validated against the published R LOESS output shipped in the
ipop_aging repo -- see validate_against_R.py.
"""
from __future__ import annotations
import numpy as np


def _tricube(u):
    w = np.clip(1.0 - np.abs(u) ** 3, 0.0, None)
    return w ** 3


def loess_fit(x, y, xout, span=0.75, degree=2):
    """Local polynomial regression, R `loess(..., surface="direct")` semantics.

    span > 1 inflates the neighbourhood bandwidth by span**(1/degree_of_freedom)
    exactly as R does (bandwidth = max distance * span**(1/1) in 1-D).
    """
    x = np.asarray(x, float)
    y = np.asarray(y, float)
    ok = np.isfinite(x) & np.isfinite(y)
    x, y = x[ok], y[ok]
    order = np.argsort(x)
    x, y = x[order], y[order]
    n = x.size
    xout = np.atleast_1d(np.asarray(xout, float))

    q = int(np.floor(n * span)) if span <= 1 else n
    q = max(q, degree + 1)

    out = np.empty(xout.size)
    for i, x0 in enumerate(xout):
        d = np.abs(x - x0)
        idx = np.argpartition(d, q - 1)[:q] if q < n else np.arange(n)
        dq = d[idx]
        h = dq.max()
        if span > 1:                       # R inflates the bandwidth beyond span=1
            h *= span
        w = _tricube(dq / h) if h > 0 else np.ones_like(dq)
        keep = w > 0
        if keep.sum() < degree + 1:        # degenerate neighbourhood
            out[i] = np.nan
            continue
        xi, yi, wi = x[idx][keep], y[idx][keep], w[keep]
        X = np.vander(xi - x0, degree + 1, increasing=True)
        sw = np.sqrt(wi)
        beta, *_ = np.linalg.lstsq(X * sw[:, None], yi * sw, rcond=None)
        out[i] = beta[0]                   # value at x0 (centred design)
    return out


def optimize_loess_span(x, y, span_range=(0.3, 0.4, 0.5, 0.6), degree=2):
    """Leave-one-out RMSE span selection, mirroring `optimize_loess_span` in
    ipop_aging/artifactual-waves-of-aging (interior points only: indices 2..n-1 in R)."""
    x = np.asarray(x, float)
    y = np.asarray(y, float)
    order = np.argsort(x)
    x, y = x[order], y[order]
    n = x.size
    best, best_rmse = span_range[0], np.inf
    rmses = {}
    for span in span_range:
        preds, reals = [], []
        for idx in range(1, n - 1):        # R's 2:(nrow-1), 1-based -> interior
            m = np.ones(n, bool)
            m[idx] = False
            try:
                p = loess_fit(x[m], y[m], x[idx], span=span, degree=degree)[0]
            except Exception:
                p = np.nan
            if np.isfinite(p):
                preds.append(p)
                reals.append(y[idx])
        rmse = np.sqrt(np.mean((np.array(reals) - np.array(preds)) ** 2)) if preds else np.nan
        rmses[span] = rmse
        if np.isfinite(rmse) and rmse < best_rmse:
            best, best_rmse = span, rmse
    return best, rmses


def run_loess(expression, ages, grid=None, span_range=(0.3, 0.4, 0.5, 0.6),
              suppress_negative=False, degree=2):
    """LOESS-interpolate a (variables x samples) matrix onto `grid`.

    Mirrors ipop_aging/1-code/.../loess.R: per-variable CV span selection then
    prediction on a half-year grid. Returns (variables x len(grid)) ndarray.
    """
    expression = np.asarray(expression, float)
    ages = np.asarray(ages, float)
    if grid is None:
        grid = np.arange(26, 75.5, 0.5)
    out = np.empty((expression.shape[0], len(grid)))
    spans = np.empty(expression.shape[0])
    for i in range(expression.shape[0]):
        span, _ = optimize_loess_span(ages, expression[i], span_range, degree)
        spans[i] = span
        out[i] = loess_fit(ages, expression[i], grid, span=span, degree=degree)
    if suppress_negative:
        out[out < 0] = 0.0
    return out, np.asarray(grid, float), spans
