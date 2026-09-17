"""Check the Python LOESS against the R LOESS output published with Shen et al. 2024."""
import numpy as np, pandas as pd, sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from rloess import run_loess

ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', 'datasets')
LAYER = sys.argv[1] if len(sys.argv) > 1 else 'plasma_proteomics'
N = int(sys.argv[2]) if len(sys.argv) > 2 else 25

expr = pd.read_csv(f'{ROOT}/ipop_cross_section/{LAYER}__expression.csv', index_col=0)
meta = pd.read_csv(f'{ROOT}/ipop_cross_section/{LAYER}__sample_info.csv')
ref = pd.read_csv(f'{ROOT}/ipop_loess_R_csv/{LAYER}_object_cross_section_loess__expression.csv', index_col=0)

expr = expr[meta.subject_id.astype(str).tolist()]
common = [v for v in ref.index if v in expr.index][:N]
grid = ref.columns.astype(float).values

py, _, spans = run_loess(expr.loc[common].values, meta.adjusted_age.values, grid=grid)
R = ref.loc[common].values

good = np.isfinite(py) & np.isfinite(R)
rel = np.abs(py - R)[good] / (np.abs(R)[good] + 1e-9)
cors = [np.corrcoef(py[i][good[i]], R[i][good[i]])[0, 1] for i in range(len(common))]
print(f'{LAYER}: {len(common)} variables x {len(grid)} grid points')
print(f'  median |rel err| {np.median(rel):.2e}   p90 {np.quantile(rel, .9):.2e}')
print(f'  per-variable corr(python, R): min {np.min(cors):.6f}  median {np.median(cors):.6f}')
print(f'  spans chosen: {dict(zip(*np.unique(spans, return_counts=True)))}')
