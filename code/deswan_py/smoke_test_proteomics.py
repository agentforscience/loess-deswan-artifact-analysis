"""End-to-end smoke test on iPOP proteomics: LOESS+DE-SWAN vs DE-SWAN alone.

Expected (Carbonneau et al. 2026, Fig. 6): LOESS+DE-SWAN produces large counts
with crests; DE-SWAN on the raw per-subject data produces ~none at FDR 0.05.
"""
import numpy as np, pandas as pd, sys, os
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE); sys.path.insert(0, os.path.join(HERE, '..', 'loess_py'))
from rloess import run_loess
from deswan import deswan_wilcoxon, count_significant

D = os.path.join(HERE, '..', '..', 'datasets', 'ipop_cross_section')
expr = pd.read_csv(f'{D}/plasma_proteomics__expression.csv', index_col=0)
meta = pd.read_csv(f'{D}/plasma_proteomics__sample_info.csv')
expr = expr[meta.subject_id.astype(str).tolist()]
age = meta.adjusted_age.values
mids = np.arange(40, 65)

_, Qr, sizes = deswan_wilcoxon(expr.values, age, mids, half_width=10)
raw = count_significant(Qr)

L, grid, spans = run_loess(expr.values, age, grid=np.arange(26, 75.5, 0.5))
_, Ql, _ = deswan_wilcoxon(L, grid, mids, half_width=10)
loe = count_significant(Ql)

out = pd.DataFrame({'midpoint': mids, 'n_young': sizes[:, 0], 'n_old': sizes[:, 1],
                    'n_sig_deswan_only': raw, 'n_sig_loess_deswan': loe})
print(out.to_string(index=False))
print(f'\nDE-SWAN alone total significant across windows: {raw.sum()}  (max {raw.max()})')
print(f'LOESS+DE-SWAN peaks at ages: {mids[np.argsort(loe)[::-1][:4]]}  (max {loe.max()} of {expr.shape[0]})')
os.makedirs(os.path.join(HERE, 'smoke_out'), exist_ok=True)
out.to_csv(os.path.join(HERE, 'smoke_out', 'proteomics_smoke.csv'), index=False)
