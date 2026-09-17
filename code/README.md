# Code

Three cloned upstream repositories plus two modules written here to remove the R dependency.

Cloned repos are gitignored (`ipop_aging` alone is ~2.9 GB) — re-clone with the commands below.

---

## Cloned

### `artifactual-waves-of-aging/` — the critique's own code

- https://github.com/QuackenbushLab/artifactual-waves-of-aging (Carbonneau, Shutta et al. 2026)
- `git clone https://github.com/QuackenbushLab/artifactual-waves-of-aging`

| path | what it does |
|---|---|
| `r-notebooks/preprint_figure_hierarchical.Rmd` | Artifact 1 — smooth-cluster-interpret on null data |
| `r-notebooks/preprint_null_single_gene_example.Rmd` | Artifact 2 — single-gene LOESS+DE-SWAN inflation |
| `r-notebooks/preprint_four_dist_loess_deswan.Rmd` | Artifact 2 — four age distributions |
| `r-notebooks/preprint_transcriptomic_null_example.Rmd` | Artifact 2 — iPOP age distribution, null expression |
| `ipop_data_analysis/transcriptomics_replication.ipynb` | Reproduces Shen et al. Fig. 4d, then reruns without LOESS |
| `ipop_data_analysis/proteomics_replication.ipynb` | Same for proteomics |
| `ipop_data_analysis/loess.R` | LOESS fitting copied from `ipop_aging` |
| `ipop_data_analysis/permutation_rslts/*.csv` | **Precomputed permutation results** — 80 transcriptomic, plus proteomic. Saves days of compute |
| `simulations/DE-SWAN/{uniform,normal,bimodal,bimodal_2}_dist/` | Artifact 3 — age-distribution effects (Fig. 8) |
| `simulations/DE-SWAN/{middle,ends}_var/`, `{exp_inc,exp_dec}/` | Artifact 3 — heteroskedasticity (Fig. 10) |
| `simulations/DE-SWAN/outliers/` | Artifact 3 — outlier clusters (Fig. S3) |
| `simulations/DE-SWAN/de_SWAN_helpers.py` | Shared simulation helpers |

**Runnable as-is?** The Python simulations, yes — `requirements.txt` matches what is installed
here (numpy/scipy/statsmodels/pandas/matplotlib). Each setting has `XXX_configs.py`,
`run_one_sim.py <SIM_NUM>`, `summarize_de_swan.py`, `plot_sample_exp.py`; run from the repo root.
The R notebooks and the two Jupyter replications need R + `rpy2` and will **not** run here —
`../loess_py/` exists to remove that blocker.

**Note on data**: the replication notebooks expect CSVs in `ipop_data_analysis/nonlinear_aging_data/`
with filenames left blank for the user. `../../datasets/ipop_cross_section/` supplies exactly
these three tables per layer (expression, sample_info, variable_info).

### `ipop_aging/` — the original analysis code (Shen et al. 2024)

- https://github.com/jaspershen-lab/ipop_aging
- `1-code/N_<omic>/` — one folder per layer: `1_data_preparation_cross_sectional.R`,
  `2_loess_fit_cross_section.R`, `3_DEG.R`, `3_DEG_permutation_pvalue.R`,
  `3_fuzzy_c_means_clustering_cross_section_loess.R`, `4_hclust_section_loess.R`
- `1-code/11_combined_omics/` — the combined-omics DE-SWAN that produces the 44/60 crests
- `3-data_analysis/combined_omics/DE_SWAN/` — crest molecule lists (`transcript_crest1.csv`,
  `transcript_crest2.csv`), `changed_molecules2.csv`, and `temp_data*` for bucket widths 15/20/25/30
- `3-data_analysis/combined_omics/data_preparation/cross_section_loess/expression_data.xlsx` —
  12 MB, the **LOESS-interpolated combined matrix that DE-SWAN was actually run on**
- Requires R + `tidymass`/`massdataset`/`Mfuzz`/`DEswan`. See `../../datasets/README.md` for how
  the raw inputs were recovered from git history.

### `DEswan/` — the reference DE-SWAN R package (Lehallier)

- https://github.com/lehallib/DEswan
- `R/DEswan.R` (core), `R/q.DEswan.R` (BH per window), `R/nsignif.DEswan.R` (counts),
  `R/reshape.DEswan.R`; `data/agingplasmaproteome.rda` (de-identified 171 × 1,305 fixture);
  vignette in `docs/articles/DEswan.html`
- Original test: `Feature ~ β₁·below/above + covariates`, type-II ANOVA via `car::Anova` —
  **not** the Wilcoxon variant Shen et al. used. Both are implemented in `../deswan_py/`.
- Read `R/DEswan.R` for the default window/quantile behaviour: Carbonneau et al. §4.1.1 argue the
  quantile-based window centres do **not** equalise group sizes, because window length is fixed.

---

## Written here (no R required)

### `loess_py/` — R-compatible LOESS in Python

`rloess.py` reimplements the 1-D case R's `stats::loess` is used for in this pipeline:
`family="gaussian"`, `degree=2`, tricube weights, exact local fit (equivalent to
`surface="direct"`), plus `optimize_loess_span` — the leave-one-out RMSE span search over
{0.3, 0.4, 0.5, 0.6} that both upstream repos use. `run_loess` reproduces `loess.R`'s
per-variable pipeline onto the half-year grid `seq(26, 75, 0.5)`.

**Validated** against the published R output (`../../datasets/ipop_loess_R_csv/`):

```
$ python code/loess_py/validate_against_R.py plasma_proteomics 25
median |rel err| 1.78e-03   p90 1.45e-02
per-variable corr(python, R): min 0.771636  median 0.998859
```

Fixing the span to R's choice tightens this to median relative error 1.6e-3, max 1.2e-2. So the
two sources of disagreement are (a) CV span tie-breaking, which flips the chosen span for a
minority of variables, and (b) R's default `surface="interpolate"` kd-tree approximation vs. our
exact fit. **Treat the Python fit as faithful but not bit-identical**, and re-run
`validate_against_R.py` on any layer before relying on it there.

### `deswan_py/` — DE-SWAN in Python

`deswan.py`:
- `deswan_wilcoxon(...)` — the "modified DE-SWAN" of Shen et al.: Mann-Whitney U between the
  younger and older half of each window, BH-adjusted **within** each window centre.
- `deswan_linear(...)` — the original Lehallier formulation: OLS on a below/above indicator,
  optional covariates.
- Both return `(p, q, sizes)` where `sizes` is the per-window `(n_young, n_old)` — report these
  next to every curve; unequal and varying group sizes are the mechanism behind Artifact 3.
- `count_significant(Q, threshold)` — the "number of significant molecules" that gets plotted.

Window convention: `half_width` is explicit, because the literature is ambiguous. Shen et al.
describe a "20-year window"; the replication code uses `window=10` as the half-width. So
`half_width=10` ⇒ a 20-year total window. State this in any writeup.

`smoke_test_proteomics.py` — end-to-end check on iPOP proteomics (LOESS ≈ 1 min for 302 proteins):

```
$ python code/deswan_py/smoke_test_proteomics.py
DE-SWAN alone total significant across windows: 1  (max 1)
LOESS+DE-SWAN peaks at ages: [44 45 43 42]  (max 245 of 302)
```

LOESS+DE-SWAN calls 174–245 of 302 proteins significant at FDR 0.05 in *every* window, peaking at
**44** — the published crest. DE-SWAN on the same raw per-subject data calls **one** test
significant across all 25 windows. This reproduces Carbonneau et al.'s Fig. 6 result on an
independent reconstruction of the data, and is the baseline every Phase 2 experiment builds on.
Output: `code/deswan_py/smoke_out/proteomics_smoke.csv`.

**Performance note.** `run_loess` is a pure-Python double loop: ~1 min for 302 variables × 102
subjects, so the 8,556-gene transcriptome will take ~20–30 min single-threaded. Vectorise or
parallelise across variables before running permutation studies (80 permutations × transcriptome
is otherwise ~40 hours).
