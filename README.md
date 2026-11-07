# UNIONS-3500 Paper II: B-mode validation, as a Lightcone analysis

Daley et al. 2026, *UNIONS-3500 Weak Lensing: II. B-mode validation for cosmic
shear* (A&A, accepted; arXiv [2604.03227](https://arxiv.org/abs/2604.03227)),
specified in ASTRA and executed with `lc`.

**Read it:** https://myst-viewer.vercel.app/cailmdaley/unions-bmodes-mystra.
The body follows the accepted paper, with quantitative macros and generated
figures linked to analysis outputs. Start at **Results** for the PTE table,
**Catalog comparison** for the three representations, or **Conclusions** for the
paper's takeaways. **Supporting Documents → ASTRA Inventory** lists the root
analysis's inputs, decisions, outputs, findings and prior insights, and a
figure's **Provenance** control shows the decisions behind it.

## How the analysis is organised

`astra.yaml` holds the shared inputs and decisions, the paper's findings and its
figures and tables. The tree under `analyses/` reads in the order of the argument:
the covariances, the three data vectors they serve, the null tests on those data
vectors, and the validation behind both.

| Sub-analysis | What it does |
|---|---|
| `covariances` | `CosmoCov` ξ± and NaMaster Cℓ covariances |
| `data_vectors` | The three representations of the B-mode signal |
| `data_vectors/real_space` | TreeCorr ξ± on the reporting and integration grids, and the pure E/B decomposition, localised in angular separation |
| `data_vectors/cosebis` | COSEBI modes from ξ± |
| `data_vectors/harmonic_space` | Catalog-based Cℓ^BB and Cℓ^EB |
| `null_tests` | PTE matrices over scale cuts, the primary null test, and the numbers the paper quotes |
| `validation` | The evidence that the null tests can be trusted |
| `validation/glass_mocks` | Validation of estimators and covariances on 350 GLASS mocks, and the PTE threshold |
| `validation/estimator_checks` | The adopted estimators against a 10,000-bin grid, COSEBIs from harmonic-space bandpowers against those from ξ±, and the B-mode covariances across n(z) blinds |

`data_vectors` and `validation` hold no computation of their own: each declares
what its children read and re-exports what the rest of the tree reads. Data flows
one way, from covariances through data vectors and validation to the null tests.

Each scan repeats one measurement across catalog versions, mock realisations or
scale cuts, and gathers the members into a single product. A scan sits in
`analyses/` beside the sub-analysis it serves, declared in a `*.sweep.yaml` file
from which the sibling `astra.yaml` is generated.

A catalog's n(z) file, mask power spectrum, area, n_eff and σ_e are facts from
Paper I, not steps of this analysis. They live in the option tables of
`src/workflow/scripts/decision_values.py`, selected by `fiducial_catalog` and
`fiducial_blind`, and recipes read them through that table.

Code in `src/` is vendored from `sp_validation`; the recipe on each output names
the exact command.

## Universes

| Universe | What it is |
|---|---|
| `baseline` | The accepted A&A paper: fixed-quadrature pure E/B operator, analytic covariance validated on GLASS mocks, mock-calibrated minimum-PTE threshold (0.01). |
| `submitted` | The paper as submitted (arXiv v1; v2 reports the same results): adaptive-quadrature pure modes, Monte Carlo covariance, per-statistic 0.05 threshold. |

The decisions that differ between the two are the referee round's methodological
changes. `universes/<id>.yaml` names each selection, and each sub-analysis's own
`universes/` holds its local ones.

Some decisions act more narrowly than their labels suggest, because the paper
reports both options side by side. `fiducial_scale_cuts` and `mode_subsets`
select what `primary_null_test_report` records; the report's tables and figures
always show both ranges and all six tests. `scale_cut_convention` moves COSEBIs
only, and `pte_threshold` sets the heatmap colour scale, the marked line and the
quoted value. `mc_samples` acts only under the Monte Carlo covariance.

## Reproducing

`lc status` shows what state each output is in, and `lc materialize <cluster>`
produces the stale outputs on a cluster allocated with `lc compute`. The catalogs
are UNIONS-internal, at `/n17data/UNIONS/WL/` on candide, and are not public.

The results under `results/` were made by `lc` on candide, last on 4 November
2026, and each manifest records the commit and recipe that made them. Recipes and
scripts have since been rewritten to take decisions as typed flags, without a
cluster rerun. The 32 outputs that can run without the catalogs (the tables,
macros and figures drawn from recorded evidence) reproduce their records: the
same numbers to 2.5e-14 and the same pixels. The measurement scripts have not been
run in this form, so a full rerun is still owed.

Recorded results sit in the present tree with their manifests as written, so a
manifest's output ids and paths are those of the run that made it, and `lc` treats
the outputs as changed. The members of the recorded `cl_sweep` bundle still carry
that run's `pseudo_cl_*` file names, so its gather reruns before anything reads it.

Code is not a declared input, so editing a script does not mark outputs stale;
each output's manifest records the commit that made it. `CosmoCov` and its GSL
library live outside the project at the paths in
`src/workflow/scripts/run_cosmocov_chain.sh`, and a rerun has to make them visible
to the covariance recipes.

## Checks

Run from the project root, without the scientific stack:

```bash
uvx --from astra-tools astra validate
uv run --no-project --with pyyaml python tools/check_decisions.py
uv run --no-project --with pyyaml python tools/check_links.py
uv run --no-project --with pyyaml python tools/check_reachable.py
uv run --no-project --with pyyaml python tools/check_recipe_flags.py
uv run --no-project python tools/check_report_links.py
uv run --no-project python tools/test_report_links.py
myst build
uv run --no-project python tools/check_report_build.py
```

The report pins a scoped [article theme](https://github.com/cailmdaley/paper-ii-article-theme)
that renders record links, inherited decision provenance and excluded-option
reasons. Its source is on [a fork of astra-theme](https://github.com/cailmdaley/astra-theme/tree/cail/paper-ii-dialog-polish).
