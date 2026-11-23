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

`astra.yaml` holds the inputs, the shared decisions and the paper's findings, and
re-exports each figure, table and supporting data product the paper reports from
the sub-analysis that computes it, and computes nothing itself. The inputs are
what came from outside this analysis: the shear catalogs of Paper I, the n(z),
and the Paper I survey areas (`data/survey_areas.yaml`). Everything else is an
output, the GLASS mocks included: `mocks` generates each realization, measures
it with the data's TreeCorr and NaMaster settings, and keeps the measurements
and the catalog summary of realization 00001. The tree under `analyses/` lists
the conclusions first: the significance of the B modes, the validation that
licenses it, then the covariances, the data vectors and the mocks it rests on.

| Sub-analysis | What it does |
|---|---|
| `significance` | Where data meet covariance: data-vector figures with their error bars, PTE matrices over scale cuts and catalog versions, the mock calibration of the minimum PTE, the PTE tables and heatmaps, the primary null test, and the numbers the paper quotes |
| `validation` | The evidence that the null tests can be trusted: the estimators against a 10,000-bin grid, COSEBIs from harmonic-space bandpowers against those from ξ±, and the B-mode covariances across n(z) blinds |
| `validation/glass_mocks` | Validation of estimators and covariances on 350 GLASS mocks, and the per-mock statistics |
| `covariances` | Every covariance: `CosmoCov` ξ±, NaMaster Cℓ, and the pure E/B (M C Mᵀ) and COSEBI (Tᵀ C T) covariances; the reference covariances the mocks are held against; and the scan over n(z) blinds |
| `catalog_processing` | Everything read from the catalogs: ξ± and Cℓ, including the 10,000-bin ξ± of the convergence check, and the catalog summary (count maps, footprints, weight and shape-noise sums) the covariances and mocks take; plus pure E/B modes and COSEBIs as linear transforms of ξ± |
| `mocks` | The 350 GLASS realizations of the selected `mock_suite`, generated on the footprint in the catalog summary and measured as ξ± and Cℓ stacks |

Only `catalog_processing` and the GLASS realizations read a catalog.
`catalog_processing` reads no covariance. `covariances` reads the ξ± bin centres and
pair counts from `catalog_processing`, which the pure E/B and COSEBI transforms are
built on, and none of the shear, and the NaMaster covariances take their window and
shot noise from the catalog summary. Data flows one way: from `catalog_processing` and `mocks` through
`covariances` and `validation` to `significance`.

`mock_suite` selects the GLASS suite and `reference_covariance` the covariances
its scatter is held against. A data-matched suite, at the fiducial catalog's
footprint, count and σ_e, is held against the data's own covariances
(`data_covariance`). The published suite carries 17% more galaxies than the data,
so it is held against the fiducial recipes rebuilt at its noise (`mockmatched`):
CosmoCov at that noise, and NaMaster on the summary of its realization 00001. Under
`data_covariance` the reference covariances are copies of the data covariances; the
null tests read the data covariances under either selection.

Each scan repeats one measurement across catalog versions, n(z) blinds or scale
cuts, and gathers the members into a single product. A scan sits in
`analyses/` beneath the sub-analysis it serves, declared in a `*.sweep.yaml` file
from which the sibling `astra.yaml` is generated.

The CosmoCov covariances take each catalog's area, n_eff and σ_e from Paper I, as
facts like the n(z): they live in the option tables of
`src/workflow/scripts/decision_values.py`, selected by `fiducial_catalog`, and
recipes read them through that table, as they do the noise of each
`reference_covariance` option. The catalog summary reproduces Paper I's
weight-effective counts and σ_e; its occupied-pixel areas are 3–5% below Paper I's
mask areas. Every catalog version shares one n(z), the fiducial
catalog's at the selected `fiducial_blind`: the n(z) never enters the measured data
vectors, and the B-mode covariance barely depends on it.

Code in `src/` is vendored from `sp_validation`; the recipe on each output names
the exact command.

## Universes

| Universe | What it is |
|---|---|
| `baseline` | The accepted paper's selection on the data-matched GLASS suite, its scatter held against the data's own covariances. |
| `published` | The accepted A&A paper: fixed-quadrature pure E/B operator, analytic covariance validated on GLASS mocks, mock-calibrated minimum-PTE threshold (0.01). |
| `preprint_v2` | The paper before the referee report (arXiv 2604.03227v2): adaptive-quadrature pure modes, Monte Carlo covariance, per-statistic 0.05 threshold. |

`published` and `preprint_v2` validate on the published GLASS suite against
mock-matched covariances; the decisions that differ between them are the referee
round's methodological changes. `baseline` differs from `published` in `mock_suite`
and `reference_covariance` alone. `universes/<id>.yaml` names each selection, and each sub-analysis's own
`universes/` holds its local ones.

Some decisions act more narrowly than their labels suggest, because the paper
reports both options side by side. `fiducial_scale_cuts` and `mode_subsets`
select what `primary_null_test_report` records; the report's tables and figures
always show both ranges and all six tests. `pte_threshold` sets the heatmap
colour scale, the marked line and the quoted value.

## Reproducing

`lc status` shows the state of every output across `baseline`, `published`, and
`preprint_v2`, and `lc materialize <cluster>` produces stale outputs on a
cluster allocated with `lc compute`.

Products derived directly from a shear catalog — everything under
`results/*/catalog_processing/`, the 350 per-realization GLASS measurements per
universe (`results/*/mocks/glass_suite/realization_*.tar.gz`), and realization
00001's summary — are tracked in plain Git (`.gitattributes`), so a clone from
GitHub can regenerate every covariance, validation check, significance product,
table, and figure without the UNIONS-internal FITS catalogs. The catalogs and
the downstream products stay annexed locally; every downstream output carries a
manifest recording the commit, environment, and recipe that built it.

`CosmoCov` and its GSL library live outside the project at the paths in
`src/workflow/scripts/run_cosmocov_chain.sh`, and covariance recipes run inside
the container image built by `lc`. Code is not a declared input, so editing a
script does not mark outputs stale automatically; each manifest records the
`git_sha` that produced it.

## Checks

Run from the project root, without the scientific stack:

```bash
uvx --from astra-tools astra validate
uv run --no-project --with pyyaml python tools/check_decisions.py
uv run --no-project --with pyyaml python tools/check_links.py
uv run --no-project --with pyyaml python tools/check_reachable.py
uv run --no-project --with pyyaml python tools/check_recipe_flags.py
uv run --no-project --with pyyaml python tools/check_report_links.py
uv run --no-project --with pyyaml --with numpy python -m unittest discover -s tools -p 'test_*.py'
myst build
uv run --no-project python tools/check_report_build.py
```

The report pins a scoped [article theme](https://github.com/cailmdaley/paper-ii-article-theme)
that renders record links, inherited decision provenance and excluded-option
reasons. Its source is on [a fork of astra-theme](https://github.com/cailmdaley/astra-theme/tree/cail/paper-ii-dialog-polish).
