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
re-exports each figure, table and metric the paper reports from the sub-analysis
that computes it, and computes nothing itself. The inputs are what came from
outside this analysis: the shear catalogs of Paper I, the n(z) and the GLASS suite
of Paper III. Everything the analysis made is an output. The tree under `analyses/`
reads in the order of the argument: the covariances, the data vectors they belong
to, the validation behind both, and the null tests where they meet.

| Sub-analysis | What it does |
|---|---|
| `covariances` | Every covariance: `CosmoCov` ξ±, NaMaster Cℓ, and the pure E/B (M C Mᵀ) and COSEBI (Tᵀ C T) covariances; the reference covariances the mocks are held against; the scan over n(z) blinds; and the footprint power spectrum |
| `data_vectors` | ξ± and Cℓ measured from the catalogs, including the 10,000-bin ξ± of the convergence check, and pure E/B modes and COSEBIs as linear transforms of ξ± |
| `validation` | The evidence that the null tests can be trusted: the estimators against a 10,000-bin grid, COSEBIs from harmonic-space bandpowers against those from ξ±, and the B-mode covariances across n(z) blinds |
| `validation/glass_mocks` | Validation of estimators and covariances on 350 GLASS mocks, and the per-mock statistics |
| `null_tests` | Where data meet covariance: data-vector figures with their error bars, PTE matrices over scale cuts and catalog versions, the mock calibration of the minimum PTE, the PTE tables and heatmaps, the primary null test, and the numbers the paper quotes |

`data_vectors` reads no covariance. `covariances` reads the ξ± bin centres and pair
counts from `data_vectors`, which the pure E/B and COSEBI transforms are built on,
and none of the shear. Data flows one way: from `data_vectors` through
`covariances` and `validation` to `null_tests`.

`reference_covariance` governs the covariances the GLASS mocks are held against.
They are the fiducial recipes under one decision: CosmoCov at the noise an option
sets, and NaMaster on the catalog an option names, realisation 00001 of the suite
under `mockmatched`. They stay apart from the data covariances, which the null
tests read under either selection.

Each scan repeats one measurement across catalog versions, mock realisations or
scale cuts, and gathers the members into a single product. A scan sits in
`analyses/` beneath the sub-analysis it serves, declared in a `*.sweep.yaml` file
from which the sibling `astra.yaml` is generated.

A catalog's n(z) file, area, n_eff and σ_e are facts from Paper I, not steps of
this analysis. They live in the option tables of
`src/workflow/scripts/decision_values.py`, selected by `fiducial_catalog` and
`fiducial_blind`, and recipes read them through that table, as they do the noise
and catalog of each `reference_covariance` option.

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

One output stays outside `lc`. `xi_10000_bins` is a three-node MPI job that has to
run on the host, because the container's OpenMPI cannot start under candide's
Slurm: its recipe submits `producers/run_10k_xi_barehost.sbatch` from the project
root.

The results under `results/` were made by `lc` on candide, last on 4 November
2026, and each manifest records the commit and recipe that made them. Recipes and
scripts have since been rewritten to take decisions as typed flags, without a
cluster rerun. The outputs that can run without the catalogs (the tables, macros
and figures drawn from recorded evidence) reproduce their records: the same
numbers to 7.8e-14 and the same pixels. The measurement scripts have not been run
in this form, so a full rerun is still owed.

Recorded results sit in the present tree with their manifests as written, so a
manifest's output ids and paths are those of the run that made it, and `lc` treats
the outputs as changed. The recorded `pure_eb_data`, `pure_eb_sweep` and
`cosebis_modes_data` still hold their covariances alongside the data vectors, the
members of the recorded `cl_sweep` bundle still carry that run's `pseudo_cl_*`
file names and hold the covariances too, and `pure_eb_covariance`,
`pure_eb_covariance_sweep`, `cosebis_covariance`, `cosebis_covariance_sweep`,
`cosebis_modes_sweep` and `cl_cov_sweep` have no record yet.
`tools/split_recorded_covariances.py` makes the first three pairs from recorded
files; the report's copies were made that way.

`cov_integration_reference`, `cl_cov_reference`, `footprint_power_spectrum`,
`mock_footprint_audit` and `blind_cov_sweep` have no record in the tree either:
their recorded files sit on candide, and the results that read them were made from
those. `footprint_power_spectrum.py` and `mock_footprint_audit.py` have not
reproduced them, and `run_2pcf_highres.py` has not reproduced the recorded
`xi_10000_bins`.

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
uv run --no-project --with pyyaml python tools/check_report_links.py
uv run --no-project --with pyyaml --with numpy python -m unittest discover -s tools -p 'test_*.py'
myst build
uv run --no-project python tools/check_report_build.py
```

The report pins a scoped [article theme](https://github.com/cailmdaley/paper-ii-article-theme)
that renders record links, inherited decision provenance and excluded-option
reasons. Its source is on [a fork of astra-theme](https://github.com/cailmdaley/astra-theme/tree/cail/paper-ii-dialog-polish).
