# UNIONS-3500 Paper II — B-mode validation, as a Lightcone analysis

Daley et al. 2026, *UNIONS-3500 Weak Lensing: II. B-mode validation for cosmic
shear* (A&A, accepted; arXiv [2604.03227](https://arxiv.org/abs/2604.03227)),
specified in ASTRA and executed with `lc`.

- **Read it:** https://myst-viewer.vercel.app/cailmdaley/unions-bmodes-mystra.
  The body follows the accepted paper, with quantitative macros and generated
  figures linked to analysis outputs. Start at **Results** for the PTE table,
  **Catalog comparison** for the three representations, or **Conclusions** for
  the paper's takeaways. **Supporting Documents → ASTRA Inventory** lists the
  root analysis's inputs, decisions, outputs, findings, and prior insights;
  output dependencies lead into the nested analyses. A figure's **Provenance**
  control exposes its decision dependencies.
- **Spec:** `astra.yaml` (root decisions, findings, figure outputs) and
  `analyses/<sub>/astra.yaml` for the nine sub-analyses: catalog, two_point,
  pure_eb, cosebis, cl_bb, covariance, glass_mocks, convergence, null_tests.
- **Code:** `src/` (vendored from sp_validation); recipes on each output name
  the exact command.

## Universes

| Universe | What it is |
|---|---|
| `baseline` | The accepted A&A paper: fixed-quadrature pure E/B operator, analytic covariance validated on GLASS mocks, mock-calibrated minimum-PTE threshold (0.01). |
| `submitted` | The paper as submitted (arXiv v2): adaptive-quadrature pure modes, Monte Carlo covariance, per-statistic 0.05 threshold. |

`universes/<id>.yaml` names each decision's selection and its source of truth.
The decisions that differ between the two universes are the referee round's
methodological changes.
The live report displays the accepted-paper reference comparison; it does not
switch universes interactively.
The `selected_config` output records the raw or leakage-corrected shear
selection, and `primary_null_test_report` records the selected scale range and
COSEBI subsets, recalibrated from stored mock PTEs.
Those selections do not turn the report's complete two-range, six-test tables
and figures into a different paper.

## Presentation checks

The report pins a scoped [article theme](https://github.com/cailmdaley/paper-ii-article-theme)
that renders record links, inherited decision provenance, and excluded-option reasons.
Its source is on [Cail's astra-theme fork](https://github.com/cailmdaley/astra-theme/tree/cail/paper-ii-dialog-polish);
it doesn't require a shared Lightcone theme release.
The `astra-output--no-caption` class suppresses generated captions where the
report supplies its own caption, keeping the numerical value bindings intact.

Run these checks from the project root without installing the scientific stack:

```bash
uv run --no-project --with pyyaml python tools/check_prose_only.py <base-revision>
uv run --no-project --with pyyaml python tools/check_links.py
uv run --no-project --with pyyaml python tools/check_decisions.py
uv run --no-project python tools/check_report_links.py
uv run --no-project python tools/test_report_links.py
myst build
uv run --no-project python tools/check_report_build.py
```

The build check verifies that every authored `astra:value` resolves in the built
report, including values in handwritten captions.

## Reproducing

`lc status` shows what state each output is in. `lc materialize <cluster>`
produces the stale outputs on a cluster allocated with `lc compute`. The
catalogs are UNIONS-internal: they come from `/n17data/UNIONS/WL/` on candide
and are not public.
