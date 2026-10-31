# UNIONS-3500 Paper II — B-mode validation, as a Lightcone analysis

Daley et al. 2026, *UNIONS-3500 Weak Lensing: II. B-mode validation for cosmic
shear* (A&A, accepted; arXiv [2604.03227](https://arxiv.org/abs/2604.03227)),
specified in ASTRA and executed with `lc`.

- **Read it:** https://myst-viewer.vercel.app/cailmdaley/unions-bmodes-mystra.
  The body is the accepted paper's text. Every inline number and figure
  resolves from this analysis's outputs. *Supporting documents → ASTRA
  Inventory*, at the bottom, lists every input, decision, output, finding and
  prior insight.
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

## Reproducing

`lc status` shows what state each output is in. `lc materialize <cluster>`
produces the stale outputs on a cluster allocated with `lc compute`. The
catalogs are UNIONS-internal: they come from `/n17data/UNIONS/WL/` on candide
and are not public.
