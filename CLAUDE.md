# UNIONS-3500 Paper II as an ASTRA analysis

Daley et al. 2026, *UNIONS-3500 Weak Lensing: II. B-mode validation for cosmic
shear*, accepted by A&A. The paper shows which catalog and scale cuts pass
pure-mode, COSEBI and harmonic-space B-mode null tests. Cosmological inference
belongs to Papers III and IV.

## Layout

- `astra.yaml`: root decisions, findings, prior insights and re-exports of the
  report's outputs.
- `analyses/**/astra.yaml`: the sub-analyses, four at the top, flat apart from
  `validation/glass_mocks` and the scans, which sit in `analyses/<sub>/analyses/`.
  README explains the split: every covariance lives in `covariances`,
  `data_vectors` reads no covariance, and data meet covariance in `null_tests`.
- `<scan>.sweep.yaml`: scan declarations, in each scan's own directory. The
  sibling `astra.yaml` is generated from each; edit the declaration and regenerate
  with `sweep_expand.py` from the `astra-sweep` practice.
- `universes/`: `baseline` is the accepted paper, `submitted` the paper as
  submitted. Each sub-analysis has a matching pair under its own `universes/`,
  and a parent binds its children by name.
- `index.md`: the report, which is the accepted paper's text with values and
  figures bound to outputs.
- `src/`: analysis code, vendored from `sp_validation`.
- `producers/`: the host MPI job behind `xi_10000_bins`, which `lc` cannot run.
- `tools/`: repository checks and the mirror sync.

## Conventions

- `index.md` is the authority for numbers, scope and wording. Preserve the
  paper's hedges.
- The root re-exports and computes nothing.
- The spec describes the science. Inputs are data, never code. Prose says what is
  true now, with no history or tooling talk.
- Inputs are what came from outside this analysis, and sit at the root; anything
  it made is an output with a recipe.
- Each decision reads as a timeless account of a scientific choice: the
  rationale says what is at stake and why the default is preferred, each option
  says what it does on its own terms, and a label or option id names the
  method, never its status. Version talk
  belongs in `universes/` alone.
- A decision with several non-excluded options must reach the code it governs.
  A decision with one is a documented fixed choice. Excluded options carry a
  reason a scientist would accept.
- `src/workflow/scripts/decision_values.py` maps options to the values scripts
  compute with, including a catalog's n(z) file and survey properties and the
  noise and catalog of each reference covariance, which recipes read through it
  and never as spec inputs; `numerics.py` holds the shared null-test numerics.

## Checks

```bash
uvx --from astra-tools astra validate
uv run --no-project --with pyyaml python tools/check_decisions.py
uv run --no-project --with pyyaml python tools/check_links.py
uv run --no-project --with pyyaml python tools/check_reachable.py
uv run --no-project --with pyyaml python tools/check_recipe_flags.py
uv run --no-project --with pyyaml python tools/check_report_links.py
uv run --no-project --with pyyaml --with numpy python -m unittest discover -s tools -p 'test_*.py'
myst build && uv run --no-project python tools/check_report_build.py
```

`tools/check_prose_only.py <rev>` proves a spec diff touches only prose fields.

## Report deployment

The report builds from the separate repository `cailmdaley/unions-bmodes-mystra`
at `https://myst-viewer.vercel.app/cailmdaley/unions-bmodes-mystra`. Changes here
deploy once the mirror is staged and pushed: run `tools/sync_mirror.sh
<mirror-clone>` from a checkout with the materialized results, review it, run
`myst build` and `HOST=127.0.0.1 myst build --html`, then commit and push the
mirror. Report figures are 300-dpi PNGs, because browsers cannot reliably show
PDFs and the viewer has no converter.
