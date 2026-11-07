# UNIONS-3500 Paper II as an ASTRA analysis

Daley et al. 2026, *UNIONS-3500 Weak Lensing: II. B-mode validation for cosmic
shear*, accepted by A&A. The paper shows which catalog and scale cuts pass
pure-mode, COSEBI and harmonic-space B-mode null tests. Cosmological inference
belongs to Papers III and IV.

## Layout

- `astra.yaml`: root decisions, findings, prior insights and re-exports of the
  report's outputs.
- `analyses/**/astra.yaml`: the sub-analyses, four at the top with `data_vectors`
  and `validation` holding their children in `analyses/<sub>/analyses/`. README
  explains the split.
- `<scan>.sweep.yaml`: scan declarations, in each scan's own directory. The
  sibling `astra.yaml` is generated from each; edit the declaration and regenerate
  with `sweep_expand.py` from the `astra-sweep` practice.
- `universes/`: `baseline` is the accepted paper, `submitted` the paper as
  submitted. Each sub-analysis has a matching pair under its own `universes/`,
  and a parent binds its children by name.
- `index.md`: the report, which is the accepted paper's text with values and
  figures bound to outputs.
- `src/`: analysis code, vendored from `sp_validation`.
- `tools/`: repository checks and the mirror sync.

## Conventions

- `index.md` is the authority for numbers, scope and wording. Preserve the
  paper's hedges.
- The root re-exports and computes nothing.
- The spec describes the science. Inputs are data and analysis products, never
  code. Prose says what is true now, with no history or tooling talk.
- Each decision reads as a timeless account of a scientific choice: the
  rationale says what is at stake and why the default is preferred, each option
  says what it does on its own terms, and a label or option id names the
  method, never its status. Version talk
  belongs in `universes/` alone.
- A decision with several non-excluded options must reach the code it governs.
  A decision with one is a documented fixed choice. Excluded options carry a
  reason a scientist would accept.
- `src/workflow/scripts/decision_values.py` maps options to the values scripts
  compute with, including a catalog's n(z) file, mask spectrum and survey
  properties, which recipes read through it and never as spec inputs;
  `numerics.py` holds the shared null-test numerics.

## Checks

```bash
uvx --from astra-tools astra validate
uv run --no-project --with pyyaml python tools/check_decisions.py
uv run --no-project --with pyyaml python tools/check_links.py
uv run --no-project --with pyyaml python tools/check_reachable.py
uv run --no-project --with pyyaml python tools/check_recipe_flags.py
uv run --no-project python tools/check_report_links.py
uv run --no-project python tools/test_report_links.py
uv run --no-project --with pyyaml --with numpy python tools/test_decision_values.py
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
