# Changing a decision on this laptop

This branch carries the measured ξ± and their covariances for every catalog variant
(`data/products/`), so decisions can be changed and figures remade here. A request such as
"switch the catalog to the initial one" is a new *universe*, one selection of every decision;
never edit `baseline`, `published` or `preprint_v2`. Answer it in these four commands, without
surveying the repository first:

```bash
python3 tools/new_universe.py initial fiducial_catalog=initial      # name, then decision=option ...
git add -A && git commit -qm "Add the initial-catalog universe"     # lc needs a clean tree
lc compute status local >/dev/null 2>&1 || lc compute launch --wait
lc materialize local initial/significance.pure_eb_data_vector_figure
```

Then read both figures, `results/baseline/significance/pure_eb_data_vector_figure.png` and
`results/initial/…`, and say in two or three sentences what changed, quoting the B-mode PTEs the
run prints. Only the outputs on the figure's path are made; the run takes about 35 seconds.

| Decision | Options that run here |
|---|---|
| `fiducial_catalog` | `size_cut` (the paper's), `initial`, `masked`, `relaxed_flags` |
| `leakage_correction` | `applied` (the paper's), `uncorrected` |

Several can be changed at once (`name a=x b=y`).

A change to a plot rather than a decision (colors, an axis range, a marker) is an edit to the
script its recipe names; commit it. Code is not an input, so delete the figure and its manifest
(`results/<universe>/<analysis>/<id>.png` and `.<id>.manifest.json`), commit that too, and
materialize it again.

---

# UNIONS-3500 Paper II as an ASTRA analysis

Daley et al. 2026, *UNIONS-3500 Weak Lensing: II. B-mode validation for cosmic
shear*, accepted by A&A. The paper shows which catalog and scale cuts pass
pure-mode, COSEBI and harmonic-space B-mode null tests. Cosmological inference
belongs to Papers III and IV.

## Layout

- `astra.yaml`: root decisions, findings, prior insights and re-exports of the
  report's outputs.
- `analyses/**/astra.yaml`: the sub-analyses, five at the top, flat apart from
  `validation/glass_mocks`, `mocks/generation`, `mocks/glass_suite` and the scans,
  which sit in `analyses/<sub>/analyses/`. README explains the split: every
  covariance lives in `covariances`, `catalog_processing` is the only
  sub-analysis that reads the catalogs and reads no covariance, `mocks`
  generates and measures the GLASS realizations, and data meet covariance in
  `significance`.
- What is stored: products derived directly from a catalog (everything under
  `results/*/catalog_processing/`, the GLASS members and realization 00001's
  summary) are plain git, so a clone reruns the rest without the catalogs;
  everything else is annexed on fir only and regenerates (`.gitattributes`).
- `<scan>.sweep.yaml`: scan declarations, in each scan's own directory. The
  sibling `astra.yaml` is generated from each; edit the declaration and regenerate
  with `sweep_expand.py` from the `astra-sweep` practice.
- `universes/`: `baseline` is the accepted selection on the data-matched
  GLASS suite, held against the data's own covariances; `published` is the
  accepted paper and `preprint_v2` the paper before the referee report, both
  validated on the published GLASS suite. The report renders the
  alphabetically first universe (MySTRA 0.0.8), so `baseline` must sort first. Each sub-analysis has the same three under its own `universes/`,
  and a parent binds its children by name.
- `index.md`: the report, which is the accepted paper's text with values and
  figures bound to outputs.
- `src/`: analysis code, vendored from `sp_validation`.
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
- A rationale argues from principle and quotes no configuration or results:
  numbers belong in option labels and descriptions, results in findings and
  output descriptions.
- A decision with several non-excluded options must reach the code it governs.
  A decision with one is a documented fixed choice. Excluded options carry a
  reason a scientist would accept.
- `src/workflow/scripts/decision_values.py` maps options to the values scripts
  compute with, including a catalog's n(z) file and survey properties and the
  noise of each mock-matched reference covariance, which recipes read through it
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
