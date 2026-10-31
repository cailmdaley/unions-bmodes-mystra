# lc 0.5.0 migration contract (shared conventions for all unit editors)

Read this whole file before editing anything. The survey inventory (per-unit
plans) is at
`/home/cdaley/loom/.felt/science/unions/shear_2d/bmodes_2d/pure-eb-lc050-migration/survey.json`.
Follow your unit's `plan` entries unless they contradict this contract; the
contract wins.

## The lc 0.5.0 output contract

- One output = ONE FILE, at `results/<universe>/<output_id>.<format>` — lc
  composes the path; the recipe writes to `{output}`, opened exactly as given.
  Never `mkdir -p {output}`, never join a basename onto it.
- Every executable output carries `format:` (extension token, no dot: `txt`,
  `fits`, `npz`, `json`, `png`, `csv`, `tex`, `tar.gz`).
- Recipe placeholders `{inputs.<id>}`, `{decisions.<id>}`, `{output}` only;
  every placeholder must be declared in the output's `inputs:` / `decisions:`.

## Container & environment (the big rewrite)

- lc now OWNS the container. Strip every
  `apptainer exec --bind ... /n17data/cdaley/containers/containers/ /usr/local/bin/python`
  wrapper down to bare `python`. Drop `--env PYTHONPATH=...` — the
  `sp_validation` package is installed in the project environment
  (sha-pinned git dep in pyproject.toml, same sha as the worktree).
- Scripts are NOT importable from the package: they are referenced via a
  declared input. The ROOT astra.yaml declares:

  ```yaml
  - id: spv_scripts
    type: data
    description: >
      Analysis scripts from sp_validation branch repro/paper-ii-astra
      (sha 1804f63e, same as the locked package dep). Declared as an input
      so the sandbox can read them and script changes can cascade.
    source: "/n17data/cdaley/unions/code/sp_validation.worktrees/repro-paper-ii-astra"
  ```

  Each sub-analysis that runs scripts aliases it (`- {id: spv_scripts,
  from: ../spv_scripts}`) and recipes use e.g.
  `python {inputs.spv_scripts}/workflow/scripts/run_2pcf.py ...`.
  (ROOT editor adds the declaration; unit editors add the alias — do not
  assume someone else added your alias.)

## Cross-unit data flow

- An upstream output's file is consumed via `{inputs.<id>}` where the input
  aliases the upstream output (`from: ../two_point.xi_integration` etc.).
  NEVER reference `results/...` by literal path in a recipe.
- Upstream artifacts land at `results/baseline/<output_id>.<format>`; the
  formats are frozen in the survey (`proposed_format`) — treat that as the
  interface. Notable: all two_point xi dumps are `txt`; `xi_sweep` is
  `tar.gz` (see below); catalog outputs: `shear_catalog.fits`,
  `redshift_distribution.txt`, `mask_power_spectrum.txt`,
  `survey_properties.json`.
- `xi_sweep` is a `tar.gz` bundle. Consumers (cosebis version comparison,
  pure_eb sweep) untar it into a temp dir first:
  `d=$(mktemp -d) && tar xzf {inputs.xi_sweep} -C $d && python ... --xi-sweep-dir $d`.
  (`$TMPDIR` is writable in the sandbox; `{{`/`}}` escape literal braces if
  you need them.)

## Side artifacts

- Drop unconsumed side artifacts (xi FITS dumps, `variants/`, duplicate
  PNGs, `.snakemake_timestamp`) — do not carry them into the new contract.
- `evidence.json` files: if the mystra `index.md` (`{astra:value}` refs) or a
  finding's `evidence.artifact` consumes the values, promote it to its own
  declared output (`type: metric`, `format: json`) produced by the same
  script with a flag, or keep the metric file as THE output where the
  directory's primary artifact was the evidence itself. If nothing consumes
  it, drop it. Record what you dropped in your unit's report.
- Where a recipe produced data + figure in one pass, split into two declared
  outputs (figure output re-runs the script with a `--plot-only`-style flag
  or a small plotting recipe consuming the data output). Prefer
  figure-consumes-data over recompute.

## Script CLI convention

Several scripts are invoked from more than one unit (`generate_paper_macros.py`,
`config_space_pte_matrices.py`, `harmonic_space_pte_matrices.py`,
`harmonic_config_cosebis_comparison.py`). ONE flag family per script, and it is
the one the `null_tests` unit established. Reconciliation ruling, binding:

- **Compute once, name the artifact with `--emit`.**
  `--emit <artifact> --out-file <path>` runs the heavy pass and writes exactly
  one file. The `<artifact>` vocabulary is per script, e.g.
  `--emit {evidence,macros,table}` for `generate_paper_macros.py`,
  `--emit {evidence,matrices}` for the two `*_pte_matrices.py`.
  No script-specific `--<thing>-out-file` variants (`--matrices-out-file` is
  retired in favour of `--emit matrices --out-file`).
- **Render from what was computed, never recompute.**
  `--from-evidence {inputs.<evidence_output>} --figure <which> --out-file <path>`
  for anything rendering off an evidence JSON, and the parallel
  `--<x>-from-modes <npz> --out-file <path>` for the NPZ-backed passes of
  `harmonic_config_cosebis_comparison.py` (`--metrics-from-modes`,
  `--figure-from-modes`) and `cosebis_data_vector.py` (`--plot-from-modes`).
  `<which>` is `fiducial` for the sub-analysis panel and `appendix` for the
  root's per-version composite.
- **Consequence for the root unit.** The root no longer re-runs any PTE sweep or
  COSEBI comparison: `pte_table_appendix` gathers with `--emit table`;
  `fig_pte_heatmaps` / `fig_cl_pte_heatmap` render `--from-evidence` off the
  `null_tests` evidence outputs; `fig_harmonic_config_cosebis{,_full}` render
  `--figure-from-modes` off the `cosebis` modes NPZs. Root figures that ARE a
  sub-analysis figure are aliases (`fig_cosebis_fiducial` included), never a
  second recipe against the same script.
- **A shared invocation is argument-identical.** If two units call the same
  script for the same pass, the flags match (`--blind A` is passed everywhere or
  nowhere). Divergence here means one of the two is wrong.

The consolidated, deduplicated list of script edits this implies is in
`work/notes/lc050-script-changes.md` — that file, not the individual unit
reports, is the input to the coordinated script-editing pass.

## What NOT to touch

- `results/` — the July baseline stays in place for regression comparison.
  Do not delete, rename, or write anything there.
- The mystra `index.md` / `myst.yml` — report re-wiring is a later phase.
- `universes/`, decisions, findings, prior_insights — unless your unit's
  recipe rewiring requires a new input declaration.
- Scripts in the sp_validation worktree — flag needed script changes in your
  unit report instead of editing (script edits are a separate coordinated
  pass; many scripts write parameter-encoded filenames into a directory and
  will need an `--out <file>` mode — say exactly which).

## Validation

After editing: `uvx astra-tools@0.2.17 validate` from the project root must
pass (warnings about missing `format:` on OTHER units are fine; your unit
must be warning-free). Do not commit — the orchestrator commits.
