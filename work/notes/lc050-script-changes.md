# lc 0.5.0 — consolidated script changes

Merged and deduplicated from the seven unit-editor reports of the lc 0.5.0 spec
migration, reconciled against the **Script CLI convention** section of
`lc050-migration-contract.md` (which is binding where a unit report disagrees).

Scope: `/n17data/cdaley/unions/code/sp_validation.worktrees/repro-paper-ii-astra`.
No script has been edited yet — the spec is written against these interfaces, so
the spec is the requirement and this file is the work order.

Two changes recur everywhere and are stated once here rather than repeated:

- **`--out <dir>` → `--out-file <path>`.** Every script must be able to write
  exactly one file to a caller-given path, creating no directory and no
  siblings (no `evidence.json`, no `variants/`, no parameter-encoded filename,
  no `.pdf` twin). Where a directory mode is still needed by a sweep driver,
  keep it alongside — noted per script below.
- **Figures are PNG at 300 dpi**, not PDF, wherever a figure is the output file.

---

## Shared scripts (multi-unit — settle these first)

### `papers/bmodes/scripts/config_space_pte_matrices.py`
Consumers: `null_tests.config_space_pte_data`, `null_tests.config_space_pte_evidence`,
root `fig_pte_heatmaps`.

- Add `--emit {evidence,figure}` + `--out-file <path>`, replacing the `--out
  <dir>` multi-artifact write (two PDFs + `evidence.json`).
- Add `--from-evidence <json> --figure {fiducial,appendix} --out-file <png>`:
  a pure render that recomputes nothing. **The evidence JSON must therefore
  carry the full per-version (θ_min, θ_max) PTE matrices**, not just the
  adopted-cut scalars, or the appendix composite cannot be rendered from it.
  This is the load-bearing requirement — the root's 4×3 appendix composite now
  renders from the same payload as the fiducial 1×3 panel.
- `--pte-intermediate-dir` / `--cosebis-pte-dir` must tolerate an untarred
  scratch dir (it globs, so likely a no-op — confirm).

### `papers/bmodes/scripts/harmonic_space_pte_matrices.py`
Consumers: `null_tests.harmonic_space_pte_data`, `null_tests.harmonic_space_pte_evidence`,
`cl_bb.cl_pte_per_cut`, root `fig_cl_pte_heatmap`.

- Add `--emit {evidence,matrices,figure}` + `--out-file <path>`.
  `matrices` writes the per-version Cℓ^BB PTE grid NPZ (today
  `cl_pte_matrices_A.npz` inside an output dir). **`--matrices-out-file` from
  the cl_bb draft is retired** — one flag family per the contract.
- Add `--from-evidence <json> --figure {fiducial,appendix} --out-file <png>`.
  Evidence payload must carry the full (ℓ_min, ℓ_max) PTE grid per version.
- `--results-dir` must accept a flat untarred `cl_sweep` layout.

### `papers/bmodes/scripts/generate_paper_macros.py`
Consumers: `null_tests.pte_summary_evidence`, `null_tests.pte_summary_macros`,
root `pte_table_appendix`.

- Add `--emit {evidence,macros,table}` + `--out-file <path>`; today `--out
  <dir>` writes `claims_macros.tex` + `evidence.json` + three `.tex` + the CSV
  with filenames composed inside.
  - `evidence` → aggregated PTE JSON; `macros` → `claims_macros.tex`;
    `table` → the N×6 PTE CSV (the root's table output).
- Drop the appendix/results `.tex` side writes (`pte_table_appendix.tex` is
  superseded by the CSV, `pte_table_results.tex` is dead per its own comment).
- Make `--claims-dir` optional — no recipe passes it; any non-table
  `claims_macros` entries must come from the two evidence payloads or be
  dropped.

### `papers/bmodes/scripts/harmonic_config_cosebis_comparison.py`
Consumers: `cosebis.cosebis_harmonic_modes{,_full}`,
`cosebis.cosebis_harmonic_metrics{,_full}`, root
`fig_harmonic_config_cosebis{,_full}`.

- `--out-file <path.npz>` writing only the En/Bn/covariance NPZ (both paths,
  both cuts) for the compute pass.
- `--metrics-from-modes <npz> --out-file <path.json>`: emit the
  `harmonic_b_mode_ptes` / `config_b_mode_ptes` PTE+χ² dict (currently buried
  in `evidence.json`).
- `--figure-from-modes <npz> --out-file <path.png>`: render the paper Figure 8
  composite for the requested `--angular-range` from the same NPZ. **New**, and
  required — the root figures no longer recompute.
- Make `--cosmo-val-dir` / `--covariance-dir` optional now that the fiducial
  pseudo-Cℓ / ξ / covariance arrive as explicit paths; the two render modes pass
  neither.
- Drop the per-version PNG / evidence.json side writes.

### `papers/bmodes/scripts/cosebis_data_vector.py`
Consumers: `cosebis.cosebis_modes_data`, `cosebis.fig_cosebis_data_vector`
(the root's `fig_cosebis_fiducial` is an alias of the latter, not a second call).

- `--out-file <path.npz>` writing the modes NPZ to exactly that path (today
  `cosebis_modes_<ver>.npz` plus figure plus `evidence.json` inside `--out`).
- `--plot-from-modes <npz> --out-file <png>` plot-only mode. **The root's
  proposed `--plot-only --cosebis-modes-npz` spelling is dropped** in favour of
  this one.

---

## Single-consumer scripts

### two_point
- `workflow/scripts/run_2pcf.py` — add `--out-file <path>` writing the TreeCorr
  text dump to exactly that path, no directory, no
  `{ver}_xi_minsep=..._nbins=..._npatch=...txt` name composition, no FITS /
  patches / leakage siblings. `--no-fits` already exists. The leakage-stats
  side effect (a `leakage_<ver>/` subdir relative to the output dir) must be
  suppressed or redirected to scratch under `--out-file`.
- `papers/bmodes/scripts/run_xi_sweep.py` — no change required; `--out <dir>`
  semantics preserved (the recipe passes a mktemp dir and tars it). An optional
  `--no-fits` passthrough would shrink the tarball; not needed for correctness.

### pure_eb
- `papers/bmodes/scripts/run_pure_eb_semianalytic.sh` — add `--out-file
  <path>`: run the 20 chunks in a `mktemp -d` scratch (currently `$OUT/chunks`)
  and have the gather step write the single gathered NPZ directly to that path.
  Keep `--out` for the sweep driver, which still needs dir mode.
- `papers/bmodes/scripts/pure_eb_data_vector.py` — `--out-file <path>`: one
  300-dpi PNG, no `pure_eb_data_vector.pdf` / `evidence.json`.
- `papers/bmodes/scripts/pure_eb_covariance.py` — same.
- `papers/bmodes/scripts/pure_eb_version_comparison.py` — same;
  `--results-dir` semantics unchanged (now receives an untarred temp dir).
- `papers/bmodes/scripts/calculate_pure_eb_ptes.py` — add `--out-file <path>`;
  keep the `--out` dir mode, on which `run_pure_eb_ptes_sweep.sh` relies for the
  canonical per-version filenames.
- `run_pure_eb_sweep.sh`, `run_pure_eb_ptes_sweep.sh` — no change (recipes tar
  the dir output).

### cosebis
- `papers/bmodes/scripts/compute_cosebis_pte_single.py` — `--out-file <path>`
  writing the single gathered NPZ (today `cosebis_ptes_<ver>_<blind>.npz`
  inside `--out`).
- `papers/bmodes/scripts/run_cosebis_ptes_sweep.py` — keep the directory
  `--out` (a genuine 7-file sweep), but guarantee it writes ONLY the seven
  `cosebis_ptes_{ver}_{blind}.npz` (no manifest, no timestamp) so
  `tar czf {output} -C $o .` is clean. Confirm `--xi-sweep-dir` /
  `--cov-sweep-dir` tolerate an untarred layout.
- `papers/bmodes/scripts/cosebis_binning_comparison.py` — `--out-file
  <path.png>`, single figure.
- `papers/bmodes/scripts/cosebis_version_comparison.py` — `--out-file
  <path.png>` writing only the stacked PNG (today `figure_stacked.png` +
  `cosebis_bmode_stacked.pdf` + `evidence.json`). The root's
  `fig_versions_cosebis` is now an alias of this output, not a PDF copy.

### cl_bb
- `workflow/scripts/generate_pseudo_cl.py` — `--out-file <path>` for the single
  decoupled-bandpower FITS (today `--out <dir>` + native `pseudo_cl_<ver>.fits`,
  with skip-if-exists caching keyed on the native name — that cache must not
  silently defeat `--out-file`).
- `papers/bmodes/scripts/run_cl_sweep.py` — `--out-file <path.tar.gz>`: run the
  7-version × 2-artifact sweep in a temp dir under the existing canonical names
  (`pseudo_cl_{ver}_blind=A_powspace_nbins=32.fits`, `pseudo_cl_cov_...`) and
  tar.gz the 14 files flat to that path.
- `papers/bmodes/scripts/cl_data_vector.py` — `--out-file <path.png>`, only the
  fiducial paper Fig. 4 PNG (no evidence.json, no variants, no per-version loop
  in that mode). Must accept a `--results-dir` of flat canonically-named FITS.
- `papers/bmodes/scripts/cl_version_comparison.py` — `--out-file <path.png>`,
  only the 4-version overlay.

### covariance
- `workflow/scripts/run_cosmocov_chain.sh` — `--emit
  {blocks,processed,processed_g,plot}` + `--out-file <file>`; run the chain in a
  private scratch (`mktemp -d` / `$TMPDIR`) and copy exactly the selected
  product out: `blocks` → `covariance.txt`, `processed` →
  `covariance_processed.txt`, `processed_g` → `covariance_processed_g.txt`,
  `plot` → the correlation-structure figure **as 300-dpi PNG**. The chain must
  leave no `cov_tmp*` / `order_cov_tmp*` / `.ini` / `.log` outside scratch.
- `workflow/scripts/generate_pseudo_cl_cov.py` — add `--out-file <file>`
  (single FITS) alongside/replacing `--out` (currently `required=True`,
  help "Output directory (lc {output})", writes `pseudo_cl_cov_{ver}.fits`).
- `papers/bmodes/scripts/run_cov_sweep.sh` — `--out <dir>` → `--out-file
  <file.tar.gz>`: run the 7 per-version chains in scratch and tar the seven
  `covariance_{ver}_A_g_minsep=0.5_maxsep=300.0_nbins=1000_masked/{base}_processed.txt`
  preserving that directory structure inside the archive (consumers untar and
  pass the root as `--cov-dir`); exclude all scratch.
- **Open, not yet wired:** `run_cov_sweep.sh --mask-base <dir>` is still an
  absolute path. Unlike the five single-file `--mask-cls` uses (now
  `{inputs.mask_power_spectrum}`, the catalog unit's output), the sweep selects
  a mask *per version* and so needs the directory holding both footprints. It
  wants either a second declared input for the star-halo mask, or a `--mask-cls`
  / `--mask-cls-starhalo` pair replacing `--mask-base`. The latter is preferable
  — it keeps every mask an lc-declared file.
- `workflow/scripts/compute_glass_mock_covariance.py` — make all seven output
  args optional (currently `required=True`): `--xi-covariance`,
  `--cl-covariance`, `--combined-covariance`, `--combined-correlation-plot`,
  `--xi-mean`, `--cl-mean`, `--combined-mean`; write only those supplied (each
  GLASS recipe passes exactly one).

### null_tests
- `papers/bmodes/scripts/bb_covariance_blind_independence.py` — `--out <dir>` →
  `--out-file <path>` writing the single PNG (the `evidence.json` side write is
  dropped). Its three `--*-dir` arguments now arrive as declared read-only
  archive inputs; no path-composition change needed.

### catalog
- No script changes (its recipes are `cp` / inline `python -c` over
  `cat_config.yaml`; the mask power spectrum now resolves through the declared
  `footprint_masks` input).
