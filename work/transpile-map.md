# Transpile map — UNIONS Paper II B-modes → MyST/MySTRA report

Maps the paper TeX (`paper/unions_bmodes.tex`, 649 lines) and its 89 claim
macros (`paper/claims_macros.tex`) onto the ASTRA reproduction spec
(`astra.yaml` + 7 sub-analyses) and the materialized `results/baseline/`
tree. Active universe: `baseline` (only file in `universes/`).

Path grammar per `mystra-syntax.md`: paths resolve from the **root
analysis** on every page; the `analyses.` prefix is implied, so a
sub-analysis output is `null_tests.outputs.config_space_pte_evidence`.

---

## 0. Headline finding — the value-binding gap

**None of the 89 macros is `{astra:value}`-queryable as the spec stands.**
Every numeric value lives in a materialized artifact, but not in a form the
`{astra:value}` role can read:

- The two `type: table` outputs `outputs.pte_table_results` /
  `outputs.pte_table_appendix` materialize only a LaTeX `tabular` + a
  **metadata-only** `evidence.json` (`{"output": {"macros":
  "claims_macros.tex"}}`) — no columnar CSV/JSON. `{astra:value col= where=}`
  has no table to read.
- The structured per-version PTE numbers live in **`type: figure`** outputs
  (`config_space_pte_evidence`, `harmonic_space_pte_evidence`,
  `cosebis_harmonic_modes`, `pure_eb.pure_eb_covariance_figure`) as a
  **nested dict** (`versions.<ver>.<stat>_stats.pte_at_{fiducial,full_range}`),
  not a flat table. MySTRA treats a figure output as an image; it will not
  interpolate a scalar from it.

**Bridge required before drafting numeric prose:** add `type: table` outputs
that emit a flat CSV — columns `version, statistic, cut, pte, chi2, dof` for
the PTE grid; `block, condition_number` for the covariance — so
`{astra:value col=pte where="version=size_cut statistic=combined cut=fiducial"}`
resolves. Until then, numbers reach the page only by **embedding the finding**
that quotes them (four findings carry the headline PTEs in their `claim`/`notes`
prose) — which renders them but gives no per-macro inline substitution.

The bound/gap split below is therefore **"does a materialized astra artifact
carry this value at all"** (the precondition for the bridge), *not* "is it
queryable today" (none are). "GAP" = the number is absent from every local
reproduction artifact (only in the paper's own `claims_macros.tex`).

---

## 1. Section inventory

`\section` / `\subsection` in document order (no `\appendix` — the "appendix"
PTE table and composite figures are `\input`/`\includegraphics` inside the
Methods and Results floats). Front matter: title/authors L70–152, abstract
L136–150.

| id | title | tex_lines |
|----|-------|-----------|
| introduction | Introduction | 156–213 |
| data | Data | 215–232 |
| methods | Methods | 234–248 (intro + Fig 1 float) |
| methods.two_point | Two-point correlation measurements | 249–273 |
| methods.pure_eb | Pure E/B decomposition | 274–300 |
| methods.cosebis | COSEBIs | 301–352 |
| methods.catalog_cls | Catalog-based harmonic-space power spectra | 354–371 |
| methods.covariance | Covariance estimation | 373–420 |
| methods.scale_cuts | B-mode significance | 421–438 |
| methods.mock_validation | Validation on simulations | 439–451 |
| results | Results | 452–524 |
| discussion | Discussion | 525–599 |
| conclusions | Conclusions | 601–619 |
| data_availability | Data availability (`\section*`) | 621–624 |
| acknowledgements | Acknowledgements | 626–641 |

---

## 2. Macro binding table (all 89)

**Version token → catalog variant / evidence key:**
`Five` = SP_v1.4.5 (initial) · `SixThree` = SP_v1.4.6.3 (size_cut / fiducial)
· `Eight` = SP_v1.4.8 (masked) · `ElevenThree` = SP_v1.4.11.3 (relaxed_flags).
Evidence keys carry the `_leak_corr` suffix.
**Suffix:** no `Full` → `pte_at_fiducial`; `Full` → `pte_at_full_range`.

### 2a. `configPte*` (40 macros) — BOUND
Config-space six-statistic PTE grid.
**Carrier:** `null_tests.outputs.config_space_pte_evidence`
(`evidence.json → versions.<ver>.<stat>_stats.pte_at_{fiducial,full_range}`).
Also duplicated at `outputs.fig_pte_heatmaps` (same producer/evidence).
Statistic token → key: `Xip`→`xip_stats`, `Xim`→`xim_stats`,
`Combined`→`combined_stats`, `Cosebis`→`cosebis_stats` (n≤6),
`CosebisTwenty`→`cosebis_20_stats` (n≤20).

| macro group | count | tree-path |
|---|---|---|
| `configPte{Five,SixThree,Eight,ElevenThree}{Xip,Xim,Combined,Cosebis,CosebisTwenty}` (fiducial) | 20 | `null_tests.outputs.config_space_pte_evidence` → `versions.<ver>.<stat>_stats.pte_at_fiducial` |
| same × `…Full` (full range) | 20 | `null_tests.outputs.config_space_pte_evidence` → `versions.<ver>.<stat>_stats.pte_at_full_range` |

Spot-check: `configPteSixThreeCombined`=0.18 ↔ `combined_stats.pte_at_fiducial`=0.1758 ✓;
`configPteSixThreeCosebisFull`=1.37e-5 ↔ `cosebis_stats.pte_at_full_range`=1.366e-5 ✓.

### 2b. `clPte*` (8 macros) — BOUND
Harmonic-space Cℓ^BB PTE per version.
**Carrier:** `null_tests.outputs.harmonic_space_pte_evidence`
(`versions.<ver>.pte_at_{fiducial,full_range}`); dup at `outputs.fig_cl_pte_heatmap`.

| macro group | count | tree-path |
|---|---|---|
| `clPte{Five,SixThree,Eight,ElevenThree}Fid` | 4 | `null_tests.outputs.harmonic_space_pte_evidence` → `versions.<ver>.pte_at_fiducial` |
| `clPte{…}Full` | 4 | `null_tests.outputs.harmonic_space_pte_evidence` → `versions.<ver>.pte_at_full_range` |

### 2c. `cfgCosebis*` (16 macros) — 12 BOUND / 4 GAP
Config-path COSEBI n≤6 B-mode, PTE + χ², per version, fiducial + full range.
**Carrier (fiducial):** `cosebis.outputs.cosebis_harmonic_modes`
(`config_b_mode_ptes.<ver>.{pte,chi2}`, dof=6). The full-range **PTE** also
sits in `config_space_pte_evidence.<ver>.cosebis_stats.pte_at_full_range`.

| macro group | count | tree-path / verdict |
|---|---|---|
| `cfgCosebisPte{...}Fid` | 4 | `cosebis.outputs.cosebis_harmonic_modes` → `config_b_mode_ptes.<ver>.pte` — BOUND |
| `cfgCosebisChisq{...}Fid` | 4 | `cosebis.outputs.cosebis_harmonic_modes` → `config_b_mode_ptes.<ver>.chi2` — BOUND |
| `cfgCosebisPte{...}Full` | 4 | `null_tests.outputs.config_space_pte_evidence` → `versions.<ver>.cosebis_stats.pte_at_full_range` — BOUND |
| `cfgCosebisChisq{...}Full` | 4 | `cosebis.outputs.cosebis_harmonic_modes_full` → `config_b_mode_ptes.<ver>.chi2` — BOUND (materialized 2026-07-20) |

Spot-check: `cfgCosebisPteSixThreeFid`=0.78 ↔ `config_b_mode_ptes.…6.3.pte`=0.7766 ✓;
`cfgCosebisChisqSixThreeFid`=3.25 ↔ `…chi2`=3.252 ✓.

### 2d. `harmCosebis*` (16 macros) — BOUND
Harmonic-path COSEBI n≤6 B-mode, PTE + χ², per version, fiducial + full.
**Carrier (fiducial):** `cosebis.outputs.cosebis_harmonic_modes`;
**(full range):** `cosebis.outputs.cosebis_harmonic_modes_full` (same script,
`--angular-range full`) — both `harmonic_b_mode_ptes.<ver>.{pte,chi2}`, dof=6.

| macro group | count | tree-path / verdict |
|---|---|---|
| `harmCosebisPte{...}Fid` | 4 | `cosebis.outputs.cosebis_harmonic_modes` → `harmonic_b_mode_ptes.<ver>.pte` — BOUND |
| `harmCosebisChisq{...}Fid` | 4 | `cosebis.outputs.cosebis_harmonic_modes` → `harmonic_b_mode_ptes.<ver>.chi2` — BOUND |
| `harmCosebisPte{...}Full` | 4 | `cosebis.outputs.cosebis_harmonic_modes_full` → `harmonic_b_mode_ptes.<ver>.pte` — BOUND (materialized 2026-07-20) |
| `harmCosebisChisq{...}Full` | 4 | `cosebis.outputs.cosebis_harmonic_modes_full` → `harmonic_b_mode_ptes.<ver>.chi2` — BOUND (materialized 2026-07-20) |

Spot-check: `harmCosebisPteSixThreeFid`=0.60 ↔ `harmonic_b_mode_ptes.…6.3.pte`=0.5979 ✓.
**Full-range reproduction (2026-07-20) is exact to the paper's quoted precision:**
`harmCosebisPteSixThreeFull` 1.608e-5 ↔ paper 1.61e-5, χ² 32.03 ↔ 32.03;
all four versions' PTE+χ² and all four `cfgCosebisChisq*Full` (32.40/57.15/54.00/34.69)
land on the macro values. `harmCosebisPteSixThreeFull` also promoted into the flat
`pte_promoted_values` table (6th row) for `{astra:value}` queryability.

### 2e. `ebcov*` (5 macros) — BOUND
Pure E/B covariance block condition numbers + bin count.
**Carrier:** `pure_eb.outputs.pure_eb_covariance_figure` (`evidence.json`).

| macro | tree-path |
|---|---|
| `ebcovCondE` (1.5e5) | `pure_eb.outputs.pure_eb_covariance_figure` → `block_analysis.xi_E.condition_number` |
| `ebcovCondB` (2.0e5) | `…block_analysis.xi_B.condition_number` |
| `ebcovCondAmb` (1.8e10) | `…block_analysis.xi_amb.condition_number` |
| `ebcovCondFull` (2.7e10) | `…condition_number` (full 120×120) |
| `ebcovNbins` (20) | `…n_bins` |

### 2f. `ebtheta*` (4 macros) — BOUND (decision-carried)
Adopted config-space scale-cut boundaries (12′, 83′).
**Carrier:** `decisions.fiducial_scale_cuts` option `adopted` ("12′–83′").
Also in `config_space_pte_evidence` / `cosebis_harmonic_modes.scale_cut`
=[11.997, 82.86].

| macro | tree-path |
|---|---|
| `ebthetaXipMin`, `ebthetaXimMin` (12) | `decisions.fiducial_scale_cuts` (option `adopted`) |
| `ebthetaXipMax`, `ebthetaXimMax` (83) | `decisions.fiducial_scale_cuts` (option `adopted`) |

*Queryability nuance:* `{astra:value}`decisions.fiducial_scale_cuts`` renders the
option **label**, not the bare integer 12/83 — a metric or table bridge is
needed to emit `\num{12}`.

**Macro tally: 89 total → 89 BOUND, 0 GAP** (2026-07-20).
The former 12-macro gap (`cfgCosebisChisq*Full` ×4 + `harmCosebis*Full` ×8) all
traced to one missing computation — the harmonic-path COSEBIs at the full angular
range. Materializing `cosebis.cosebis_harmonic_modes_full` (the `--angular-range
full` twin of `cosebis_harmonic_modes`) carries every one of them, and the
reproduced values match the paper macros exactly.

---

## 3. Figure table

8 figure floats (Fig 5 is a 3-panel composite) + 1 in-text table.
Root results at `results/baseline/<id>.pdf` and `results/baseline/<id>/`.

| # | TeX file (label) | tex_line | tree-path | verdict |
|---|---|---|---|---|
| 1 | `pure_eb_data_vector.pdf` (`fig:pure_eb_decomposition`) | 239 | `outputs.fig_pure_eb_decomposition` | BOUND |
| 2 | `cosebis_data_vector.pdf` (`fig:cosebis_fiducial`) | 321 | `outputs.fig_cosebis_fiducial` | BOUND |
| 3 | `cl_data_vector.pdf` (`fig:cl_fiducial`) | 365 | `outputs.fig_cl_fiducial` | BOUND |
| 4 | `eb_covariance.pdf` (`fig:eb_covariance`) | 404 | `pure_eb.outputs.pure_eb_covariance_figure` | BOUND — **not promoted to a root output** (no `fig_eb_covariance`); reference the sub-analysis path |
| — | Table 1 `\input{pte_table_appendix}` | 418 | `outputs.pte_table_appendix` | BOUND (embed); value cells not `{astra:value}`-queryable (§0) |
| 5 | `pure_eb_versions.pdf` + `cl_versions.pdf` + `cosebis_bmode_stacked.pdf` (`fig:pure_eb_versions`) | 457–459 | `outputs.fig_pure_eb_versions` (materializes all 3 PDFs) | BOUND |
| 6 | `config_space_pte_composite_appendix.pdf` (`fig:pte_heatmaps`) | 475 | `outputs.fig_pte_heatmaps` | BOUND |
| 7 | `cl_pte_composite_appendix.pdf` (`fig:pte_cl`) | 486 | `outputs.fig_cl_pte_heatmap` | BOUND |
| 8 | `harmonic_config_cosebis_full.pdf` (`fig:harmonic_config_cosebis_full`) | 540 | `outputs.fig_harmonic_config_cosebis_full` | BOUND (materialized 2026-07-20) — the paper's `…_full.pdf` variant now has its own root output; `fig_harmonic_config_cosebis` (`…_fiducial.pdf`) remains as the non-paper fiducial panel |

**Figure tally: 8 figures → 8 BOUND, 0 GAP** (1 caveat: Fig 4 not root-promoted — reference the sub-analysis path. Fig 8 full variant materialized 2026-07-20).
`pte_table_results` (2-row size-cut summary) exists as `outputs.pte_table_results`
but is **not `\input` anywhere** in the current TeX (only `pte_table_appendix` is).

---

## 4. Proposed report page structure

Multi-page, mirroring the 7-sub-analysis DAG decomposition (paper-reproduction
mode structures by the DAG, harvesting the TeX for content). Root page is the
end-to-end view; details telescope into sub-analysis pages.

| page file | scope | content |
|---|---|---|
| `index.md` | root | Intro (question + headline: size-cut uniquely passes) · end-to-end Methods overview in DAG order · Results synthesis · Discussion + Conclusions. Embeds Table 1 (`outputs.pte_table_appendix`), version-comparison Fig 5, the two PTE composites (Figs 6–7), the 4 findings, and sub-analysis nav cards (`:::{astra} <sub>:::`). Drafted last. |
| `catalog.md` | `catalog` | Four v1.4 variants, masks, n(z); `decisions.fiducial_catalog`, `leakage_correction`. Thin — foldable into index if too light. |
| `two_point.md` | `two_point` | TreeCorr ξ± (20-bin reporting, 1000-bin integration); `decisions.reporting_binning`, `integration_grid`. |
| `pure_eb.md` | `pure_eb` | Pure E/B decomposition (§3.2); Fig 1, Fig 4 (`pure_eb_covariance_figure`); `decisions.mc_samples`, `hartlap_application`. |
| `cosebis.md` | `cosebis` | COSEBIs (§3.3); Fig 2; harmonic-vs-config cross-check Fig 8; `findings.cosebis_full_range_failure`, `filter_sets_sensitivity`. |
| `cl_bb.md` | `cl_bb` | Catalog-based pseudo-Cℓ (§3.4); Fig 3; `decisions.estimator_choice`, `bandpower_binning`. |
| `covariance.md` | `covariance` | CosmoCov + MC propagation + GLASS validation (§3.5, §3.7); `decisions.covariance_engine`, `glass_validation_strategy`. |
| `null_tests.md` | `null_tests` | PTE framework (§4.1), Figs 6–7, Table 1 detail; `findings.size_cut_unique_pass`, `multi_statistic_consensus`; `decisions.pte_definition`, `healthy_pte_range`, `fiducial_scale_cuts`. |

Planned page files: `index.md`, `catalog.md`, `two_point.md`, `pure_eb.md`,
`cosebis.md`, `cl_bb.md`, `covariance.md`, `null_tests.md` (8 pages).
Register all in `myst.yml` `toc:`. `myst.yml` does not yet exist in the
project root — `lc init` scaffold or minimal template needed before first build.
