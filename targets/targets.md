# UNIONS B-Modes Paper — Replication Targets

## Headline Claims

1. **All three B-mode statistics pass null tests at fiducial scale cuts for the size-cut catalog.** (Minimum PTE = 0.18 across ξ₊ᴮ, ξ₋ᴮ, ξ_totᴮ combined; Abstract, §4.1, Eq. 141) — *This is the core validation result: the adopted cuts eliminate B-mode contamination simultaneously across pure-mode, COSEBI, and harmonic statistics.*

2. **Pure E/B decomposition ξ_totᴮ (combined) PTE = 0.40 at full range; 0.18 at fiducial cuts (12–83 arcmin, size-cut catalog).** (§4.1, Table 1, claims_macros.tex configPteSixThreeCombined/Combined Full) — *The fiducial cut reduces PTE failure from near-significance to robust acceptance.*

3. **Configuration-space COSEBI Bₙ (n≤6) PTE = 0.78 at fiducial; full-range oscillatory pattern drops to PTE = 1.37×10⁻⁵ at n≤6, consistent with repeating CCD-scale additive bias.** (§4.1, Table 1, Abstract line 142, claims_macros.tex) — *Multi-mode COSEBI signature traces detector-level systematics across multiple Stage-III surveys.*

4. **Harmonic-space Cℓᴮᴮ PTE = 0.30 at fiducial (300 < ℓ < 1600, size-cut catalog).** (§4.1, Table 1, claims_macros.tex clPteSixThreeFid) — *Full-range harmonic space passes (PTE=0.13) but at lower sensitivity than fiducial cuts due to filter-function differences.*

5. **Only the size-cut catalog passes all three statistics simultaneously at fiducial cuts; masked, initial, and relaxed-flags versions each fail in at least one framework.** (Abstract, §4.1, §5, Table 1) — *This validates that multi-framework consistency is essential for systematic contamination assessment; single-statistic tests are insufficient.*

## Figure Targets

| ID | Caption Summary | What it Shows | Panels/Details | Source File | Sub-Analysis |
|:---|:---|:---|:---|:---|:---|
| fig_pure_eb_decomposition | Pure E/B-mode decomposition of measured shear correlations, fiducial catalog | ξ₊ and ξ₋ decomposed into E-mode (lensing), B-mode (systematics), ambiguous (boundary), and total; error bars from semi-analytical covariance propagation | Two panels (ξ₊, ξ₋); each shows total (black circles), E-mode (teal squares), B-mode (crimson ×), ambiguous (purple triangles); shaded bands mark fiducial scale cuts | pure_eb_data_vector.pdf | pure_eb |
| fig_eb_covariance | Pure E/B-mode correlation-function covariance matrix from 2000 MC realizations, fiducial catalog | Covariance structure showing block correlations within E, B, ambiguous modes | 20×20 bin matrix; heatmap with log scale | eb_covariance.pdf | pure_eb covariance |
| fig_cosebis_fiducial | COSEBI B-mode amplitudes (Bₙ/σₙ) for first 20 modes, fiducial catalog | Full angular range (1–250 arcmin) shows >4σ oscillatory pattern in first mode; fiducial cuts (12–83 arcmin) suppress to <1σ | Orange squares (full range), blue circles (fiducial); 20 modes; error bars shown | cosebis_data_vector.pdf | cosebis |
| fig_cl_fiducial | Harmonic-space power spectra Cℓᴮᴮ and Cℓᴱᴮ normalized by uncertainty, fiducial catalog | Low-level positive offset in Cℓᴮᴮ across fiducial range (300 < ℓ < 1600); outliers near ℓ≈125, 250, >1600; Cℓᴱᴮ consistent with zero | Two panels (Cℓᴮᴮ, Cℓᴱᴮ) vs ℓ; error bars; shaded fiducial range | cl_data_vector.pdf | cl_bb |
| fig_pure_eb_versions | B-mode data vectors across four catalog versions (initial, fiducial, masked, relaxed-flags) | Overlaid ξ₊ᴮ, ξ₋ᴮ, Cℓᴮᴮ, Cℓᴱᴮ, COSEBI Bₙ across versions; rectangles show spread, horizontal line marks fiducial | Three rows (pure-mode ξ, harmonic Cℓ, COSEBI Bₙ); each with full range and fiducial cuts marked; four versions color-coded | pure_eb_versions.pdf, cl_versions.pdf, cosebis_bmode_stacked.pdf | all sub-analyses |
| fig_pte_heatmaps | Configuration-space PTE maps for all four catalog versions vs angular scale cuts | 2D grid of PTE values for ξ₊ᴮ, ξ₋ᴮ, COSEBI Bₙ as function of lower/upper θ cuts | Four rows (catalog versions); three columns (ξ₊ᴮ, ξ₋ᴮ, COSEBI); blue=failing (PTE<0.05), white/red=consistent; black square marks adopted cuts (12–83 arcmin) | config_space_pte_composite_appendix.pdf | pure_eb, cosebis |
| fig_pte_cl | Harmonic-space Cℓᴮᴮ PTE maps for all four catalog versions | 2D grid of PTE values vs ℓ_min, ℓ_max cuts for each catalog version | Four rows (catalog versions); blue=failing, white/red=passing; adopted cuts (300 < ℓ < 1600) marked | cl_pte_composite_appendix.pdf | cl_bb |
| fig_harmonic_config_cosebis_full | COSEBI Eₙ and Bₙ/σₙ computed from configuration-space ξ± vs harmonic-space Cℓ bandpowers | Direct comparison: both computation paths yield identical E_n and B_n; validates filter-function explanation for disagreements | Two rows (Eₙ top, Bₙ/σₙ bottom); filled circles=config-space, open squares=harmonic-space; gray shading marks n≤6; 20 modes | harmonic_config_cosebis_full.pdf | pure_eb, cosebis, cl_bb |

## Table Targets

### Table 1: B-Mode PTE Summary (Main Results, §4.1)

**Caption:** "B-mode PTE values across catalog versions at fiducial and full-range scale cuts."

**Columns:** 
- Scale cuts (fiducial / full range)
- COSEBI Bₙ (n≤6)
- COSEBI Bₙ (n≤20)
- Pure E/B ξ₊ᴮ
- Pure E/B ξ₋ᴮ
- Pure E/B ξ_totᴮ (combined)
- Harmonic Cℓᴮᴮ

**Rows:**
- Fiducial (size-cut catalog): 0.78 | 0.94 | 0.31 | 0.26 | 0.18 | 0.30
- Full range: **1.37e-05** | **1.08e-04** | 0.45 | 0.16 | 0.40 | 0.13

**Interpretation:** Fiducial cuts transform from low-PTE COSEBI (detector-scale contamination) to stable acceptance across all statistics. Full-range harmonic space passes but configuration-space statistics fail.

**Source:** pte_table_results.tex (auto-generated from claims evidence)

---

### Appendix Table: PTE Comparison Across Catalog Versions

**Caption:** "B-mode PTE comparison for initial, size-cut (fiducial), masked, and relaxed-flags catalog versions at fiducial and full-range scale cuts."

**Grid Structure:**
```
Version              Fiducial                                   Full Range
                 COSEBI6 COSEBI20 ξ+B   ξ-B   ξ_tot  CℓBB | COSEBI6 COSEBI20 ξ+B   ξ-B   ξ_tot  CℓBB
Initial          0.94    0.61    0.51  0.005 0.030 0.003  | 7.37e-10 4.15e-07 0.026 7.71e-04 0.004 2.53e-04
Size-cut/fiducial 0.78    0.94    0.31  0.26  0.18  0.30  | 1.37e-05 1.08e-04 0.45  0.16     0.40  0.13
Masked           0.60    0.81    0.020 0.047 0.014 0.43  | 4.94e-06 2.50e-05 0.039 0.020    0.018 0.11
Relaxed-flags    0.82    0.83    0.78  9.4e-04 0.014 5.74e-10 | 1.70e-10 1.16e-08 0.10  2.94e-06 9.07e-04 1.78e-12
```

**Key Results:** 
- Initial catalog fails in ξ₋ᴮ (PTE=0.005) and full-range harmonic space
- Fiducial (size-cut) is only version passing all statistics at fiducial cuts
- Masked catalog paradoxically fails pure-mode tests (pure-mode failures despite harmonic/COSEBI n≤6 passes) — suggests mask-induced ambiguous-mode sensitivity
- Relaxed-flags fails in ξ₋ᴮ at fiducial (approaches threshold at PTE=0.78 vs 0.26 fiducial)

**Source:** pte_table_appendix.tex

---

## Quantitative Results Table

| Quantity | Value | Uncertainty / Notes | Paper Location |
|:---|:---|:---|:---|
| **Survey Properties** | | | |
| Survey area | 2894 | deg² | Abstract |
| Effective galaxy density (n_eff) | 5.0 | arcmin⁻² | Abstract |
| Catalog version (fiducial) | SP_v1.4.6 (size-cut) | + PSF leakage correction | §2, §4 |
| Galaxy size cut threshold | r_h,gal / r_h,psf > 0.707 | (vs initial 0.5) | §2 |
| Stellar masking | r<12 (bright), 12<r<18 (faint) | 13% area loss (masked version only) | §2 |
| **Configuration-Space Fiducial Scale Cuts** | | | |
| Angular range (fiducial, ξ±) | 12–83 | arcmin | §4.1, claims_macros |
| Angular bins | 20 | logarithmic | §3.1 |
| Angular range (integration) | 0.5–300 | arcmin | §3.2 |
| Integration bins | 1000 | logarithmic | §3.2 |
| **Harmonic-Space Fiducial Scale Cuts** | | | |
| Multipole range (fiducial, Cℓᴮᴮ) | 300 < ℓ < 1600 | — | §4.1, claims_macros |
| Full multipole range | 8 ≤ ℓ ≤ 2048 | — | §4.1, claims_macros |
| **Covariance Matrix Condition Numbers (Pure E/B)** | | | |
| E-mode block (full cov) | 1.5e+05 | — | claims_macros |
| B-mode block (full cov) | 2.0e+05 | — | claims_macros |
| Ambiguous-mode block | 1.8e+10 | — | claims_macros |
| Full covariance matrix | 2.7e+10 | — | claims_macros |
| **Fiducial-Cut PTEs (Size-Cut Catalog)** | | | |
| ξ₊ᴮ (pure-mode) | 0.31 | — | Table 1 |
| ξ₋ᴮ (pure-mode) | 0.26 | — | Table 1 |
| ξ_totᴮ (combined, fiducial cuts) | **0.18** | Minimum PTE, fiducial | Table 1, Abstract |
| COSEBI Bₙ (n≤6) | 0.78 | — | Table 1 |
| COSEBI Bₙ (n≤20) | 0.94 | — | Table 1 |
| Cℓᴮᴮ (harmonic) | 0.30 | — | Table 1 |
| **Full-Range PTEs (Size-Cut Catalog)** | | | |
| ξ₊ᴮ (full range) | 0.45 | — | Table 1 |
| ξ₋ᴮ (full range) | 0.16 | — | Table 1 |
| ξ_totᴮ (full range) | 0.40 | — | Table 1 |
| COSEBI Bₙ (n≤6, full range) | **1.37e-05** | Strong oscillatory pattern | Table 1 |
| COSEBI Bₙ (n≤20, full range) | **1.08e-04** | — | Table 1 |
| Cℓᴮᴮ (full range) | 0.13 | — | Table 1 |
| **COSEBI Harmonic vs Config Agreement (Fiducial, n≤6)** | | | |
| Harmonic-space path (config-computed COSEBIs) | 0.60 (fid) / 1.61e-05 (full) | χ²=4.59 (fid) / 32.03 (full) | §5, claims_macros |
| Config-space path | 0.78 (fid) / 1.37e-05 (full) | χ²=3.25 (fid) / 32.40 (full) | §5, claims_macros |
| **Initial Catalog Fiducial-Cut Results** | | | |
| ξ₋ᴮ PTE | 0.005 | **FAILS** null test | Appendix Table |
| ξ_totᴮ PTE | 0.030 | **FAILS** null test | Appendix Table |
| Cℓᴮᴮ PTE | 0.003 | **FAILS** null test | Appendix Table |
| **Masked Catalog Fiducial-Cut Results** | | | |
| ξ₊ᴮ PTE | 0.020 | **FAILS** pure-mode | Appendix Table |
| ξ₋ᴮ PTE | 0.047 | **FAILS** pure-mode (margin) | Appendix Table |
| ξ_totᴮ PTE | 0.014 | **FAILS** pure-mode | Appendix Table |
| Cℓᴮᴮ PTE | 0.43 | **PASSES** harmonic | Appendix Table |
| **Relaxed-Flags Catalog Fiducial-Cut Results** | | | |
| ξ₊ᴮ PTE | 0.78 | **PASSES** pure-mode plus statistic | Appendix Table, §4.1 |
| ξ₋ᴮ PTE | 9.40e-04 | **FAILS** pure-mode minus statistic | Appendix Table, §4.1 |
| Cℓᴮᴮ PTE | 5.74e-10 | **FAILS** harmonic fiducial cut | Appendix Table |
| Cℓᴮᴮ PTE (full range) | 1.78e-12 | **FAILS** harmonic full range | Appendix Table |
| **Simulation Validation** | | | |
| Simulation ensemble | 350 | lognormal GLASS realizations | §3.4 |
| Simulated B-mode recovery | consistent with zero | across all three statistics | §3.4 |
| COSEBI mock mean bias | <0.3σ | per mode over 100-realisation subset | §3.4 |
| Mock scatter | comparable | to fiducial error bars for each statistic | §3.4 |
| Cosmology (simulations + analysis) | Planck 2018 | fiducial | §3.4 |
| **MegaCam Detector Geometry** | | | |
| CCD size in field | 6–14 | arcmin | §5, discussion |
| UNIONS dither size | ~1/3 | field of view | §5 |
| COSEBI feature suppression threshold | ≥12 | arcmin (adopted lower cut) | §5 |

---

## Out-of-Scope

Explicitly NOT part of reproduction targets:

- **Cosmological constraints (S₈, Ωm, σ₈):** These appear only in context (e.g., KiDS-Legacy S₈ = 0.815±0.016, DES Y6 S₈ = 0.798±0.015) as Stage-III benchmarks for tension with Planck. Cosmological inference is in Papers III and IV; we validate data-quality here only.

- **Leakage correction details and PSF diagnostics:** The Paper I leakage model and PSF diagnostic machinery are out of scope. Paper II's limited B-mode statement — PTE shifts of <0.05 with no pass/fail changes — is represented by `decisions.leakage_correction`; no Paper I diagnostic products are reproduced.

- **Full simulation validation pipeline:** Papers III and IV perform separate cosmological-inference validation on GLASS simulations. This paper validates only that B-mode estimators recover zero on simulations (§3.4).

- **Companion paper results:** References to Papers III (Goh et al., config-space cosmology), IV (Guerrini et al., harmonic-space cosmology), V (Hervas-Peters et al., simulations) are context only. We validate the scale cuts they adopt, not their results.

- **Multi-field or tomographic extensions:** This paper uses a single tomographic bin (all galaxies combined). Multi-bin analyses are future work.

- **High-order astrophysical B-modes:** Lens-lens coupling, source clustering, intrinsic alignment—mentioned as theoretical possibilities but orders of magnitude below current sensitivity; not target of validation.

---

## Notes on Metrics

1. **PTE interpretation:** PTE < 0.05 conventional threshold for "failure." This paper emphasizes that no single statistic is sufficient; multi-framework consensus is required. Fiducial cuts achieve PTE > 0.18 across all statistics in size-cut catalog.

2. **Scale-cut notation:** Fiducial cuts (12–83 arcmin, 300 < ℓ < 1600) are chosen conservatively based on broad stable PTE regions, not driven by PTE thresholds alone. Papers III & IV apply identical cuts for blinded cosmological inference.

3. **Covariance treatment:** All null tests use Gaussian-only covariance (conservative, underestimating errors). Non-Gaussian terms present in data but omitted from reported PTEs to match Paper IV's harmonic-space methodology.

4. **Oscillatory COSEBI pattern:** Traced to CCD-scale repeating additive bias (detector-level effect fixed in focal-plane coordinates, imprinted across multiple MegaCam exposures). Dither strategy and potential WCS recentering issues discussed in §5 and Paper V.

5. **Fiducial adoption logic:** Size-cut catalog chosen not purely on lowest PTE, but on simultaneous passage of all three statistics. Masked catalog paradox (passes harmonic/COSEBI n≤6 but fails pure-mode tests) shows that filter-function sensitivity varies; multi-framework requirement is essential.

---

## Discharge Ledger

Where each target above is represented in the spec. Use this to verify that a row's evidence is cited by some declared structure (input, output, decision, finding, prior insight, or sub-analysis) and reachable from a narrative anchor. Rows discharged as "Out of scope" are deliberately not represented; they are listed so the omission itself is auditable.

### Headline Claims

| # | Discharge |
|:--|:--|
| 1. All three statistics pass at fiducial cuts | `findings.size_cut_unique_pass` (root) — claim + notes line up with paper Abstract / §4.1 / Eq. 141; cited in `narrative.findings`. Numerical evidence anchored to `outputs.pte_table_appendix`. |
| 2. ξ_tot^B PTE 0.40 (full) → 0.18 (fiducial) | `findings.size_cut_unique_pass.notes` includes 0.18; full-range 0.40 traces through `outputs.pte_table_results` (Table 1) which is generated by `analyses.null_tests.outputs.pte_summary_evidence`. Paper macros: `configPteSixThreeCombined`, `configPteSixThreeCombinedFull`. |
| 3. Full-range COSEBIs PTE = 1.37×10⁻⁵ → repeating additive bias | `findings.cosebis_full_range_failure` (root) — quotes the macro `cfgCosebisPteSixThreeFull`; cited in `narrative.findings`; the systematic interpretation links to `prior_insights.repeating_additive_bias_signature` via the same narrative paragraph. |
| 4. Cℓ^BB PTE = 0.30 fiducial / 0.13 full | `outputs.pte_table_results` (root) row; `analyses.cl_bb` discharges via `outputs.cl_pte_per_cut`; macros `clPteSixThreeFid`, `clPteSixThreeFull`. |
| 5. Only size-cut passes simultaneously | `findings.size_cut_unique_pass` + `findings.multi_statistic_consensus` (root); both cited in `narrative.findings`. Per-version evidence in `outputs.pte_table_appendix`. |

### Figure Targets

| Target ID | Discharge |
|:--|:--|
| fig_pure_eb_decomposition | `outputs.fig_pure_eb_decomposition` (root). Recipe: `snakemake pure_eb_data_vector`. Source: `analyses.pure_eb.outputs.pure_eb_data_vector_figure`. |
| fig_eb_covariance | `analyses.pure_eb.outputs.pure_eb_covariance_figure` (recipe `pure_eb_covariance`); narrative cites it in pure_eb summary + outputs. The root does not re-export this paper Figure 2 since it lives in the sub-analysis where it is generated; it is reachable from root narrative through `analyses.pure_eb`. |
| fig_cosebis_fiducial | `outputs.fig_cosebis_fiducial` (root) ↔ `analyses.cosebis.outputs.cosebis_modes_data` (data) → recipe: `snakemake cosebis_data_vector`. Cited in `findings.cosebis_full_range_failure.evidence`. |
| fig_cl_fiducial | `outputs.fig_cl_fiducial` (root). Recipe: `snakemake cl_data_vector`. Source: `analyses.cl_bb.outputs.cl_data_vector_figure`. |
| fig_pure_eb_versions (paper Figure 5 stack) | `outputs.fig_pure_eb_versions` (root). Recipe: `snakemake pure_eb_version_comparison cl_version_comparison cosebis_version_comparison`. Sources: `analyses.cosebis.outputs.cosebis_version_comparison` + per-statistic siblings. |
| fig_pte_heatmaps (config-space PTE composite) | `outputs.fig_pte_heatmaps` (root). Recipe: `snakemake config_space_pte_matrices`. Source: `analyses.null_tests.outputs.config_space_pte_evidence`. |
| fig_pte_cl (harmonic-space PTE composite) | `outputs.fig_cl_pte_heatmap` (root). Recipe: `snakemake harmonic_space_pte_matrices`. Source: `analyses.null_tests.outputs.harmonic_space_pte_evidence`. ID divergence: targets.md uses `fig_pte_cl`, spec uses `fig_cl_pte_heatmap`; equivalent. |
| fig_harmonic_config_cosebis_full | `outputs.fig_harmonic_config_cosebis` (root). Recipe: `snakemake harmonic_config_cosebis_comparison`. Source: `analyses.cosebis.outputs.cosebis_harmonic_modes`. ID divergence: targets.md uses `fig_harmonic_config_cosebis_full`, spec uses `fig_harmonic_config_cosebis`; equivalent. |

### Table Targets

| Target | Discharge |
|:--|:--|
| Table 1 (Main PTE Summary) | `outputs.pte_table_results` (root). Recipe: `snakemake paper_macros`. Source evidence: `analyses.null_tests.outputs.pte_summary_evidence`. |
| Appendix PTE Table (all four versions) | `outputs.pte_table_appendix` (root). Same recipe; same evidence input. Anchored from `findings.size_cut_unique_pass.evidence` and `findings.multi_statistic_consensus.evidence`. |

### Quantitative Results

| Group / Quantity | Discharge |
|:--|:--|
| Survey area 2894 deg² | `inputs.footprint_masks` (root) — described as "standard 2894 deg² footprint"; baseline mask selection in `universes/baseline.yaml`. |
| n_eff = 5.0 arcmin⁻² | Carried in `analyses.catalog.outputs` (catalog packaging exposes survey properties); paper macros render this from Paper I tables. Not surfaced as its own root input/decision since it is derived. |
| Fiducial catalog SP_v1.4.6.3 | `decisions.fiducial_catalog.options.size_cut` (root) + baseline universe selection. |
| Galaxy size cut threshold r_h,gal/r_h,psf > 0.707 | `decisions.fiducial_catalog.options.size_cut.description` (root) — encodes the threshold; carried by Paper I rather than re-derived here. |
| Stellar masking (r<12, 12<r<18) | `decisions.fiducial_catalog.options.masked.description`. Out of scope for fiducial spec since `masked` is not the default. |
| Angular range (fiducial ξ±) 12–83 arcmin | `decisions.fiducial_scale_cuts.options.adopted` (root). |
| Angular bins = 20 (logarithmic reporting) | `analyses.two_point.decisions.reporting_binning` — captures both bin count and log spacing. |
| Integration grid 0.5–300 arcmin / 1000 bins | `decisions.integration_grid.options.bins_1000` (root); see `implementation-notes.md` for the paper's "0.5–500" prose discrepancy (operational parity with code is 300). |
| Multipole range 300 < ℓ < 1600 | `decisions.fiducial_scale_cuts.options.adopted` (root) — carries both θ and ℓ bounds. |
| Full multipole range 8 ≤ ℓ ≤ 2048 | Implicit in the harmonic full-range PTE row of `outputs.pte_table_results`; `analyses.cl_bb` describes the production pseudo-Cℓ binning. Not surfaced as its own decision. |
| Covariance condition numbers (E, B, amb, full) | Operational diagnostics on the MC pure-mode covariance; not represented as paper findings. Documented in `analyses.pure_eb.outputs.pure_eb_semianalytic_data.description` as the four `\ebcovCond{E,B,Amb,Full}` paper macros; rendered by `paper_macros` and consumed only by paper prose. |
| Fiducial-cut PTEs (six values) | `outputs.pte_table_results` (root) — generated row-by-row by `analyses.null_tests.outputs.pte_summary_evidence`; macro names are `configPteSixThree*`, `cfgCosebisPteSixThreeFid`, `clPteSixThreeFid`. |
| Full-range PTEs (six values) | Same, with macro names `configPteSixThree*Full`, `cfgCosebisPteSixThreeFull`, `clPteSixThreeFull`. |
| Harmonic-vs-config COSEBI agreement | `findings.filter_sets_sensitivity` (root) carries the synthesis claim that filter functions, not basis, set null-test sensitivity. `analyses.cosebis.outputs.cosebis_harmonic_modes` carries the harmonic path; `analyses.cosebis.outputs.cosebis_modes_data` carries the config path. Paper Figure 8 = `outputs.fig_harmonic_config_cosebis` (root). PTE values at fid / full appear via `paper_macros` (`harmCosebisPteSixThreeFid` / `Full`). |
| Initial / masked / relaxed-flags catalog PTEs | Per-row in `outputs.pte_table_appendix`; the underlying multiverse choices are `decisions.fiducial_catalog.options.{initial,masked,relaxed_flags}`. |
| Simulation ensemble = 350 GLASS realisations | `inputs.glass_mock_realisations` (root); covariance validation via `analyses.covariance.decisions.glass_validation_strategy`. |
| Simulated B-mode recovery consistent with zero | Discharged as a covariance-validation cross-check; not a paper finding. Anchored in `analyses.covariance` narrative, including the <0.3σ COSEBI mock-mean bias and comparable-scatter claims from paper §3.4. |
| MegaCam CCD geometry 6–14 arcmin | Mentioned in `findings.cosebis_full_range_failure.claim` and in the root narrative findings paragraph (just above the adopted-cut motivation) as the rationale for the 12' lower cut. Not stored as a paper macro — appears as inline numeric prose in `paper/unions_bmodes.tex` §5. |
| COSEBI feature suppression threshold ≥ 12' | Equal to the lower bound of `decisions.fiducial_scale_cuts.options.adopted`; same anchor. |
| B-mode bias subtraction rejected | `decisions.fiducial_scale_cuts.options.bmode_bias_subtraction.excluded_reason` records the §5 argument: B-mode amplitudes diagnose contamination but do not determine the E-mode bias without a mechanism-specific E/B relation. |

### Out-of-Scope (acknowledged absences)

| Item | Reason absent |
|:--|:--|
| Cosmological constraints (S₈, Ω_m, σ_8) | Paper III (config-space) and Paper IV (harmonic-space) carry these; explicitly out of scope per the constitution and the spec's `narrative.summary`. |
| Leakage correction & PSF diagnostics | Paper I details remain absent; Paper II's <0.05/no-pass-fail B-mode statement is surfaced in `decisions.leakage_correction.rationale`, not as a finding. |
| Full simulation validation pipeline | Papers III/IV; this spec keeps only the GLASS-mock B-mode recovery cross-check used by `analyses.covariance`. |
| Companion paper results | Papers III/IV/V; cited in narrative summary as scope context, not as findings. |
| Multi-field / tomographic extensions | Single tomographic bin is the paper's scope; no decision exposed for multi-bin. |
| High-order astrophysical B-modes | Below current sensitivity; mentioned only in narrative discussion if at all. |

### Completeness Status

- Headline claims: 5/5 represented.
- Figure targets: 9/9 represented (id divergence on `fig_pte_cl` and `fig_harmonic_config_cosebis_full` recorded above; equivalent structures present).
- Table targets: 2/2 represented.
- Quantitative results: 20/20 represented (codex iter-4 review added the §3.4 mock-validation rows: COSEBI mean B_n bias < 0.3σ and mock scatter comparable to fiducial error bars).
- Out-of-scope: 6/6 acknowledged.
