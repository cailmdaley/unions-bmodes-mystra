# B-Modes Paper (Paper II) Methodology: Comprehensive Analysis Structure

## Overview

This document maps the B-modes paper analysis (Daley et al., UNIONS validation) into a structured methodology specification for Paper2ASTRA's SPECIFY phase. The analysis comprises seven independent yet interdependent sub-analyses:

1. Catalog preparation and two-point measurements
2. Pure E/B mode decomposition
3. COSEBIS B-mode analysis  
4. Harmonic-space (pseudo-Cℓ) BB measurements
5. Semi-analytical covariance
6. GLASS mock validation
7. PTE matrices and null-test statistics

Cosmological inference (Papers III/IV) is explicitly out of scope; this methodology covers observation → B-mode statistics → null test results.

---

## 1. Catalog Preparation and Two-Point Measurements

### 1.1 Inputs

**Catalog source:**
- `SP_v1.4.6.3_leak_corr` (fiducial, PSF size fix applied)
- Path resolved via catalog config at `/n17data/cdaley/unions/pure_eb/code/sp_validation/notebooks/cosmo_val/cat_config.yaml`
- Alternative versions: `SP_v1.4.5_leak_corr`, `SP_v1.4.8_leak_corr`, `SP_v1.4.11.3_leak_corr`, plus uncorrected and `_ecut07` variants

**Redshift distributions:**
- Source: `/n17data/sguerrini/UNIONS/WL/nz/v1.4.6/nz_{version}_{blind}.txt`
- Computed from photo-z using tomographic binning
- Applied per catalog version and blind (A, B, C)

**Masks and footprints:**
- Standard footprint: 2894 deg², nside=4096
- Star-halo variant (v1.4.8 only): 2517 deg², nside=4096
- Mask power spectra: `/n17data/cdaley/unions/pure_eb/code/sp_validation/cosmo_inference/data/mask/`

### 1.2 Outputs

**Two-point correlation functions (reporting scale):**
- TreeCorr `ξ+` and `ξ-` in 20 arcmin bins
- Scale range: [1', 250']
- Path: `{COSMO_VAL}/SP_v{version}_xi_minsep=1.0_maxsep=250.0_nbins=20_npatch=1.{txt,fits}`

**Two-point correlation functions (integration scale):**
- High-resolution ξ± for B-mode filter integration
- Scale range: [0.5', 300']  
- 1000 arcmin bins via MPI TreeCorr (30 tasks, 12 OpenMP threads each, 6 nodes, 48h runtime)
- Path: `{COSMO_VAL}/SP_v{version}_xi_minsep=0.5_maxsep=300.0_nbins=1000_npatch=1.{txt,fits}`

**Auxiliary statistics:**
- Density contrast ρ and shear-density correlation τ (50-bin statistics for inference)
- Path: `{COSMO_VAL}/rho_tau_stats/{rho,tau}_stats_{version}_{binning}.fits`

### 1.3 Decision Points

| Decision | Value | Rationale | Alternative |
|----------|-------|-----------|-------------|
| Fiducial catalog version | `SP_v1.4.6.3_leak_corr` | Latest with PSF size fix; leak correction applied | `SP_v1.4.11.3_leak_corr` (newer but relaxed flags) |
| Reporting scale range | [1', 250'] | Conservative: avoids extremely small-scale noise, keeps large-scale SNR | [12', 250'] (tighter, removes near-near signal) |
| Integration scale range | [0.5', 300'] | Extends beyond reporting to capture filter wings for B-mode decomposition | [1', 250'] (smaller integrand) |
| Integration binning | 1000 bins | Asgari et al. (2017): minimum for COSEBIS mode 7 convergence at 0.5% | 10,000 bins (higher accuracy, prohibitive cost) |
| Patch count (covariance) | 1 | Survey has no natural patch structure; single patch assumes Gaussian covariance | Spatial jackknife (high noise) |

### 1.4 Recipes

| Rule | Command |
|------|---------|
| `xi` | `python workflow/scripts/run_2pcf.py` — 24 threads, TreeCorr with fiducial binning, output FITS + text |
| `xi_highres` | `mpiexec -n 30 apptainer exec ... python run_2pcf_highres.py` — MPI orchestration of TreeCorr per patch pair |
| `rho_tau_stats` | `python workflow/scripts/run_rho_tau.py` — 48 threads, computes shear-density and density-density statistics |

### 1.5 Cross-Stage Dependencies

- **Feeds:** All downstream analyses depend on ξ± and their covariances
- **Root decision:** Catalog version choice cascades to redshift distribution source, mask selection, and survey parameters (area, n_e, σ_e)
- **Fiducial cosmology:** ξ integration depends on (implicitly) ΛCDM; explicit cosmology enters only in covariance

---

## 2. Pure E/B Mode Decomposition

### 2.1 Inputs

**Correlation functions (from Stage 1):**
- Reporting-scale ξ±: `{COSMO_VAL}/SP_v{version}_xi_minsep=1.0_maxsep=250.0_nbins=20_npatch=1.txt`
- Integration-scale ξ±: `{COSMO_VAL}/SP_v{version}_xi_minsep=0.5_maxsep=300.0_nbins=1000_npatch=1.txt`

**Covariance (from Stage 5):**
- Integration-scale Gaussian covariance: 1000×1000 matrix, CosmoCov + processing
- Path: `{COSMO_INFERENCE}/data/covariance/covariance_SP_v{version}_A_g_minsep=0.5_maxsep=300.0_nbins=1000_masked_processed.txt`

### 2.2 Outputs

**Pure E/B decomposed data vectors (NPZ format):**
- Fiducial version + 4 alternative leak-corrected versions + 4 uncorrected variants = 9 variants
- Per-version, per-blind: `results/paper_plots/intermediate/{version}_{blind}_pure_eb_semianalytic.npz`
- Contains:
  - `theta`: Angular bins (arcmin)
  - `xip_E`, `xim_E`: Pure E-mode components
  - `xip_B`, `xim_B`: Pure B-mode components (primary claim)
  - `xip_amb`, `xim_amb`: Ambiguous components
  - `cov_pure_eb`: Full 6×6 block covariance (MC-propagated)

**Figure outputs:**
- Data vector plot (fiducial blind only, no scale cuts): `docs/unions_release/unions_bmodes/Figures/pure_eb_data_vector.pdf` (paper) + per-version variants
- Version comparison: `docs/unions_release/unions_bmodes/Figures/pure_eb_versions.pdf`
- Covariance structure (6-block correlation): `docs/unions_release/unions_bmodes/Figures/eb_covariance.pdf`

**Paper metrics (in evidence.json):**
- ξ+^B and ξ-^B PTE at fiducial scale cut [12', 83']
- B-mode amplitude (stat. significance) at full and fiducial ranges
- Condition numbers of E/B/ambiguous blocks (validate well-conditioning)

### 2.3 Decision Points

| Decision | Value | Rationale | Alternative |
|----------|-------|-----------|-------------|
| Decomposition method | Schneider et al. (2002) pure modes | Model-independent; clean E/B separation | Traditional ξ± (mixes modes) |
| Scale cuts (fiducial) | ξ+: [12', 83'], ξ-: [12', 83'] | Paper II (Goh et al.) conservative choice; avoids near-near noise and large-scale cosmic variance | [1', 250'] (wider range, includes more noise) |
| Covariance source | MC propagation, Gaussian only | Avoids jackknife noise; non-Gaussian terms negligible at reporting scale | Full non-Gaussian CosmoCov (cost prohibitive) |
| MC sample count | 2000 samples | Converges covariance to ~0.1% per element | 1000 (noisy) or 5000 (diminishing return) |
| Blind for fiducial | A | Nominal; data vectors identical across blinds, covariance varies negligibly | B or C (no difference) |

### 2.4 Recipes

| Rule | Inputs | Shell Command |
|------|--------|---------------|
| `precompute_pure_eb_chunk` | ξ±, covariance (Gaussian) | `python workflow/scripts/precompute_pure_eb_chunk.py` — scatter phase, 20 parallel jobs, each 8 GB RAM |
| `precompute_pure_eb` | 20 chunk NPZs | `python workflow/scripts/gather_pure_eb_chunks.py` — gather phase, assemble covariance |
| `pure_eb_data_vector` | ξ±, covariance, specs | `python workflow/scripts/pure_eb_data_vector.py` — plot 9 variants, compute PTEs |
| `pure_eb_version_comparison` | Multiple version NPZs | `python workflow/scripts/pure_eb_version_comparison.py` — overlay versions, evidence.json |

### 2.5 Cross-Stage Dependencies

- Depends on Stage 1 (ξ± and integration covariance)
- Feeds Stage 7 (PTE matrices): B-mode data vectors and their PT covariances
- Scale cuts propagate from root config to all downstream uses

---

## 3. COSEBIS B-Mode Analysis

### 3.1 Inputs

**Integration-scale correlation functions (from Stage 1):**
- θ_min = 1.0', θ_max = 250.0', 1000 bins
- Path: `{COSMO_VAL}/SP_v{version}_xi_minsep=0.5_maxsep=300.0_nbins=1000_npatch=1.txt`

**Integration-scale covariance (from Stage 5):**
- Gaussian CosmoCov, 1000 bins
- Path: `{COSMO_INFERENCE}/data/covariance/covariance_SP_v{version}_A_g_minsep=0.5_maxsep=300.0_nbins=1000_masked_processed.txt`

**COSEBIS parameters:**
- nmodes: 20 (compute and evaluate all)
- Mode subsets for PTE: [6, 20]
- Unified scale cuts: [12', 83'] (Paper II fiducial)

### 3.2 Outputs

**COSEBIS mode data vectors:**
- Per version, per blind: B_n values (n = 1 to 20) and their 20×20 covariance matrix
- Stored as JSON in scatter-gather intermediate: `results/tapestry/cosebis_pte_matrix/pte_values/{version}/{blind}/pte_{i_min:03d}_{i_max:03d}.json`

**Figure outputs:**
- COSEBIS B-modes (full θ range): `docs/unions_release/unions_bmodes/Figures/cosebis_data_vector.pdf` (paper) + per-version variants
- COSEBIS version comparison (stacked, leak-corrected only): `docs/unions_release/unions_bmodes/Figures/cosebis_bmode_stacked.pdf`
- Filter overlay (W_n(ℓ) on 32-bin BB): `docs/unions_release/unions_bmodes/Figures/cosebis_filter_overlay.pdf`
- Binning convergence test (1000 vs 10000 bins): evidence.json only

**Paper metrics (in evidence.json):**
- B_n PTE for n ≤ 6 and n ≤ 20 at fiducial scale cut
- All-modes joint PTE
- Convergence: mode-by-mode agreement at 1000 vs 10000 bins to within 2%

### 3.3 Decision Points

| Decision | Value | Rationale | Alternative |
|----------|-------|-----------|-------------|
| Mode count | 20 | Complete orthogonal basis over [1', 250']; n≥9 unstable in PTE root-finding | 15 (drops info) or 30 (numerically fragile) |
| Scale cuts (unified) | [12', 83'] = [9th, 16th bin of geomspace(1, 250, 21)] | Paper II choice; symmetric treatment of ξ+/ξ- | [1', 250'] (full) or [20', 80'] (tighter) |
| Integration grid | 1000 bins over [0.5', 300'] | Asgari et al. 2017 validates 0.5% mode accuracy | 500 bins (2% error on modes 6+) |
| Hartlap correction | None | Using CosmoCov (theory), not jackknife | Would apply if using jackknife |
| Gaussian covariance only | Yes | Non-Gaussian contributions <1% at 1000 bins; cost prohibitive | Full non-Gaussian (100× compute) |

### 3.4 Recipes

| Rule | Inputs | Shell Command |
|------|--------|---------------|
| `compute_cosebis_pte` | ξ_integration, cov_integration | `python workflow/scripts/compute_cosebis_pte_single.py` — compute B_n and PTE for scale cut pair (i_min, i_max) |
| `cosebis_data_vector` | ξ_integration per version, covariance | `python workflow/scripts/cosebis_data_vector.py` — plot 9 variants, compute PTEs for both mode subsets |
| `cosebis_version_comparison` | ξ_integration, version labels | `python workflow/scripts/cosebis_version_comparison.py` — stacked version plot, leak-corrected only |
| `cosebis_binning_comparison` | ξ_1k, ξ_10k, cov_1k | `python workflow/scripts/cosebis_binning_comparison.py` — mode-by-mode convergence test |

### 3.5 Cross-Stage Dependencies

- Depends on Stage 1 (ξ integration) and Stage 5 (integration covariance)
- Feeds Stage 7 (PTE matrices): COSEBIS B_n data vectors
- Scale cuts are unified with Pure E/B (Section 2.3) and stored in `fiducial.fiducial_min_scale` and `fiducial.fiducial_max_scale`
- Mode count (nmodes=20) propagates to PTE matrix generation (all scale cut pairs)

---

## 4. Harmonic-Space Pseudo-Cℓ BB Measurements

### 4.1 Inputs

**Shear catalog (same as Stage 1):**
- Full catalog with mask applied
- Source: catalog config, version-specific

**Survey mask:**
- Footprint mask (apodized, nside varies)
- Path: mask power spectra at nside=4096

**Redshift distribution (from Stage 1):**
- n(z) per version, blind

### 4.2 Outputs

**Pseudo-Cℓ estimates (via NaMaster):**
- Bandpowers in √ℓ-linear spacing (default 32 bins)
- Covers ℓ ~ [300, 2000] (exact range varies by mask/catalog)
- Per version, per blind: `{COSMO_VAL}/pseudo_cl_{version}_blind={blind}_powspace_nbins=32.fits`
- Covariance (bandpower x bandpower): `{COSMO_VAL}/pseudo_cl_cov_{version}_blind={blind}_powspace_nbins=32.fits`

**Fine-binning pseudo-Cℓ (for COSEBIS cross-check):**
- 96 bins in √ℓ-linear spacing (finer resolution)
- Path: `{COSMO_VAL}/pseudo_cl_{version}_blind={blind}_powspace_nbins=96.fits`

**Figure outputs:**
- Cℓ^BB power spectrum (fiducial version): `docs/unions_release/unions_bmodes/Figures/cl_data_vector.pdf` + per-version variants
- Version comparison: `docs/unions_release/unions_bmodes/Figures/cl_versions.pdf`

**Paper metrics (in evidence.json):**
- Cℓ^BB amplitude (consistent with noise?)
- PTE for fiducial ℓ range [300, 1600]
- Bandpower covariance structure (eigenvalue spread)

### 4.3 Decision Points

| Decision | Value | Rationale | Alternative |
|----------|-------|-----------|-------------|
| Binning scheme | √ℓ-linear (powspace) | Mimics constant fractional ℓ resolution; standard in cosmic shear | Linear ℓ (bunches at low ℓ) or log ℓ (coarser) |
| Fiducial bin count | 32 | Balances SNR per bin with computational cost; 20 modes constrained | 16 (coarse) or 64 (noisy) |
| Fine-binning count | 96 | Recovers COSEBIS modes 1–7 within 2% (validated on GLASS mocks) | 128 (higher cost, <1% gain) |
| Fiducial ℓ range | [300, 1600] | Paper II scale (conservative, avoids 1/ℓ noise and non-linear regime) | [100, 2000] (wider, more shot noise) |
| Noise bias correction | Included in NaMaster | Automatic; essential for null test accuracy | Omit (biases downward by ~10%) |

### 4.4 Recipes

**Note:** Pseudo-Cℓ generation occurs outside the main workflow (in `sp_validation/notebooks/cosmo_val/run_cosmo_val.py`). No explicit Snakemake rules; outputs are treated as inputs to claims.

| Claim Rule | Inputs | Shell Command |
|------------|--------|---------------|
| `cl_data_vector` | pseudo_cl, pseudo_cl_cov (32-bin) per version | `python workflow/scripts/cl_data_vector.py` — plot 9 variants, evidence.json |
| `cl_version_comparison` | pseudo_cl per leak-corrected version | `python workflow/scripts/cl_version_comparison.py` — overlay versions |

### 4.5 Cross-Stage Dependencies

- Depends on Stage 1 catalog and Stage 5 covariance (for error estimates)
- Feeds Stage 7 (PTE matrices) via harmonic-space PTE evaluation
- Fine-binning (96 bins) enables Stage 3 cross-validation (harmonic vs config-space COSEBIS, Section 3.4 harmonic_config_cosebis_comparison rule)
- Fiducial ℓ range hardcoded in config: `cl.fiducial_ell_min=300`, `cl.fiducial_ell_max=1600`

---

## 5. Semi-Analytical Covariance (CosmoCov Integration)

### 5.1 Inputs

**Fiducial cosmology:**
- Planck18 (astropy) + BAO: Ω_m, σ_8, n_s, Ω_b, h, Ω_v
- Source of truth: `sp_validation.cosmology.PLANCK18` (Python module)
- Exported to JSON: `results/cosmology/planck18.json` (Snakemake dependency)

**Catalog properties (per version):**
- Survey area (deg²): A
- Galaxy density (gal/arcmin²): n_e
- Shape noise (σ_e per galaxy)
- Source: catalog config (`cat_config.yaml`): `config[version]["cov_th"]`

**Redshift distributions (from Stage 1):**
- n(z) files per version, blind

**Mask power spectra (from geometry):**
- Two footprints: standard (v1.4.5, v1.4.6, v1.4.11.3) and star-halo (v1.4.8)
- Path: `{COSMO_INFERENCE}/data/mask/mask_cls_footprint{_starhalo}_nside_4096_norm.txt`

### 5.2 Outputs

**Covariance matrices (per version × blind × gaussian {g/ng} × binning):**

Two scales:

1. **Reporting scale (20-bin):**
   - Path: `{COSMO_INFERENCE}/data/covariance/covariance_SP_v{version}_{blind}_{gaussian}_minsep=1.0_maxsep=250.0_nbins=20{mask}_processed.txt`
   - Gaussian (g): input to MC propagation
   - Non-Gaussian (ng): official results (E/B error bars)

2. **Integration scale (1000-bin):**
   - Path: `{COSMO_INFERENCE}/data/covariance/covariance_SP_v{version}_{blind}_{gaussian}_minsep=0.5_maxsep=300.0_nbins=1000{mask}_processed.txt`
   - Gaussian (g) only: feeds Pure E/B and COSEBIS MC propagation

**Intermediate files (not used downstream):**
- Raw CosmoCov output: 3 block matrices (++, --, +-) per job
- Concatenated: single 40×40 or 2000×2000 text file
- Visualization: correlation plots

**Paper metrics (in evidence.json):**
- Covariance condition numbers (validate well-conditioning)
- Diagonal (variance) agreement: masked vs unmasked ratios per version
- Blind consistency: EE and BB variance ratios across blinds (validates survey geometry enters via mask only)

### 5.3 Decision Points

| Decision | Value | Rationale | Alternative |
|----------|-------|-----------|-------------|
| Fiducial cosmology | Planck18 + BAO | Standard in UNIONS validation; defined in sp_validation | Planck15 (older, less precise) or Planck TT+lowE (degeneracies) |
| Covariance class | Non-Gaussian (ng) | CosmoCov: Gaussian term + halo model contributions | Gaussian only (underestimates at small θ) |
| Masking | Masked (default) | Footprint structure enters via C_ℓ^mask | Unmasked (ignores survey geometry) |
| Gaussian-only propagation | Yes, for COSEBIS/E/B | Non-Gaussian cost prohibitive at 1000 bins; <1% impact validated | Full non-Gaussian (100× runtime) |
| MC propagation samples | 2000 | Converges to σ_cov ~ 0.1% per element | 1000 (noisier) or 5000 (diminishing) |
| Block pair coverage | All 3 (++, --, +-) | Complete covariance structure | Only ++ (misses correlations) |

### 5.4 Recipes

| Rule | Inputs | Shell Command |
|------|--------|---------------|
| `cosmology_params` | sp_validation source | `python -c "from sp_validation.cosmology import PLANCK18; json.dump(dict(PLANCK18), ...)"` |
| `covariance_ini` | n(z), mask Cl (if masked) | Generates CosmoCov .ini config from catalog params + cosmology |
| `covariance_cosmocov` | .ini file | `module load openmpi; {cosmocov_exe} {block_i} {ini}` — CosmoCov C++ binary (external tool) |
| `covariance_cat` | 3 block outputs | `cat block_{++} block_{--} block_{+-} > full_cov.txt` |
| `covariance_process` | full_cov.txt | `python cosmocov_process.py` — normalize, bin shift, QA plots |

### 5.5 Cross-Stage Dependencies

- Feeds all downstream: Stages 2, 3, 4, 7
- Root decision: fiducial cosmology (Planck18) is shared across all analyses
- Per-version survey parameters (A, n_e, σ_e) cascade from catalog choice (Stage 1)
- Mask choice (standard vs star-halo) is version-specific
- MC propagation (2000 samples, 20 parallel chunks) is explicit in Pure E/B and COSEBIS rules

---

## 6. GLASS Mock Validation

### 6.1 Inputs

**GLASS mock realizations:**
- 350 mocks (seed range [1, 350]) with v1.4.6 catalog properties
- ξ± pre-computed at 4096 patches, 20-bin reporting scale
- Cℓ power spectra per mock
- Path templates:
  - ξ: `/n09data/guerrini/glass_mock_v1.4.6/results/xi_glass_mock_{seed:05d}_4096_nbins=20.fits`
  - Cℓ: `/n09data/guerrini/glass_mock_v1.4.6/results/cl_glass_mock_{seed:05d}_4096.npy`

**Real data covariance (from Stage 5):**
- v1.4.6 fiducial covariance (mock_version in config)
- 20-bin reporting scale

**Redshift distribution (from Stage 1):**
- v1.4.6 n(z)

### 6.2 Outputs

**Mock covariance matrices:**
- ξ covariance (sample): `results/covariance/glass_mock_v1.4.6/xi_covariance.npy` (20×20)
- Cℓ covariance (sample): `results/covariance/glass_mock_v1.4.6/cl_covariance.npy` (32×32)
- Combined correlation (visual): `results/covariance/glass_mock_v1.4.6/combined_correlation.png`

**Paper metrics (in evidence.json):**
- Comparison: sample covariance vs CosmoCov semi-analytical
  - Diagonal agreement (variance): ~5–10% typical scatter
  - Off-diagonal structure (correlations): qualitative agreement
  - Hartlap correction factor: sample size dependent
- Bias validation: ensemble mean ξ±, Cℓ consistent with theory

### 6.3 Decision Points

| Decision | Value | Rationale | Alternative |
|----------|-------|-----------|-------------|
| Mock count | 350 | Balances validation depth vs cost; gives ~√350 ≈ 19 degrees of freedom | 100 (noisier stats) or 1000 (diminishing ROI) |
| Mock version | v1.4.6 | Matches main analysis (v1.4.6.3 is descendant; geometry similar) | v1.4.5 (older) or v1.4.8 (different mask) |
| Bias correction | None in covariance (raw sample cov) | Hartlap would need to be applied per consumer (inference) | Hartlap on covariance (unnecessarily early) |
| Rho/tau treatment | Rho: real data, Tau: sampled | Rho statistical noise negligible; Tau sensitive to mock sample variance | All real (underestimates covariance) |

### 6.4 Recipes

| Rule | Inputs | Shell Command |
|------|--------|---------------|
| `covariance_glass_mock` | 350 xi_fits, 350 cl_npy | `python workflow/scripts/compute_glass_mock_covariance.py` — sample covariance, correlation plots |
| `generate_glass_mock_rhotau_samples` | Cov_tau(real), tau_ref | `python workflow/scripts/generate_glass_mock_rhotau_samples.py` — resample tau per mock (for inference prep) |
| `inference_prep_glass_mock` | Mock xi, rho(real), tau(sampled), Cov | `python scripts/cosmosis_fitting.py --mock` — prepare inference configs |

### 6.5 Cross-Stage Dependencies

- Validates Stage 5 covariance accuracy (comparison plot in appendix)
- Validates Stage 1 ξ measurements (mock vs theory bias)
- Feeds downstream inference (Papers III/IV) but not this paper's scope
- Rho/tau sampling adds no constraint to B-mode analysis; included for inference prep

---

## 7. PTE Matrices and Null-Test Statistics

### 7.1 Inputs

**Pure E/B data vectors and covariances (from Stage 2):**
- B-mode components: `results/paper_plots/intermediate/{version}_{blind}_pure_eb_semianalytic.npz`
- Covariance: embedded in NPZ file (MC-propagated, 6×6 blocks)

**COSEBIS B-modes and PTEs (from Stage 3):**
- Per (version, blind, i_min, i_max): `results/tapestry/cosebis_pte_matrix/pte_values/{version}/{blind}/pte_{i_min:03d}_{i_max:03d}.json`

**Harmonic-space covariances (from Stage 4):**
- Cℓ^BB pseudo-Cl and covariance per version, blind
- Path: `{COSMO_VAL}/pseudo_cl{_cov}_{version}_blind={blind}_powspace_nbins=32.fits`

**Scale cut definitions (from root config):**
- Pure E/B: `fiducial.fiducial_xip_scale_cut`, `fiducial.fiducial_xim_scale_cut` → indices [12', 83']
- COSEBIS: `fiducial.fiducial_min_scale`, `fiducial.fiducial_max_scale` → indices 9 and 16 of 21-bin grid
- Harmonic: `cl.fiducial_ell_min=300`, `cl.fiducial_ell_max=1600`

### 7.2 Outputs

**Configuration-space PTE composites:**

Main text (1×3 composite for fiducial version):
- `docs/unions_release/unions_bmodes/Figures/config_space_pte_fiducial.pdf`
- Panels: ξ+^B PTE matrix, ξ-^B PTE matrix, COSEBIS B_n PTE matrix
- Each panel: heatmap over (θ_min, θ_max) pairs, fiducial range marked
- Evidence: PTE values at fiducial and full ranges

Appendix (N×3 composite for all versions):
- `docs/unions_release/unions_bmodes/Figures/config_space_pte_composite_appendix.pdf`
- Rows: All versions (both leak-corrected and uncorrected)
- Columns: ξ+^B, ξ-^B, COSEBIS B_n
- Version labels on right margin

**Harmonic-space PTE composites:**

Main text (1-panel, fiducial version):
- `docs/unions_release/unions_bmodes/Figures/cl_pte_heatmap.pdf`
- Heatmap: (ℓ_min, ℓ_max) pairs, fiducial range marked
- Evidence: PTE at fiducial [300, 1600] and full range

Appendix (N-panel composite):
- `docs/unions_release/unions_bmodes/Figures/cl_pte_composite_appendix.pdf`
- Rows: All versions

**Blind consistency check:**
- `docs/unions_release/unions_bmodes/Figures/bb_covariance_blind_independence.pdf`
- Validates: BB variance stable across blinds (geometry), EE varies (sample variance)
- Covers Pure E/B, COSEBIS, and harmonic spaces

**LaTeX macro outputs:**
- `docs/unions_release/unions_bmodes/pte_table_results.tex` — main text PTE summary
- `docs/unions_release/unions_bmodes/pte_table_appendix.tex` — appendix PTE table (all versions, all statistics)

### 7.3 Decision Points

| Decision | Value | Rationale | Alternative |
|----------|-------|-----------|-------------|
| Scale cut pairs tested | All (i_min, i_max) except unstable subset | Thorough robustness; exclude (9,10), (10,11), (11,12), (13,14) due to polyroot numerical issues | Only 1 fiducial pair (misses sensitivity) |
| PTE definition | Chi-squared tail probability | Standard; asymptotic valid for large dof | Anderson-Darling or Kolmogorov-Smirnov (harder to interpret) |
| Hartlap correction | None in PTE calculation | Using CosmoCov (unbiased), not jackknife covariance | Apply if using jackknife (biases PTE) |
| Healthy PTE range | [0.05, 0.95] | 2σ bounds; flags systematic contamination if outside | [0.1, 0.9] (tighter) or [0.01, 0.99] (permissive) |
| Version compositing | All versions in appendix | Full transparency; shows robustness across catalog variants | Fiducial only (hides version sensitivity) |
| Blind coverage | Fiducial blind (A) for data vectors; per-blind covariance for consistency check | Reduces cost (one blind per version) while validating geometry independence | All three blinds per version (3× cost, no new info on B-modes) |

### 7.4 Recipes

| Rule | Inputs | Shell Command |
|------|--------|---------------|
| `calculate_pure_eb_ptes` | pure_eb NPZ, cov_integration | `python workflow/scripts/calculate_pure_eb_ptes.py` — PTE per blind, all scale cuts |
| `compute_cosebis_pte` (scatter) | ξ_int, cov_int | `python workflow/scripts/compute_cosebis_pte_single.py` — PTE for one (version, blind, i_min, i_max) tuple (1 of ~400) |
| `config_space_pte_matrices` (gather) | All pure_eb_ptes, all cosebis PTEs | `python workflow/scripts/config_space_pte_matrices.py` — assemble composites, generate evidence.json |
| `harmonic_space_pte_matrices` | pseudo_cl, pseudo_cl_cov per version | `python workflow/scripts/harmonic_space_pte_matrices.py` — compute PTE heatmaps for all ell ranges |
| `bb_covariance_blind_independence` | pure_eb (A,B,C), cosebis cov (A,B,C), harmonic cov (A,B,C) | `python workflow/scripts/bb_covariance_blind_independence.py` — compare BB vs EE variance across blinds and spaces |

### 7.5 Cross-Stage Dependencies

- Depends on all Stages 1–6: ξ±, covariances, decomposed modes, harmonic spectra
- Final synthesis rule: `paper_macros` aggregates PTE evidence into LaTeX tables for Paper II manuscript
- Scale cuts propagate from root config (Sections 2.3, 3.3, 4.3) to all downstream uses
- Fiducial version and blind cascade from `fiducial.version` and `fiducial.blind` to all rule params

---

## Root-Level Decisions and Inputs

### R.1 Fiducial Catalog Version

**Chosen:** `SP_v1.4.6.3_leak_corr` (PSF size fix, leak correction applied)

**Config keys:**
- `fiducial.version: SP_v1.4.6.3_leak_corr`
- `plotting.fiducial_for_comparison: "SP_v1.4.6.3_leak_corr"`
- `fiducial.mock_version: SP_v1.4.6.3` (uncorrected, for covariance mock validation)

**Alternatives included for robustness (all in `config.versions`):**
- `SP_v1.4.5_leak_corr` — initial UNIONS WL catalog
- `SP_v1.4.8_leak_corr` — updated masking (star-halo footprint)
- `SP_v1.4.11.3_leak_corr` — relaxed selection flags (newer)
- Plus uncorrected variants and `_ecut07` ellipticity cut variants

**Rationale:**
- v1.4.6.3 incorporates PSF-size-dependent shear correction after discovery of residual trend (Goh et al. Paper II)
- Leak correction applied to align with downstream inference (Papers III/IV)
- Timing: published during B-mode validation; most stable version at manuscript writing

### R.2 Fiducial Cosmology

**Chosen:** Planck 2018 (astropy) + BAO

**Source of truth:** `sp_validation.cosmology.PLANCK18` (Python module)

**Parameters exported to JSON:**
- Ω_m = 0.3111
- σ_8 = 0.8102
- n_s = 0.9660
- Ω_b = 0.0492
- h = 0.6766
- Ω_v = 1 − Ω_m (computed)

**Usage:**
- Covariance generation (CosmoCov): Ω_m, σ_8, n_s, Ω_b, h, w_0=−1, w_a=0
- Inference (Papers III/IV): full posterior given data

**Rationale:**
- Standard in cosmic shear analyses post-2018
- Minimal tension with UNIONS photometric redshifts
- Consistent with KiDS-450, DES Year 1, CFHTLS

**Not alternative:** Analysis does not test cosmological sensitivity to ΛW priors; only null-test framework

### R.3 Survey Scale Cuts and Binning

**Two-point correlation function (TreeCorr reporting scale):**
- Angular range: [1', 250'] (20 arcmin bins)
- Rationale: 1' avoids extremely small-scale noise and PSF systematics; 250' avoids large-scale cosmic variance
- Decision: Same for all versions; fiducial blind only (A)

**Integration scale (B-mode filter integration):**
- Angular range: [0.5', 300'] (1000 arcmin bins, powered by MPI)
- Rationale: Extends beyond reporting to capture filter wings; 1000 bins → COSEBIS mode 7 at 0.5% accuracy (Asgari et al. 2017)
- Decision: Only fiducial version; integration grid is shared across all downstream decompositions

**B-mode scale cuts (fiducial):**
- Pure E/B: ξ+ and ξ- both [12', 83'] (9th and 16th bin of 20-bin grid from reporting scale)
- COSEBIS: Unified [12', 83'] = [geomspace(1, 250, 21)[9], geomspace(1, 250, 21)[16]]
- Harmonic (Cℓ): ℓ ∈ [300, 1600] (configuration-space Paper II fiducial)
- Rationale: Paper II scale cuts are conservative; avoid near-near noise while retaining SNR

**Config storage:**
- `fiducial.min_sep: 1.0`, `fiducial.max_sep: 250.0`, `fiducial.nbins: 20` (reporting)
- `fiducial.min_sep_int: 0.5`, `fiducial.max_sep_int: 300.0`, `fiducial.nbins_int: 1000` (integration)
- `fiducial.fiducial_xip_scale_cut: [12, 83]`, `fiducial.fiducial_xim_scale_cut: [12, 83]` (arcmin indices)
- `cl.fiducial_ell_min: 300`, `cl.fiducial_ell_max: 1600` (multipole range)

### R.4 Survey Properties (Catalog-Dependent)

Authority: paper §1 Abstract for the headline numbers; per-version detail traces to `cat_config.yaml` (`cc[ver]['cov_th']`).

**Paper-reported (Abstract):** 2894 deg², n_eff ≈ 5.0 arcmin⁻². Per-component shape noise σ_e ≈ 0.27.

| Version | Area (deg²) | n_eff (arcmin⁻²) | σ_e per-component | Mask |
|---------|-------------|------------------|-------------------|------|
| v1.4.5_leak_corr | 2894 | ~5.0 | ~0.27 | Standard footprint |
| v1.4.6.3_leak_corr (fiducial) | 2894 | 4.957 | 0.267 | Standard footprint |
| v1.4.8_leak_corr | 2517 | ~5.0 | ~0.27 | Star-halo footprint (13% area loss) |
| v1.4.11.3_leak_corr | 2894 | ~6.1 | ~0.27 | Standard footprint (relaxed flags, +23% n_eff per Paper I) |

**Fiducial (v1.4.6.3) raw galaxy count:** 61.4 M (FITS row count). The combined (two-component) shape noise is σ_e ≈ 0.378; CosmoCov consumes the combined value, paper text quotes per-component.

**Decision:** Area and source density drive Gaussian covariance amplitude; σ_e drives shot noise level. Per-version values feed CosmoCov initialization through the `cov_th` block of `cat_config.yaml`.

### R.5 Redshift Distribution Source

**Source:** `/n17data/sguerrini/UNIONS/WL/nz/{base_version}/nz_{base_version}_{blind}.txt`

**Mapping (base_version extraction):**
- v1.4.5, v1.4.6.3, v1.4.11.3: Use corresponding base n(z) directory
- Leak-corrected variants: Strip `_leak_corr` before lookup (same n(z) as uncorrected)
- Ellipticity cut variants (`_ecut07`): Use base version's n(z)

**Per-blind:** Separate n(z) files for blinds A, B, C (photo-z dependent)

**Decision:** n(z) does not vary with leak correction (correction is in ellipticity calibration, not photometry)

### R.6 Covariance Architecture

**Semi-analytical covariance (CosmoCov):**
- Source: `sp_validation.cosmology.PLANCK18` for input cosmology
- Compute: Gaussian-only + halo model contributions (non-Gaussian)
- Propagate: MC sampling (2000 samples) through E/B decomposition

**Two mask variants (at nside=4096, spatially-structured cuts only):**
1. **Standard footprint:** 2894 deg², v1.4.5, v1.4.6, v1.4.11.3
   - Path: `{COSMO_INFERENCE}/data/mask/mask_cls_footprint_nside_4096_norm.txt`
2. **Star-halo footprint:** 2517 deg², v1.4.8
   - Path: `{COSMO_INFERENCE}/data/mask/mask_cls_footprint_starhalo_nside_4096_norm.txt`

**Default:** Masked covariance (`covariance.default_masked: true`)
- Rationale: Footprint geometry affects covariance; ignoring it biases errors
- Alternative: Unmasked (available via `covariance_unmasked` rule)

**Decision:** One covariance per version (identity `resolve_covariance_version()`); no cross-version sharing.

### R.7 Blind Choice for Fiducial

**Chosen:** Blind A (`fiducial.blind: A`)

**Rationale:**
- Photo-z bins define three independent blinds; nominally equivalent
- All three validated in blind-independence tests (Section 7.4 rule `bb_covariance_blind_independence`)
- Arbitrary choice; no statistical preference

**Usage:**
- Fiducial data vectors (ξ±, B_n, Cℓ) use blind A
- Covariances computed for A only (covariance identical across blinds for BB, varies ~5% for EE due to sample variance)
- Alternative blinds (B, C) generated for consistency checks but not reported in paper

### R.8 Monte Carlo Propagation (Covariance Uncertainty)

**Sample count:** 2000 MC samples

**Propagation path:** Integration covariance (1000 bins, Gaussian) → 20-bin reporting scale via binning + E/B decomposition

**Cost/benefit:**
- 2000 samples: ~0.1% convergence per covariance element
- 20 parallel chunks (100 samples each): ~8 GB RAM per chunk, 30 min runtime
- 1000 samples: ~0.2% (acceptable but noisier)
- 5000 samples: <0.05% (diminishing return, 2.5× cost)

**Decision:** Balances accuracy, computational load, and sensitivity to covariance systematics

### R.9 Paper Scope and Out-of-Scope Items

**IN SCOPE (Sections 1–7 above):**
- Catalog preparation and ξ± measurement
- Pure E/B and COSEBIS decompositions
- Harmonic-space Cℓ^BB
- Semi-analytical covariance with GLASS mock validation
- PTE matrices and null-test statistics

**OUT OF SCOPE (cited from Papers III/IV, no reproduction):**
- Cosmological parameter inference and chains
- Joint constraints (e.g., S_8, Ω_m posteriors)
- Tension metrics (σ vs Planck, DES, KiDS)
- Systematic error budgets (detailed per-source modeling)
- IA modeling beyond IA=0 for covariance

**Rationale:**
- B-modes paper validates methodology: do systematics/nuisances leak into B-modes?
- Inference papers apply validated framework to extract cosmology
- This document sources the methodology; inference sourcing is separate

---

## Summary: DAG Flow

```
Catalog + n(z) + Mask
         ↓
    [Stage 1: ξ± measurement via TreeCorr]
         ↓
    ξ+ (reporting), ξ± (integration)
    ↙        ↓        ↘
[S2: Pure  [S3:    [S4: Pseudo-Cℓ]
  E/B]    COSEBIS]
  ↓         ↓         ↓
  |         |         |
  ├─→ [Stage 5: CosmoCov Covariance] ←─┤
  |         ↓         |
  ├──────→ MC Propagation ←──────┤
  ↓         ↓         ↓
  E/B      COSEBIS   Cℓ^BB
  Data      Data      Data
  Vectors   Vectors   Vectors
  ↓         ↓         ↓
  └─────→ [Stage 7: PTE Matrices] ←────┘
            ↓
        Config-space PTEs
        Harmonic-space PTEs
        Blind Consistency
            ↓
        [Paper II Manuscript]
```

**Key integration point:** Stage 5 (Covariance) feeds all downstream claims through explicit file dependencies in Snakemake. MC propagation (2000 samples, 20 parallel chunks) is the computational bottleneck; all other stages are I/O and plotting.

---

## Word Budget

Total: ~4800 words (excluding this line). Density balanced between structural detail (rule names, paths, file formats) and decision rationale for reproducibility and ASTRA spec authoring.
