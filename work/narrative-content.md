# Narrative content saved from astra.yaml (removed during astra-spec 0.0.11 migration, 2026-07-20)

authors (removed from astra.yaml):
authors:
  - "C. Daley et al. (UNIONS Collaboration)"

narrative block (verbatim):

narrative:
  summary: |
    Reproduction of the UNIONS-3500 B-mode validation paper
    (arXiv to be assigned, "Paper II" of the five-paper UNIONS-3500
    cosmic-shear release). The paper validates that the
    [size-cut ShapePipe v1.4 catalog with PSF-leakage correction
    ](#decisions.fiducial_catalog) is internally consistent with
    zero B-mode signal across three E/B-separable estimators —
    pure-mode correlation functions, COSEBIs, and harmonic-space
    pseudo-Cℓ^BB — under the
    [adopted scale cuts](#decisions.fiducial_scale_cuts) of
    12–83 arcmin in configuration space and 300 < ℓ < 1600 in
    harmonic space. [The size-cut catalog is the only one of the
    four ShapePipe variants to pass all three statistics
    simultaneously at fiducial cuts](#findings.size_cut_unique_pass),
    and [its full-range COSEBI failure traces to detector-scale
    repeating additive bias](#findings.cosebis_full_range_failure).
    Cosmological inference is out of scope, and is the subject of
    the configuration-space companion (Goh et al. 2026, Paper III)
    and the harmonic-space companion (Guerrini et al. 2026,
    Paper IV).

    Reproduction status: this is a real, executed reproduction.
    Every paper-facing figure, table, and finding under
    `results/baseline/` and `analyses/<sub>/results/baseline/` is
    materialised by `lc run` executing the recipe's actual compute
    inside the sp_validation apptainer container — no imported or
    spoofed artifacts. The fiducial catalogue
    (`SP_v1.4.6.3_leak_corr`) chain traces end-to-end through lc; the
    version-comparison figures read lc outputs for the fiducial row
    and the canonical CosmoStat tree for the other seven catalog
    variants, which are not re-run through lc. The headline COSEBI
    B-mode PTEs reproduce byte-exact (0.78 at the adopted 12'–83'
    cut, 1.37×10⁻⁵ at full range), and the PTE table matches the
    paper's LaTeX macros to the last displayed digit for 47 of 48
    values — the exception being the full-range Cℓ^BB PTE
    (reproduced 0.14 vs the paper's 0.13), which differs by a single
    display-rounding step traced to a ~1% difference in the CosmoCov
    covariance vintage used for that full-range diagnostic; the
    fiducial-cut Cℓ^BB PTE (0.30) matches exactly. Deferred as
    out-of-scope diagnostics, not re-run:
    `two_point.xi_integration_convergence` and
    `cosebis.cosebis_binning_convergence` (10,000-bin grids),
    `covariance.cov_glass_*` (GLASS validation), and the
    `catalog.*` path-resolution helpers.
  inputs: |
    Three input classes. The
    [four ShapePipe v1.4 catalog variants](#inputs.shapepipe_v14_catalogs)
    (initial / size-cut / masked / relaxed-flags, each with a
    leakage-corrected and an uncorrected variant) drive every
    downstream measurement; the spec exposes these as a
    [single decision](#decisions.fiducial_catalog) so the
    multiverse explores all four. The
    [SOM-calibrated photo-z distributions](#inputs.redshift_distributions)
    are per-blind (A, B, C) and per base-version, but invariant
    under leakage correction. Two
    [survey footprint masks](#inputs.footprint_masks) at
    HEALPix nside=4096 supply the geometry that enters the
    Gaussian covariance: the standard 2894 deg² footprint
    (v1.4.5 / v1.4.6.3 / v1.4.11.3) and a star-halo-masked
    2517 deg² variant (v1.4.8 only). The
    [Planck-2018 ΛCDM cosmology](#inputs.planck18_cosmology)
    enters via `sp_validation.cosmology.PLANCK18` and is the
    sole cosmology assumed for covariance theory.
  methods: |
    The pipeline walks the snakemake DAG bottom-up.
    [Catalog preparation](#analyses.catalog) materialises the
    four variants under their per-version masks and packages
    them with their photo-z and survey properties for downstream
    consumption.
    [Two-point measurement](#analyses.two_point) runs TreeCorr in
    two configurations sharing the same catalog: a 20-bin
    reporting grid over [1', 250'] for pure-E/B and Cℓ data
    vectors, and a 1000-bin integration grid over
    [0.5', 300'] feeding COSEBI mode integration and Monte-Carlo
    propagation of the covariance.

    Three E/B-separable estimators decompose the same ξ± into
    B-mode statistics: [pure-mode ξ_±^B](#analyses.pure_eb)
    via the [finite-interval pure-mode construction of Schneider
    et al. (2022)](#prior_insights.pure_mode_decomposition_unique_on_finite_interval),
    [COSEBI Bₙ amplitudes](#analyses.cosebis) via the
    [logarithmic complete-orthogonal basis of Schneider, Eifler
    & Krause (2010)](#prior_insights.cosebis_complete_eb_separation)
    integrated to mode 20 on the
    [1000-bin grid](#decisions.integration_grid) whose
    resolution is cross-checked against a one-shot 10,000-bin
    run (paper §3.3), and
    [pseudo-Cℓ^BB bandpowers](#analyses.cl_bb) via the
    [catalog-based pseudo-Cℓ estimator of Wolz, Alonso &
    Nicola 2025](#prior_insights.wolz25_catalog_based_pseudo_cl)
    operating on the same catalog — evaluating spin-2
    transforms at source positions avoids the pixelization
    artefacts that map-based pseudo-Cℓ develops at the UNIONS
    footprint angular resolution, and the
    [analytic shot-noise correction](#cl_bb.decisions.noise_bias_correction)
    is part of that estimator choice. A
    [Fourier band-power family-mate (HybridEB; Becker & Rozo
    2016)](#prior_insights.hybrideb_filter_optimisation_alternative)
    is not adopted here — Jefferson et al. (2025) report that
    DES-Y3 and HSC-Y3 fail B-mode null tests under HybridEB
    while passing under their original estimators, motivating
    the multi-statistic strategy this paper applies to UNIONS. The shared
    [covariance sub-analysis](#analyses.covariance) supplies
    CosmoCov matrices on the 20-bin reporting grid (in
    Gaussian and Gaussian + non-Gaussian flavours) and the
    1000-bin integration grid (Gaussian-only). All B-mode
    null tests use
    [Gaussian-only covariance](#covariance.decisions.covariance_class) —
    non-Gaussian connected-trispectrum and super-sample terms
    are retained only as a reporting-grid diagnostic — and
    the integration-grid matrix is
    [propagated through the pure-mode and COSEBI filters by 2000
    Monte-Carlo samples](#decisions.mc_samples). Estimator
    recovery and the propagated scatter are cross-checked
    against
    [350 GLASS lognormal mock realisations matched to v1.4.6
    catalog properties](#inputs.glass_mock_realisations)
    (Hervas-Peters et al. 2026, "Paper V").

    The decisive movement is in the choice of the
    [fiducial catalog](#decisions.fiducial_catalog): blinded
    null tests on the initial v1.4.5 catalog showed PTE failures
    in ξ_-^B and Cℓ^BB; tightening the size cut from
    r_h,gal/r_h,psf > 0.5 to > 0.707 (the size-cut variant) was
    the change that eliminated those failures while keeping the
    per-component shape noise stable, even though the masked
    variant fails in all three pure-mode statistics. The
    [adopted scale cuts](#decisions.fiducial_scale_cuts) follow
    the same logic: they were chosen for stable PTE regions
    across statistics, not by minimising any single PTE.

    [Null tests](#analyses.null_tests) are the spec's product.
    PTE matrices over (θ_min, θ_max) and (ℓ_min, ℓ_max) for each
    catalog variant make the scale-cut sensitivity inspectable;
    the [PTE summary table](#outputs.pte_table_results) reports
    fiducial and full-range values for the size-cut catalog,
    and the [appendix table](#outputs.pte_table_appendix)
    tabulates all four variants. The
    [healthy PTE band](#decisions.healthy_pte_range) is the 2σ
    interval [0.05, 0.95]; the [chi-squared
    definition](#decisions.pte_definition) is the standard tail
    probability. The
    [theory-covariance Hartlap correction](#decisions.hartlap_correction)
    is `none` for COSEBIs and Cℓ^BB because their inputs are
    semi-analytical CosmoCov and iNKA matrices; the
    [pure-mode χ² stage](#pure_eb.decisions.hartlap_application)
    carries its own Hartlap factor on the MC-propagated
    covariance (≤ 2 % at N_samples = 2000, N_obs ≤ 40).
  findings: |
    The headline result is structural rather than numerical:
    [the size-cut catalog is the only ShapePipe v1.4 variant to
    pass all three E/B-separable statistics at the adopted scale
    cuts](#findings.size_cut_unique_pass), with a minimum PTE of
    0.18 across pure-mode, COSEBI, and harmonic-space tests.
    Both a tighter (masked) and a looser (relaxed-flags)
    selection fail at least one of the three, demonstrating
    that [multi-statistic concurrence is the operational test
    of catalog systematics](#findings.multi_statistic_consensus)
    rather than any single estimator's PTE.

    The most informative failure mode is COSEBIs over the full
    angular range: [ξ_+^B and ξ_-^B integrated against the
    T_±,n filters at θ ∈ [1', 250'] show a PTE of 1.37×10⁻⁵ in
    n ≤ 6 modes for the size-cut catalog](#findings.cosebis_full_range_failure),
    persisting even after restriction to the size-cut sample.
    The pattern is consistent with a [repeating additive shear
    bias on instrument-fixed CCD scales
    ](#prior_insights.repeating_additive_bias_signature) —
    Asgari et al. (2019a) traced a similar feature in CFHTLenS
    to detector-fixed offsets, and the lower bound of the
    adopted [12', 83'] cut is set just above the 6–14 arcmin
    MegaCam CCD footprint for that reason.
    The harmonic-space and full-θ pure-mode PTEs do not see
    this feature with comparable significance because their
    filter functions weight the affected scales differently,
    sharpening the case that
    [no single estimator suffices](#findings.multi_statistic_consensus).
    [Computing COSEBIs from harmonic-space pseudo-Cℓ
    bandpowers](#findings.filter_sets_sensitivity) demonstrates
    this directly: the same Cℓ^BB that passes its own null test
    fails when re-projected through the W_n(ℓ) filters, so the
    differing significance is set by filter functions rather
    than by the angular-versus-harmonic basis (paper Figure 8,
    §5).
  outputs: |
    The paper-facing products are two PTE tables and a small
    set of composite figures. The
    [main-text PTE summary](#outputs.pte_table_results) is
    Table 1 — six PTE values for the size-cut catalog at
    fiducial and full-range cuts, generated from the null-test
    sub-analysis. The
    [appendix PTE table](#outputs.pte_table_appendix)
    tabulates the same six statistics across all four catalog
    variants.

    Three single-version data-vector figures display the
    fiducial (size-cut) measurement in each statistic family:
    the [pure E/B decomposition of ξ_±
    ](#outputs.fig_pure_eb_decomposition) (paper Figure 1), the
    [COSEBI B_n / σ_n spectrum
    ](#outputs.fig_cosebis_fiducial) (paper Figure 3) which
    carries the visible evidence behind the
    [full-range COSEBI failure](#findings.cosebis_full_range_failure),
    and the [harmonic-space Cℓ^BB and Cℓ^EB
    ](#outputs.fig_cl_fiducial) (paper Figure 4). The
    [B-mode versions composite](#outputs.fig_pure_eb_versions)
    overlays the four catalog variants for each statistic
    family. PTE-matrix figures expose scale-cut sensitivity:
    the [config-space PTE composite
    ](#outputs.fig_pte_heatmaps) and the
    [harmonic-space Cℓ^BB PTE composite
    ](#outputs.fig_cl_pte_heatmap), the latter showing that
    the source workflow uses
    [adopted blind A](#decisions.fiducial_blind) for the
    reported data vectors. The
    [configuration-vs-harmonic COSEBI cross-check
    ](#outputs.fig_harmonic_config_cosebis) (paper Figure 8)
    establishes that the two computation paths yield identical
    mode amplitudes within filter-function differences. All
    figures and tables consume sub-analysis outputs by anchored
    reference; paper macros in `paper/claims_macros.tex` carry
    the rendered numerical values.


# Sub-analysis narrative blocks (removed 2026-07-31)

The seven `analyses/*/astra.yaml` files carried `narrative:` blocks long after the root
narrative was removed. astra-spec 0.0.11 rejects the field outright, but root validation
never recursed into sub-analyses, so they survived unvalidated. Archived here, verbatim.


## analyses/catalog/astra.yaml

```yaml
narrative:
  summary: |
    Materialises one of the four
    [ShapePipe v1.4 catalog variants](#../decisions.fiducial_catalog)
    along with its photo-z distribution per blind, footprint
    mask, and survey properties (area, n_eff, σ_e). All
    downstream sub-analyses consume these products via the
    [packaged shear catalog](#outputs.shear_catalog), the
    [redshift distribution](#outputs.redshift_distribution),
    and the [survey-property summary](#outputs.survey_properties);
    no sub-analysis re-resolves catalog paths from `cat_config.yaml`.
  inputs: |
    The four variants enter via the parent's
    [ShapePipe v1.4 catalogs](#inputs.shapepipe_v14_catalogs);
    their photo-z files via the parent's
    [redshift distributions](#inputs.redshift_distributions); the
    standard / star-halo footprint masks via the parent's
    [footprint masks](#inputs.footprint_masks).
  methods: |
    The variant is selected by the parent's
    [fiducial_catalog](#../decisions.fiducial_catalog) decision
    and combined with the
    [PSF-leakage correction](#../decisions.leakage_correction) toggle
    inherited from the parent. Path resolution goes through
    `cat_config.yaml` in the source repo's
    `code/sp_validation/notebooks/cosmo_val/` directory:
    `cc[ver]['shear']['path']` is the FITS catalog,
    `cc[ver]['shear']['redshift_path']` is the n(z), and
    `cc[ver]['mask']` is the appropriate footprint mask.
    Each base version has a single n(z) regardless of leakage
    correction, so leak-corrected variants strip `_leak_corr`
    before lookup. Loading is read-only — no shape recalibration
    or reweighting happens here.
  outputs: |
    The packaged
    [shear catalog](#outputs.shear_catalog),
    [redshift distribution](#outputs.redshift_distribution), and
    [footprint mask power spectrum](#outputs.mask_power_spectrum)
    are passed to
    [two_point](#../analyses.two_point) and
    [covariance](#../analyses.covariance); the
    [survey-property summary](#outputs.survey_properties)
    feeds CosmoCov's `cov_th` block (n_eff, σ_e, area).
```


## analyses/cl_bb/astra.yaml

```yaml
narrative:
  summary: |
    Pseudo-Cℓ^EB and Cℓ^BB bandpowers measured directly on the
    discrete shear catalog via the
    [Wolz et al. 2025 catalog-based
    estimator](#decisions.estimator_choice) — spin-2 spherical
    harmonics evaluated at source positions through the `ducc`
    library, with mask-induced mode coupling deconvolved by
    NaMaster (MASTER algorithm; Hivon et al. 2002; Alonso et
    al. 2019). The
    [paper data vector](#outputs.cl_bandpowers_reporting) is on
    the [32-bandpower √ℓ-linear reporting grid
    ](#decisions.bandpower_binning) over 8 ≤ ℓ ≤ 2048; the
    [adopted scale cut](#../decisions.fiducial_scale_cuts) is
    300 < ℓ < 1600. A
    [96-bandpower fine variant](#outputs.cl_bandpowers_fine)
    is materialised for the [COSEBI harmonic-space validation
    path](#../cosebis.outputs.cosebis_harmonic_modes); both
    share the same NaMaster invocation pipeline.
  inputs: |
    The harmonic-space estimator consumes the
    [shear catalog](#inputs.shear_catalog),
    [redshift distribution](#inputs.redshift_distribution),
    and [footprint mask power spectrum
    ](#inputs.mask_power_spectrum) materialised by the
    [catalog sub-analysis](#../analyses.catalog). The
    bandpower covariance is supplied by the
    [covariance sub-analysis](#../analyses.covariance) via its
    [iNKA Gaussian Cℓ covariance
    ](#../analyses.covariance) — the same NaMaster invocation
    that produces the data-vector bandpowers also produces
    its iNKA covariance, but the spec splits the production:
    the bandpower estimator lives here, and the bandpower
    covariance lives in the covariance sub-analysis (where it
    co-locates with the configuration-space CosmoCov
    matrices).
  methods: |
    The Wolz et al. (2025) catalog-based estimator represents
    the survey mask as delta functions at source positions and
    computes the spin-2 transforms on this irregular grid using
    the `ducc` library (Reinecke & Seljebotn 2023), avoiding
    the pixelization artifacts that complicate map-based
    pseudo-Cℓ implementations on the UNIONS footprint. The
    same metacalibration weights and ellipticities that enter
    the [TreeCorr ξ_± measurement](#../analyses.two_point) are
    used here (paper Eq. 9). NaMaster decouples the measured
    bandpowers via the mode-coupling matrix induced by the
    [footprint mask power spectrum
    ](#inputs.mask_power_spectrum) and applies an
    [analytic shot-noise correction
    ](#decisions.noise_bias_correction) following the
    implementation in Paper IV; no Monte-Carlo bias subtraction
    is needed because the Wolz et al. construction yields
    unbiased Cℓ̂^{αβ} estimators directly.

    Two binning schemes coexist on the same NaMaster call. The
    [32-bandpower reporting grid
    ](#decisions.bandpower_binning) is the paper data-vector
    binning and the input to the
    [Cℓ^BB null test](#outputs.cl_pte_per_cut); the
    [96-bandpower fine grid](#outputs.cl_bandpowers_fine) is
    the COSEBI cross-check binning, dense enough to resolve
    the W_n(ℓ) filter oscillations through n ≤ 6 and chosen
    by [cosebis](#../analyses.cosebis) for the harmonic-path
    Eq. 19. Both grids span the same 8 ≤ ℓ ≤ 2048 native
    range; the
    [√ℓ-linear spacing](#decisions.bandpower_binning) (`power
    = 0.5` in the source `pseudo_cl` rule) gives constant
    fractional ℓ-resolution, the standard cosmic-shear
    convention.

    Cℓ^BB null tests use the iNKA Gaussian
    bandpower covariance from the
    [covariance sub-analysis](#../analyses.covariance), which
    is theory-derived (no Hartlap correction needed, paralleling
    the [root hartlap_correction decision
    ](#../decisions.hartlap_correction)) and validated by Paper
    IV against 10,000 Gaussian simulations. The
    [PTE-per-(ℓ_min, ℓ_max)
    matrix](#outputs.cl_pte_per_cut) drives the
    harmonic-space PTE composite figure that
    [null_tests](#../analyses.null_tests) gathers,
    and follows the same
    [blind-A paper strategy](#decisions.pte_blind_strategy) as
    [pure_eb](#../analyses.pure_eb): paper-facing PTEs use
    blind A only, with per-blind generation reserved for the
    BB-covariance blind-independence cross-check.
  outputs: |
    [cl_bandpowers_reporting](#outputs.cl_bandpowers_reporting)
    is the primary data vector consumed by the
    [paper Cℓ figure rule](#outputs.cl_data_vector_figure)
    that produces
    [paper Figure 4](#../outputs.fig_cl_fiducial) and by the
    [version comparison](#outputs.cl_version_comparison_figure)
    that feeds the
    [paper Figure 5 stack](#../outputs.fig_pure_eb_versions).
    [cl_bandpowers_fine](#outputs.cl_bandpowers_fine) is the
    sole consumer of the 96-bin grid: the
    [cosebis sub-analysis](#../analyses.cosebis) reads it for
    the harmonic-path Eq. 19 evaluation. The
    [PTE matrix](#outputs.cl_pte_per_cut) feeds the
    [harmonic Cℓ^BB PTE composite
    ](#../outputs.fig_cl_pte_heatmap) and the corresponding
    appendix table row when
    [null_tests](#../analyses.null_tests) gathers it.
```


## analyses/cosebis/astra.yaml

```yaml
narrative:
  summary: |
    COSEBI E_n / B_n mode amplitudes for n = 1 … 20 computed
    from ξ_±(θ) via the T_{±n}(θ) filter integration of Eq. 16
    in `paper/unions_bmodes.tex` §3.2, evaluated on the
    [1000-bin integration grid](#../decisions.integration_grid)
    over θ ∈ [θ_min, θ_max]. The
    [load-bearing PTE is the n ≤ 6 subset
    ](#decisions.mode_subsets), where the cosmological E-mode
    information saturates and where the
    [full-range first-mode excess](#../findings.cosebis_full_range_failure)
    of >4σ collapses below 1σ once the
    [adopted 12'–83' cut is applied
    ](#../decisions.fiducial_scale_cuts). The [10,000-bin
    convergence variant](#outputs.cosebis_binning_convergence)
    confirms B_n stability against integration-grid resolution
    (paper §3.3 convergence check), and the
    [harmonic-space evaluation
    ](#outputs.cosebis_harmonic_modes) of Eq. 19 cross-checks
    the angular-space modes against the same filters expressed
    in ℓ.
  inputs: |
    The COSEBI integrands are the
    [1000-bin ξ_±](#inputs.xi_integration) materialised by
    [two_point](#../analyses.two_point) and the corresponding
    [Gaussian covariance](#inputs.cov_integration_g) from the
    [covariance sub-analysis](#../analyses.covariance). The
    one-shot
    [10,000-bin convergence run
    ](#inputs.xi_integration_convergence) feeds the
    binning-comparison output; it is not used by paper PTEs.
    The harmonic cross-check consumes the
    [fine-binned pseudo-Cℓ](#inputs.cl_bandpowers_fine) at 96
    √ℓ-linear bandpowers from the cl_bb sub-analysis.
  methods: |
    The angular-space mode amplitudes follow paper Eq. 16 with
    logarithmic T_{±n}(θ) filter functions whose orthogonality
    on [θ_min, θ_max] removes the
    [ambiguous-mode content](#../analyses.pure_eb)
    that pure-mode ξ_±^{E/B} carry as a separate component. We
    evaluate the integral on the
    [1000-bin grid](#../decisions.integration_grid) for paper
    PTEs and on the 10,000-bin grid for the §3.3
    [convergence cross-check](#outputs.cosebis_binning_convergence).
    Mode amplitudes are computed for all n up to the
    [paper's mode-count cap of 20](#decisions.mode_count);
    PTEs are reported for the
    [two mode subsets](#decisions.mode_subsets) that the paper
    table carries (n ≤ 6 cosmologically informative, n ≤ 20
    full transparency).

    The B_n covariance is the linear transform
    `C_{B_n B_m} = T^T C_ξ T` of the
    [Gaussian integration-grid covariance
    ](#../analyses.covariance), exploiting the closed-form
    derivative of the COSEBI mode integral with respect to ξ_±.
    Equivalently — and as a self-consistency check — the same
    covariance is also produced by passing the
    [2000 Monte-Carlo samples](#../decisions.mc_samples) of
    ξ_±(θ) through the COSEBI integrals and computing the
    sample covariance of the resulting B_n vectors. The
    [PTE scatter-gather](#outputs.cosebis_pte_per_cut)
    materialises one JSON per (version, blind, i_min, i_max)
    tuple from the
    [stable-pair scale-cut grid](#decisions.scale_cut_pair_grid)
    surrounding the
    [adopted 12'–83' fiducial cut](#../decisions.fiducial_scale_cuts);
    the gather step folds them into the
    [config-space PTE composite](#../outputs.fig_pte_heatmaps).

    The
    [harmonic-space cross-check](#outputs.cosebis_harmonic_modes)
    realises Eq. 19, summing the same FFT-log Fourier-dual
    filters W_n(ℓ) against the
    [96 √ℓ-linear pseudo-Cℓ bandpowers
    ](#decisions.harmonic_validation_binning); this bin count
    balances W_n(ℓ) oscillation resolution against
    mode-coupling matrix conditioning, and the GLASS
    validation in the
    [covariance sub-analysis](#../analyses.covariance)
    confirms agreement to within 2% for n ≤ 6, the modes
    carrying nearly all cosmological information (paper §3.2,
    Figure 8).
  outputs: |
    [cosebis_modes_data](#outputs.cosebis_modes_data) is the
    primary B_n / E_n data vector consumed by the paper's
    [single-version COSEBI plot](#../outputs.fig_cosebis_fiducial)
    and by the
    [version-comparison composite](#outputs.cosebis_version_comparison)
    feeding the
    [paper Figure 5 stack](#../outputs.fig_pure_eb_versions).
    [cosebis_pte_per_cut](#outputs.cosebis_pte_per_cut)
    discharges the COSEBI rows of the
    [paper PTE table](#../outputs.pte_table_appendix)
    after the [null_tests sub-analysis](#../analyses.null_tests)
    composes the multi-statistic PTE matrices. The
    [10,000-bin binning convergence
    artifact](#outputs.cosebis_binning_convergence) backs the
    paper's §3.3 convergence statement, and the
    [harmonic-config comparison
    ](#outputs.cosebis_harmonic_modes) feeds the cross-check
    figure [(paper Figure 8)](#../outputs.fig_harmonic_config_cosebis).
```


## analyses/covariance/astra.yaml

```yaml
narrative:
  summary: |
    Builds the [CosmoCov](#decisions.covariance_engine) covariance
    matrices that all three B-mode statistics consume, on two
    grids and in [Gaussian-only and Gaussian + non-Gaussian
    flavours](#decisions.covariance_class). The reporting-grid
    20-bin matrix is produced both with and without the
    non-Gaussian connected-trispectrum / super-sample terms; the
    integration-grid 1000-bin matrix is Gaussian-only because
    [the paper's null-test PTEs use Gaussian-only covariance
    throughout](#decisions.covariance_class) for tractability
    on the fine grid and for parity with Paper IV's iNKA
    covariance in harmonic space. Validation against
    [350 GLASS lognormal mock realisations](#inputs.glass_mock_realisations)
    is computed in parallel and is the evidence behind the
    paper's mock-recovery statement; the validation sample
    covariances do not feed paper PTEs.
  inputs: |
    The covariance build consumes the
    [redshift distribution](#inputs.redshift_distribution),
    [footprint mask power spectrum](#inputs.mask_power_spectrum),
    and per-version
    [survey-property summary](#inputs.survey_properties)
    materialised by the [catalog sub-analysis](#../analyses.catalog),
    together with the parent's
    [Planck-2018 fiducial cosmology](#../inputs.planck18_cosmology)
    serialised once via the
    [cosmology snapshot](#outputs.cosmology_snapshot). The
    GLASS validation arm consumes the parent's
    [350 GLASS mock realisations](#inputs.glass_mock_realisations)
    pre-computed at v1.4.6 catalog properties; its outputs are
    sample covariances on the same 20-bin reporting grid the
    [two-point sub-analysis](#../analyses.two_point) uses, not
    new mock generations.
  methods: |
    CosmoCov is invoked once per (catalog version, blind, grid,
    [Gaussian/non-Gaussian flavour](#decisions.covariance_class),
    [masked / unmasked](#decisions.mask_application)) tuple. Each
    invocation begins with a generated `.ini` config file emitted by
    `rule covariance_ini`, which serialises the
    [Planck-2018](#../decisions.fiducial_cosmology) cosmology,
    the per-version `(area, n_eff, σ_e)` triple from
    [survey_properties](#inputs.survey_properties), the grid
    `(min_sep, max_sep, n_bins)`, the IA-off settings
    (A_IA = 0), and either the
    [standard 2894 deg² or star-halo 2517 deg² mask power
    spectrum](#inputs.mask_power_spectrum) for `_masked` runs.
    The CosmoCov C++ binary is then invoked by
    `rule covariance_cosmocov` for each of the three
    [shear-shear block pairs](#decisions.block_pairs), producing
    [block files](#outputs.cov_blocks) `++`, `--`, `+-`
    independently (one block per Slurm task); the three
    [blocks](#outputs.cov_blocks) are concatenated by `rule
    covariance_cat` and post-processed by
    `cosmocov_process.py` (`rule covariance_process`) into the
    `_processed.txt` matrix consumed downstream, the
    `_processed_g.txt` Gaussian-only extract used by the
    [pure-mode and COSEBI Monte-Carlo
    propagation](#../decisions.mc_samples), and a
    [QA structure plot](#outputs.cov_processed_qa) used to
    inspect diagonal variance and condition-number behaviour
    during covariance review.

    Two sub-trees of paper consumption hang off the processed
    matrices. The
    [reporting-grid 20-bin matrices](#outputs.cov_reporting_g)
    feed the harmonic-grid pseudo-Cℓ pipeline directly (as
    Knox-formula bandpower variances) and feed the
    [pure E/B](#../analyses.pure_eb) and
    [COSEBI](#../analyses.cosebis) data-vector covariances
    when reporting-scale error bars are needed at the same
    binning as the data vector itself. The
    [integration-grid 1000-bin Gaussian matrix
    ](#outputs.cov_integration_g) is the Gaussian draw source
    for the [2000 Monte-Carlo
    realisations](#../decisions.mc_samples) that propagate
    through the pure-mode and COSEBI filters into the
    reporting-grid B-mode covariance. The pure-mode case
    requires this MC propagation because no closed-form linear
    transform exists from ξ_±-space to ξ_±^{E,B,amb}-space on a
    finite interval (Schneider et al. 2010 / 2022); the COSEBI
    case has the closed-form `C_{B_n B_m} = T^T C_ξ T` linear
    transform but is also computed via the same MC pipeline
    for self-consistency and to recover the cross-statistic
    covariances used by the joint ξ_tot^B test.

    The GLASS validation arm computes the
    [ξ_± sample covariance](#outputs.cov_glass_xi_reporting)
    and the
    [pseudo-Cℓ sample covariance](#outputs.cov_glass_cl_reporting)
    from the 350 v1.4.6-matched GLASS realisations on the
    20-bin reporting grid; this is the
    [validation-only role](#decisions.glass_validation_strategy)
    of the mocks. Because the mocks contain no PSF-related
    systematics, they test whether the estimators generate
    spurious B-modes on a noise-dominated input; the paper
    reports B-modes consistent with zero across all three
    statistics, COSEBI mean B_n bias below 0.3σ per mode over a
    100-realisation subset, and mock scatter comparable to the
    adopted fiducial error bars. These sample covariances are
    not used to weight the paper PTEs — the paper PTEs are
    CosmoCov-weighted, so Hartlap correction does not enter on
    the covariance side.
  outputs: |
    Downstream consumers anchor here.
    [pure_eb](#../analyses.pure_eb) consumes
    [cov_integration_g](#outputs.cov_integration_g) for MC
    propagation and [cov_reporting_ng](#outputs.cov_reporting_ng)
    for cross-checks against the propagated matrix.
    [cosebis](#../analyses.cosebis) consumes
    [cov_integration_g](#outputs.cov_integration_g) for the
    same MC pipeline, and the closed-form linear-transform
    covariance is recomputed there from
    [cov_integration_g](#outputs.cov_integration_g) when
    needed. [cl_bb](#../analyses.cl_bb) is supplied via the
    [pseudo-Cℓ covariance](#outputs.cl_cov_reporting) which is
    produced inline with the harmonic-space pseudo-Cℓ
    estimator in the same NaMaster invocation as the data
    vector itself; the harmonic CosmoCov entries are exposed
    here as a cross-check on the iNKA covariance Paper IV
    uses for its own cosmological inference. The
    [null_tests](#../analyses.null_tests) sub-analysis is
    where the GLASS sample covariance is consumed
    qualitatively as the mock-pipeline check.
```


## analyses/null_tests/astra.yaml

```yaml
narrative:
  summary: |
    Gather stage. Composes per-statistic per-(scale-cut) PTEs
    from [pure_eb](#../analyses.pure_eb),
    [cosebis](#../analyses.cosebis), and
    [cl_bb](#../analyses.cl_bb) into the
    [paper PTE table](#../outputs.pte_table_appendix) (all
    four catalog variants × six statistics, at
    [adopted](#../decisions.fiducial_scale_cuts) and
    [full-range](#../decisions.fiducial_scale_cuts) cuts), and
    two PTE-matrix composite figures (one for configuration
    space, one for harmonic space). Also runs the
    [BB-covariance blind-independence cross-check
    ](#outputs.bb_covariance_blind_independence_evidence) that
    underwrites the source workflow's
    [blind-A fiducial choice](#../decisions.fiducial_blind) for the
    reported data vectors. Cosmological inference is out of
    scope for this paper; this sub-analysis exists to discharge
    the paper's measurement-side findings.
  inputs: |
    The PTE composition consumes scale-cut grids of B-mode PTEs
    from the three statistic families: the
    [pure-mode PTE matrix](#inputs.pure_eb_pte_per_cut), the
    [COSEBI PTE scatter](#inputs.cosebis_pte_per_cut), and the
    [Cℓ^BB PTE matrix](#inputs.cl_pte_per_cut), each
    materialised by the corresponding sub-analysis. The
    blind-independence cross-check pulls
    per-(blind ∈ {A,B,C}) MC-propagated pure-mode covariances
    from [pure_eb](#inputs.pure_eb_per_blind_covariances), the
    [integration covariance](#inputs.cov_integration_per_blind)
    from [covariance](#../analyses.covariance), and per-blind
    [pseudo-Cℓ bandpower covariances](#inputs.cl_cov_per_blind);
    these are produced at the v1.4.6
    [mock_version](#decisions.bb_blind_check_version) since
    blind independence is a property of survey geometry, not
    catalog version.
  methods: |
    The PTE composition is a thin gather: it reads per-statistic
    PTEs over all scale-cut pairs and arranges them into the
    paper's [main-text](#outputs.pte_summary_evidence) and
    [appendix](#outputs.pte_summary_evidence) tables. The
    [healthy PTE band](#../decisions.healthy_pte_range) of
    [0.05, 0.95] is the colour-coding threshold;
    [PTE definition](#../decisions.pte_definition) is the
    χ²-tail probability, with no Hartlap correction at the χ²
    stage for COSEBIs and Cℓ^BB (the inputs are theory
    covariances) and the
    [pure_eb Hartlap correction](#../analyses.pure_eb)
    applied internally before the PTE values land here.

    The [config-space PTE composite
    ](#outputs.config_space_pte_evidence) renders three columns
    (ξ_+^B, ξ_-^B, COSEBI B_n) over (i_min, i_max) bandpower
    pairs of `geomspace(1', 250', 21)`, with rows for each
    catalog variant and the
    [adopted 12'–83' fiducial cut](#../decisions.fiducial_scale_cuts)
    marked as a black square; the three statistic panels
    share that θ-grid so the adopted-cut cell sits at the
    same (i, j) position in every panel
    ([scale-cut grid alignment
    ](#decisions.scale_cut_pair_grid_alignment)). The
    [harmonic-space PTE composite
    ](#outputs.harmonic_space_pte_evidence) renders the same
    grid over (ℓ_min, ℓ_max) bandpower pairs, with the
    [adopted 300 < ℓ < 1600 cut
    ](#../decisions.fiducial_scale_cuts) marked. Both
    composites are produced as one paper-facing PDF (fiducial
    catalog only) plus an appendix PDF that tiles all four
    catalog variants.

    The
    [BB-covariance blind-independence
    check](#outputs.bb_covariance_blind_independence_evidence)
    is a separate one-shot analysis. It generates per-blind
    pure-mode, COSEBI, and harmonic covariances at the
    [v1.4.6 mock version](#decisions.bb_blind_check_version)
    and compares the diagonal BB-block variance across blinds
    A, B, C against the analogous EE-block variance. Geometry
    determines the BB variance under the null hypothesis, so
    the check predicts B-mode covariance ratios closer to unity
    than E-mode ratios, with MC sampling noise limiting that
    statement in the pure-mode branch. The figure is evidence
    for the source workflow's
    [blind-A fiducial](#../decisions.fiducial_blind) being
    representative; the paper does not include it in the
    main-text figures, but it underwrites the validity of the
    blind-A PTE table.

    Three of the paper's measurement-side findings are
    evidenced from the same JSON:
    the [size-cut catalog's unique simultaneous pass
    ](#../findings.size_cut_unique_pass) and the
    [multi-statistic concurrence requirement
    ](#../findings.multi_statistic_consensus)
    both read the [appendix PTE table evidence
    ](#outputs.pte_summary_evidence), while the
    [full-range COSEBI failure
    ](#../findings.cosebis_full_range_failure) reads its COSEBI
    rows at the full-range cut. The fourth paper finding,
    [filter-function sensitivity
    ](#../findings.filter_sets_sensitivity), is evidenced
    upstream in the [COSEBI sub-analysis](#../analyses.cosebis)
    from the
    [harmonic-vs-config COSEBI cross-check
    ](#../cosebis.outputs.cosebis_harmonic_modes), not from this
    PTE gather.
  outputs: |
    The [PTE summary JSON](#outputs.pte_summary_evidence) is
    the evidence backing the
    [paper PTE table](#../outputs.pte_table_appendix); the
    `paper_macros` rule reads it to render the LaTeX `\num{}`
    macros in `paper/claims_macros.tex`.
    The [configuration-space PTE evidence
    ](#outputs.config_space_pte_evidence) feeds the
    [paper config-space PTE
    composite](#../outputs.fig_pte_heatmaps), and
    the [harmonic-space PTE evidence
    ](#outputs.harmonic_space_pte_evidence) feeds the
    [paper Cℓ^BB PTE
    composite](#../outputs.fig_cl_pte_heatmap). The
    [BB blind-independence figure
    ](#outputs.bb_covariance_blind_independence_evidence)
    is not in the paper but is referenced by the
    [workflow blind-A justification](#../decisions.fiducial_blind).
```


## analyses/pure_eb/astra.yaml

```yaml
narrative:
  summary: |
    Schneider et al. (2022) decomposition of ξ_±(θ) into pure
    E-mode, B-mode, and ambiguous components on the
    [1000-bin integration grid](#../decisions.integration_grid),
    extended to 0.5'–300' to suppress edge effects in the
    boundary functions S_±(θ) and V_±(θ) (paper §3.1, Eq. 13).
    The pure-mode B-mode vectors are weighted by a covariance
    [Monte-Carlo-propagated](#decisions.mc_propagation) from the
    [Gaussian integration-grid CosmoCov matrix
    ](#../analyses.covariance) through 2000 draws of the same
    integral transforms; the resulting ξ_±^B covariance is
    then [Hartlap-corrected](#decisions.hartlap_application) at
    the χ² stage, the operationally important Hartlap
    application in this paper.
  inputs: |
    The decomposition consumes the
    [reporting-grid ξ_±](#inputs.xi_reporting) for the 20-bin
    output binning the paper figure carries, the
    [1000-bin integration ξ_±](#inputs.xi_integration) for the
    integral kernels of paper Eq. 13, and the
    [Gaussian integration-grid covariance
    ](#inputs.cov_integration_g) materialised by the
    [covariance sub-analysis](#../analyses.covariance) as the
    Gaussian draw source.
  methods: |
    The pipeline is a scatter-gather around the Monte-Carlo
    propagation. `precompute_pure_eb_chunk` runs 20 parallel
    Slurm tasks per the
    [chunked propagation strategy](#decisions.mc_propagation),
    each drawing 100 of the 2000
    [MC samples](#../decisions.mc_samples) of ξ_±(θ) from the
    [Gaussian integration covariance
    ](#inputs.cov_integration_g) and applying the
    [Schneider 2022 kernels](#decisions.decomposition_method)
    to each draw via the
    [cosmo_numba implementation](#decisions.implementation)
    (Guinot et al. 2026), producing
    per-chunk MC artifacts;
    `precompute_pure_eb` gathers the 20 chunks
    into the per-(version, blind) NPZ holding ξ_±^E, ξ_±^B,
    ξ_±^amb data vectors, the diagonal σ used for the
    [paper data-vector figure](#../outputs.fig_pure_eb_decomposition),
    and the 6-block (E+, E−, B+, B−, amb+, amb−) covariance
    visualised in the
    [covariance correlation matrix](#outputs.pure_eb_covariance_figure).

    The B-mode null tests apply Eq. 22 of the paper:
    χ²_B = d_B^T C_BB^−1 d_B with the
    [Hartlap finite-sample correction](#decisions.hartlap_application)
    factor (N_samples − N_obs − 2)/(N_samples − 1), which
    reduces the inverse covariance by at most 2% for
    N_samples = 2000 and N_obs ≤ 40, because the MC-propagated
    pure-mode covariance is a sample-derived estimator (paper
    §3.4). This is the operationally meaningful Hartlap
    application in this paper — the
    [root hartlap_correction decision](#../decisions.hartlap_correction)
    refers to the CosmoCov input matrix itself, which is
    theory-derived and does not need correction; the MC
    output covariance does. The
    [PTE matrix](#outputs.pure_eb_pte_per_cut) is computed
    over the same scale-cut grid that
    [cosebis](#../analyses.cosebis) uses; the joint ξ_tot^B
    test concatenates the ξ_+^B and ξ_-^B data vectors and
    weights against the full block matrix from the same NPZ.
    Per the
    [PTE blind strategy](#decisions.pte_blind_strategy),
    paper-facing PTEs are evaluated at blind A only; the
    `{blind}` wildcard is reserved for the BB-covariance
    blind-independence cross-check downstream.
  outputs: |
    [pure_eb_semianalytic_data](#outputs.pure_eb_semianalytic_data)
    is the per-(version, blind) primary artifact consumed by
    every downstream pure-mode plot or PTE: it carries the
    decomposed data vectors and the MC covariance. The
    [data-vector figure](#outputs.pure_eb_data_vector_figure)
    materialises the
    [paper Figure 1](#../outputs.fig_pure_eb_decomposition);
    the [covariance correlation
    matrix](#outputs.pure_eb_covariance_figure) materialises
    the paper Figure 2 (`eb_covariance.pdf`); the
    [version-comparison figure
    ](#outputs.pure_eb_version_comparison_figure) feeds the
    pure-mode row of the
    [paper Figure 5 stack](#../outputs.fig_pure_eb_versions);
    [pure_eb_pte_per_cut](#outputs.pure_eb_pte_per_cut)
    discharges the ξ_+^B / ξ_-^B / ξ_tot^B rows of the
    [paper PTE table](#../outputs.pte_table_appendix) once the
    [null_tests sub-analysis](#../analyses.null_tests) gathers
    them.
```


## analyses/two_point/astra.yaml

```yaml
narrative:
  summary: |
    TreeCorr ξ_±(θ) measurement of the
    [packaged shear catalog](#../catalog.outputs.shear_catalog) on
    two grids: a 20-bin
    [reporting grid](#decisions.reporting_binning) over [1', 250']
    used by the pure-mode and Cℓ data vectors, and a
    [1000-bin integration grid](#../decisions.integration_grid)
    over [0.5', 300'] feeding COSEBI mode integration and
    Monte-Carlo propagation of the covariance. Both grids run on
    the same catalog; the integration run is single-threaded
    sequential at 1000 bins or MPI-distributed at 10,000 for the
    one-shot convergence test.
  inputs: |
    The two-point measurement consumes the
    [packaged shear catalog](#inputs.shear_catalog) materialised
    by the [catalog sub-analysis](#../analyses.catalog), which
    has already resolved per-version PSF-leakage correction and
    masking. No other input enters at this stage; the redshift
    distribution and footprint mask flow downstream to
    [covariance](#../analyses.covariance) directly from
    [catalog](#../analyses.catalog).
  methods: |
    Both grids invoke the same TreeCorr `GGCorrelation`
    estimator weighted by the `metacalibration` shear weights
    `w_i` from Paper I, accumulating tangential and cross
    products in the angular bins defined by
    [reporting_binning](#decisions.reporting_binning) and the
    parent's
    [integration_grid](#../decisions.integration_grid)
    decisions. The reporting run uses a single TreeCorr process
    on 24 OpenMP threads; the 1000-bin integration run uses the
    same process with a longer wall time, while the 10,000-bin
    convergence variant distributes patch-pair correlations
    over 30 MPI tasks (paper §3.3 convergence check).

    The patch decomposition is held at
    [npatch = 1](#decisions.patch_count): the survey has no
    natural patch structure for spatial jackknife, and the
    semi-analytical CosmoCov covariance described in
    [covariance](#../analyses.covariance) replaces empirical
    jackknife. The [120-patch
    alternative](#decisions.patch_count.options.npatch_120) was
    considered but does not enter any paper PTE; it is exposed
    for cross-check sensitivities only.
  outputs: |
    [Reporting-scale ξ_±](#outputs.xi_reporting) is the input to
    the [pure E/B sub-analysis](#../analyses.pure_eb) and the
    [Cℓ^BB sub-analysis](#../analyses.cl_bb).
    [Integration-scale ξ_±](#outputs.xi_integration) is the
    input to the
    [pure E/B sub-analysis](#../analyses.pure_eb) again (for
    the integration-grid Gaussian draws via the
    [covariance sub-analysis](#../analyses.covariance)) and to
    the [COSEBI sub-analysis](#../analyses.cosebis). The
    [10,000-bin convergence variant
    ](#outputs.xi_integration_convergence) is materialised once
    for the size-cut catalog and is the evidence behind the
    paper's §3.3 convergence statement.
```
