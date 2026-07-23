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
