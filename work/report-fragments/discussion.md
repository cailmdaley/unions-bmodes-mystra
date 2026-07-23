<!-- Fragment: discussion (TeX lines 525–599) — \section{Discussion} \label{sec:discussion} -->
<!-- Macro bindings per work/transpile-map.md. Paper-series macros (\paperconfig=III,
     \paperharmonic=IV, \papercatalog=I, \papersims=V) rendered as literal paper numbers. -->

## Discussion

The pure-mode correlation functions, COSEBIs, and harmonic-space power spectra presented in this paper weight angular scales and treat ambiguous modes differently, so they do not respond identically to the same contamination.
When they disagree, we treat that as evidence of contamination in the sample even if each statistic is not equally sensitive to it, and require the adopted catalog and cuts to pass in all three frameworks.

At full range, all four catalogs show a COSEBI first-mode excess exceeding $4\,\sigma$ and oscillatory structure across all twenty modes.
Because COSEBIs are orthogonal over a finite angular range ({astra}`cosebis`), a feature compact in angle spreads across the mode spectrum; the full-range oscillations are qualitatively consistent with contamination at CCD angular scales ([Results](#results), {astra:ref}`outputs.fig_cosebis_fiducial`).
In the pure-mode and harmonic-space data vectors, no comparably localized feature is apparent: $\xi_+^{\mathrm{B}}$ shows mild excess at the smallest and largest separations, and $C_\ell^{BB}$ has a few ${\gtrsim}2\,\sigma$ outliers at low $\ell$ and elevated power in the highest multipole bins, but neither displays a sharp angular signature that would point unambiguously to CCD-scale contamination.
<!-- ebthetaXipMin → decisions.fiducial_scale_cuts (map §2f); note: renders the adopted-option label, not the bare integer 12 — metric/table bridge needed to emit \num{12} -->
MegaCam's 40 CCDs span $6$--$14$ arcmin individually, and our adopted lower cut of {astra:value}`decisions.fiducial_scale_cuts` arcmin excludes most pair separations at these scales; at the adopted scale cuts, the COSEBI features are suppressed and all four catalogs pass.
Pure-mode correlation functions and $C_\ell^{BB}$ on the other hand pass at the adopted cuts only for the size-cut catalog and for the size-cut and masked versions, respectively.
These equivocal results lead us to argue against relying on any single null test; taken together, we choose the size-cut catalog and our adopted scale cuts as the combination least affected by systematic contamination for cosmological inference in Papers III and IV.

<!-- Fig 8 (fig:harmonic_config_cosebis_full): outputs.fig_harmonic_config_cosebis — BOUND, but note angular-range mismatch: root recipe renders the *fiducial* variant while the paper includes the *full* variant (map §3) -->
:::{astra} outputs.fig_harmonic_config_cosebis
:caption: COSEBI $E_n$ (top row) and $B_n/\sigma_n$ (bottom row) for the fiducial catalog, computed from configuration-space $\xi_\pm$ (filled circles) and harmonic-space bandpowers (open squares). Gray shading marks modes $n \leq 6$, where the cosmological information is concentrated and we confirm consistency between the two paths on GLASS simulations.
:::

Although the masked catalog passes in harmonic space and COSEBIs ($n \leq 6$), this ostensibly more conservative choice fails all three pure-mode tests.
Excising stellar halos introduces irregular gaps in the survey footprint, and because $B$-modes are non-local, a more complex mask geometry may increase the fraction of modes classified as ambiguous, leaving the remaining pure-mode estimates less stable; alternatively, removing regions of elevated shape noise near bright stars could increase the significance of contamination elsewhere.
We cannot distinguish these explanations with the present data.
The initial and relaxed-flags catalogs fail in $\xi_-^{\mathrm{B}}$ and harmonic space, consistent with poorly resolved or deblended objects contributing to the contamination; the relaxed-flags catalog also shifts the total $\xi_+$ systematically higher than the other versions beyond 16 arcmin ({astra:ref}`outputs.fig_pure_eb_versions`), suggesting that admitting these objects affects the signal itself, not only the $B$-mode diagnostics.

Configuration-space and harmonic-space estimators respond differently to the same contamination.
Pure-mode correlation functions are localized in angular separation, so cuts that remove CCD-scale separations directly reduce the contamination.
In multipole space, a feature narrow in angle maps onto a broad range of $\ell$, so it cannot be isolated by cutting a localized multipole range.
COSEBIs apply filter functions to the two-point statistics that concentrate cosmological information into the first few modes; contamination at a characteristic angular scale spreads across the mode spectrum rather than isolating in a single mode.
Pure-mode and COSEBI filters also discard ambiguous modes ({astra}`pure_eb`) that cannot be uniquely classified as $E$ or $B$ on a finite angular interval.
The harmonic-space estimator used here does not deproject ambiguous modes; it assigns all power to $C_\ell^{EE}$ or $C_\ell^{BB}$, and any residual mask-induced $E$-to-$B$ leakage after mode-coupling deconvolution is absorbed into $C_\ell^{BB}$.

We show directly that $C_\ell^{BB}$ and COSEBIs have different sensitivity to systematic features, even when computed from the same data, by transforming the harmonic-space bandpowers into COSEBIs ([](#eq-cosebi-harmonic)).
The $C_\ell^{BB}$ power spectrum passes the null test at our fiducial multipole range, yet COSEBIs computed from the same $C_\ell^{BB}$ fail at full angular range---because the COSEBI filter functions $W_n(\ell)$ weight the power spectrum differently, concentrating sensitivity on the scales where the contamination resides.
This provides empirical support for the argument of {cite:t}`asgari.etal19a` that a null test in one representation cannot validate another; {cite:t}`jefferson.etal25` reached the same conclusion across DES-Y3, KiDS-1000, and HSC-Y3, finding that tomographic bins passing pseudo-$C_\ell$ $B$-mode tests can fail when reanalyzed with HybridEB.
The converse also holds: a $B$-mode failure in any statistic should raise concern for cosmological inference in all of them, since the underlying contamination is present in the data regardless of how it is projected.

{astra:ref}`Fig. %s <outputs.fig_harmonic_config_cosebis>` compares COSEBIs computed from configuration-space $\xi_\pm$ and from harmonic-space bandpowers.
Both $E_n$ and $B_n$ agree across the two computation paths for all twenty modes, with the tightest correspondence in the first six modes validated on GLASS simulations.
<!-- harmCosebisPteSixThreeFid → cosebis.outputs.cosebis_harmonic_modes → evidence.json harmonic_b_mode_ptes.SP_v1.4.6.3_leak_corr.pte (0.60); cfgCosebisPteSixThreeFid → same output → config_b_mode_ptes.….pte (0.78). Nested-dict figure evidence — value bridge pending (map §0). -->
At fiducial scale cuts, both paths yield $B$-modes consistent with zero (PTE $=$ {astra:value}`cosebis.outputs.cosebis_harmonic_modes` harmonic, {astra:value}`cosebis.outputs.cosebis_harmonic_modes` configuration; $n \leq 6$).
<!-- TODO-GAP: harmCosebisPteSixThreeFull = \num{1.61e-05} — harmonic full-range COSEBI not materialized (root fig_harmonic_config_cosebis recipe uses --angular-range fiducial; map §2d). TeX value kept below. -->
<!-- cfgCosebisPteSixThreeFull → null_tests.outputs.config_space_pte_evidence → versions.SP_v1.4.6.3_leak_corr.cosebis_stats.pte_at_full_range (1.37e-5) -->
At full range, both detect the same structure (PTE $= 1.61\times10^{-5}$ <!-- TODO-GAP harmCosebisPteSixThreeFull --> and {astra:value}`null_tests.outputs.config_space_pte_evidence`).
This confirms that the sensitivity to systematic contamination is set by the filter functions, not by the basis in which the two-point function is measured.
Throughout this analysis, we use Gaussian-only covariance, which underestimates variances for the UNIONS geometry ({astra}`covariance`); the reported PTEs are therefore conservative.

{cite:t}`asgari.etal19a` identified an oscillatory COSEBI $B$-mode pattern in multiple Stage-III surveys---CFHTLenS, KiDS-450, and DES-SV---tracing it to a repeating additive shear bias (i.e. a spatially varying $c$-term) from detector-level effects fixed in focal-plane coordinates.
Its appearance across three Stage-III surveys points to a detector-level phenomenon rather than a pipeline-specific artifact, a point especially relevant for UNIONS because CFHTLenS also used MegaCam on CFHT {cite:p}`heymans.etal12; guinot.etal22`.
The dither strategy, however, differs substantially.
CFHTLenS and KiDS used small dithers designed to bridge chip gaps, so the focal-plane pattern maps almost directly onto the sky.
UNIONS dithers by ${\sim}1/3$ of the MegaCam field of view {cite:p}`gwyn.etal25`, larger than a single CCD, so each sky position averages the additive bias from several different CCDs across exposures.
{cite:t}`asgari.etal19a` note that the $B$-mode signature of a repeating additive pattern depends on the dithering strategy, pointing to DES-SV's half-field dithers as one reason its pattern differs from the small-dither surveys.
The larger UNIONS dithers may similarly attenuate the amplitude of the repeating pattern, but are unlikely to fully eliminate it: the effective pattern is the per-CCD bias convolved with the dither geometry, potentially spreading the characteristic angular scale beyond the $6$--$14$ arcmin CCD footprint.

Several mechanisms can create additive shear structure at CCD angular scales.
The $\rho$-statistics and $\xi_{\mathrm{sys}}$ diagnostics in Paper I show scale-dependent PSF residuals with a feature near the ${\sim}10$ arcmin MegaCam CCD scale, consistent with per-CCD PSF contributions; astrometric residuals may also play a role.
Leakage correction has negligible impact on $B$-mode PTEs at the adopted scales ([Results](#results)).
However, the adopted cuts themselves exclude the large scales where PSF leakage is strongest, so this does not rule out leakage as a contributor to the full-range signal.
At the fiducial scales, residual contamination from higher-order PSF moments or other effects beyond a second-order leakage model may still contribute {cite:p}`zhang.etal23`.

Paper V identifies a pipeline-specific source arising from the improper propagation of exposure offsets within the shape measurement pipeline.
In the current multi-exposure fitting procedure, objects are effectively assumed to be centered on the pixel grid of each individual exposure.
However, a consistent treatment would require recentering the model using the World Coordinate System (WCS), anchored to the detection position defined at the tile level.
In areas where CCDs overlap, the relative pixel offsets between exposures can induce a preferred orientation in the inferred galaxy shapes, corresponding to the vector displacements between exposure pixel grids.
As a result, a coherent additive shear pattern may be imprinted across the CCD overlap regions, on angular scales between the per-CCD footprint ($\sim 6$--$14$ arcmin) and the UNIONS dither spacing ($1/3$ MegaCam FoV, or $\sim 19$ arcmin).
The observed $B$-mode contamination concentrates at $\theta_\mathrm{min} \lesssim 12$ arcmin ({astra:ref}`outputs.fig_pte_heatmaps`), close to the per-CCD scale but below the dither offsets.
Whether this mechanism is responsible for the observed $B$ modes remains an open question.
Even so, a contamination localized near the CCD scale would be cut away in configuration space, while its Hankel transform to harmonic space would spread across multipoles, beyond any single $\ell$ cut.
If a comparable additive contribution were present in the $E$ modes, it would lift the harmonic-space $S_8$ relative to configuration space, and could contribute to the mild offset seen between Papers III and IV.
Including configuration-space scales down to $5$ arcmin does move its $S_8$ toward the harmonic value (Paper III).
While a fraction of these effects may be absorbed by standard multiplicative bias corrections, we mitigate any remaining impact on cosmological constraints by cutting the majority of CCD scales; the adopted cuts pass without fine-tuning.

An alternative to excluding contaminated scales would be to model the $B$-mode signal as a bias and subtract it from the total measurement.
However, the relationship between $E$- and $B$-mode contamination depends on the physical mechanism; {cite:t}`asgari.etal19a` showed that different systematics produce markedly different $E$/$B$ ratios, so the $B$-mode amplitude alone does not determine the $E$-mode contamination.
We therefore treat $B$-modes as a diagnostic rather than a correction: scale cuts that remove $B$-mode contamination also remove any corresponding contamination of the total signal, without assuming a model for the relationship between them.
