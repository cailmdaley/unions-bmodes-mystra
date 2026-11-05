---
title: "UNIONS-3500 Weak Lensing: II. B-mode validation for cosmic shear"
---

+++ {"part": "abstract"}

At Stage-III sensitivities, cosmic shear $B$ modes unambiguously indicate systematic contamination and are often used to inform data selection and scale cuts for cosmological inference. We validate $B$ modes for the Ultraviolet Near-Infrared Optical Northern Survey (UNIONS)-3500 using three $E$/$B$-separable statistics: pure-mode correlation functions $\xi_\pm^{\mathrm{B}}(\theta)$, Complete Orthogonal Sets of $E$/$B$-mode Integrals (COSEBIs) $B_n$, and harmonic-space power spectra $C_\ell^{BB}$. We perform null tests by computing probability-to-exceed (PTE) values for each statistic and comparing the minimum PTE to a threshold of {astra:value col=display where="macro=mockMinPteThreshold"}`outputs.paper_claims`, calibrated on systematics-free simulations. $B$-mode detections on an initial catalog led us to repeat these tests for different sample selections and scale ranges before settling on fiducial cuts of 0.707 in galaxy-to-PSF size ratio, {astra:value col=display where="macro=ebthetaXipMin"}`outputs.paper_claims`--{astra:value col=display where="macro=ebthetaXipMax"}`outputs.paper_claims` arcmin in configuration space, and $300 < \ell < 1600$ in harmonic space. At these cuts the fiducial catalog passes with a minimum PTE of {astra:value col=display where="macro=fiducialMinPte"}`outputs.paper_claims`; {astra:value col=display where="macro=mockPGlobal"}`outputs.paper_claims` of the systematics-free simulations have a lower PTE. When using the full range of angular scales between 1 and 250 arcmin, we measure an oscillatory COSEBI $B$-mode pattern that could be explained by a repeating additive shear bias, a detector-level effect seen across multiple Stage-III surveys including CFHTLenS, which used the same MegaCam camera; scale cuts that exclude the charge-coupled device (CCD) angular scale suppress this pattern. Although the three $B$-mode statistics probe the same two-point shear field, scale cuts applied in one basis do not map exactly onto cuts in another, because their filter functions weight angular scales differently. We therefore require scale and sample selections that pass null tests across all frameworks simultaneously, an approach that applies directly to Stage-IV surveys where systematic errors dominate.

+++


(reading-the-analysis)=
:::{note} Reading this reproduction
This report pairs the accepted paper with its analysis records.
Start with the [results and PTE table](#sec:results), the [comparison across catalog versions](#sec:catalog_comparison), or the [conclusions](#sec:conclusions).
Open the [ASTRA Inventory](#astra-inventory) to inspect inputs, decisions, outputs, and findings; a figure's **Provenance** control traces the choices behind it.

The [repository guide](https://github.com/cailmdaley/unions-bmodes-mystra#readme) describes the analysis structure and the accepted and submitted universes.
The report shows the accepted-paper reference comparison, not an interactive universe switcher.
The tables and figures retain the complete two-range, six-test reference comparison; the repository guide explains the configurable primary-report selections.
:::

# Introduction

Weak gravitational lensing distorts the observed shapes of distant galaxies, encoding information about the intervening matter distribution along the line of sight {cite:p}`blandford.etal91,miraldaescude91,kaiser92`. First detections came in 2000 {cite:p}`bacon.etal00,kaiser.etal00,vanwaerbeke.etal00,wittman.etal00`. The scalar gravitational potential produces a shear field that is curl-free to leading order and can be decomposed into gradient-like $E$ modes and curl-like $B$ modes; under the Born approximation, scalar density perturbations source only $E$-mode power {cite:p}`stebbins96,kamionkowski.etal97`. Higher-order effects (lens-lens coupling, source clustering, and intrinsic alignments) can produce $B$ modes {cite:p}`hilbert.etal09,schneider.etal02b,crittenden.etal02`, but these are orders of magnitude below current sensitivity {cite:p}`cooray.hu02`. They are likely only accessible to future ultra-high-precision measurements, potentially through cross-correlation techniques involving CMB lensing rotation reconstructions and large-scale-structure tracers {cite:p}`robertson.etal25`. For galaxy shear measurements at current sensitivity, detected $B$-mode power indicates residual observational systematic effects rather than cosmological signal, though the absence of $B$ modes does not rule out systematic contamination.

Stage-III cosmic shear surveys have measured the structure growth parameter $S_8 \equiv \sigma_8 \sqrt{\Omega_{\mathrm{m}}/0.3}$ at percent-level precision, with some reporting up to ${\sim}3\,\sigma$ tension with *Planck* {cite:p}`asgari.etal21,abbott.etal22,li.etal23,planck20`. The Kilo-Degree Survey (KiDS)-Legacy {cite:p}`wright.etal25` found $S_8 = 0.815^{+0.016}_{-0.021}$ (0.73 $\sigma$ from *Planck*) and the Dark Energy Survey (DES) Y6 {cite:p}`abbott.etal26` found $S_8 = 0.798^{+0.014}_{-0.015}$ ($1.1\,\sigma$); both achieved $2\%$ precision on $S_8$, comparable to CMB constraints from *Planck* alone. Stage-IV surveys will reduce statistical uncertainties further, shifting the error budget toward systematic effects; current data offer the opportunity to validate these methods before that transition.

Several methods separate $E$- and $B$-mode contributions in shear two-point statistics. Aperture mass dispersion {cite:p}`schneider96,schneider.etal98` is the foundational real-space method for $E$/$B$ separation. In practice, however, its evaluation requires knowledge of the shear correlation functions down to zero separation, which is unavailable on a finite observed interval and leads to $E$/$B$ mixing {cite:p}`schneider.eifler.krause10,schneider.etal22`. Complete Orthogonal Sets of $E$/$B$-mode Integrals (COSEBIs; {cite:t}`schneider.eifler.krause10`, {cite:t}`asgari.schneider.simon12`) generalize this approach, compressing all $E$/$B$-separable information from $\xi_\pm(\theta)$ on a finite angular interval into discrete orthogonal modes. A small number of modes captures the cosmological signal; by construction, the COSEBI basis excludes the ambiguous modes that cannot be assigned pure $E$ or $B$ character on a finite interval. Pure-mode correlation functions $\xi_\pm^{\mathrm{E/B}}(\theta)$ {cite:p}`schneider.etal22` are an alternative real-space representation that decomposes $\xi_\pm$ into $E$-mode, $B$-mode, and ambiguous components through integral transforms derived from the COSEBI basis functions. Harmonic-space estimators measure $C_\ell^{EE}$ and $C_\ell^{BB}$ directly but require careful treatment of mask-induced $E$/$B$ mixing; catalog-based methods that evaluate spherical harmonic transforms at source positions {cite:p}`wolz.etal25` avoid the pixelization artifacts that complicate standard implementations. The **HybridEB** Fourier band-power estimator {cite:p}`becker.rozo16` constructs $\ell$-space band-powers from linear combinations of binned $\xi_\pm$, projecting out ambiguous modes to minimize $E$/$B$ mixing. DES Y6 adopted this estimator {cite:p}`abbott.etal26`, although COSEBIs are more widely used in the Stage-III literature, where they have traced systematic signatures across multiple surveys {cite:p}`asgari.etal19a`. These statistics represent $B$-mode power in complementary bases, and the mapping between them is not one-to-one, so contamination excluded by scale cuts in one basis may remain visible in another. Comparing across bases can expose systematics that any single framework would absorb.

Stage-III analyses have found $B$-mode contamination from a range of systematic effects. COSEBI analyses of the Canada-France-Hawaii Telescope Lensing Survey (CFHTLenS), KiDS-450, and DES-SV detected significant $B$ modes at 2--5 $\sigma$, linked to repeating additive shear bias, PSF leakage, and photometric selection effects {cite:p}`asgari.etal19a`. Repeating additive shear biases can arise from detector-level defects that imprint a constant ellipticity offset on each CCD, varying across the focal plane; because this pattern is fixed in instrument coordinates, it repeats with every pointing. The bias is constant within each CCD but jumps at chip boundaries; these discontinuities carry curl, generating $B$-mode power at CCD angular scales. Hyper Suprime-Cam (HSC) Y3 used harmonic-space power spectra and found significant large-scale $B$ modes {cite:p}`li.etal23,dalal.etal23`. A dedicated PSF-systematics study identified PSF fourth-moment leakage as a leading additive contaminant in the affected region {cite:p}`zhang.etal23`. KiDS-Legacy tested with COSEBIs, traced an initial failure to astrometric systematics, and applied conservative masking with negligible impact on $S_8$ {cite:p}`wright.etal25`. DES Y3 reported null detections using pseudo-$C_\ell$ and COSEBIs {cite:p}`gatti.etal21`. However, {cite:t}`jefferson.etal25` applied the **HybridEB** estimator uniformly to DES-Y3 and HSC-Y3 data and found $B$-mode failures in specific tomographic bins not seen with the estimators originally used. DES Y6 subsequently adopted the **HybridEB** estimator and passed {cite:p}`abbott.etal26`. Because each estimator weights angular scales differently and deprojects ambiguous modes differently, the same data can pass one null test and fail another. The most informative validation requires consistent passage across multiple statistics.

The Ultraviolet Near-Infrared Optical Northern Survey (UNIONS; {cite:t}`gwyn.etal25`) is surveying 6250 deg² of northern sky with multi-band imaging from telescopes in Hawai'i: the Canada-France Imaging Survey (CFIS) contributes $u$- and $r$-band imaging from CFHT, Pan-STARRS contributes $i$- and $z$-band data, and Subaru adds $z$-band imaging through the Wide Imaging with Subaru HSC of the *Euclid* Sky (WISHES) and $g$-band imaging through the Waterloo-Hawai'i Institute for Astronomy $g$-band Survey (WHIGS). Nearly all Stage-III cosmic shear survey area is in the southern sky; UNIONS provides the first deep, wide-area shape catalog in the north, with $r$-band median seeing of $0.7^{\prime\prime}$, enabling high-quality shape measurements. {cite:t}`guinot.etal22` presented the first UNIONS weak-lensing analysis over 1700 deg² using the `ShapePipe` pipeline {cite:p}`farrens.etal22`; the present analysis expands this to 3500 deg² with an updated catalog ([section](#sec:data)).

We present $B$-mode validation tests for the first UNIONS-3500 cosmic shear cosmology analysis, using pure-mode correlation functions $\xi_\pm^{\mathrm{E/B}}(\theta)$, COSEBIs, and harmonic-space power spectra $C_\ell^{BB}$ to define scale cuts in both angular and multipole space within a blinded analysis framework. This paper is one of five coordinated UNIONS-3500 weak-lensing publications: catalog construction (Paper I; {cite:t}`hervaspaters.etal26a`); the present $B$-mode validation (Paper II); configuration-space cosmological constraints (Paper III; {cite:t}`goh.etal26`); harmonic-space cosmological constraints (Paper IV; {cite:t}`guerrini.etal26`); and image simulations (Paper V; Hervas-Peters et al., in prep.). We validate the fiducial scale cuts adopted for cosmological inference in Papers III and IV, both of which analyze the shear field without tomographic binning; a UNIONS tomographic analysis is now in progress. [section](#sec:data) describes the UNIONS shear catalogs, their evolution, and the simulations we use for validation. [section](#sec:methods) details the three $E$/$B$-separable statistics and covariance framework. [section](#sec:results) presents the measurements and null-test results. [section](#sec:discussion) discusses systematic trends and methodological implications.

(sec:data)=
# Data

UNIONS is the largest deep weak-lensing survey of the northern sky, with $r$-band depth of $r \approx 24.2$ at $10\,\sigma$ {cite:p}`gwyn.etal25`. We analyze four variants of the `ShapePipe` v1.4 catalog, summarized in [table](#tab:survey_properties); all survey-property statistics quoted in this article are derived from that table. The initial catalog uses a loose size cut ($r_{\mathrm{h, gal}}/r_{\mathrm{h, psf}} > 0.5$). For the size-cut catalog, we tighten this threshold to $r_{\mathrm{h, gal}}/r_{\mathrm{h, psf}} > 0.707$, removing galaxies where shape measurement and leakage correction are least reliable. For the masked catalog, we retain the tightened size cut and add bright-star ($r < 12$) and faint-star ($12 < r < 18$) halo masking, reducing sky coverage by 13%. [figure](#fig:footprint_zoom) shows the holes this masking cuts into the footprint, typically ${\sim}8'$ across, comparable to the short side of a CCD. The relaxed-flags catalog uses the same sample selection as the size-cut catalog but relaxes the `SExtractor` flags criterion to include deblended objects, increasing the effective galaxy density by 26%. We introduced these refinements after blinded validation tests showed early signs of systematic contamination (Paper I). We adopt the size-cut catalog as the fiducial for Papers II, III, and IV.

For each catalog version, we produce shear measurements both with and without the object-wise PSF-leakage correction described in Paper I, which models PSF contamination via binned regression in signal-to-noise and galaxy-PSF size ratio; {cite:t}`guerrini.etal25` provide the broader PSF-leakage diagnostic framework used to assess its impact. We report results from the leakage-corrected catalogs, which are used for cosmological inference in Papers III and IV.

To validate our estimators and covariances, we use {astra:value col=display where="macro=mockNRealizations"}`outputs.paper_claims` lognormal simulated catalogs from the suite described in Paper III, generated with the `GLASS` package {cite:p}`tessore.etal23`. Each realization samples galaxy positions from a full-sky lognormal density field, within the survey footprint and following the fiducial redshift distribution. The galaxies carry unit weights and Gaussian intrinsic ellipticities matching the data's $\sigma_e$. The realizations thus contain sample variance (including from modes larger than the footprint), lognormal non-Gaussianity, and shape noise, but no PSF leakage or other imaging systematics. Because of an error in the catalog generation, their galaxy number density is {astra:value col=display where="macro=mockDensityExcess"}`outputs.paper_claims` higher than that of the data, which we account for when comparing covariances ([section](#sec:mock_validation)). The companion cosmology papers validate their pipelines and covariances on the same suite (Paper III for `CosmoCov`, Paper IV for `NaMaster`).

(tab:survey_properties)=

Survey properties for UNIONS shear catalog versions.

| Catalog | Area (deg²) | $n_{\mathrm{eff}}$ (arcmin⁻²) | $\sigma_e$ |
|---|---:|---:|---:|
| Initial | 2894 | 6.48 | 0.28 |
| Size cuts (fiducial) | {astra:value col=display where="macro=surveyArea"}`outputs.paper_claims` | {astra:value col=display where="macro=surveyNeff"}`outputs.paper_claims` | {astra:value col=display where="macro=surveySigmaE"}`outputs.paper_claims` |
| Masked | 2517 | 4.98 | 0.27 |
| Relaxed flags | 2894 | 6.26 | 0.27 |

$n_{\mathrm{eff}}$ follows {cite:t}`heymans.etal12` and $\sigma_e$ is the per-component ellipticity dispersion, both from metacalibration-weighted samples. All versions share the same footprint except the masked catalog.



(fig:footprint_zoom)=
:::{note} Footprint illustration
The footprint illustration is Figure 1 of the [paper](https://arxiv.org/abs/2604.03227).
It isn't included among this reproduction's generated outputs.
:::



(sec:methods)=
# Methods

(fig:pure_eb_decomposition)=
:::{astra} outputs.fig_pure_eb_decomposition
:::

Pure $E$/$B$-mode decomposition of the measured shear correlation functions for the fiducial catalog. Panels show $\xi_+$ (left) and $\xi_-$ (right) decomposed into total (black circles), $E$-mode (teal squares), $B$-mode (crimson $\times$), and ambiguous (purple triangles) components following {cite:t}`schneider.etal22`. The $E$ modes trace the lensing signal, while $B$ modes---which should vanish for pure gravitational lensing---diagnose systematic contamination. Ambiguous modes capture boundary effects from the finite angular range. Shaded regions indicate fiducial scale cuts.



## Two-point correlation measurements

All $B$-mode diagnostics in this analysis use the discrete shear catalog $\mathcal{C} = \{(\gamma_i, w_i, \hat{\mathbf{n}}_i)\}_{i=1}^N$, comprising complex ellipticities $\gamma_i$, shear weights $w_i$ (Paper I), and sky positions $\hat{\mathbf{n}}_i$ for $N$ source galaxies. The ellipticities $\gamma_i$ are shear estimates calibrated via metacalibration {cite:p}`huff.mandelbaum17,sheldon.huff17`: galaxy images are artificially sheared to measure the $2 \times 2$ shear response matrix $\mathbf{R}$, yielding $\boldsymbol{\gamma}_i = \mathbf{R}^{-1}(\mathbf{e}_i - \mathbf{c}_i)$ from the observed ellipticity $\mathbf{e}_i$ and additive bias $\mathbf{c}_i$. The two-point shear correlation functions $\xi_\pm(\theta)$ are formed from products of tangential and cross-component ellipticities for galaxy pairs at separation $\theta$. We estimate $\xi_\pm$ by weighting each pair by metacalibration weights and accumulating contributions in angular bins: $$\begin{equation}
{\xi}_\pm(\theta_k) = \frac{\sum_{i \neq j} w_i w_j \, \Pi(\theta_{ij}, \theta_k) \left[\gamma_{\mathrm{t},i} \gamma_{\mathrm{t},j} \pm \gamma_{\times,i} \gamma_{\times,j}\right]}{\sum_{i \neq j} w_i w_j \, \Pi(\theta_{ij}, \theta_k)}.
\label{eq:xi_estimator}
\end{equation}$$ Here $\gamma_{\mathrm{t}}$ and $\gamma_{\times}$ are the tangential and cross components of the ellipticity in the coordinate frame defined by the pair-separation vector, $\Pi(\theta_{ij}, \theta_k)$ is the binning kernel selecting pairs with separations in bin $k$, and $\theta_{ij} = \arccos(\hat{\mathbf{n}}_i \cdot \hat{\mathbf{n}}_j)$. We evaluate [equation](#eq:xi_estimator) using the tree-based correlation package `TreeCorr` {cite:p}`jarvis04`. We bin into 20 logarithmic bins over $\theta = 1$--$250$ arcmin, matching the configuration-space data-vector binning adopted in Paper III.

In harmonic space, the analogous statistics are the angular power spectra $C_\ell^{EE}$ and $C_\ell^{BB}$. We estimate these using the MASTER algorithm {cite:p}`hivon.etal02`, which recovers the true $C_\ell$ by deconvolving the mode-coupling matrix induced by the survey geometry, as implemented in `NaMaster` {cite:p}`alonso.etal19`. Traditional implementations operate on pixelized maps, but {cite:t}`wolz.etal25` introduced a catalog-based estimator that evaluates spherical harmonic transforms directly at source positions, avoiding pixelization artifacts. The same weighted ellipticities from [equation](#eq:xi_estimator) enter this estimator: $$\begin{equation}
\tilde{C}_\ell^{\alpha\beta} = \sum_{i,j} w_i w_j \, \gamma_i^a \, \gamma_j^b \; \frac{1}{2\ell+1} \sum_m {}_2 Y_{\ell m}^{a\alpha *}(\hat{\mathbf{n}}_i) \; {}_2 Y_{\ell m}^{b\beta}(\hat{\mathbf{n}}_j),
\label{eq:cl_definition}
\end{equation}$$ where $a, b \in \{1, 2\}$ label the two ellipticity components, ${}_2 Y_{\ell m}^{a\alpha}$ are spin-2 spherical harmonics, and $\alpha, \beta \in \{E, B\}$. We then decouple the measured bandpowers and correct for noise bias following {cite:t}`wolz.etal25`, yielding unbiased harmonic-space bandpower estimates $\hat{C}_\ell^{\alpha\beta}$.

(sec:pure_eb)=
## Pure $E$/$B$ decomposition

We decompose the measured correlation functions $\xi_\pm(\theta)$ into pure $E$-mode, $B$-mode, and ambiguous components following {cite:t}`schneider.etal22`: $$\begin{align}
\xi_+(\theta) &= \xi_+^{\mathrm{E}}(\theta) + \xi_+^{\mathrm{B}}(\theta) + \xi_+^{\mathrm{amb}}(\theta); \label{eq:decomp_plus} \\
\xi_-(\theta) &= \xi_-^{\mathrm{E}}(\theta) - \xi_-^{\mathrm{B}}(\theta) + \xi_-^{\mathrm{amb}}(\theta). \label{eq:decomp_minus}
\end{align}$$ Ambiguous modes $\xi_\pm^{\mathrm{amb}}(\theta)$ arise from shear patterns that cannot be assigned unique $E$- or $B$-mode origin on a finite angular interval; they take the functional forms $\xi_+^{\mathrm{amb}}(\theta) = a + b\theta^2$ and $\xi_-^{\mathrm{amb}}(\theta) = c\theta^{-2} + d\theta^{-4}$. They therefore appear on large scales in $\xi_+$ and on small scales in $\xi_-$, as seen in [figure](#fig:pure_eb_decomposition). The pure-mode correlation functions are computed from integral transforms over the measured $\xi_\pm(\theta)$:

$$
\begin{aligned}
\xi_+^{\mathrm{E,B}}(\theta) &= \frac{1}{2}\left[\xi_+(\theta) \pm \xi_-(\theta) \pm \int_\theta^{\theta_{\max}} \frac{\mathrm{d}\theta'}{\theta'}\,
\xi_-(\theta')\left(4-\frac{12\theta^2}{\theta'^2}\right)\right] \\
&\quad - \frac{1}{2}\left[S_+(\theta) \pm S_-(\theta)\right];
\end{aligned}
$$ (eq:xip_EB)

$$
\begin{aligned}
\xi_-^{\mathrm{E,B}}(\theta) &= \frac{1}{2}\left[\xi_+(\theta) \pm \xi_-(\theta) + \int_{\theta_{\min}}^\theta \frac{\mathrm{d}\theta'\,\theta'}{\theta^2}\,
\xi_+(\theta')\left(4-\frac{12\theta'^2}{\theta^2}\right)\right] \\
&\quad - \frac{1}{2}\left[V_+(\theta) \pm V_-(\theta)\right].
\end{aligned}
$$ (eq:xim_EB)

Here, the $+$ sign gives $E$-modes and the $-$ sign gives $B$-modes. The auxiliary functions $S_\pm(\theta)$ and $V_\pm(\theta)$ enforce boundary conditions that remove ambiguous modes: $\xi_+^{\mathrm{amb}} = S_+$ and $\xi_-^{\mathrm{amb}} = V_-$ {cite:p}`schneider.etal22`.

The integral transforms require finely binned $\xi_\pm$ over an extended angular range. We evaluate these transforms on a separate 1000-bin logarithmic grid spanning $\theta = 0.5$--$300$ arcmin, and the limits $\theta_{\min}$ and $\theta_{\max}$ in [equation](#eq:xip_EB) and [equation](#eq:xim_EB) are this full extent; the extended range reduces edge effects in the boundary functions $S_\pm$ and $V_\pm$. We apply the {cite:t}`schneider.etal22` filter kernels to this grid as a fixed linear operator, using the `cosmo_numba`[^1] implementation {cite:p}`guinot.etal26` and integrating each grid cell with a Gauss--Legendre rule. We then average the resulting pure modes over each of the 20 reporting bins, weighting by pair counts. [figure](#fig:pure_eb_decomposition) shows the pure $E$/$B$ decomposition for the fiducial catalog.

(sec:cosebis)=
## COSEBIs

(fig:cosebis_fiducial)=
:::{astra} outputs.fig_cosebis_fiducial
:::

COSEBI $B$-mode amplitudes for the fiducial catalog, normalized by uncertainty ($B_n/\sigma_n$), at the fiducial scale cuts of {astra:value col=display where="macro=ebthetaXipMin"}`outputs.paper_claims`--{astra:value col=display where="macro=ebthetaXipMax"}`outputs.paper_claims` arcmin (filled circles) and over the full angular range of $1$--$250$ arcmin (open diamonds). Modes $n > 6$ are shaded; the $n \leq 6$ subset captures nearly all cosmological information.



COSEBIs compress $\xi_\pm(\theta)$ into discrete orthogonal modes $E_n$ and $B_n$, with clean $E$/$B$ separation on finite angular intervals {cite:p}`schneider.eifler.krause10`. The filter functions $T_{\pm n}(\theta)$ are constructed to satisfy boundary conditions that remove the ambiguous modes described in [section](#sec:pure_eb). Because the $T_{\pm n}$ are orthonormal, different modes capture largely independent information. Each $T_{\pm n}(\theta)$ oscillates $n+1$ times across the angular interval; higher-$n$ modes are therefore sensitive to increasingly oscillatory features in $\theta$-space. Conversely, a feature localized in angle projects onto many modes at once, producing an oscillatory signature across the basis.

We compute the COSEBI mode amplitudes $X_n \in \{E_n, B_n\}$ over a given angular range $[\theta_{\mathrm{min}}, \theta_{\mathrm{max}}]$ arcmin using logarithmic filter functions: $$\begin{equation}
X_n = \frac{1}{2} \int_{\theta_{\mathrm{min}}}^{\theta_{\mathrm{max}}} \mathrm{d}\theta \, \theta \left[ T_{+n}(\theta) \xi_+(\theta) \pm T_{-n}(\theta) \xi_-(\theta) \right],
\label{eq:cosebi_xn}
\end{equation}$$ where the upper (lower) sign gives $E_n$ ($B_n$). We evaluate these integrals on the same 1000-bin grid used for the pure-mode transforms ([section](#sec:pure_eb)). As a convergence check, we repeat the calculation on a 10 000-bin grid for the fiducial catalog and find the $B_n$ amplitudes and their significance to be stable.

[figure](#fig:cosebis_fiducial) shows COSEBI $B$-mode amplitudes for the first $n_{\mathrm{max}} = 20$ modes of the fiducial catalog. For the UNIONS data vector, the $E$-mode signal saturates by $n \approx 5$--$6$ (see [section](#sec:discussion)), so we take $n \leq 6$ as the set of modes carrying most of the cosmological information, while also testing the $B$ modes for all $n \leq 20$.

Pure-mode correlation functions ([section](#sec:pure_eb)) can be reconstructed from the COSEBI coefficients: $$\begin{equation}
\xi_\pm^{X}(\theta) = \frac{\bar{\theta}^2}{\beta_\theta} \sum_{n=1}^{n_{\mathrm{max}}} X_n \, T_{\pm n}(\theta),
\label{eq:xi_XEB_reconstruction}
\end{equation}$$ where $\bar{\theta} = (\theta_{\mathrm{min}} + \theta_{\mathrm{max}})/2$ and $\beta_\theta = (\theta_{\mathrm{max}} - \theta_{\mathrm{min}})/(\theta_{\mathrm{max}} + \theta_{\mathrm{min}})$. In practice, we compute the pure-mode correlation functions directly from the finely binned $\xi_\pm$ measurements rather than reconstructing them from the COSEBI coefficients. COSEBIs require consistent angular ranges for $\xi_+$ and $\xi_-$, while the pure-mode integral transforms [equation](#eq:xip_EB) and [equation](#eq:xim_EB) permit independent scale cuts for each.

The COSEBI modes can equivalently be computed from harmonic-space power spectra: $$\begin{equation}
X_n = \sum_i \frac{\Delta\ell_i \, \ell_i}{2\pi} \, W_n(\ell_i) \, C_{\ell_i}^{XX},
\label{eq:cosebi_harmonic}
\end{equation}$$ where $W_n(\ell)$ are the harmonic-space COSEBI filter functions (Fourier duals of $T_{\pm n}(\theta)$, computed via FFT-log) and the sum runs over bandpower bins. Each $W_n(\ell)$ oscillates $n+1$ times, mirroring its real-space counterpart. A fine binning of the bandpowers is needed to resolve these oscillations, so for this transform we use 96 square-root-spaced bins, distinct from the 32 bandpowers of the direct $C_\ell^{BB}$ null tests ([section](#sec:catalog_cls)). The number of bins balances two requirements: too few cannot resolve the filter oscillations at high $n$, while too many lead to an unstable mode-coupling matrix deconvolution. We validated the transform on the GLASS simulations, confirming that modes $n \leq 6$---the range carrying nearly all cosmological $E$-mode information---agree between the two paths.

(fig:cl_fiducial)=
:::{astra} outputs.fig_cl_fiducial
:::

Harmonic-space power spectra for the fiducial catalog, normalized by uncertainty. Top: $B$-mode auto-power $C_\ell^{BB}/\sigma$. Bottom: $E$-$B$ cross-power $C_\ell^{EB}/\sigma$. $C_\ell^{BB}$ is predominantly positive across the fiducial range. $C_\ell^{EB}$ scatters symmetrically around zero, suggesting that any $B$-mode contamination is not strongly correlated with the lensing signal.



(sec:catalog_cls)=
## Catalog-based harmonic-space power spectra

We estimate $C_\ell^{EB}$ and $C_\ell^{BB}$ [equation](#eq:cl_definition) directly from the discrete catalog, providing $\ell$-dependent $B$-mode diagnostics that complement the angular-space statistics. The {cite:t}`wolz.etal25` estimator represents the survey mask as delta functions at source positions and computes the spin-2 transforms on this irregular grid using the `ducc` library {cite:p}`reinecke.seljebotn23`. We subtract an analytic shot-noise contribution following the implementation in Paper IV. We bin into 32 bandpowers with square-root spacing ($\ell^{0.5}$) over $8 \leq \ell \leq 2048$. [figure](#fig:cl_fiducial) shows the measured $C_\ell^{BB}$ and $C_\ell^{EB}$ for the fiducial catalog.

(fig:eb_covariance)=
:::{astra} analyses.pure_eb.outputs.pure_eb_covariance_figure
:::

Pure $E$/$B$-mode correlation-function covariance matrix for the fiducial catalog, propagated linearly from the `CosmoCov` $\xi_\pm$ covariance ([section](#sec:covariance)). The matrix is organized into six blocks: $\xi_+^{\mathrm{E}}$, $\xi_-^{\mathrm{E}}$, $\xi_+^{\mathrm{B}}$, $\xi_-^{\mathrm{B}}$, $\xi_+^{\mathrm{amb}}$, $\xi_-^{\mathrm{amb}}$, with 20 angular bins per block. Off-diagonal correlations arise from integral coupling in the pure-mode filters.



(sec:covariance)=
## Covariance estimation

(fig:mock_validation)=
:::{astra} analyses.glass_mocks.outputs.mock_validation_figure
:::

Validation of the $B$-mode estimators and covariances on {astra:value col=display where="macro=mockNRealizations"}`outputs.paper_claims` `GLASS` mock realizations. Top: mean $B$-mode data vector across realizations, in units of the analytic $\sigma$. Bottom: ratio of the analytic $\sigma$ used in the PTE calculations to the scatter across realizations. Left: COSEBIs at the fiducial cuts (filled) and over all scales (open); center: $\xi_+^{\mathrm{B}}$ (filled) and $\xi_-^{\mathrm{B}}$ (open); right: $C_\ell^{BB}$. Shaded regions mark scales outside the fiducial cuts.



For each $B$-mode statistic, we start from the same covariance prescriptions adopted in the companion cosmology papers. In configuration space, Paper III uses `CosmoCov` {cite:p}`krause.eifler17,fang.eifler.krause20` to model the covariance of the shear correlation functions $\xi_\pm$; we compute this covariance on a 1000-bin logarithmic grid spanning $\theta = 0.5$--$300$ arcmin, using the survey properties of Paper I and the Planck 2018 fiducial cosmology {cite:p}`planck20`. For COSEBIs and the pure-mode correlation functions, we compute the covariance as the linear transformation $\mathbf{C} = \mathbf{T}^{\mathsf{T}} \mathbf{C}_{\xi} \mathbf{T}$, where $\mathbf{T}$ encodes the COSEBI filter functions or the pure-mode transforms of [section](#sec:pure_eb); we show the pure-mode covariance in [figure](#fig:eb_covariance). In harmonic space, we use the Gaussian covariance from `NaMaster`'s iNKA framework for $C_\ell^{BB}$ and $C_\ell^{EB}$, as Paper IV does for its $B$-mode null tests, adding non-Gaussian terms only for the $EE$ cosmology analysis.

For all $B$-mode null tests, we retain only the Gaussian part of the covariance. In configuration space, this keeps the calculation tractable on the finely binned integration grid. This approximation underestimates the true variances because it omits connected trispectrum and super-sample covariance terms, which can also add off-diagonal structure. In Paper IV we find that adding these terms, computed with `OneCovariance` {cite:p}`reischke.etal25`, increases the diagonal errors by about 5% at $\ell = 300$. In configuration space, the non-Gaussian terms contribute about {astra:value col=display where="macro=ngFractionDiagonalFiducial"}`outputs.paper_claims` of the diagonal variance of $\xi_\pm$ within the fiducial scale cuts. The null tests are therefore conservative: a test that passes with the Gaussian-only covariance would generally pass with a more complete covariance estimate, as the simulations of [section](#sec:mock_validation) confirm.

Two independent checks support these covariances. Paper IV validates the harmonic-space $B$-mode covariance against Gaussian simulations (Appendix C therein), and its $E$-mode covariance against the `OneCovariance` code {cite:p}`reischke.etal25` and the GLASS mocks of [section](#sec:mock_validation). During development we also compared the semi-analytic $\xi_\pm$ covariance against a jackknife estimate from the data and found ten-percent-level agreement over the fiducial angular range.

(sec:mock_validation)=
## Validation on simulations


No systematics are injected into the `GLASS` realizations ([section](#sec:data)), so they allow us to test directly whether our estimators generate spurious $B$ modes from a noise-dominated signal. For these comparisons we recompute the analytic covariances of [section](#sec:covariance) at the galaxy density of the realizations. Across the ensemble, the pure-mode decomposition, COSEBIs, and harmonic-space power spectra all recover $B$ modes consistent with zero (top row of [figure](#fig:mock_validation)); the mean COSEBI $B_n$ deviates from zero by at most {astra:value col=display where="macro=mockMeanBiasCosebis"}`outputs.paper_claims` $\sigma$ in any mode.

We also use the simulations to test the covariances, comparing the analytic $\sigma$ of each statistic with the scatter across realizations (bottom row of [figure](#fig:mock_validation)). At the fiducial cuts, the median ratio of analytic to mock $\sigma$ over bins is {astra:value col=display where="macro=mockSigmaRatioPureMin"}`outputs.paper_claims`--{astra:value col=display where="macro=mockSigmaRatioCosebis"}`outputs.paper_claims` in configuration space and {astra:value col=display where="macro=mockSigmaRatioCl"}`outputs.paper_claims` in harmonic space. The agreement degrades on the large scales outside the fiducial cuts. There the analytic $\sigma$ is about {astra:value col=display where="macro=mockPureLargeThetaTighter"}`outputs.paper_claims` too small for the first COSEBI mode and for $\xi_\pm^{\mathrm{B}}$ above $100'$, while for $C_\ell^{BB}$ below $\ell \approx 100$ it is too large because $E$-to-$B$ leakage is treated only approximately {cite:p}`nicola.etal21`. PSF leakage already excludes these scales (Paper I), so the miscalibration does not affect our choice of scale cuts.

(sec:scale_cuts)=
## $B$-mode significance

We assess $B$-mode significance across a grid of angular scale-cut combinations. For each scale range, we test the $B$-mode data vector $\vec{d}_B$ against zero with $\chi^2_B = \vec{d}_B^{\mathsf{T}} \mathbf{C}_{BB}^{-1} \vec{d}_B$, where $\mathbf{C}_{BB}$ is the $B$-mode covariance submatrix for the selected scale range ([section](#sec:covariance)). The probability-to-exceed, $\mathrm{PTE} = P(\chi^2 > \chi^2_B \,|\, \nu)$, is the probability that a $\chi^2$ variable with $\nu$ degrees of freedom exceeds the measured value; we set $\nu$ equal to the number of data points in the range. For the pure-mode decomposition, a scale cut selects reporting bins of the full-range decomposition, and we compute PTEs separately for $\xi_+^{\mathrm{B}}$ and $\xi_-^{\mathrm{B}}$, as well as a joint test $\xi_{\mathrm{tot}}^{\mathrm{B}}$ using the concatenated data vector and full cross-covariance. For COSEBIs, which we recompute for each angular range, we evaluate PTEs using both $n \leq 6$ and $n \leq 20$ ([section](#sec:cosebis)). For $C_\ell^{BB}$, we apply the same framework across multipole ranges. We compute a PTE for each of the six statistics at each scale range and take the smallest $\mathrm{PTE}_{\min}$ as our test statistic. Although each PTE is expected to be uniformly distributed in the absence of $B$-modes, their minimum is skewed toward low values because the six statistics are not perfectly correlated. We therefore set the passing threshold $p_{\mathrm{thr}}$ on the systematics-free mock realizations of [section](#sec:mock_validation) such that a clean sample fails 5% of the time, as would be the case for a single test at PTE $=0.05$. We adopt a representative single value, $p_{\mathrm{thr}} = ${astra:value col=display where="macro=mockMinPteThreshold"}`outputs.paper_claims`, which is conservative for nearly all scale-cut combinations; at the fiducial cuts the calibrated value is {astra:value col=display where="macro=mockMinPteFloor"}`outputs.paper_claims` ([figure](#fig:mock_min_pte)). We do not apply a look-elsewhere correction for the many scale-cut combinations tested, because the $B$ modes were only one of four factors that informed the scale cuts, alongside PSF leakage, modeling systematics, and the scale-cut dependence of the blinded cosmological constraints (Papers III and IV).[^2]

We use the same realizations to test the $\chi^2_\nu$ sampling distribution assumed for the PTE. The $C_\ell^{BB}$ PTEs are uniform, while the configuration-space PTEs skew mildly low, with $\langle\chi^2\rangle/\nu = ${astra:value col=display where="macro=mockChiTwoNuConfigMin"}`outputs.paper_claims`--{astra:value col=display where="macro=mockChiTwoNuConfigMax"}`outputs.paper_claims`. This excess matches the covariance comparison above. In configuration space the analytic $\sigma$ is slightly smaller than the mock scatter, which raises $\chi^2$ by the inverse square of their ratio; in harmonic space the two agree. Calibrating $p_{\mathrm{thr}}$ on the same realizations partially absorbs this skew.

(fig:mock_min_pte)=
:::{astra} analyses.glass_mocks.outputs.mock_min_pte_figure
:::

Cumulative distribution of the minimum PTE over the six statistics at the fiducial cuts, across the {astra:value col=display where="macro=mockNRealizations"}`outputs.paper_claims` `GLASS` realizations. Gray dashed lines mark the 5% quantile, which calibrates the adopted threshold $p_{\mathrm{thr}} = ${astra:value col=display where="macro=mockMinPteThreshold"}`outputs.paper_claims` (solid gray lines), chosen to be conservative across nearly all scale cuts. Red lines mark the minimum PTE of the data, which {astra:value col=display where="macro=mockPGlobal"}`outputs.paper_claims` of realizations fall at or below. The dotted curve shows the expectation for six independent tests. The mocks generally follow this expectation because of two effects that largely cancel: correlations between the statistics raise the minimum PTE, while the slightly underestimated configuration-space covariance ([section](#sec:mock_validation)) lowers it.



These tests were carried out within the blinded analysis framework shared with Papers III and IV. Blinding was performed on the source redshift distribution. An external coordinator generated three $n(z)$ variants, two with random shifts, and all analysis choices were made before the true variant was revealed (Paper III). The $n(z)$ does not enter the $B$-mode data vectors, and its effect on their covariances is negligible. We computed the covariances for each of the three variants. In configuration space they are identical, since the Gaussian $B$-mode covariance of a null signal is shape noise alone; in harmonic space they agree to better than 0.5% over the fiducial multipole range, the residual coming from $E$-to-$B$ leakage from the mask. We therefore use a single variant throughout.

(sec:results)=
# Results

(tab:pte_results)=
:::{astra} outputs.pte_table_appendix
:::

$B$-mode PTE values across catalog versions at fiducial and full-range scale cuts.

Fiducial scale cuts are {astra:value col=display where="macro=ebthetaXipMin"}`outputs.paper_claims`--{astra:value col=display where="macro=ebthetaXipMax"}`outputs.paper_claims` arcmin in configuration space and $300 < \ell < 1600$ in harmonic space; the full range is $1$--$250$ arcmin / $8 \leq \ell \leq 2048$.
A row passes the null test when its minimum PTE is at least $p_{\mathrm{thr}} = ${astra:value col=display where="macro=mockMinPteThreshold"}`outputs.paper_claims`; only the fiducial catalog at the adopted cuts passes.

For the fiducial catalog before scale cuts, COSEBIs fail with a $4\,\sigma$ first mode and both pure-mode correlation functions fail, while $C_\ell^{BB}$ passes but shows low-level $B$-mode structure ([table](#tab:pte_results)). At the adopted scale cuts, the fiducial catalog passes with a minimum PTE of {astra:value col=display where="macro=fiducialMinPte"}`outputs.paper_claims`; {astra:value col=display where="macro=mockPGlobal"}`outputs.paper_claims` of the systematics-free realizations have a lower $\mathrm{PTE}_{\min}$ ([figure](#fig:mock_min_pte)). We chose those cuts using the $B$-mode tests together with additional systematics checks and blinded inference-stability checks in Papers III and IV. The object-wise PSF-leakage correction (Paper I) shifts individual PTEs by at most {astra:value col=display where="macro=leakCorrMaxPteShift"}`outputs.paper_claims` across all catalog versions and statistics; no pass/fail conclusion changes. We report leakage-corrected results throughout.

The low-level features that drive the full-range behavior leave different signatures in each statistic. In the pure-mode decomposition ([figure](#fig:pure_eb_decomposition)), $\xi_+^{\mathrm{B}}/\sigma$ shows a broad positive excess on the smallest angular scales that tapers toward zero by $\theta \sim 10$ arcmin, apart from a ${\sim}3\,\sigma$ outlier in the 9--12 arcmin bin, just below the fiducial lower cut. This outlier appears in all four catalog versions ([figure](#fig:pure_eb_versions)) and coincides with a jump in the total $\xi_+$, at roughly the angular size of a MegaCam CCD. The corresponding $\xi_-^{\mathrm{B}}$ is quieter on the smallest scales but remains broadly elevated from roughly 7 arcmin to the largest scales. The COSEBI data vector ([figure](#fig:cosebis_fiducial)) compresses the same structure into a smaller number of modes: on the full angular range, the first mode reaches $4\,\sigma$ and the higher modes show coherent oscillations, whose possible origins we discuss in [section](#sec:discussion); after applying the fiducial angular cuts, the first-mode excess drops below 1 $\sigma$ and the oscillatory pattern largely disappears. In harmonic space ([figure](#fig:cl_fiducial)), $C_\ell^{BB}$ shows a low-level positive offset across much of the fiducial range, with the most conspicuous outliers near $\ell \approx 125$ and $250$, and above $\ell \approx 1600$. Even so, the fiducial range $300 < \ell < 1600$ passes the null test, while $C_\ell^{EB}$ remains consistent with zero, suggesting that any $B$-mode contamination is not strongly correlated with the lensing signal.

[figure](#fig:pure_eb_versions) compares the $B$-mode measurements across catalog versions, and [table](#tab:pte_results) summarizes the PTEs. At the fiducial cuts, $\xi_+^{\mathrm{B}}$ and COSEBIs pass for all four versions, whereas $\xi_-^{\mathrm{B}}$ shows a broad positive excess from roughly 7 arcmin to the largest scales in every version, driving failures in all but the fiducial catalog. The fiducial catalog passes, but at the edge of its $\xi_-^{\mathrm{B}}$ acceptance region ([figure](#fig:pte_heatmaps)). In harmonic space ([figure](#fig:pte_cl)), the initial and relaxed-flags catalogs fail across most multipole combinations, whereas the fiducial and masked catalogs show broad acceptance regions. Only the fiducial catalog passes all statistics at the adopted cuts.

At full range, all versions show low-order oscillatory structure consistent with the repeating additive pattern discussed below, driving COSEBI ($n \leq 6$) PTEs to {astra:value col=display where="macro=configPteSixThreeCosebisFull"}`outputs.paper_claims` or below; at fiducial cuts, all versions pass, although the initial catalog retains a mild oscillatory pattern.

[figure](#fig:pte_cl) and [figure](#fig:pte_heatmaps) map null-test performance over the full grid of scale-cut combinations. In the configuration-space maps, $\xi_+^{\mathrm{B}}$ and COSEBI failures concentrate at small $\theta_{\min}$, below the fiducial cut. For COSEBIs, whose data vector is indexed by mode number rather than angle, the PTE heatmap translates null-test performance back into angular space.

(sec:catalog_comparison)=
## Catalog comparison

The following three figures compare $B$-mode data vectors across catalog versions: initial (red circles), size-cut/fiducial (gold squares), masked (teal diamonds), and relaxed-flags (purple triangles).
Rectangles delineate each bin in the horizontal direction and span the range of values across versions in the vertical direction; horizontal lines mark the fiducial value.
Only the size-cut catalog passes all statistics at the fiducial scale cuts ([table](#tab:pte_results)).

(fig:pure_eb_versions)=
:::{astra} outputs.fig_versions_pure_eb
:::

Pure $E$/$B$-mode correlation functions across catalog versions.
Upper panels show $\xi_+$ and $\xi_-$ scaled by $\theta$; lower panels show $B$-mode components $\xi_\pm^{\mathrm{B}}/\sigma$.
Shaded bands mark fiducial scale cuts.

(fig:cl_versions)=
:::{astra} outputs.fig_versions_cl
:::

Harmonic-space power spectra $C_\ell^{BB}/\sigma$ and $C_\ell^{EB}/\sigma$ across catalog versions.
Shaded bands mark fiducial scale cuts.

(fig:cosebis_versions)=
:::{astra} outputs.fig_versions_cosebis
:::

COSEBI $B$-mode amplitudes $B_n/\sigma_n$ across catalog versions at full range (upper) and fiducial cuts (lower).

(sec:scale_cut_maps)=
## Scale-cut maps



(fig:pte_heatmaps)=
:::{astra} outputs.fig_pte_heatmaps
:::

Configuration-space PTE maps for all four catalog versions ([table](#tab:survey_properties)). Columns show $\xi_+^{\mathrm{B}}$, $\xi_-^{\mathrm{B}}$, and COSEBI $B_n$ ($n \leq 6$) PTEs as a function of angular scale cuts; rows show different catalog versions. Solid blue cells indicate PTE $< ${astra:value col=display where="macro=mockMinPteThreshold"}`outputs.paper_claims`; whites and reds indicate consistency with zero. Black squares mark the adopted cuts ({astra:value col=display where="macro=ebthetaXipMin"}`outputs.paper_claims`--{astra:value col=display where="macro=ebthetaXipMax"}`outputs.paper_claims` arcmin). The size cut (fiducial versus initial) expands the acceptance region, while stellar masking contracts it.



(fig:pte_cl)=
:::{astra} outputs.fig_cl_pte_heatmap
:::

Harmonic-space $C_\ell^{BB}$ PTE maps for all four catalog versions ([table](#tab:survey_properties)). Solid blue cells indicate failing PTEs; whites and reds indicate consistency with zero. The initial and relaxed-flags catalogs show widespread failures across multipole ranges; the fiducial and masked catalogs pass across nearly all multipole combinations.



(sec:discussion)=
# Discussion

The pure-mode correlation functions, COSEBIs, and harmonic-space power spectra presented in this paper weight angular scales and treat ambiguous modes differently, so they do not respond identically to the same contamination. When they disagree, we treat that as evidence of contamination in the sample even if each statistic is not equally sensitive to it, and require the adopted catalog and cuts to pass in all three frameworks.

At full range, all four catalogs show a COSEBI first-mode excess exceeding $4\,\sigma$ and oscillatory structure across all twenty modes. Because COSEBIs are orthogonal over a finite angular range ([section](#sec:cosebis)), a feature compact in angle spreads across the mode spectrum; the full-range oscillations are qualitatively consistent with contamination at CCD angular scales ([section](#sec:results), [figure](#fig:cosebis_fiducial)). In the pure-mode data vectors, $\xi_+^{\mathrm{B}}$ shows a ${\sim}3\,\sigma$ outlier in the 9--12 arcmin bin, at CCD scales, where the total $\xi_\pm$ also show noticeable breaks ([figure](#fig:pure_eb_decomposition)). In harmonic space, $C_\ell^{BB}$ has a few ${\gtrsim}2\,\sigma$ outliers at low $\ell$ and elevated power in the highest multipole bins, without a sharp angular signature that would point unambiguously to CCD-scale contamination. MegaCam's 40 CCDs ([figure](#fig:footprint_zoom)) span $6$--$14$ arcmin individually, and our adopted lower cut of {astra:value col=display where="macro=ebthetaXipMin"}`outputs.paper_claims` arcmin excludes most pair separations at these scales; at the adopted scale cuts, the COSEBI features are suppressed and the COSEBI tests pass for all four catalogs. Pure-mode correlation functions and $C_\ell^{BB}$ on the other hand pass at the adopted cuts only for the size-cut catalog and for the size-cut and masked versions, respectively. These equivocal results lead us to argue against relying on any single null test; taken together, we choose the size-cut catalog and our adopted scale cuts as the combination least affected by systematic contamination for cosmological inference in Papers III and IV.

(fig:harmonic_config_cosebis_full)=
:::{astra} outputs.fig_harmonic_config_cosebis_full
:::

COSEBI $E_n$ (top row) and $B_n/\sigma_n$ (bottom row) for the fiducial catalog, computed from configuration-space $\xi_\pm$ (filled circles) and harmonic-space bandpowers (open squares). Modes $n > 6$ are shaded; for $n \leq 6$, where the cosmological information is concentrated, we confirm consistency between the two paths on GLASS simulations.



Although the masked catalog passes in $\xi_+^{\mathrm{B}}$, $\xi_{\mathrm{tot}}^{\mathrm{B}}$, harmonic space, and COSEBIs ($n \leq 6$), this ostensibly more conservative choice fails in $\xi_-^{\mathrm{B}}$. Excising stellar halos introduces irregular gaps in the survey footprint ([figure](#fig:footprint_zoom)), and because $B$-modes are non-local, a more complex mask geometry may increase the fraction of modes classified as ambiguous, leaving the remaining pure-mode estimates less stable; alternatively, removing regions of elevated shape noise near bright stars could increase the significance of contamination elsewhere. We cannot distinguish these explanations with the present data. Our footprints are reconstructed by binning the catalog's per-object mask flags onto a HEALPix grid (Paper I) rather than taken directly from the underlying pixelized masks; the UNIONS tomographic analysis now in progress defines the footprint from the pixelized masks themselves, which will allow a better understanding of the impact of the masks on the data. The initial and relaxed-flags catalogs fail in $\xi_-^{\mathrm{B}}$ and harmonic space, consistent with poorly resolved or deblended objects contributing to the contamination; the relaxed-flags catalog also shifts the total $\xi_+$ systematically higher than the other versions beyond 16 arcmin ([figure](#fig:pure_eb_versions)), suggesting that admitting these objects affects the signal itself, not only the $B$-mode diagnostics.

Configuration-space and harmonic-space estimators respond differently to the same contamination. Pure-mode correlation functions are localized in angular separation, so cuts that remove CCD-scale separations directly reduce the contamination. In multipole space, a feature narrow in angle maps onto a broad range of $\ell$, so it cannot be isolated by cutting a localized multipole range. COSEBIs apply filter functions to the two-point statistics that concentrate cosmological information into the first few modes; contamination at a characteristic angular scale spreads across the mode spectrum rather than isolating in a single mode. Pure-mode and COSEBI filters also discard ambiguous modes ([section](#sec:pure_eb)) that cannot be uniquely classified as $E$ or $B$ on a finite angular interval. The harmonic-space estimator used here does not deproject ambiguous modes; it assigns all power to $C_\ell^{EE}$ or $C_\ell^{BB}$, and any residual mask-induced $E$-to-$B$ leakage after mode-coupling deconvolution is absorbed into $C_\ell^{BB}$.

We show directly that $C_\ell^{BB}$ and COSEBIs have different sensitivity to systematic features, even when computed from the same data, by transforming the harmonic-space bandpowers into COSEBIs ([equation](#eq:cosebi_harmonic)). The $C_\ell^{BB}$ power spectrum passes the null test at our fiducial multipole range, yet COSEBIs computed from the same harmonic-space data fail at full angular range---because the COSEBI filter functions $W_n(\ell)$ weight the power spectrum differently, concentrating sensitivity on the scales where the contamination resides. This provides empirical support for the argument of {cite:t}`asgari.etal19a` that a null test in one representation cannot validate another; {cite:t}`jefferson.etal25` reached the same conclusion across DES-Y3, KiDS-1000, and HSC-Y3, finding that tomographic bins passing pseudo-$C_\ell$ $B$-mode tests can fail when reanalyzed with **HybridEB**. The converse also holds: a $B$-mode failure in any statistic should raise concern for cosmological inference in all of them, since the underlying contamination is present in the data regardless of how it is projected.

[figure](#fig:harmonic_config_cosebis_full) compares COSEBIs computed from configuration-space $\xi_\pm$ and from harmonic-space bandpowers. Both $E_n$ and $B_n$ agree across the two computation paths for all twenty modes, with the tightest correspondence in the first six modes validated on GLASS simulations. At fiducial scale cuts, both paths yield $B$-modes consistent with zero (PTE $= ${astra:value col=display where="macro=harmCosebisPteSixThreeFid"}`outputs.paper_claims` harmonic, {astra:value col=display where="macro=cfgCosebisPteSixThreeFid"}`outputs.paper_claims` configuration; $n \leq 6$). At full range, both detect the same structure (PTE $= ${astra:value col=display where="macro=harmCosebisPteSixThreeFull"}`outputs.paper_claims` and {astra:value col=display where="macro=cfgCosebisPteSixThreeFull"}`outputs.paper_claims`). This confirms that the sensitivity to systematic contamination is set by the filter functions, not by the basis in which the two-point function is measured. Throughout this analysis, we use Gaussian-only covariance, which underestimates variances for the UNIONS geometry ([section](#sec:covariance)); the reported PTEs are therefore conservative.

{cite:t}`asgari.etal19a` identified an oscillatory COSEBI $B$-mode pattern in multiple Stage-III surveys---CFHTLenS, KiDS-450, and DES-SV---tracing it to a repeating additive shear bias (i.e. a spatially varying $c$-term) from detector-level effects fixed in focal-plane coordinates. Its appearance across three Stage-III surveys points to a detector-level phenomenon rather than a pipeline-specific artifact, a point especially relevant for UNIONS because CFHTLenS also used MegaCam on CFHT {cite:p}`heymans.etal12,guinot.etal22`. The dither strategy, however, differs substantially. CFHTLenS and KiDS used small dithers designed to bridge chip gaps, so the focal-plane pattern maps almost directly onto the sky. UNIONS dithers by ${\sim}1/3$ of the MegaCam field of view {cite:p}`gwyn.etal25`, larger than a single CCD, so each sky position averages the additive bias from several different CCDs across exposures. {cite:t}`asgari.etal19a` note that the $B$-mode signature of a repeating additive pattern depends on the dithering strategy, pointing to DES-SV's half-field dithers as one reason its pattern differs from the small-dither surveys. The larger UNIONS dithers may similarly attenuate the amplitude of the repeating pattern, but are unlikely to fully eliminate it: the effective pattern is the per-CCD bias convolved with the dither geometry, potentially spreading the characteristic angular scale beyond the $6$--$14$ arcmin CCD footprint.

A repeating additive shear bias could itself be created by PSF leakage that is tied to an instrumental or observational scale. The $\rho$-statistics and $\xi_{\mathrm{sys}}$ diagnostics of Paper I show a feature near the ${\sim}10$ arcmin size of a MegaCam CCD, consistent with CCD-scale leakage. In the simulations of {cite:t}`asgari.etal19a`, spatially varying PSF leakage mainly raises the low COSEBI modes, while a repeating additive bias produces oscillations across modes; our full-range first-mode excess and oscillations resemble these signatures respectively. However, leakage correction does not significantly change either feature ([section](#sec:results)), pointing to either a different origin for the $B$-modes or leakage beyond our second-order model from higher-order PSF moments {cite:p}`zhang.etal23`.

Paper V discusses a possible pipeline-specific source arising from the improper propagation of exposure offsets within the shape measurement pipeline. In the current multi-exposure fitting procedure, objects are effectively assumed to be centered on the pixel grid of each individual exposure. However, a consistent treatment would require recentering the model using the World Coordinate System (WCS), anchored to the detection position defined at the tile level. In areas where CCDs overlap, the relative pixel offsets between exposures can induce a preferred orientation in the inferred galaxy shapes, corresponding to the vector displacements between exposure pixel grids. As a result, a coherent additive shear pattern may be imprinted across the CCD overlap regions, on angular scales between the per-CCD footprint ($\sim 6$--$14$ arcmin) and the UNIONS dither spacing ($1/3$ MegaCam FoV, or $\sim 19$ arcmin). The observed $B$-mode contamination concentrates at $\theta_\mathrm{min} \lesssim ${astra:value col=display where="macro=ebthetaXipMin"}`outputs.paper_claims` arcmin ([figure](#fig:pte_heatmaps)), close to the per-CCD scale but below the dither offsets. Whether this mechanism is responsible for the observed $B$ modes remains an open question. Even so, a contamination localized near the CCD scale would be cut away in configuration space, while its Hankel transform to harmonic space would spread across multipoles, beyond any single $\ell$ cut. If a comparable additive contribution were present in the $E$ modes, it would lift the harmonic-space $S_8$ relative to configuration space, and could contribute to the mild offset seen between Papers III and IV. Including configuration-space scales down to $5$ arcmin does move its $S_8$ toward the harmonic value (Paper III). While a fraction of these effects may be absorbed by standard multiplicative bias corrections, we mitigate any remaining impact on cosmological constraints by cutting the majority of CCD scales.

An alternative to excluding contaminated scales would be to model the $B$-mode signal as a bias and subtract it from the total measurement. However, the relationship between $E$- and $B$-mode contamination depends on the physical mechanism; {cite:t}`asgari.etal19a` showed that different systematics produce markedly different $E$/$B$ ratios, so the $B$-mode amplitude alone does not determine the $E$-mode contamination. We therefore treat $B$-modes as a diagnostic rather than a correction: scale cuts that remove $B$-mode contamination also remove any corresponding contamination of the total signal, without assuming a model for the relationship between them.

(sec:conclusions)=
# Conclusions

UNIONS-3500 weak-lensing $B$ modes are consistent with zero at the adopted scale cuts across pure-mode correlation functions, COSEBIs, and harmonic-space power spectra. Without scale cuts, $\xi_+^{\mathrm{B}}$ and the COSEBIs show significant $B$ modes on small scales, whereas $\xi_-^{\mathrm{B}}$ and $C_\ell^{BB}$ show a broad, low-level excess with little scale dependence; this excess is significant in $\xi_-^{\mathrm{B}}$ on large scales and marginal at the fiducial cuts.

On the full angular range, all catalog versions show an oscillatory COSEBI $B$-mode pattern plausibly caused by a repeating additive shear bias at CCD angular scales {cite:p}`asgari.etal19a`. This is especially relevant for UNIONS because a similar signature also appeared in CFHTLenS, which likewise used MegaCam on CFHT {cite:p}`heymans.etal12,guinot.etal22`. Our adopted lower cut suppresses this pattern without discarding the angular scales that carry most of the cosmological signal. By computing COSEBIs from the harmonic-space bandpowers, we confirm that the disagreement is not a matter of harmonic versus real space: $C_\ell^{BB}$ passes where COSEBIs fail, because the COSEBI filter functions concentrate sensitivity on the contaminated scales.

Of the four catalog variants tested, only the fiducial passes in every representation. Each statistic responds differently: tightening the galaxy-size cut suppresses $B$-mode power in $\xi_-^{\mathrm{B}}$ and harmonic space; adding stellar-halo masks causes a failure in $\xi_-^{\mathrm{B}}$; and COSEBIs ($n \leq 6$) pass for all versions, even though the full-range measurements show strong low-order structure. These differences arise because the same contamination projects differently into angular, modal, and multipole representations; demanding that all three frameworks pass simultaneously places tighter constraints than any single null test. Papers III and IV adopt these fiducial cuts for cosmological inference.

Because the underlying contamination is present in the data regardless of how it is projected, a $B$-mode failure in any one framework warrants scrutiny across all of them. Rather than modeling $B$-modes as a bias to subtract, we use them as a diagnostic: scale cuts that remove $B$-mode contamination also remove any corresponding contamination of the cosmological signal, without assuming a model for the $E$/$B$ ratio. Whether the same CCD-scale effects arise in *Euclid* and the Legacy Survey of Space and Time (LSST), with their different detector geometries and dithering strategies, remains to be seen. As surveys become systematics-limited, disagreements between $B$-mode estimators will carry as much information as their agreements; multi-framework validation will be needed to distinguish instrumental contamination, pipeline effects, and residual astrophysical signals.

(data-availability .unnumbered)=
# Data availability

A subset of the raw data underlying this article is publicly available via the Canadian Astronomical Data Centre at <http://www.cadc-ccda.hia-iha.nrc-cnrc.gc.ca/en/megapipe/>. The remaining raw data and all processed data are available to members of the Canadian and French communities via reasonable requests to the principal investigators of the Canada-France Imaging Survey, Alan McConnachie and Jean-Charles Cuillandre. A public release of the shear catalog will follow the acceptance of Paper I; the products of this paper are available from the authors on request.

(acknowledgements .unnumbered)=
# Acknowledgements

We are honored and grateful for the opportunity of observing the Universe from Maunakea and Haleakala, which both have cultural, historical and natural significance in Hawaii. This work is based on data obtained as part of the Canada-France Imaging Survey, a CFHT large program of the National Research Council of Canada and the French Centre National de la Recherche Scientifique. Based on observations obtained with MegaPrime/MegaCam, a joint project of CFHT and CEA Saclay, at the Canada-France-Hawaii Telescope (CFHT) which is operated by the National Research Council (NRC) of Canada, the Institut National des Sciences de l'Univers (INSU) of the Centre National de la Recherche Scientifique (CNRS) of France, and the University of Hawaii. This research used the facilities of the Canadian Astronomy Data Centre operated by the National Research Council of Canada with the support of the Canadian Space Agency. This research is based in part on data collected at Subaru Telescope, which is operated by the National Astronomical Observatory of Japan. Pan-STARRS is a project of the Institute for Astronomy of the University of Hawaii, and is supported by the NASA SSO Near Earth Observation Program under grants 80NSSC18K0971, NNX14AM74G, NNX12AR65G, NNX13AQ47G, NNX08AR22G, 80NSSC21K1572 and by the State of Hawaii. CD and MK acknowledge support from the Agence Nationale de la Recherche (ANR-22-CE31-0014-01) TOSCA project. This work was made possible by utilizing the CANDIDE cluster at the Institut d'Astrophysique de Paris. The cluster was funded through grants from the PNCG, CNES, DIM-ACAV, the *Euclid* Consortium, and the Danish National Research Foundation Cosmic Dawn Center (DNRF140); it is maintained by Stephane Rouberol. The authors acknowledge the use of the Canadian Advanced Network for Astronomy Research (CANFAR) Science Platform operated by the Canadian Astronomy Data Centre (CADC) and the Digital Research Alliance of Canada (DRAC), with support from the National Research Council of Canada (NRC), the Canadian Space Agency (CSA), CANARIE, and the Canada Foundation for Innovation (CFI). LWKG thanks the University of Edinburgh School of Physics and Astronomy for a postdoctoral Fellowship. AG acknowledges the support of a grant from the Simons Foundation (Simons Investigator in Astrophysics, Award ID 620789). HH is supported by a DFG Heisenberg grant (Hi 1495/5-1), the DFG Collaborative Research Center SFB1491, an ERC Consolidator Grant (No. 770935), and the DLR project 50QE2305. MJH acknowledges support from NSERC through a Discovery Grant. LVW acknowledges support from NSERC through a Discovery Grant. We would like to thank our external blinding coordinator, Koen Kuijken. We thank Douglas Scott for a thorough reading and detailed scientific and editorial comments. This work made use of large language models (Claude, Anthropic; Codex, OpenAI) for code development, data analysis, and manuscript preparation. All AI-assisted output was reviewed and validated by the lead author.

[^1]: <https://github.com/aguinot/cosmo-numba>

[^2]: The most pessimistic alternative treats the two scale-cut boundaries as fitted parameters and reduces $\nu$ by two for every statistic. The minimum PTE at the fiducial cuts, from $\xi_-^{\mathrm{B}}$, then falls from {astra:value col=display where="macro=fiducialMinPte"}`outputs.paper_claims` to {astra:value col=display where="macro=fiducialMinPteNuPenalty"}`outputs.paper_claims`, just below the threshold of {astra:value col=display where="macro=mockMinPteFloor"}`outputs.paper_claims` calibrated at these cuts.
