<!-- Fragment: introduction (TeX lines 156–213) — no claim macros or figure floats in range -->

(introduction)=
## Introduction

Weak gravitational lensing distorts the observed shapes of distant galaxies, encoding information about the intervening matter distribution along the line of sight {cite:p}`blandford.etal91,miraldaescude91,kaiser92`.
First detections came in 2000 {cite:p}`bacon.etal00,kaiser.etal00,vanwaerbeke.etal00,wittman.etal00`.
The scalar gravitational potential produces a shear field that is curl-free to leading order and can be decomposed into gradient-like $E$ modes and curl-like $B$ modes; under the Born approximation, scalar density perturbations source only $E$-mode power {cite:p}`stebbins96,kamionkowski.etal97`.
Higher-order effects (lens-lens coupling, source clustering, and intrinsic alignments) can produce $B$ modes {cite:p}`hilbert.etal09,schneider.etal02b,crittenden.etal02`, but these are orders of magnitude below current sensitivity {cite:p}`cooray.hu02`.
They are likely only accessible to future ultra-high-precision measurements, potentially through cross-correlation techniques involving CMB lensing rotation reconstructions and large-scale-structure tracers {cite:p}`robertson.etal25`.
For galaxy shear measurements at current sensitivity, detected $B$-mode power indicates residual observational systematic effects rather than cosmological signal, though the absence of $B$ modes does not rule out systematic contamination.

Stage-III cosmic shear surveys have measured the structure growth parameter $S_8 \equiv \sigma_8 \sqrt{\Omega_{\rm m}/0.3}$ at percent-level precision, with some reporting up to ${\sim}3\,\sigma$ tension with *Planck* {cite:p}`asgari.etal21,abbott.etal22,li.etal23,planck20`.
The Kilo-Degree Survey (KiDS)-Legacy {cite:p}`wright.etal25` found $S_8 = 0.815^{+0.016}_{-0.021}$ (0.73 $\sigma$ from *Planck*) and the Dark Energy Survey (DES) Y6 {cite:p}`abbott.etal26` found $S_8 = 0.798^{+0.014}_{-0.015}$ ($1.1\,\sigma$); both achieved $2\%$ precision on $S_8$, comparable to CMB constraints from *Planck* alone.
Stage-IV surveys will reduce statistical uncertainties further, shifting the error budget toward systematic effects; current data offer the opportunity to validate these methods before that transition.

Several methods separate $E$- and $B$-mode contributions in shear two-point statistics.
Aperture mass dispersion {cite:p}`schneider96,schneider.etal98` is the foundational real-space method for $E$/$B$ separation.
In practice, however, its evaluation requires knowledge of the shear correlation functions down to zero separation, which is unavailable on a finite observed interval and leads to $E$/$B$ mixing {cite:p}`schneider.eifler.krause10,schneider.etal22`.
Complete Orthogonal Sets of $E$/$B$-mode Integrals (COSEBIs; {cite:t}`schneider.eifler.krause10,asgari.schneider.simon12`) generalize this approach, compressing all $E$/$B$-separable information from $\xi_\pm(\theta)$ on a finite angular interval into discrete orthogonal modes.
A small number of modes captures the cosmological signal; by construction, the COSEBI basis excludes the ambiguous modes that cannot be assigned pure $E$ or $B$ character on a finite interval.
Pure-mode correlation functions $\xi_\pm^{\mathrm{E/B}}(\theta)$ {cite:p}`schneider.etal22` are an alternative real-space representation that decomposes $\xi_\pm$ into $E$-mode, $B$-mode, and ambiguous components through integral transforms derived from the COSEBI basis functions.
Harmonic-space estimators measure $C_\ell^{EE}$ and $C_\ell^{BB}$ directly but require careful treatment of mask-induced $E$/$B$ mixing; catalog-based methods that evaluate spherical harmonic transforms at source positions {cite:p}`wolz.etal25` avoid the pixelization artifacts that complicate standard implementations.
The **HybridEB** Fourier band-power estimator {cite:p}`becker.rozo16` constructs $\ell$-space band-powers from linear combinations of binned $\xi_\pm$, projecting out ambiguous modes to minimize $E$/$B$ mixing.
DES Y6 adopted this estimator {cite:p}`abbott.etal26`, although COSEBIs are more widely used in the Stage-III literature, where they have traced systematic signatures across multiple surveys {cite:p}`asgari.etal19a`.
These statistics represent $B$-mode power in complementary bases, and the mapping between them is not one-to-one, so contamination excluded by scale cuts in one basis may remain visible in another.
Comparing across bases can expose systematics that any single framework would absorb.

Stage-III analyses have found $B$-mode contamination from a range of systematic effects.
COSEBI analyses of the Canada-France-Hawaii Telescope Lensing Survey (CFHTLenS), KiDS-450, and DES-SV detected significant $B$ modes at 2--5 $\sigma$, linked to repeating additive shear bias, PSF leakage, and photometric selection effects {cite:p}`asgari.etal19a`.
Repeating additive shear biases can arise from detector-level defects that imprint a constant ellipticity offset on each CCD, varying across the focal plane; because this pattern is fixed in instrument coordinates, it repeats with every pointing.
The bias is constant within each CCD but jumps at chip boundaries; these discontinuities carry curl, generating $B$-mode power at CCD angular scales.
Hyper Suprime-Cam (HSC) Y3 used harmonic-space power spectra and found significant large-scale $B$ modes {cite:p}`li.etal23,dalal.etal23`. A dedicated PSF-systematics study identified PSF fourth-moment leakage as a leading additive contaminant in the affected region {cite:p}`zhang.etal23`.
KiDS-Legacy tested with COSEBIs, traced an initial failure to astrometric systematics, and applied conservative masking with negligible impact on $S_8$ {cite:p}`wright.etal25`.
DES Y3 reported null detections using pseudo-$C_\ell$ and COSEBIs {cite:p}`gatti.etal21`.
However, {cite:t}`jefferson.etal25` applied the **HybridEB** estimator uniformly to DES-Y3 and HSC-Y3 data and found $B$-mode failures in specific tomographic bins not seen with the estimators originally used.
DES Y6 subsequently adopted the **HybridEB** estimator and passed {cite:p}`abbott.etal26`.
Because each estimator weights angular scales differently and deprojects ambiguous modes differently, the same data can pass one null test and fail another.
The most informative validation requires consistent passage across multiple statistics.

The Ultraviolet Near-Infrared Optical Northern Survey (UNIONS; {cite:t}`gwyn.etal25`) is surveying 6250 deg² of northern sky with multi-band imaging from telescopes in Hawai'i:
the Canada-France Imaging Survey (CFIS) contributes $u$- and $r$-band imaging from CFHT,
Pan-STARRS contributes $i$- and $z$-band data,
and Subaru adds $z$-band imaging through the Wide Imaging with Subaru HSC of the *Euclid* Sky (WISHES) and $g$-band imaging through the Waterloo-Hawai'i Institute for Astronomy $g$-band Survey (WHIGS).
Nearly all Stage-III cosmic shear survey area is in the southern sky; UNIONS provides the first deep, wide-area shape catalog in the north, with $r$-band median seeing of $0\overset{\prime\prime}{.}7$, enabling high-quality shape measurements.
{cite:t}`guinot.etal22` presented the first UNIONS weak-lensing analysis over 1700 deg² using the `ShapePipe` pipeline {cite:p}`farrens.etal22`;
the present analysis expands this to 3500 deg² with an updated catalog ([](#data)).

We present $B$-mode validation tests for the first UNIONS-3500 cosmic shear cosmology analysis, using pure-mode correlation functions $\xi_\pm^{\mathrm{E/B}}(\theta)$, COSEBIs, and harmonic-space power spectra $C_\ell^{BB}$ to define scale cuts in both angular and multipole space within a blinded analysis framework.
This paper is one of five coordinated UNIONS-3500 weak-lensing publications:
catalog construction (Paper I; {cite:t}`hervaspaters.etal26a`);
the present $B$-mode validation (Paper II);
configuration-space cosmological constraints (Paper III; {cite:t}`goh.etal26`);
harmonic-space cosmological constraints (Paper IV; {cite:t}`guerrini.etal26`);
and image simulations (Paper V; {cite:t}`hervaspaters.etal26b`).
We validate the fiducial scale cuts adopted for cosmological inference in Papers III and IV.
[](#data) describes the UNIONS shear catalogs and their evolution.
[](#methods) details the three $E$/$B$-separable statistics and covariance framework.
[](#results) presents measurements and PTE assessments.
[](#discussion) discusses systematic trends and methodological implications.
