<!-- Transpiled from paper/unions_bmodes.tex L215-232 (section "Data", id: data).
     No claim macros or figure embeds in range (numbers below are hard-typed
     in the TeX source, not claims_macros.tex macros).
     Paper-number macros expanded per TeX preamble: \papercatalog=I,
     \paperbmodes=II, \paperconfig=III, \paperharmonic=IV. -->

(data)=
# Data

UNIONS is the largest deep weak-lensing survey of the northern sky, with $r$-band depth of $r \approx 24.2$ at $10\,\sigma$ {cite:p}`gwyn.etal25`.
We analyze four variants of the `ShapePipe` v1.4 catalog (Paper I, Table H.1).
The initial catalog uses a loose size cut ($r_{\mathrm{h, gal}}/r_{\mathrm{h, psf}} > 0.5$).
For the size-cut catalog, we tighten this threshold to $r_{\mathrm{h, gal}}/r_{\mathrm{h, psf}} > 0.707$, removing galaxies where shape measurement and leakage correction are least reliable.
For the masked catalog, we retain the tightened size cut and add bright-star ($r < 12$) and faint-star ($12 < r < 18$) halo masking, reducing sky coverage by 13%.
The relaxed-flags catalog uses the same sample selection as the size-cut catalog but relaxes the `SExtractor` flags criterion to include deblended objects, increasing the effective galaxy density by 23% (Paper I).
We introduced these refinements after blinded validation tests showed early signs of systematic contamination (Paper I).
We adopt the size-cut catalog as the fiducial for Papers II, III, and IV.

For each catalog version, we produce shear measurements both with and without the object-wise PSF-leakage correction described in Paper I, which models PSF contamination via binned regression in signal-to-noise and galaxy-PSF size ratio; {cite:t}`guerrini.etal25` provide the broader PSF-leakage diagnostic framework used to assess its impact.
We report results from the leakage-corrected catalogs, which are used for cosmological inference in Papers III and IV.

<!-- TeX comment preserved: Table moved to Paper I (Table 2) as the canonical
     reference for survey properties. -->
