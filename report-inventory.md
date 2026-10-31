# Paper report inventory

The accepted source is `paper/unions_bmodes.tex`.
This inventory records numerical literals that the paper types directly rather than supplying through `claims_macros.tex`.
Repeated occurrences are grouped below; the line references identify the source context.
Numbers in citations, author identifiers, grant identifiers, equation subscripts, and mathematical coefficients are not measured results; they remain in the transcription as written.

## Hard-typed quantities and settings

- **Abstract** (lines 132–136): UNIONS-3500; galaxy-to-PSF size ratio 0.707; multipoles 300–1600; full angular range 1–250 arcmin.
- **Introduction** (lines 153–193): first detections in 2000; $S_8$ convention 0.3 and approximately 3σ tension; KiDS-Legacy $S_8=0.815^{+0.016}_{-0.021}$ and 0.73σ; DES Y6 $S_8=0.798^{+0.014}_{-0.015}$ and 1.1σ; 2% precision; 2–5σ earlier B-mode detections; 6250 deg² survey, 0.7 arcsec seeing, and 1700 and 3500 deg² analysis areas.
- **Data** (lines 212–218): depth 24.2 at 10σ; catalog v1.4; size cuts 0.5 and 0.707; stellar cuts 12 and 18 magnitudes; 13% area reduction, approximately 8 arcmin halo holes, and 26% density increase.
- **Survey table** (lines 243–246): initial catalog 2894 deg², 6.48 arcmin⁻², 0.28; masked 2517 deg², 4.98 arcmin⁻², 0.27; relaxed flags 2894 deg², 6.26 arcmin⁻², 0.27.
  The fiducial row instead uses claims macros.
- **Footprint caption** (lines 256–258): 2.4° × 1.2°, 40 CCDs, approximately 8 arcmin holes.
- **Methods and captions** (lines 281–431): 2 × 2 response; 20 reporting bins over 1–250 arcmin; 1000-bin integration grid over 0.5–300 arcmin; 10,000-bin convergence grid; first 20 COSEBI modes, saturation at approximately 5–6 and tests at $n\leq6$ and $n\leq20$; 96 and 32 harmonic bins; square-root spacing $\ell^{0.5}$ over 8–2048; 20 bins per covariance block; Planck 2018; approximately 5% error increase at $\ell=300$.
- **Null-test table and calibration** (lines 448–500): multipoles 300–1600, full range 1–250 arcmin and 8–2048; 100 arcmin and $\ell\approx100$ covariance-comparison boundaries; 5% clean-sample failure rate and 0.05 single-test PTE; 5% quantile; less than 0.5% blind-to-blind harmonic covariance difference.
- **Results** (lines 507–523): 4σ first COSEBI mode; approximately 10 arcmin taper, approximately 3σ outlier in the 9–12 arcmin bin, approximately 7 arcmin onset, below 1σ after cuts; multipole outliers near 125 and 250 and above 1600; fiducial multipoles 300–1600.
- **Discussion and captions** (lines 528–643): $n\leq6$ COSEBI tests; above 4σ first mode; approximately 3σ in 9–12 arcmin and at least 2σ harmonic outliers; 40 CCDs spanning 6–14 arcmin; 16 arcmin catalog-comparison scale; UNIONS dithers of approximately 1/3 field; approximately 10 arcmin CCD-scale leakage feature; approximately 19 arcmin dither spacing; 5 arcmin alternative inference cut.
- **Conclusions** (lines 653–661): UNIONS-3500 and $n\leq6$.

## Figures without a baseline output

`Figures/footprint_zoom` is omitted from the report image but retains its caption and an explanatory note.
All other `\includegraphics` targets in the accepted TeX have corresponding report output embeds.
The accepted TeX contains no appendix section or appendix `\input`; its appendix-named PTE table is embedded in the main text.
