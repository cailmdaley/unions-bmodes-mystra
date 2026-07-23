# Conclusions

UNIONS-3500 weak-lensing $B$ modes are consistent with zero at the adopted scale cuts across pure-mode correlation functions, COSEBIs, and harmonic-space power spectra.
The cuts lie in broad stable regions of the PTE maps.

On the full angular range, all catalog versions show an oscillatory COSEBI $B$-mode pattern consistent with repeating additive shear bias at CCD angular scales {cite:p}`asgari.etal19a`. This is especially relevant for UNIONS because a similar signature also appeared in CFHTLenS, which likewise used MegaCam on CFHT {cite:p}`heymans.etal12, guinot.etal22`.
Our adopted lower cut suppresses this pattern without discarding the angular scales that carry most of the cosmological signal.
By computing COSEBIs from the harmonic-space bandpowers, we confirm that the disagreement is not a matter of harmonic versus real space: $C_\ell^{BB}$ passes where COSEBIs fail, because the COSEBI filter functions concentrate sensitivity on the contaminated scales.

Of the four catalog variants tested, only the fiducial passes in every representation.
Each statistic responds differently: tightening the galaxy-size cut suppresses $B$-mode power in $\xi_-^{\mathrm{B}}$ and harmonic space; adding stellar-halo masks introduces new pure-mode failures; and COSEBIs ($n \leq 6$) pass for all versions, even though the full-range measurements show strong low-order structure.
These differences arise because the same contamination projects differently into angular, modal, and multipole representations; demanding that all three frameworks pass simultaneously places tighter constraints than any single null test.
Papers III and IV adopt these fiducial cuts for cosmological inference.

Because the underlying contamination is present in the data regardless of how it is projected, a $B$-mode failure in any one framework warrants scrutiny across all of them.
Rather than modeling $B$-modes as a bias to subtract, we use them as a diagnostic: scale cuts that remove $B$-mode contamination also remove any corresponding contamination of the cosmological signal, without assuming a model for the $E$/$B$ ratio.
Whether the same CCD-scale effects arise in *Euclid* and the Legacy Survey of Space and Time (LSST), with their different detector geometries and dithering strategies, remains to be seen.
As surveys become systematics-limited, disagreements between $B$-mode estimators will carry as much information as their agreements; multi-framework validation will be needed to distinguish instrumental contamination, pipeline effects, and residual astrophysical signals.
