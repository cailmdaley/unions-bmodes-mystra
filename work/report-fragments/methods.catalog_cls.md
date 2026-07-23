<!-- Fragment: methods.catalog_cls (TeX lines 354–371) — §3.4 Catalog-based harmonic-space power spectra -->

### Catalog-based harmonic-space power spectra

We estimate $C_\ell^{EB}$ and $C_\ell^{BB}$ [](#eq-cl_definition) directly from the discrete catalog, providing $\ell$-dependent $B$-mode diagnostics that complement the angular-space statistics.
The {cite:t}`wolz.etal25` estimator represents the survey mask as delta functions at source positions and computes the spin-2 transforms on this irregular grid using the `ducc` library {cite:p}`reinecke.seljebotn23`.
We subtract an analytic shot-noise contribution following the implementation in Paper IV.
We bin into 32 bandpowers with square-root spacing ($\ell^{0.5}$) over $8 \leq \ell \leq 2048$.
{astra:ref}`Fig. %s <outputs.fig_cl_fiducial>` shows the measured $C_\ell^{BB}$ and $C_\ell^{EB}$ for the fiducial catalog.

:::{astra} outputs.fig_cl_fiducial
:caption: Harmonic-space power spectra for the fiducial catalog, normalized by uncertainty. Top: $B$-mode auto-power $C_\ell^{BB}/\sigma$. Bottom: $E$-$B$ cross-power $C_\ell^{EB}/\sigma$. $C_\ell^{BB}$ is predominantly positive across the fiducial range. $C_\ell^{EB}$ scatters symmetrically around zero, suggesting that any $B$-mode contamination is not strongly correlated with the lensing signal.
:::
