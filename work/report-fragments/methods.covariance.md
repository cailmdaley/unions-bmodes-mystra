<!-- Fragment: methods.covariance (TeX lines 373–420) — Covariance estimation -->
<!-- Cross-fragment refs: equations eq:xip_EB / eq:xim_EB are defined in the
     pure_eb fragment (methods.pure_eb); anchors assumed as #eq-xip-eb / #eq-xim-eb. -->

### Covariance estimation

For each $B$-mode statistic, we start from the same covariance prescriptions adopted in the companion cosmology papers.
In configuration space, Paper III uses `CosmoCov` {cite:p}`krause.eifler17; fang.eifler.krause20` to model the covariance of the shear correlation functions $\xi_\pm$.
In harmonic space, Paper IV uses the Gaussian covariance from `NaMaster`'s iNKA framework for the $EB$ and $BB$ spectra, with separate non-Gaussian terms added only for the $EE$ cosmology analysis.
Here we use those same inputs to construct the covariances needed for the derived $B$-mode data vectors.

For $C_\ell^{BB}$ and $C_\ell^{EB}$, the harmonic-space covariance is used directly.
Paper IV validates this $B$-mode covariance against 10,000 Gaussian simulations and finds good agreement away from the lowest multipoles (Appendix C therein).
For the configuration-space statistics, we begin from the finely binned `CosmoCov` covariance of $\xi_\pm(\theta)$ on a 1000-bin logarithmic grid spanning $\theta = 0.5$--$300$ arcmin, using the survey properties of Paper I and the Planck 2018 fiducial cosmology {cite:p}`planck20`.

For COSEBIs, the covariance follows directly from the linear transformation

$$
\mathbf{C}_{B_n B_m} = \mathbf{T}^{\mathsf{T}} \mathbf{C}_{\xi} \mathbf{T},
$$ (eq-cosebi-cov)

where $\mathbf{T}$ encodes the COSEBI filter functions.
For the pure-mode correlation functions, no equivalent direct transformation is available.
Instead, we draw $N_{\mathrm{samples}} = 2000$ realizations of $\xi_\pm(\theta)$ from the same finely binned `CosmoCov` covariance, propagate each realization through the integral transforms in [](#eq-xip-eb) and [](#eq-xim-eb), and compute the empirical covariance of the reconstructed $\xi_\pm^{\mathrm{E/B}}$.
This captures the cross-covariance between $\xi_+$ and $\xi_-$, correlations between $E$/$B$/ambiguous mode types, and angular-bin correlations induced by the overlapping filter support.
When computing $\chi^2$ statistics for PTE assessments, we apply the correction factor $(N_{\mathrm{samples}} - N_{\mathrm{obs}} - 2)/(N_{\mathrm{samples}} - 1)$ to the inverse covariance of the Monte Carlo pure-mode estimate {cite:p}`hartlap07`;
for $N_{\mathrm{samples}} = 2000$ and the data-vector dimensions used here ($N_{\mathrm{obs}} \leq 40$), this amounts to at most a 2% correction.

For all $B$-mode null tests, we retain only the Gaussian part of the covariance.
In configuration space, this keeps the calculation tractable on the finely binned integration grid, and in harmonic space it matches the Gaussian iNKA covariance used for the $BB$ and $EB$ null tests in Paper IV.
This approximation underestimates the true variances because it omits connected trispectrum and super-sample covariance terms.
The reported PTEs are therefore conservative: a null test that passes with Gaussian-only covariance would generally pass with a more complete covariance estimate.
{astra:ref}`pure_eb.outputs.pure_eb_covariance_figure` shows the pure $E$/$B$ covariance structure.

:::{astra} pure_eb.outputs.pure_eb_covariance_figure
:caption: Pure $E$/$B$-mode correlation-function covariance matrix for the fiducial catalog from semi-analytical propagation of 2000 Monte Carlo realizations. The matrix is organized into six blocks: $\xi_+^{\mathrm{E}}$, $\xi_-^{\mathrm{E}}$, $\xi_+^{\mathrm{B}}$, $\xi_-^{\mathrm{B}}$, $\xi_+^{\mathrm{amb}}$, $\xi_-^{\mathrm{amb}}$, with 20 angular bins per block. Off-diagonal correlations arise from integral coupling in the pure-mode filters.
:::

<!-- ebthetaXipMin/ebthetaXipMax bound via decisions.fiducial_scale_cuts;
     {astra:value} renders the adopted option's label ("12′–83′"), not the
     bare integers — a metric/table bridge is needed for exact "12–83 arcmin"
     phrasing (transpile-map §2f queryability nuance). -->
:::{astra} outputs.pte_table_appendix
:caption: '$B$-mode PTE values across catalog versions at fiducial and full-range scale cuts. Fiducial scale cuts are {astra:value}`decisions.fiducial_scale_cuts` in configuration space and $300 < \ell < 1600$ in harmonic space; the full range is $1$--$250$ arcmin / $8 \leq \ell \leq 2048$. Only the fiducial catalog passes all statistics at the adopted cuts. Bold values indicate PTE $< 0.05$ (null-test failure).'
:::
