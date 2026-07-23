## Two-point correlation measurements

All $B$-mode diagnostics in this analysis use the discrete shear catalog $\mathcal{C} = \{(\gamma_i, w_i, \hat{\mathbf{n}}_i)\}_{i=1}^N$, comprising complex ellipticities $\gamma_i$, shear weights $w_i$ (Paper I), and sky positions $\hat{\mathbf{n}}_i$ for $N$ source galaxies.
The ellipticities $\gamma_i$ are shear estimates calibrated via metacalibration {cite:p}`huff.mandelbaum17, sheldon.huff17`: galaxy images are artificially sheared to measure the $2 \times 2$ shear response matrix $\mathbf{R}$, yielding $\boldsymbol{\gamma}_i = \mathbf{R}^{-1}(\mathbf{e}_i - \mathbf{c}_i)$ from the observed ellipticity $\mathbf{e}_i$ and additive bias $\mathbf{c}_i$.
The two-point shear correlation functions $\xi_\pm(\theta)$ are formed from products of tangential and cross-component ellipticities for galaxy pairs at separation $\theta$.
We estimate $\xi_\pm$ by weighting each pair by metacalibration weights and accumulating contributions in angular bins:

$$
{\xi}_\pm(\theta_k) = \frac{\sum_{i \neq j} w_i w_j \, \Pi(\theta_{ij}, \theta_k) \left[\gamma_{\mathrm{t},i} \gamma_{\mathrm{t},j} \pm \gamma_{\times,i} \gamma_{\times,j}\right]}{\sum_{i \neq j} w_i w_j \, \Pi(\theta_{ij}, \theta_k)}.
$$ (eq-xi_estimator)

Here $\gamma_{\mathrm{t}}$ and $\gamma_{\times}$ are the tangential and cross components of the ellipticity in the coordinate frame defined by the pair-separation vector, $\Pi(\theta_{ij}, \theta_k)$ is the binning kernel selecting pairs with separations in bin $k$, and $\theta_{ij} = \arccos(\hat{\mathbf{n}}_i \cdot \hat{\mathbf{n}}_j)$.
We evaluate [](#eq-xi_estimator) using the tree-based correlation package `TreeCorr` {cite:p}`jarvis04`.
We bin into 20 logarithmic bins over $\theta = 1$--$250$~arcmin, matching the configuration-space data-vector binning adopted in Paper III.

In harmonic space, the analogous statistics are the angular power spectra $C_\ell^{EE}$ and $C_\ell^{BB}$.
We estimate these using the MASTER algorithm {cite:p}`hivon.etal02`, which recovers the true $C_\ell$ by deconvolving the mode-coupling matrix induced by the survey geometry, as implemented in `NaMaster` {cite:p}`alonso.etal19`.
Traditional implementations operate on pixelized maps, but {cite:t}`wolz.etal25` introduced a catalog-based estimator that evaluates spherical harmonic transforms directly at source positions, avoiding pixelization artifacts.
The same weighted ellipticities from [](#eq-xi_estimator) enter this estimator:

$$
\tilde{C}_\ell^{\alpha\beta} = \sum_{i,j} w_i w_j \, \gamma_i^a \, \gamma_j^b \; \frac{1}{2\ell+1} \sum_m {}_2 Y_{\ell m}^{a\alpha *}(\hat{\mathbf{n}}_i) \; {}_2 Y_{\ell m}^{b\beta}(\hat{\mathbf{n}}_j),
$$ (eq-cl_definition)

where $a, b \in \{1, 2\}$ label the two ellipticity components, ${}_2 Y_{\ell m}^{a\alpha}$ are spin-2 spherical harmonics, and $\alpha, \beta \in \{E, B\}$.
We then decouple the measured bandpowers and correct for noise bias following {cite:t}`wolz.etal25`, yielding unbiased harmonic-space bandpower estimates $\hat{C}_\ell^{\alpha\beta}$.
