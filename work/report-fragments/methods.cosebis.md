<!-- Fragment: methods.cosebis (TeX lines 301–352) — §3.3 COSEBIs -->

(methods-cosebis)=
### COSEBIs

COSEBIs compress $\xi_\pm(\theta)$ into discrete orthogonal modes $E_n$ and $B_n$, with clean $E$/$B$ separation on finite angular intervals {cite:p}`schneider.eifler.krause10`.
The filter functions $T_{\pm n}(\theta)$ are constructed to satisfy boundary conditions that remove the ambiguous modes described in [](#methods-pure-eb).
Because the $T_{\pm n}$ are orthonormal, different modes capture largely independent information.
Each $T_{\pm n}(\theta)$ oscillates $n+1$ times across the angular interval; higher-$n$ modes are therefore sensitive to increasingly oscillatory features in $\theta$-space.
Conversely, a feature localized in angle projects onto many modes at once, producing an oscillatory signature across the basis, as in the full-range measurement of {astra:ref}`outputs.fig_cosebis_fiducial`.

We compute the COSEBI mode amplitudes $X_n \in \{E_n, B_n\}$ over a given angular range $[\theta_{\mathrm{min}}, \theta_{\mathrm{max}}]$ arcmin using logarithmic filter functions:

```{math}
:label: eq-cosebi-xn
X_n = \frac{1}{2} \int_{\theta_{\mathrm{min}}}^{\theta_{\mathrm{max}}} \mathrm{d}\theta \, \theta \left[ T_{+n}(\theta) \xi_+(\theta) \pm T_{-n}(\theta) \xi_-(\theta) \right],
```

where the upper (lower) sign gives $E_n$ ($B_n$).
We evaluate these integrals on the same 1000-bin grid used for the pure-mode transforms ([](#methods-pure-eb)).
As a convergence check, we repeated the calculation on a 10\,000-bin grid for the fiducial catalog and found the $B_n$ amplitudes and null-test PTEs to be stable.

:::{astra} outputs.fig_cosebis_fiducial
:caption: COSEBI $B$-mode amplitudes for the fiducial catalog, normalized by uncertainty ($B_n/\sigma_n$). Orange squares: full angular range $1$--$250$ arcmin; blue circles: fiducial scale cuts {astra:value}`decisions.fiducial_scale_cuts` arcmin. The shaded region highlights mode numbers $n \leq 6$, which capture nearly all cosmological information. On the full range, $n=1$ shows a $>4\,\sigma$ excess, with oscillatory structure through all twenty modes, consistent with repeating additive shear bias ([](#discussion)). Scale cuts suppress the first-mode excess to below $1\,\sigma$ and largely eliminate the oscillatory signature (see [](#tab-pte-results)).
:::

{astra:ref}`outputs.fig_cosebis_fiducial` shows COSEBI $B$-mode amplitudes for the first $n_{\mathrm{max}} = 20$ modes of the fiducial catalog.
For the UNIONS data vector, the $E$-mode signal saturates by $n \approx 5$--$6$ (see [](#discussion)), so we use $n \leq 6$ as the set of modes carrying most of the cosmological information.
We also report $B$-mode probability-to-exceed (PTE) statistics for $n \leq 20$; at full angular range, the oscillatory pattern extends across the full mode range, characteristic of the repeating additive shear bias identified by {cite:t}`asgari.etal19a` (see [](#discussion)).

The COSEBI coefficients can reconstruct pure-mode correlation functions ([](#methods-pure-eb)):

```{math}
:label: eq-xi-xeb-reconstruction
\xi_\pm^{X}(\theta) = \frac{\bar{\theta}^2}{\beta_\theta} \sum_{n=1}^{n_{\mathrm{max}}} X_n \, T_{\pm n}(\theta),
```

where $\bar{\theta} = (\theta_{\mathrm{min}} + \theta_{\mathrm{max}})/2$ and $\beta_\theta = (\theta_{\mathrm{max}} - \theta_{\mathrm{min}})/(\theta_{\mathrm{max}} + \theta_{\mathrm{min}})$.
In practice, we compute the pure-mode correlation functions directly from the finely binned $\xi_\pm$ measurements rather than reconstructing them from the COSEBI coefficients.
COSEBIs require consistent angular ranges for $\xi_+$ and $\xi_-$, while the pure-mode integral transforms [](#eq-xip-eb) and [](#eq-xim-eb) permit independent scale cuts for each.

The COSEBI modes can equivalently be computed from harmonic-space power spectra:

```{math}
:label: eq-cosebi-harmonic
X_n = \sum_i \frac{\Delta\ell_i \, \ell_i}{2\pi} \, W_n(\ell_i) \, C_{\ell_i}^{XX},
```

where $W_n(\ell)$ are the harmonic-space COSEBI filter functions (Fourier duals of $T_{\pm n}(\theta)$, computed via FFT-log) and the sum runs over 96 square-root-spaced bandpower bins.
Each $W_n(\ell)$ oscillates $n+1$ times, mirroring its real-space counterpart.
The number of bandpower bins balances two requirements: too few cannot resolve the filter oscillations at high $n$, while too many lead to an unstable mode-coupling matrix deconvolution.
We chose 96 bins and validated the transform on GLASS simulations, confirming that modes $n \leq 6$---the range carrying nearly all cosmological $E$-mode information---agree between the two paths.
The comparison on data is presented in [](#discussion) and {astra:ref}`outputs.fig_harmonic_config_cosebis`.
