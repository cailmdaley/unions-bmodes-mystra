## Pure $E$/$B$ decomposition

We decompose the measured correlation functions $\xi_\pm(\theta)$ into pure $E$-mode, $B$-mode, and ambiguous components following {cite:t}`schneider.etal22`:

```{math}
:label: eq-decomp-plus
\xi_+(\theta) = \xi_+^{\mathrm{E}}(\theta) + \xi_+^{\mathrm{B}}(\theta) + \xi_+^{\mathrm{amb}}(\theta);
```

```{math}
:label: eq-decomp-minus
\xi_-(\theta) = \xi_-^{\mathrm{E}}(\theta) - \xi_-^{\mathrm{B}}(\theta) + \xi_-^{\mathrm{amb}}(\theta).
```

Ambiguous modes $\xi_\pm^{\mathrm{amb}}(\theta)$ arise from shear patterns that cannot be assigned unique $E$- or $B$-mode origin on a finite angular interval; they take the functional forms $\xi_+^{\mathrm{amb}}(\theta) = a + b\theta^2$ and $\xi_-^{\mathrm{amb}}(\theta) = c\theta^{-2} + d\theta^{-4}$.
They therefore appear on large scales in $\xi_+$ and on small scales in $\xi_-$, as seen in {astra:ref}`outputs.fig_pure_eb_decomposition`.
The pure-mode correlation functions are computed from integral transforms over the measured $\xi_\pm(\theta)$:

```{math}
:label: eq-xip-EB
\xi_+^{\mathrm{E,B}}(\theta) = \frac{1}{2}\left[\xi_+(\theta) \pm \xi_-(\theta) \pm \int_\theta^{\theta_{\max}} \frac{\mathrm{d}\theta'}{\theta'}\,
\xi_-(\theta')\left(4-\frac{12\theta^2}{\theta'^2}\right)\right]
- \frac{1}{2}\left[S_+(\theta) \pm S_-(\theta)\right];
```

```{math}
:label: eq-xim-EB
\xi_-^{\mathrm{E,B}}(\theta) = \frac{1}{2}\left[\xi_+(\theta) \pm \xi_-(\theta) + \int_{\theta_{\min}}^\theta \frac{\mathrm{d}\theta'\,\theta'}{\theta^2}\,
\xi_+(\theta')\left(4-\frac{12\theta'^2}{\theta^2}\right)\right]
- \frac{1}{2}\left[V_+(\theta) \pm V_-(\theta)\right].
```

Here, the $+$ sign gives $E$-modes and the $-$ sign gives $B$-modes.
The auxiliary functions $S_\pm(\theta)$ and $V_\pm(\theta)$ enforce boundary conditions that remove ambiguous modes: $\xi_+^{\mathrm{amb}} = S_+$ and $\xi_-^{\mathrm{amb}} = V_-$ {cite:p}`schneider.etal22`.

The integral transforms require finely binned $\xi_\pm$ over an extended angular range.
We evaluate these transforms on a separate 1000-bin logarithmic grid spanning $\theta = 0.5$–$300$ arcmin; the extended range reduces edge effects in the boundary functions $S_\pm$ and $V_\pm$.
We apply the {cite:t}`schneider.etal22` filter kernels to this grid using the [`cosmo_numba`](https://github.com/aguinot/cosmo-numba) implementation {cite:p}`guinot.etal26`.
{astra:ref}`Fig. %s <outputs.fig_pure_eb_decomposition>` shows the pure $E$/$B$ decomposition for the fiducial catalog.
