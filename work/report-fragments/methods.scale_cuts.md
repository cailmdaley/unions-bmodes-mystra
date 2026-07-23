<!-- Fragment: methods.scale_cuts (TeX lines 421–438) — no claim macros or figures in range; footnote numbers (0.18, 0.09, 0.05) are hard-typed in the TeX, not macro-indirected -->

## $B$-mode significance

We assess $B$-mode significance across a grid of angular scale-cut combinations.
For each scale range, we compute a $\chi^2$ statistic testing the $B$-mode data vector $\vec{d}_B$ against zero:

```{math}
:label: eq-chi2-bmode
\chi^2_B = \vec{d}_B^{\mathsf{T}} \, \mathbf{C}_{BB}^{-1} \, \mathbf{d}_B,
```

where $\mathbf{C}_{BB}$ is the $B$-mode covariance submatrix for the selected scale range (corrected for finite-sample bias in the MC-propagated pure-mode covariance; {astra}`covariance`).
The probability-to-exceed (PTE) is then $P(\chi^2 > \chi^2_B \,|\, \nu)$ for $\nu$ degrees of freedom equal to the number of data points in the range.
Since the scale-cut boundaries are themselves selected using the PTE grid, one could treat them as two fitted parameters and reduce $\nu$ accordingly.[^dof]
For the pure-mode decomposition, we compute PTEs separately for $\xi_+^{\mathrm{B}}$ and $\xi_-^{\mathrm{B}}$, as well as a joint test $\xi_{\mathrm{tot}}^{\mathrm{B}}$ using the concatenated data vector and full cross-covariance.
For COSEBIs, we evaluate PTEs using both $n \leq 6$ and $n \leq 20$ ({astra}`cosebis`).
For $C_\ell^{BB}$, we apply the same framework across multipole ranges.
We adopt a uniform threshold of PTE $=0.05$ for all tests.
We do not attempt an explicit correction for look-elsewhere effects across the many scale-cut combinations, because the tests are strongly correlated both within and across statistics; we instead require the adopted cuts to lie in broad stable regions that pass across all three frameworks.

[^dof]: Doing so lowers the minimum PTE across all statistics from 0.18 to 0.09, still above the 0.05 threshold.
