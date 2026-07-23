(results)=
# Results

:::{astra} outputs.fig_pure_eb_versions
:caption: >
  $B$-mode data vectors across catalog versions: initial (red circles),
  size-cut/fiducial (gold squares), masked (teal diamonds), and relaxed-flags
  (purple triangles). Rectangles delineate each bin in the horizontal
  direction and span the range of values across versions in the vertical
  direction; horizontal lines mark the fiducial value. Shaded bands mark
  fiducial scale cuts for $\xi_\pm$ and $C_\ell$. *Top:* Pure $E$/$B$-mode
  correlation functions. Upper panels show $\xi_+$ and $\xi_-$ scaled by
  $\theta$; lower panels show $B$-mode components
  $\xi_\pm^{\mathrm{B}}/\sigma$. *Middle:* Harmonic-space power spectra
  $C_\ell^{BB}/\sigma$ and $C_\ell^{EB}/\sigma$. *Bottom:* COSEBI $B$-mode
  amplitudes $B_n/\sigma_n$ at full range (upper) and fiducial cuts (lower).
  Only the size-cut catalog passes all statistics at the fiducial scale cuts
  ({astra:ref}`Table %s <outputs.pte_table_appendix>`).
:::

% TeX Fig. 5 is a three-PDF composite (pure_eb_versions / cl_versions /
% cosebis_bmode_stacked); outputs.fig_pure_eb_versions materializes all three
% per the transpile map.

:::{astra} outputs.fig_pte_heatmaps
:caption: >
  Configuration-space PTE maps for all four catalog versions (Paper I,
  Table H.1). Columns show $\xi_+^{\mathrm{B}}$, $\xi_-^{\mathrm{B}}$, and
  COSEBI $B_n$ PTEs as a function of angular scale cuts; rows show different
  catalog versions. Solid blue cells indicate PTE $< 0.05$; whites and reds
  indicate consistency with zero. Black squares mark the adopted cuts
  ({astra:value}`decisions.fiducial_scale_cuts`).
  The size cut (fiducial versus initial) expands the acceptance region, while
  stellar masking contracts it.
:::

% Caption binding note: the TeX range "$\ebthetaXipMin$--$\ebthetaXipMax$~arcmin"
% (12--83 arcmin) is carried by the single decision-carried reference
% {astra:value}`decisions.fiducial_scale_cuts`, which renders the adopted
% option's label ("12'--83'") — covering both macros ebthetaXipMin and
% ebthetaXipMax per transpile-map §2f.

:::{astra} outputs.fig_cl_pte_heatmap
:caption: >
  Harmonic-space $C_\ell^{BB}$ PTE maps for all four catalog versions
  (Paper I, Table H.1). Solid blue cells indicate failing PTEs; whites and
  reds indicate consistency with zero. The initial and relaxed-flags catalogs
  show widespread failures across multipole ranges; the fiducial and masked
  catalogs pass across nearly all multipole combinations.
:::

For the fiducial catalog before scale cuts, COSEBIs fail with a $>4\,\sigma$ first mode, while the pure $B$-mode correlation functions and $C_\ell^{BB}$ pass but show low-level $B$-mode structure ({astra:ref}`Table %s <outputs.pte_table_appendix>`).
At the adopted scale cuts, all null tests pass.
We chose those cuts using the $B$-mode tests together with additional systematics checks and blinded inference-stability checks in Papers III and IV; they are more conservative than the PTEs alone would require.
The object-wise PSF-leakage correction (Paper I) shifts individual PTEs by $\lesssim 0.05$ across all catalog versions and statistics; no pass/fail conclusion changes.
We report leakage-corrected results throughout.

The low-level features that drive the full-range behavior leave different signatures in each statistic.
In the pure-mode decomposition ({astra:ref}`Fig. %s <outputs.fig_pure_eb_decomposition>`), $\xi_+^{\mathrm{B}}/\sigma$ shows a broad positive excess on the smallest angular scales that tapers toward zero by $\theta \sim 10$ arcmin;
the corresponding $\xi_-^{\mathrm{B}}$ is quieter on the smallest scales but remains broadly elevated from roughly 7 arcmin to the largest scales.
The COSEBI data vector ({astra:ref}`Fig. %s <outputs.fig_cosebis_fiducial>`) compresses the same structure into a smaller number of modes: on the full angular range, the first mode exceeds $4\,\sigma$ and the higher modes show coherent oscillations consistent with the repeating additive pattern discussed in [Sect. Discussion](#discussion);
after applying the fiducial angular cuts, the first-mode excess drops below 1 $\sigma$ and the oscillatory pattern largely disappears.
In harmonic space ({astra:ref}`Fig. %s <outputs.fig_cl_fiducial>`), $C_\ell^{BB}$ shows a low-level positive offset across much of the fiducial range, with the most conspicuous outliers near $\ell \approx 125$ and $250$, and above $\ell \approx 1600$.
Even so, the fiducial range $300 < \ell < 1600$ passes the null test, while $C_\ell^{EB}$ remains consistent with zero, suggesting that any $B$-mode contamination is not strongly correlated with the lensing signal.

{astra:ref}`Figure %s <outputs.fig_pure_eb_versions>` compares the $B$-mode measurements across catalog versions, and {astra:ref}`Table %s <outputs.pte_table_appendix>` summarizes the PTEs.
Each statistic shows a different pattern of failures.
In configuration space, the initial and relaxed-flags catalogs fail in $\xi_-^{\mathrm{B}}$ and in $\xi_{\mathrm{tot}}^{\mathrm{B}}$ at the fiducial cuts.
The relaxed-flags catalog has a similar acceptance region to the fiducial in $\xi_+^{\mathrm{B}}$, but a narrower one in $\xi_-^{\mathrm{B}}$, with the failure boundary approaching the adopted cuts (PTE $=$ {astra:value col=pte where="version=relaxed_flags statistic=xim cut=fiducial"}`null_tests.outputs.config_space_pte_evidence`).
The masked catalog fails all three pure-mode tests despite passing in harmonic space and COSEBIs ($n \leq 6$).
The PTE maps ({astra:ref}`Fig. %s <outputs.fig_pte_heatmaps>`) show that masking increases the $B$-mode significance on both the largest and smallest scales, although the $\xi_-^{\mathrm{B}}$ PTE falls only marginally below threshold ($0.047$).
In harmonic space ({astra:ref}`Fig. %s <outputs.fig_cl_pte_heatmap>`), the initial and relaxed-flags catalogs fail across most multipole combinations, whereas the fiducial and masked catalogs show broad acceptance regions.

COSEBI $B_n$ tests with $n \leq 6$ pass for all four versions at fiducial cuts, despite pure-mode and harmonic-space failures in three of them ({astra:ref}`Table %s <outputs.pte_table_appendix>`).
Only the fiducial catalog passes all statistics at the adopted cuts.
At full range, all versions show low-order oscillatory structure consistent with the repeating additive pattern discussed below, driving COSEBI ($n \leq 6$) PTEs to $10^{-5}$ or below; at fiducial cuts, all versions pass, although the initial catalog retains a mild oscillatory pattern.

{astra:ref}`Figures %s <outputs.fig_cl_pte_heatmap>` and {astra:ref}`%s <outputs.fig_pte_heatmaps>` map null-test performance over the full grid of scale-cut combinations.
In each representation, the adopted cuts sit well inside broad acceptable regions, so moving them by several bins in either direction does not change the outcome.
For COSEBIs, whose data vector is indexed by mode number rather than angle, the PTE heatmap translates null-test performance back into angular space.
