# Inputs

The analysis's external inputs, annexed (content travels with `git annex get`, never through GitHub).

- `v1.4.x/v{version}/unions_shapepipe_cut_struc_2024_v{version}.fits`: the ShapePipe v1.4 cut catalogues,
  copied unchanged from `vos:cfis/weak_lensing/ShapePipe/v1.4.x/v{version}/`. The v1.4.5 file is the released one, not
  `_DEPRECATED`.
- `nz/nz_SP_v1.4.6.3_{A,B,C}.txt`: the v1.4.6.3 SOM n(z) per blind, made from
  `arc:projects/unions/lensing/ShapePipe/v1.4.x/v1.4.6.3/blind_nz_shapepipe_v1.4.6.3.fits` by sp_validation's
  `scripts/get_nz_txt_files_from_fits.py` (SOM-weighted histogram of `Z_{blind}`, 200 bins on 0 ≤ z ≤ 5, normalised).
- `glass_mock_v1.4.6.3_v2/`: products of the GLASS mock suite of Paper III, regenerated on fir with the published
  settings (350 realisations, seeds 1000000+N; `analyses/glass-v1.4.6.3-regen`). `mock_xi_fine.npz` and
  `mock_xi_reporting.npz` stack each realisation's 1000-bin and 20-bin TreeCorr ξ±; `cl/cl_glass_mock_*_4096.npy` are
  its NaMaster Cℓ; `unions_glass_sim_00001_4096.fits` is the catalogue of realisation 00001.
