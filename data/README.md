# Inputs

The analysis's external inputs, annexed (content travels with `git annex get`, never through GitHub).

- `v1.4.x/v{version}/unions_shapepipe_cut_struc_2024_v{version}.fits`: the ShapePipe v1.4 cut catalogues,
  copied unchanged from `vos:cfis/weak_lensing/ShapePipe/v1.4.x/v{version}/`. The v1.4.5 file is the released one, not
  `_DEPRECATED`.
- `nz/nz_SP_v1.4.6.3_{A,B,C}.txt`: the v1.4.6.3 SOM n(z) per blind, made from
  `arc:projects/unions/lensing/ShapePipe/v1.4.x/v1.4.6.3/blind_nz_shapepipe_v1.4.6.3.fits` by sp_validation's
  `scripts/get_nz_txt_files_from_fits.py` (SOM-weighted histogram of `Z_{blind}`, 200 bins on 0 ≤ z ≤ 5, normalised).
- `v1.4.x/v1.4.6/unions_shapepipe_cut_struc_2024_v1.4.6.fits`: the v1.4.6 cut catalogue, copied unchanged from
  `vos:cfis/weak_lensing/ShapePipe/v1.4.x/v1.4.6/`. Its galaxies' nside-4096 pixels are the footprint of the
  published GLASS suite (2757.25 deg²); the data-matched suite uses the fiducial catalogue's.
- `survey_areas.yaml` (plain git): each catalogue's survey area in deg², the values Paper I and its `cat_config.yaml`
  use. They come from a high-resolution HEALPix footprint that combines the structural masks with a binned
  comprehensive catalogue; sp_validation's scripts could remake it, but not from data on Fir. The galaxies themselves
  occupy less (2750.65 deg² for v1.4.6.3 at nside 4096), so the area is a convention every CosmoCov covariance and
  quoted density shares.
