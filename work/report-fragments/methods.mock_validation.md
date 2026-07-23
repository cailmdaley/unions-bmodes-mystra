## Validation on simulations

We validate the analysis pipeline on {astra}`350 lognormal simulated catalogs <inputs.glass_mock_realisations>` generated with the [`GLASS`](https://github.com/glass-dev/glass) package {cite:p}`tessore.etal23`, generated with the same survey footprint, galaxy positions, weights, and Planck 2018 cosmology as the UNIONS data.
These simulated catalogs contain no PSF-related systematics, so they provide a direct check that our estimators do not generate spurious $B$ modes when applied to a noise-dominated signal.
Across the ensemble, the pure-mode decomposition, COSEBIs, and harmonic-space power spectra all recover $B$ modes consistent with zero;
for COSEBIs specifically, the mean $B_n$ over a subset of 100 realizations is biased by less than $0.3\,\sigma$ per mode.
We also find that the scatter in the simulated realizations is comparable to the fiducial error bars adopted for each statistic, indicating that our uncertainty estimates are of the right order.

The companion cosmology papers validate their pipelines with the same simulation suite: Paper III checks the configuration-space pipeline and [`CosmoCov`](https://github.com/CosmoLike/CosmoCov) covariance, while Paper IV checks the harmonic-space pipeline and [`NaMaster`](https://github.com/LSSTDESC/NaMaster) covariance.
Here we {astra}`apply the simulations <covariance.decisions.glass_validation_strategy>` to the $B$-mode statistics specifically, complementing the cosmological-inference tests in those papers.
Taken together, these tests show that the low-level $B$-mode structure discussed below is a property of the data rather than an artifact of the estimators or their numerical implementation.
