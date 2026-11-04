# Decision-wiring identity criteria and source audit

The reference is the dereferenced pre-wiring results snapshot at `_scratch/decision-wiring-before/{baseline,submitted}` (starting checkout `716120f7`).
`tools/prove_decision_identity.py` compares every changed manifest's artifact against that snapshot, including explicitly paired neutral-name replacements.
The committed CSV is `work/notes/decision-wiring-identity.csv`.
`--rerun-targets work/notes/decision-wiring-review-reruns.txt` additionally requires every review rerun to appear in the CSV, even if its manifest is unchanged; a missing target fails the proof.

## Criteria

- Byte identity passes immediately.
- Floating arrays must satisfy max(abs(new−old))/max(abs(old)) ≤ 1e-12, independently for each array.
- Only the explicitly audited unchanged numerical paths below may instead satisfy ≤ 1e-4, labelled `PASS-arch`, under the orchestrator's 3 November ruling.
- Integer/string arrays, shapes, dtypes and nonfinite values must match exactly under both criteria.
- Archive timestamps and FITS writing dates are not scientific payloads; archive members and FITS metadata otherwise match.
- Matplotlib PDF objects are compared after removing creation/modification dates and generated six-letter font-subset identifiers; compression lengths and xref offsets are serialization bookkeeping.
  Drawing commands, coordinates, labels, font programs (apart from those generated identifiers), object identities and trailer remain exact; no visual or pixel tolerance is used.
- JSON comparisons, including a scalar `provenance` JSON string inside an NPZ, ignore only the `generated` timestamp and compare floating payloads numerically; all other nonfloating values remain exact.
- Protected mock ξ measurement manifests must remain unchanged.

New outputs (`catalog.selected_config` and `cosebis.primary_null_test_report`) have no old artifact and are labelled `new_output`, not identity passes.
Their defaults were also checked directly: both selected configurations equal the original resolved size-cut profile.
The new primary report's calibration fields (apart from its descriptive wording) match the original six-test adopted-range calibration exactly in baseline and to 6.70e-19 relative in submitted.
Multiblock TreeCorr text is checked one numeric column at a time, so its ~1e15 pair counts cannot hide a change in ~1e-5 shear correlations.

## Completed review proof

The review made 94 outputs (47 per universe), all explicitly required by `work/notes/decision-wiring-review-reruns.txt`.
All 94 pass strictly against the original snapshot: 59 byte-identical and 35 exact-payload or ≤1e-12 floating-payload passes.
Their maximum relative difference is 4.888487030639471e-15, in the submitted Cℓ sweep bundle.
The complete CSV still covers 253 products: 249 existing outputs pass strictly and four are new, with zero FAIL and zero PASS-arch.
Its overall maximum remains 5.874194113632341e-13 in the prior baseline integration ξ / bundle, neither of which was rerun for this review.
No protected mock ξ measurements or external mock products were rerun.

## Unchanged numerical paths eligible for PASS-arch

The allowlist is explicit in `tools/prove_decision_identity.py`, restricted to submitted pure-mode products.
These comparisons are against original artifacts, never against a recomputed reference.
The following commands inspect the changes relative to the starting checkout:

```bash
git diff 716120f7 -- src/papers/bmodes/scripts/precompute_pure_eb_chunk.py src/papers/bmodes/scripts/gather_pure_eb_chunks.py src/papers/bmodes/scripts/pure_eb_transform.py
git diff 716120f7 -- src/papers/bmodes/scripts/run_pure_eb_semianalytic.sh src/papers/bmodes/scripts/run_pure_eb_sweep.sh src/papers/bmodes/scripts/calculate_pure_eb_ptes.py
```

The first diff is empty.
The second contains sample-count flag parsing/forwarding on covariance producers, the required-argument check, the neutral covariance filename, and removal of the dead PTE sample-count override.
PTEs read the Hartlap count from the covariance artifact, as before for every production artifact; unrecorded historical artifacts retain the shared 2000-draw convention.
At the default, it forwards the same 2000 samples through the same 20 chunks with the same seeds and order.
The selected default catalogue configuration retains the same values and upstream data.
No integration, covariance-reduction, χ² or PTE algorithm changed.

| Eligible output | Numerical implementation and wiring-only changes |
|---|---|
| `submitted/pure_eb.pure_eb_semianalytic_data` | Unchanged precompute, transform and gather; driver parses `n_2000` as the old 2000. |
| `submitted/glass_mocks.pure_cov_mockmatched` | Same driver and unchanged numerical functions, same stored mock-matched covariance and ξ inputs. |
| `submitted/pure_eb.pure_eb_sweep` | Same numerical chain; adds sample-count forwarding, selected-config input and neutral covariance archive paths. |
| `submitted/pure_eb.pure_eb_pte_per_cut` | Same PTE calculation and artifact-recorded Hartlap count; the unused sample-count override is removed. |
| `submitted/pure_eb.pure_eb_ptes` | Same PTE implementation and artifact-recorded Hartlap count applied to the version sweep; config/path plumbing and removal of the unused override only. |

The original submitted calculation ran on AMD n31–33 (allocation 1047589).
The initial reruns used Intel Xeon Gold 5220R n23 (allocation 1051354); the final wiring reruns used AMD n31 (1051406) and restored strict identity without PASS-arch.
The review-fix reruns also use AMD n31 (1051480).
The submitted-producing commit `57d82ae4` and starting checkout also have no numerical changes in the precompute/gather/transform chain.
Original and rerun manifests record the same frozen environment hash (`ab2d7cf…`) and image (`9189e621…`) for this Monte Carlo chain.
This is the specific source/architecture evidence for the allowlist; it does not relax the criterion for newly changed computational paths.

## Collected architecture diagnostic

A bounded rerun on n23 used the original snapshot inputs, seeds and sample counts with `NUMBA_CPU_NAME=znver2` and `OPENBLAS_CORETYPE=Zen`.
Its covariance differed from the original submitted covariance by **2.523980460334549e-6 relative**, essentially the same maximum drift as native Intel (**2.523980454722103e-6**).
Against the native Intel covariance, the diagnostic differed by **5.489531309944620e-7 relative**.
Every non-covariance floating array was exactly identical to the original.
An earlier BLAS-only diagnostic also failed to restore identity.
Forcing the AMD code-generation target therefore did not reproduce the AMD-host result; no further architecture experiments were pursued, per the ruling.
These diagnostics did not replace any materialized output.

## Effective mask paths

The binary's source is `/n23data1/n06data/lgoh/scratch/UNIONS/CosmoCov`.
`covs/compute_covariances_real.c` includes `cosmolike_core/theory/covariances_3D.c` and `covariances_fourier.c`.

- `covariances_3D.c:455–539`, `w_mask()`, fills its Gaussian pair-count correction with 1 when the file is shorter than 2π/θ_min: 21,600 rows at 1 arcmin or 43,200 at 0.5 arcmin, versus the supplied 12,288 rows.
- `covariances_3D.c:238–268`, `C_survey_window()`, reads the spectrum and normalizes it by C_0 without that minimum-length check.
- `covariances_3D.c:271–315`, `sum_variance_healpix()` / `survey_variance()`, uses its first 1000 multipoles for SSC when a spectrum exists, otherwise an analytic circular-area window.
- `covariances_fourier.c:38` calls `survey_variance()` in the shear SSC contribution.

Thus the Gaussian integration covariances and Gaussian reporting extract are area-only, but the reporting G+NG diagnostic really uses the footprint through SSC.
`mask_application` governs only that reporting-grid SSC window.
Gaussian+NG is excluded on the 1000-bin integration/version grid as intractable; those recipes no longer carry mask branches or inputs.
Gaussian-only recipes pass no unused spectrum; repairing Gaussian boundary corrections is post-paper work.
Both root universes inherit the same covariance selection: the default `masked` now truthfully means footprint SSC, not a Gaussian mask correction.
The fixed Gaussian mock reference is area-only.

A bounded real-binary diagnostic at 20 reporting bins compared footprint SSC
with analytic circular-area SSC, holding all other inputs fixed.
Removing the spectrum changed the G+NG matrix by **0.007735583627879504 relative**
(max absolute difference / max absolute original entry), with a largest
per-bin diagonal change of **3.061490973785128%**.
The default footprint-SSC G+NG rerun is byte-identical to the original, as are
both Gaussian reporting extracts after dropping their unused spectrum.
The Gaussian no-mask change is therefore a numerical no-op; the SSC selector is not.

The unmasked INI uses the nonexistent filename `no_footprint_mask`, not a blank
value: `covs/init.c:152` ignores sscanf's return count and reuses the preceding
n(z) filename for an empty `c_footprint_file`.
The initial empty-value diagnostic was stopped and discarded before this
corrected comparison; no such diagnostic artifact entered `results/`.
