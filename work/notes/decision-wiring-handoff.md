# Decision wiring — completed (4 November 2026)

All nine review fixes are applied on `paper-ii-reference`, with one source commit per fix and a follow-up for the direct pure-mode variant-scan caller.
The 94 review reruns, identity CSV and final gates pass.
No agent workers were spawned.
The lc cluster `pii-review-amd` is down and allocation **1051480 / AMD n31** is released.
No materializer remains running.

## Review verdicts

| # | Verdict | Current implementation |
|---|---|---|
| 1 | Applied | `covariance_class.gaussian_plus_ng` is excluded on the 1000-bin integration grid; the integration recipe and seven generated scan members have no NG/mask branch. The shell name mapping is removed. `mask_application` governs only `cov_reporting_ng`. |
| 2 | Applied | `fiducial_catalog.initial`, `masked`, and `relaxed_flags` are excluded as primary selections. Every universe still compares their fixed scan members; the scans and 10k convergence input remain centred on v1.4.6.3. Only `catalogs/size_cut.yaml` exists. |
| 3 | Applied | `fiducial_blind.blind_b` and `blind_c` are excluded. The paper reports A; `null_tests.bb_covariance_blind_independence_evidence` remains the evidence for using one blind. Hard-coded `--blind A` recipes are unchanged. |
| 4 | Applied, live option | `mc_samples.n_1000` is live alongside 2000 and 5000. The pure-mode chunk/gather implementation accepts all three counts. The rationale no longer attributes MC sampling noise to COSEBI modes or claims unmeasured percentage convergence. COSEBI covariances are propagated linearly. |
| 5 | Applied | PTEs read the Hartlap count exclusively from the covariance artifact. The dead `--n-samples` flag, sweep forwarding and PTE-output decision entries are removed. The direct caller in `pure_eb_scan_statistics.py` also drops the removed argument; its producer records the scan's own draw count. |
| 6 | Applied | The unused `mode_subsets` imports are removed from both COSEBI PTE sweep declarations and their regenerated specs. Both mode subsets remain in the exhaustive evidence products. |
| 7 | Applied within generator constraint | Dead `config=None` branches are removed from `bundle_sweep.py`. Gathers still reconstruct the profile plus leakage correction: `sweep_expand.py:337–340` discards `gather.inputs` and replaces it with the member IDs, with no supported extra-input mechanism. |
| 8 | Applied | A supplied CosmoCov spectrum now applies to SSC unless explicitly disabled with `--mask-application unmasked`. Gaussian-only runs still omit the unused spectrum. Tests cover implicit masking and both explicit choices. |
| 9 | Applied | `cl_scan` has no nominal `nbins` matrix column. Both recipe branches and bundle naming use the `grid` role; only the selected harmonic-validation decision supplies the fine-grid count. |

No material review claim was disputed.
For item 1, the external CosmoCov source (`covariances_real_binned_fullsky_nonlimber_w.c:159–169`) nests the angular-bin pair accumulation inside the NG multipole-pair loop.
The roughly 2500-fold increase from 20 to 1000 angular bins is a cost extrapolation, not a weeks-long benchmark; no such run was attempted.

## Current decision space

| Live decision | Governs |
|---|---|
| Root `fiducial_scale_cuts` | `cosebis.primary_null_test_report`, adopted versus full-range primary report. Exhaustive evidence retains both ranges. |
| Root `mc_samples` | Pure-mode MC covariance producers, including the mock-matched covariance: 1000 / 2000 / 5000 draws. Exact-linear covariance uses no draws; PTEs inherit the artifact's recorded count. |
| Root `leakage_correction` | Materialized selected config, catalogue products, data measurements/transforms/scans and downstream summaries. Raw/corrected scan slots stay complete. |
| `cosebis.mode_subsets` | Primary report and recalibration from stored mock PTE columns: both / six-only / twenty-only. |
| Root `harmonic_validation_binning` | Fine harmonic measurements/covariances, version scans and harmonic COSEBI validation: 96 / 128 bins, both comparing six reliable modes. |
| `covariance.mask_application` | Reporting-grid G+NG SSC window only: footprint spectrum versus analytic circular-area approximation. |

Integration/version covariance is Gaussian-only and area-only.
The primary catalogue is size-cut v1.4.6.3 and the reporting blind is A.
Previously excluded alternatives remain excluded: 10k production integration grid, mock-based production Hartlap factor, AD/KS PTE definitions, map-based Cℓ estimator, MC noise subtraction, non-32 reporting bandpowers, all-blinds PTE strategies, 15/30 COSEBI modes, single-pair COSEBI grid, alternative configuration-space covariance engines, independent null-test grids, COSEBI-reconstructed pure modes, and 40-bin reporting ξ.
The corrected-only 10k convergence comparison remains explicitly unavailable for raw shear.

## Reruns and scope

The full explicit scope is **`work/notes/decision-wiring-targets.txt`**, 272 baseline/submitted producer targets.
Both the initial and final checks used that list, with explicit universe prefixes and no `--refresh`.
The pre-edit check was up to date.
After the edits, the check found 60 direct recipe changes and 140 possible downstream reruns, with no failures, blocked outputs or protected mock measurements in the plan.
The largest historical walltime × declared CPU estimate among those candidates was **21.16 core-hours**, below the 500-core-hour stop threshold.

`lc materialize` made **94 outputs**, **47 in each universe**; content identity stopped the remaining planned cascades.
The exact actual rerun list is **`work/notes/decision-wiring-review-reruns.txt`**.
Per universe, it contains:

- 20 `cl_bb.cl_scan` members and their bundle;
- seven `covariance.cov_scan` members and `covariance.cov_integration`;
- `pure_eb.pure_eb_pte_per_cut` and `pure_eb.pure_eb_ptes`;
- three harmonic-COSEBI scan members and their version-PTE gather;
- `cl_bb.cl_pte_per_cut` and its version-comparison figure;
- the mock minimum-PTE calibration and figure;
- three null-test summaries and five root reports/figures.

No pure-mode MC covariance, data ξ measurement, protected GLASS ξ remeasurement, or external mock spectrum/covariance was rerun.
The largest observed rerun walltime × declared CPUs was **0.408 core-hours**; the sum was **6.58 core-hours**.
The final scoped check reports **622 current**, including **364 environment-only behind**, and **zero planned / failed / blocked**.
Behind outputs were deliberately not refreshed.

## Identity proof

`tools/prove_decision_identity.py` compares against the original dereferenced snapshot **`_scratch/decision-wiring-before/{baseline,submitted}`**, from checkout `716120f7`, never against a recomputed reference.
Its new `--rerun-targets` option requires every explicitly rerun output to appear in the CSV, even if its manifest is unchanged; a missing target fails the proof.
Protected mock ξ manifests are asserted unchanged, and requesting one as a rerun fails too.

The regenerated **`work/notes/decision-wiring-identity.csv`** has **253 rows**: 126 baseline and 127 submitted.
All **249 existing outputs pass strictly**, with four genuinely new outputs (two selected configs and two primary reports), **zero FAIL** and **zero PASS-arch**.

| Verdict | Rows |
|---|---:|
| byte_identical | 155 |
| float_arrays_within_1e-12 | 63 |
| archive_float_arrays_within_1e-12 | 4 |
| json_float_payload_within_1e-12 | 9 |
| archive_payload_identical | 8 |
| json_payload_identical | 6 |
| array_payload_identical | 2 |
| pdf_payload_identical | 2 |
| new_output | 4 |

All 94 review reruns are included: **59 byte-identical**, 24 floating-array passes, two floating-archive passes, five floating-JSON passes and four exact archive/JSON payload passes.
Their largest relative difference is **4.888487030639471e-15**, in the submitted Cℓ sweep bundle.
The overall maximum remains **5.874194113632341e-13**, in the unchanged baseline v1.4.5 corrected integration ξ / bundle from the original wiring campaign.
The proof criteria, exact source audit and collected Intel/AMD diagnostic are in **`work/notes/decision-wiring-proof.md`**.
No architecture allowance was needed for either the completed wiring or review reruns.

The original direct checks of the four new products also stand: selected configs equal the original resolved size-cut profile; primary-report calibration fields match the original adopted-range six-test calibration exactly in baseline and to 6.70e-19 relative in submitted, apart from descriptive wording.

## Mask result and architecture lesson

The footprint spectrum is too short for CosmoCov's Gaussian `w_mask()` correction, but SSC reads its first 1000 multipoles without that minimum-length check.
The existing bounded 20-bin real-binary diagnostic found **0.773558%** maximum-relative matrix change and **3.061491%** largest per-bin diagonal change between footprint SSC and the circular-area SSC approximation.
That diagnostic was not rerun: the reporting recipe's explicit flags are unchanged.
Gaussian integration covariances remain area-only; the mask choice is live only on reporting G+NG.
The unmasked INI must name the absent file `no_footprint_mask`, not a blank value, because CosmoCov otherwise reuses the preceding n(z) filename.

Use AMD n31–n33 for bit-stable MC/adaptive-quadrature comparisons with the originals.
Forcing AMD code generation on Intel did not recover the original MC covariance; the source/architecture diagnostic is already collected, not an open investigation.
On this run, Slurm rejected `comp -w n31`: n31–n33 are in **pscomp** (also listed in pscompl).
Allocation 1051480 used `-p pscomp -w n31 -c 32 --mem=200G --time=24:00:00 --exclude=n09,n17,n36`.

## Gates and reproducible commands

All **13 synthetic tests pass**, including 1000/2000/5000 draw forwarding, artifact-owned 500/1000/2000/exact Hartlap handling, implicit/explicit masking, bundle grid naming, and required-rerun identity coverage.
Schema/semantic validation, links, reachability and `git diff --check` pass.
All changed sweep specs were regenerated with the existing expander.
No alternative universes were materialized merely to test option plumbing.

Use installed `lc`, not `uv run --project` on a CLI checkout that changes branches.
The validator is `/n17data/cdaley/.local/share/uv/tools/lightcone-cli/bin/astra`.
The frozen scientific Python is **`.lightcone/venv/bin/python` inside the image**, not `/app/.venv/bin/python`.
After allocating a new job, the proof command is:

```bash
srun --jobid=<jobid> --overlap -n 1 -c 1 env TMPDIR=/tmp \
  OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 \
  apptainer exec --bind /automnt,/n17data,/n23data1,/n09data,/home --pwd "$PWD" \
  .lightcone/images/9189e621d803ff89888c08b8947857d518a9239de1a064de22d38506f0795eab.sif \
  .lightcone/venv/bin/python tools/prove_decision_identity.py \
  --rerun-targets work/notes/decision-wiring-review-reruns.txt
```

Replace the final Python arguments with `tools/test_decision_wiring.py` for the 13 tests.
For lc materialization, launch a local cluster inside the Slurm allocation and keep that launching step alive; otherwise Slurm kills its background cluster processes when the step exits.
End the cluster before releasing the allocation.

## Friction

- **Sweep generator:** `gather.inputs` is overwritten with member IDs, so gathers cannot consume `selected_config`; profile-plus-correction reconstruction remains the supported route.
- **Validation hook:** PATH-only edit hooks cannot find `astra` and print empty validation failures; the absolute-path validator passes.
- **Compute:** AMD nodes require `pscomp`, and an lc local cluster needs its launching Slurm step held open; all scientific commands used `srun --jobid=1051480 --overlap -n 1 -c K`.
- **Freshness:** the 364 environment-only behind outputs are not stale recipes; refreshing them would needlessly touch protected mocks.
