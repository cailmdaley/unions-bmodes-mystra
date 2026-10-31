# Implementation notes

Active paper/code discrepancies, schema limitations, and prompt/tooling issues that affect SPECIFY work. This is a thin discrepancy layer — not a second narrative system. Each entry should resolve into the spec, the baseline universe, or an explicit out-of-scope discharge as the loop progresses; resolved entries get removed.

## Active

### Recipe `inputs:` semantics — declared outputs of the same analysis only

Empirically (against astra-tools 0.2.3 / astra-spec 0.0.6), `recipe.inputs` accepts **only** declared outputs of the same analysis. Cross-sub-analysis qualified syntax (`<sub>.outputs.<id>` and `<sub>.<id>`) is rejected with `[INVALID_RECIPE_INPUT]`. The aspirational wiring described in the original note (e.g. `pte_table_results.recipe.inputs: [null_tests.outputs.config_space_pte_evidence]`) is therefore not currently expressible. Workaround used throughout the drafted tree: drop `recipe.inputs` for cross-sub-analysis dependencies; the analysis-level `inputs:` block with `from: <sub>.<output>` carries the declared dependency. Local recipe-level `inputs:` continue to work (e.g. `pure_eb_data_vector_figure.recipe.inputs: [pure_eb_semianalytic_data]` within the same analysis). The root-level pte_table outputs continue to carry empty `inputs:` lists; if cross-sub-analysis recipe wiring is later supported by the validator, they should point at `null_tests.outputs.config_space_pte_evidence` etc. Upstream issue against `astra-spec` to be drafted next turn.

### RULED 2026-08-30 — `fixed_operator` adopts the single-grid transform-averaged geometry

The two-grid limitation below is resolved by ruling, from the evidence, that
`pure_eb_integration_method: fixed_operator` is a **compound option**: it
selects the fixed-quadrature operator *and* moves every consuming chain to the
single-grid transform-averaged geometry (evaluate the transform on the
fine-grid nodes, average the pure modes into the reporting bins, tmin/tmax at
the fine-grid extent). Grounds: (1) this is the exact
geometry the operator validation and the rectification MC ran — extending the
operator to off-grid evaluation would be new unvalidated numerical code
reopening the gap the validation just closed; (2) the campaign findings prose
already describes the transform-averaged geometry; (3) pointwise evaluation at
reporting `meanr` is itself an approximation of the bin-averaged estimand.
Encoded in the root `fixed_operator` option description. `adaptive_dqags`
keeps the published pointwise geometry, byte-identical to the paper.
Implementation: geometry resolves from the method inside the shared
`pure_eb_transform.py` layer / consuming scripts.
Implemented in commit b060964 (shared `AveragedPureEBTransform` in
`pure_eb_transform.py`; smokes: operator vs dqags agree to 7.7e-5 of max on
smooth synthetic input, guard exception intact).

**Amended 2026-08-30** — the option is compound in geometry only. The *weights*
of the fine-to-reporting average were welded into it and are now their own root
decision, `pure_eb_rebinning_weights` (`npairs`, the Asgari et al. 2019 Appendix A
pair-count weighting, adopted 2026-09-07 after the mock decision scan; `uniform`
the unweighted mean it replaced; `none` the pointwise chain), reaching every consuming script as a required
`--rebinning-weights` flag. `pure_eb_transform.check_rebinning_weights` enforces
the one coupling that remains: `none` in the pointwise geometry, `uniform` or
`npairs` in the averaged one.
## Notes for the IMPLEMENT phase (downstream)

_(empty — populated as the SPECIFY loop encounters issues that the implement phase will need to know)_
