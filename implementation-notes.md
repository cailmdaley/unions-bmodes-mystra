# Implementation notes

Active paper/code discrepancies, schema limitations, and prompt/tooling issues that affect SPECIFY work. This is a thin discrepancy layer — not a second narrative system. Each entry should resolve into the spec, the baseline universe, or an explicit out-of-scope discharge as the loop progresses; resolved entries get removed.

## Active

### Recipe `inputs:` semantics — declared outputs of the same analysis only

Empirically (against astra-tools 0.2.3 / astra-spec 0.0.6), `recipe.inputs` accepts **only** declared outputs of the same analysis. Cross-sub-analysis qualified syntax (`<sub>.outputs.<id>` and `<sub>.<id>`) is rejected with `[INVALID_RECIPE_INPUT]`. The aspirational wiring described in the original note (e.g. `pte_table_results.recipe.inputs: [null_tests.outputs.config_space_pte_evidence]`) is therefore not currently expressible. Workaround used throughout the drafted tree: drop `recipe.inputs` for cross-sub-analysis dependencies; the analysis-level `inputs:` block with `from: <sub>.<output>` carries the declared dependency. Local recipe-level `inputs:` continue to work (e.g. `pure_eb_data_vector_figure.recipe.inputs: [pure_eb_semianalytic_data]` within the same analysis). The root-level pte_table outputs continue to carry empty `inputs:` lists; if cross-sub-analysis recipe wiring is later supported by the validator, they should point at `null_tests.outputs.config_space_pte_evidence` etc. Upstream issue against `astra-spec` to be drafted next turn.

## Notes for the IMPLEMENT phase (downstream)

_(empty — populated as the SPECIFY loop encounters issues that the implement phase will need to know)_
