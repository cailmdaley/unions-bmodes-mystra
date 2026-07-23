# Paper2ASTRA Prompt Audit

**Scope:** Alignment between Paper2ASTRA phase prompts (summarize, extract_targets, literature, specify, review) and current ASTRA schema (v0.0.5) plus the narrative skill (paper-reproduction mode).

**Focus areas:** Schema drift, narrative skill integration, stale references, coverage gaps for rich input bundles.

**Date:** 2026-04-25  
**Auditor:** Claude Code  
**Context:** pure_eb case — paper exists, code exists, fibers known, baseline cosmology known.

---

## Resolution status (2026-04-28)

The friction-driven subset shipped as **Paper2ASTRA PR #9** "Quality-of-life prompt fixes from three reproduction loops" (5 commits, 3 files, +133/−43). Below is the final disposition of every audit item; further detail in fiber `paper-ii-astra-conversion/cross-reproduction-consolidation/audit-ledger`.

| # | Item | Disposition |
|---|---|---|
| 1.1 | W3C TextQuoteSelector prefix/suffix | Landed in PR #9 (`literature.py`) |
| 1.2 | Reserved entity names | Landed in PR #9 (`specify.py`) |
| 1.3 | Conditional narrative requirement | **Obsoleted** by `03764c3` (PR #9): narrative is delegated to the `/narrative` skill, which teaches this rule. Closed. |
| 1.4 | Decision constraint format | Filed as Paper2ASTRA #6 §1 |
| 1.5 | Anchor grammar | Landed in PR #9 (`specify.py`) + ASTRA #72/#73 (validator-side diagnostic) |
| 2.1 | Paper-reproduction mode unknown to SPECIFY | Resolved by `03764c3` (PR #9): SPECIFY delegates to `/narrative` skill |
| 2.2 | Five-key narrative structure | Resolved by `03764c3` (PR #9): same delegation |
| 2.3 | Paper-reproduction fidelity discipline | Landed in PR #9 (`specify.py`) |
| 2.4 | Reserved-name validation as SPECIFY step | Resolved by `03764c3` (PR #9): rule delegated to skill |
| 3.1 / 5.6 | Stale `astra paper verify-quotes` CLI | Landed in PR #9 (`literature.py`) |
| 3.2 | `work/` directory paths | No-action (confirmed correct) |
| 3.3 | `astra validate` CLI syntax | No-action (confirmed correct) |
| 4.1 | Code+paper+fibers reconciliation | **Out of scope.** UNIONS-specific simulation/fiber input class isn't a generalizable Paper2ASTRA concern. The baseline-direction sub-question of 4.1 filed as Paper2ASTRA #8. |
| 4.2 | Snakemake DAG awareness | **Out of scope.** UNIONS-specific optimization, audit itself flagged as nice-to-have. |

### Filed as Paper2ASTRA issues during/after the loops

- **#6** — SPECIFY/SUMMARIZE prompt gaps (covers 1.4 + 4.3 explicitly)
- **#7** — Prefer arXiv TeX source over PDF for evidence verification (architectural proposal; gates the eq/section-number-from-PDF rule's eventual shape)
- **#8** — When code and paper conflict, what takes priority for `baseline.yaml`? (worldview question; the MIRoRS commit's baseline.yaml direction reversal was excluded from PR #9 pending this resolution)

### Beyond the audit

PR #9 also folded in two MIRoRS-loop additions that aren't in the original audit:
- `literature.py` quote granularity rule (full sentences carry the claim; fragments fail the reader even when they pass fuzzy-match)
- `literature.py` cross-section attribution rule (when methodology is split across paper sections, file separate insights, one per piece)

---

## Landed in `paper-ii-conversion`

The ralph loop has applied the friction-driven subset of this audit to the
Paper2ASTRA `paper-ii-conversion` branch (commits `d33e3c8` + `cd65cd8`,
the second is codex-review fixes). Each landed change names the loop
fiber that evidenced its friction.

| Audit item | Status | Friction source |
|---|---|---|
| 1.1 Prefix/suffix W3C TextQuoteSelector convention | **Landed** (`literature.py`) | `iter6-evidence-verification` (5 prefix realignments) |
| 1.2 Reserved entity names | **Landed** (`specify.py`) | `narrative-anchor-cant-drill` |
| 1.5 Anchor grammar | **Landed** (`specify.py`) | `narrative-anchor-cant-drill` |
| 2.3 Paper-reproduction fidelity discipline | **Landed** (`specify.py`) | `iter6-evidence-verification` (fabricated `asgari17` insight) |
| 5.6 Stale `astra paper verify-quotes` CLI | **Landed** (`literature.py`) | `iter6-evidence-verification` (venv'd validator workflow) |

The remaining items (1.3 conditional narrative; 1.4 `when` / `incompatible_with`
/ `requires` expansion; 2.2 five-key checklist; 4.1 rich input bundle guidance;
4.2 Snakemake DAG awareness; 4.3 baseline-parameter discovery) are not yet
landed because they have not surfaced as friction in this loop. Per the
constitution's friction-only discipline, prompt edits require concrete
evidence from the loop's just-completed work; speculative gardening is out
of scope.

---

## 1. Schema Drift

### 1.1 Evidence quote anchoring (REQUIRED, NOT MENTIONED)

**Status:** Critical gap.

**Finding:** The literature phase (`literature.py`) specifies extraction of TextQuoteSelector evidence with `exact`, `prefix`, and `suffix` fields. The schema (insight.yaml) mandates `prefix` and `suffix` on every TextQuoteSelector.

**In prompts:**
- `literature.py:31–32` explicitly requires prefix/suffix: "Prefix and suffix context (~20-100 chars each) for robust matching"
- `literature.py:77–80` YAML output shows all three fields present

**In schema:**
- `insight.yaml:36–39` shows TextQuoteSelector with `prefix` and `suffix` as optional (`description` only, no `required: true`)

**Status:** The prompts are correctly stricter than the schema. No drift detected, but schema should enforce these fields as required for verification robustness. (Minor schema issue, not a prompt bug.)

---

### 1.2 Reserved entity names (CRITICAL, NOT MENTIONED)

**Status:** Major gap — prompts blind to constraint.

**Finding:** The narrative skill (`SKILL.md:241–248`) defines reserved names that cannot be used as entity IDs: `inputs`, `outputs`, `decisions`, `findings`, `prior_insights`, `analyses`, `options`, `content`, `narrative`.

**In prompts:**
- **None of the five prompts mention reserved names.**
- `specify.py` (lines 19–44) instructs agents to build decision maps and define Inputs/Outputs/Decisions but never warns against ID collisions with reserved words.
- `literature.py` (lines 162–182) builds extraction output with insight IDs and decision links but provides no validation rule.

**Consequence:** An agent could create an entity `id: findings` or `id: decisions`, which would silently shadow narrative anchors. The spec would serialize but break rendering and validation.

**Recommend:** Add to `specify.py` system prompt (before Outputs section):
```
Reserved names (cannot be used as entity IDs): inputs, outputs, decisions, 
findings, prior_insights, analyses, options, content, narrative. If you find 
an entity using one of these, rename it immediately.
```

---

### 1.3 Narrative conditional requirement (MAJOR, PARTIALLY MENTIONED)

**Status:** Underspecified in prompts.

**Finding:** The schema (analysis.yaml:100–114) enforces a conditional requirement: `findings`, `methods`, `inputs`, `outputs` narrative keys are required **only when** their structural counterparts exist on the Analysis node.

**In prompts:**
- `specify.py` (lines 14–17) says "produce the structured ASTRA spec and implementation notes" but does not explain the conditional narrative rule.
- No prompt mentions that a stub analysis with `summary` only is valid, or that `findings` prose should not appear if `Analysis.findings` is empty.
- The SPECIFY prompt delegates narrative to `.claude/guides/decision-guide.md` (line 15), which is project-specific and not version-controlled in these repos.

**Consequence:** Agents may over-fill narratives (writing `findings` on analyses without declared findings, violating schema) or under-fill (omitting `methods` when decisions exist).

**Recommend:** Add to `specify.py` system prompt (before narrative discussion):
```
## Conditional narrative requirement

Narrative sections are schema-optional, but validation enforces a condition:
- `findings` prose required only when Analysis.findings has entries
- `methods` prose required only when Analysis.decisions or Analysis.analyses has entries  
- `inputs` prose required only when Analysis.inputs has entries
- `outputs` prose required only when Analysis.outputs has entries
- `summary` is always optional (no structured counterpart)

If a key's structural list is empty, the narrative key should not appear. 
Write what exists; don't invent.
```

---

### 1.4 Decision `when` / `incompatible_with` / `requires` constraints (UNMENTIONED)

**Status:** Incomplete guidance.

**Finding:** The schema (analysis.yaml:314–346) defines Option attributes:
- `incompatible_with`: multivalued list of "Decision.option pairs that cannot be selected together"
- `requires`: multivalued list of "Decision.option pairs that must also be selected"

And Decision:
- `when`: conditions for when a decision is active

**In prompts:**
- `specify.py:42` mentions `when`, `incompatible_with`, and `requires` exist: "Use when, incompatible_with, and requires constraints for non-independent decisions."
- No guidance on **format** or **when to use each**.

**Schema detail:** 
- `when` format (analysis.yaml:45–49): 'decision_id.option_id' or '~decision_id.option_id', ANDed together
- `incompatible_with` / `requires` format: not explicitly documented in schema; implied to be same format

**Consequence:** Agents may write constraints in wrong format or omit them when the paper hints at dependencies.

**Recommend:** Expand `specify.py:42` with an example:
```
Use when, incompatible_with, and requires constraints for non-independent decisions.

- when: ["parent_decision.selected_option"] marks a decision active only under 
  a specific parent choice
- incompatible_with: ["other_decision.option_that_conflicts"] on an option 
  marks choices that cannot be selected together
- requires: ["other_decision.required_option"] on an option marks a forced pairing

Example: if "use_prior" requires "set_priors_from: literature", add 
requires: ["set_priors_from.literature"] to the use_prior option.
```

---

### 1.5 Anchor grammar (UNMENTIONED, CRITICAL FOR NARRATIVE SKILL)

**Status:** Major gap — prompts do not prepare for narrative authoring.

**Finding:** The narrative skill (`SKILL.md:207–227`) defines tree-path-first anchor syntax, including sub-analysis references like `#reconstruction.decisions.algorithm` (sub-analysis first, then category), NOT `#decisions.reconstruction.algorithm`.

**In prompts:**
- **No prompt mentions anchor syntax.**
- `specify.py` instructs building the ASTRA spec but does not prepare the agent for downstream narrative authoring.
- The SPECIFY prompt delegates to `.claude/guides/decision-guide.md` (project-specific, not auditable here).

**Consequence:** When the narrative skill tries to author anchors, it must re-teach the grammar. Errors cascade when agents structure sub-analyses without naming conventions that support tree-path anchoring.

**Recommend:** Add to `specify.py` before the sub-analysis discussion:
```
## Sub-analysis naming and anchor structure

Sub-analyses use tree-path-first anchor grammar:
- Reference a sub-analysis: #analyses.<sub_id>
- Reference an element inside a sub-analysis: #<sub_id>.<category>.<id>
  Examples: #reconstruction.decisions.algorithm, #preprocessing.outputs.features
  NOT #decisions.reconstruction.algorithm

Sub-analysis IDs must be lowercase_with_underscores and distinct from 
reserved names (inputs, outputs, decisions, findings, prior_insights, 
analyses, options, content, narrative).
```

---

## 2. Narrative Skill Alignment

### 2.1 SPECIFY prompt delegates narrative, creating visibility gap

**Status:** Confirmed gap — narrative discipline is invisible to SPECIFY.

**Finding:** 
- `specify.py:14–15` says: "Read CLAUDE.md for the full ASTRA spec format, decision patterns, and conventions. Read `.claude/guides/decision-guide.md` before pruning decisions."
- No mention of the narrative skill, paper-reproduction mode, or the five-key structure.
- The SPECIFY prompt treats narrative as post-spec (line 33: "Do NOT add recipes — that is the implement phase's job"), but does not clarify whether SPECIFY should draft narrative or leave it for a separate step.

**Paper-reproduction mode is entirely unknown to SPECIFY:**
- The narrative skill reference (`paper-reproduction.md`) teaches fidelity to source, voice seams, scope-limited reproductions, and harvest-don't-invent discipline.
- None of this is mentioned in SPECIFY.

**Consequence:** When SPECIFY runs on a paper reproduction case, it builds the spec structure but:
1. Does not prepare narrative anchors (uses non-canonical IDs, misses sub-analysis naming)
2. Does not mark where reproducer content diverges from the paper (no voice seams)
3. Does not scope the reproduction explicitly
4. Leaves narrative authoring to a later phase without anchoring the spec for it

**Recommend:** Add to `specify.py` system prompt (before Outputs section):

```
## Narrative preparation (for paper reproductions)

If this is a paper reproduction (published paper + ASTRA spec):
- Use sub-analysis IDs that are noun phrases matching the paper's section 
  structure where possible (e.g., `preprocessing`, `calibration`, `fitting`).
- Build decision IDs that are human-readable concepts, not acronyms 
  (e.g., `fiducial_cosmology` not `fid_cosmo`).
- Use consistent snake_case throughout.

The narrative skill will later author the five-key structure (summary, 
methods, findings, inputs, outputs). Structure the spec so every decision, 
output, and finding is independently anchored and distinct from reserved 
names. Narrative will cite these by anchor; missing or ambiguous IDs break 
rendering and validation.

Defer narrative writing itself to the narrative skill, but ensure the spec 
is ready for it: complete, validated, and anchored.
```

---

### 2.2 Five-key narrative structure unknown to SPECIFY

**Status:** Gap — SPECIFY does not prime the spec for narrative authoring.

**Finding:** The narrative skill (`SKILL.md:190–200`) maps the five keys to phases:
- `summary` — question, scope, headline
- `inputs` — provenance
- `methods` — movement of learning, decisions, sub-analyses
- `findings` — implications, findings synthesis
- `outputs` — what was promoted and why

And (`paper-reproduction.md:105–127`) specifies draft order: `inputs` → `methods` → `findings` → `outputs` → `summary`.

**In SPECIFY:**
- No mention of these five keys or their semantics.
- The spec is built as an object graph (inputs, outputs, decisions, findings, prior_insights, analyses) without priming how it will be narrated.

**Consequence:** The spec structure and the narrative structure diverge. Example: SPECIFY may create a decision that looks complete to it but is not ready for narrative authoring — missing label, lacking rationale source citations, unclear why it was included.

**Recommend:** Add to `specify.py` as a validation checklist (before the Rules section):

```
## Narrative readiness checklist (before building narrative)

Before passing the spec to narrative authoring:

- [ ] Every Decision has a label (human-readable, not an acronym)
- [ ] Every Option has a label  
- [ ] Decisions are grouped logically — if >10 at a single level, consider 
      sub-analyses
- [ ] Sub-analyses are named as noun phrases (e.g., preprocessing, fitting)
- [ ] Outputs each have a type and description
- [ ] Inputs each have a source (URL, API, package name)
- [ ] No entity IDs collide with reserved names: inputs, outputs, decisions, 
      findings, prior_insights, analyses, options, content, narrative

These are not optional — they block downstream narrative authoring.
```

---

### 2.3 Paper-reproduction mode constraints not mentioned

**Status:** Gap — SPECIFY does not prepare for reproduction fidelity.

**Finding:** The narrative skill (`paper-reproduction.md:141–170`) specifies:
- **Fidelity to source confidence:** Don't sharpen or soften claims
- **Harvest, don't invent:** Paraphrase the paper, don't originate
- **Voice seams:** Mark when reproducer-specific content enters
- **Published = done:** Present tense matching the paper
- **Scope-limited reproductions:** Name what's in and out

**In SPECIFY:**
- No guidance on fidelity, voice seams, or scope handling.
- Line 53 says "work primarily from `work/notes/` — the extract phase has already distilled the paper" but does not explain the extraction was a **distillation**, not a simplification (fidelity matters).

**Consequence:** SPECIFY may unknowingly alter the paper's claims when deciding which options to include, or may fail to note when reproducer code diverges from the published method.

**Recommend:** Add to `specify.py` when building decisions (around line 35):

```
## Paper reproductions: fidelity discipline

If the source is a published paper:
- Extract decision rationales directly from the paper or from prior_insights.
  Don't invent justifications.
- When code diverges from the paper's stated method, document both options 
  in the decision; mark the codebase choice in universes/baseline.yaml and 
  flag in implementation-notes.md.
- Don't sharpen claims ("we detect" vs "we strongly detect"); preserve the 
  paper's confidence level in all prose and decision descriptions.
- If you add an option the paper doesn't mention, mark it as a reproducer 
  extension (can use excluded: false + rationale mentioning the extension).

The narrative phase will later check fidelity; prepare rationales and options 
that can survive that audit.
```

---

### 2.4 Reserved entity names propagate to narrative anchors

**Status:** Critical cross-tool issue.

**Finding:** The reserved names constraint (`SKILL.md:241–248`) exists to prevent entity IDs from shadowing narrative anchors. An entity named `findings` would break `#findings.xyz` references.

**Impact on SPECIFY:**
- SPECIFY has no way to know about this constraint (not mentioned in prompts).
- An agent could create `decisions: { findings: { ... } }`, which passes schema validation but breaks rendering when the narrative skill tries to use `#findings.xyz` anchors.

**Recommend:** This is a spec-level validation issue. Until the schema enforces the reserved-names pattern at the ID level, SPECIFY must enforce it:

Add to `specify.py` validation section:
```
# Check for reserved names
reserved = {'inputs', 'outputs', 'decisions', 'findings', 'prior_insights', 
            'analyses', 'options', 'content', 'narrative'}
for entity_id in all_entity_ids():
  if entity_id in reserved:
    FAIL: f"Entity ID '{entity_id}' is reserved (used in narrative anchors). Rename it."
```

---

## 3. Stale References

### 3.1 CLI command: `astra paper verify-quotes` (OUTDATED OR UNAVAILABLE)

**Status:** Unverified — likely stale.

**Finding:** The literature phase (`literature.py:49`) instructs:
```bash
echo '<json>' | astra paper verify-quotes "<DOI>"
```

**Current status:**
- This command does not appear in the current Paper2ASTRA or ASTRA CLI public interfaces (checked astra-spec README, Paper2ASTRA CLI entry points).
- The ASTRA CLI likely changed; this may be a draft command that was never implemented or has been renamed.

**Consequence:** When the literature phase runs, the batch verification step will fail. Agents will retry manual verification without understanding why the command disappeared.

**Recommend:** 
1. Verify whether `astra paper verify-quotes` exists in the current ASTRA CLI.
2. If it exists, document it in ASTRA's README.
3. If it doesn't, replace with a real command (possibly `astra validate --verify-evidence` or a standalone Python script).
4. Update `literature.py:49` with the correct command and any schema changes to the verification response format.

---

### 3.2 File paths assume `work/` directory structure

**Status:** Confirmed — correct but assumes specific workspace layout.

**Finding:** All five prompts reference paths like:
- `work/notes/methodology.md`, `work/reference/document.md`, `work/notes/literature.yaml`, etc.
- `targets/` and `targets/targets.md` (extract_targets.py:47, 73)
- `universes/baseline.yaml` (specify.py:30)

**Current implementation:** Paper2ASTRA enforces this layout (workdir.py defines these paths).

**Status:** No drift. Paths are consistent and enforced by the framework.

---

### 3.3 ASTRA validation command syntax

**Status:** Confirmed — correct.

**Finding:** Multiple prompts reference:
- `astra validate astra.yaml` (specify.py:76, review.py:49)
- `astra validate astra.yaml --verify-evidence` (review.py:42)

**Current status:** These are valid ASTRA CLI commands (astra-spec schema v0.0.5 includes validation tooling).

**Status:** No drift detected.

---

## 4. Coverage Gaps for Rich Input Bundles

### 4.1 Code + paper + fibers present

**Status:** Gaps in guidance when all three signals are available.

**Finding:** Our case (pure_eb) has:
- Paper (Planck 2018, Cosmic-Microwave Background analysis)
- Reference code (Snakemake DAG known, dependencies declared)
- Extensive fibers (simulations with known cosmology, baseline parameters documented)

**Prompts assume partial information:**
- `summarize.py:118–137` includes a CODE_SYSTEM_PROMPT for code exploration "if reference code exists" (line 148).
- But there is no guidance on **reconciling** code + paper + fibers when all three disagree.

**Specific gaps:**

1. **Code vs. paper discrepancy:** `specify.py:52–55` says "If notes mention paper-vs-code discrepancies: explore both choices in the multiverse" but does not explain how to weight them or mark which is canonical.

2. **Fibers / simulation data:** Neither the SUMMARIZE nor SPECIFY prompts mention how to handle simulation inputs, baseline parameters, or fiber selection. `summarize.py:42–50` lists "Data sources" but treats them all equally — no distinction between real observational data and simulation outputs.

3. **Parameter provenance:** When code has a default (e.g., `sigma8 = 0.807`) that differs from the paper (e.g., "we use Planck 2018 cosmology"), which is the decision point? The prompts don't clarify.

**Consequence:** SUMMARIZE will extract the paper's stated method, SPECIFY will build the spec, but SPECIFY won't know whether to use paper values or code values as the baseline, and will miss the opportunity to declare fibers as a distinct input class.

**Recommend:**

Add to `summarize.py:CODE_SYSTEM_PROMPT` (around line 118):

```
## When code exists

Also document:
- **Parameter defaults** — every numeric value, threshold, or string literal 
  defined in code (logging, configuration, magic numbers). These are often 
  unstated in the paper and become decision points.
- **Baseline choices** — if code has a default cosmology, template, or set of 
  fibers, record them explicitly.
- **Code-paper reconciliation** — when code diverges from the paper's stated 
  method, note both. Flag in the summary for the specify phase to address.

And if the repository contains simulation data (mocks, fibers, synthetic 
catalogs):
- **Simulation metadata** — cosmology used, number of realizations, random seed.
- **Fiber selection** — which fibers are used, which are held out, how they 
  relate to observations.
```

Add to `specify.py` (before "Prior Insights from Literature", around line 45):

```
## When input bundles include code + fibers + paper

If notes mention simulation data (fibers, mocks, baseline cosmology):

1. **Baseline cosmology:** Create a decision for which cosmology is used. 
   If code has a default, use code value in universes/baseline.yaml even if 
   the paper states a different one. Rationale: the baseline must match what 
   the code actually runs.

2. **Fibers / simulation inputs:** Treat fibers and baseline simulations as 
   explicit inputs (type: data or analysis). Document:
   - Which fibers are used in the primary analysis
   - Baseline cosmology they assume
   - How they relate to the observational data
   - Selection criteria (if some fibers are excluded)

3. **Paper vs. code discrepancy:** If code uses different parameters or 
   methods than the paper states:
   - Create decision options for both paper-stated and code-implemented choices
   - Mark code choice as default in universes/baseline.yaml
   - Explain the discrepancy in implementation-notes.md
   - Do not guess which is more "correct" — document both and let the 
     multiverse explore them

This applies especially to cosmological parameters, prior specifications, 
and selection criteria.
```

---

### 4.2 Snakemake DAG awareness

**Status:** Not mentioned, but could improve implementation guidance.

**Finding:** Our code is organized as a Snakemake DAG with explicit rule dependencies. Neither the SUMMARIZE nor REVIEW prompts mention Snakemake or DAG structure.

**Current prompts:**
- `summarize.py:123–135` describes code exploration generically (architecture, execution flow, outputs).
- `review.py` checks data obtainability but doesn't verify that data download rules are present in the DAG.

**Opportunity:** Snakemake DAGs make dependencies explicit and traceable. If the code uses Snakemake, SUMMARIZE could extract the DAG structure and use it to validate that every input declared in `astra.yaml` has a rule that obtains it.

**Recommend:** This is a nice-to-have, not a must-fix. If adding Snakemake support:

Add to `summarize.py:CODE_SYSTEM_PROMPT` (optional):

```
## If the code uses Snakemake

Identify the main Snakefile and list:
- **Input rules** — rules that download or generate data (e.g., download_fibers, 
  generate_baseline_mocks)
- **Processing rules** — rules that transform data
- **Output rules** — rules that produce the final outputs
- **Dependencies** — how data flows from inputs to outputs (rule order matters 
  for reproducibility)

Extract the DAG's input-to-output chain. The implement phase will need these 
rules to be executable (or to know which are missing).
```

---

### 4.3 Baseline parameters and their discovery

**Status:** Underspecified — hard to extract from code without explicit guidance.

**Finding:** Our case has baseline cosmology (Planck 2018) and baseline parameters (σ₈, Ωₘ, etc.). These may be:
- Hardcoded in code
- Read from a config file
- Derived from a prior
- Stated in the paper

**Current prompts:** `summarize.py:125–135` (code exploration) asks for "key variables and parameters" but doesn't distinguish baseline parameters (decision points) from derived quantities.

**Consequence:** SUMMARIZE produces "sigma8 = 0.807" without context: is this a decision, a derived value, a prior, or a hard-coded default? SPECIFY won't know whether to include it in the multiverse or lock it.

**Recommend:** Add to `summarize.py:CODE_SYSTEM_PROMPT`:

```
## Key variables vs. baseline parameters

Distinguish:
- **Baseline parameters** — cosmological parameters, priors, or selection 
  thresholds that are fixed for the baseline run but could be varied in 
  extensions. These become decision points.
- **Derived quantities** — values computed from baseline parameters 
  (e.g., derived from σ₈ and Ωₘ).
- **Hard-coded constants** — values that are implementation details, not 
  scientific choices (e.g., array allocation size, precision thresholds 
  for numerical stability).

List baseline parameters explicitly. For each, note:
- Where it's defined (code line, config file, paper section)
- Current value
- Whether alternatives are mentioned or explored in the code/paper
```

---

## 5. Recommended Edits

### Priority: MUST-FIX FOR PAPER-II

#### 5.1 Reserve names validation in specify.py

**File:** `/home/cdaley/code/LightconeResearch/Paper2ASTRA/src/paper2astra/prompts/specify.py`

**Before (lines 10–44):**
```python
SPECIFY_SYSTEM_PROMPT = """\
You are an ASTRA specification agent. Read the paper and accumulated notes, then \
produce the structured ASTRA spec and implementation notes.

Read CLAUDE.md for the full ASTRA spec format, decision patterns, and conventions. \
Read `.claude/guides/decision-guide.md` before pruning decisions.
...
```

**After:** Insert after "Read `.claude/guides/decision-guide.md` before pruning decisions." and before "## Inputs":

```python
## Reserved Names

The following names cannot be used as entity IDs (they shadow narrative anchors):
inputs, outputs, decisions, findings, prior_insights, analyses, options, 
content, narrative. If you encounter an entity using one of these names, 
rename it immediately. The spec will serialize but rendering and validation 
will fail silently if IDs collide with reserved names.
```

---

#### 5.2 Conditional narrative requirement in specify.py

**File:** `/home/cdaley/code/LightconeResearch/Paper2ASTRA/src/paper2astra/prompts/specify.py`

**Before (lines 35–44):** Decisions section ends with decision count guidance.

**After:** Insert before "## Prior Insights from Literature":

```python
## Narrative readiness (if this is a paper reproduction)

If you are reproducing a published paper:
- Use sub-analysis IDs matching the paper's section structure where possible.
- Ensure every Decision, Option, Input, and Output has a clear label 
  (not an acronym).
- Use consistent lowercase_with_underscores for all IDs.
- Do NOT write narrative sections (summary, methods, findings, inputs, outputs) 
  yourself — defer to the narrative skill. But structure the spec so every 
  decision, output, and finding is independently anchored.

Narrative will later author the five-key structure. It requires:
- `findings` prose only if Analysis.findings has entries
- `methods` prose only if Analysis.decisions or Analysis.analyses has entries
- `inputs` prose only if Analysis.inputs has entries
- `outputs` prose only if Analysis.outputs has entries
- `summary` is always optional

Do NOT add narrative prose yourself. Ensure the spec is structurally complete 
and ready for it.
```

---

#### 5.3 Paper-vs-code reconciliation in specify.py

**File:** `/home/cdaley/code/LightconeResearch/Paper2ASTRA/src/paper2astra/prompts/specify.py`

**Before (lines 52–55):**
```python
## Paper vs Code Discrepancies

If notes mention paper-vs-code discrepancies: explore both choices in the multiverse \
and create an insight explaining the discrepancy.
```

**After:** Expand to:

```python
## Paper vs Code Discrepancies

If notes mention paper-vs-code discrepancies (different methods, parameters, or 
data selections):

1. **Create two options for the decision:** one matching the paper, one matching 
   the code. Both are valid; the multiverse explores both.
2. **Mark the code choice as default in universes/baseline.yaml.** The baseline 
   must reflect what the code actually runs, not what the paper claims.
3. **Document the discrepancy in implementation-notes.md** so the review and 
   implement phases understand why both exist.
4. **Create a prior_insight** explaining the discrepancy if the paper or code 
   provides context (e.g., "later versions of the method changed the default").

Example: if the paper says "we use Planck 2018 cosmology" but code hardcodes 
"Planck 2015", create a cosmology_version decision with options planck_2015 
and planck_2018, set planck_2015 as default, and note the discrepancy in 
implementation-notes.

This applies to all parameters, methods, and selection criteria where paper 
and code diverge.
```

---

#### 5.4 Rich input bundles guidance in specify.py

**File:** `/home/cdaley/code/LightconeResearch/Paper2ASTRA/src/paper2astra/prompts/specify.py`

**Before (lines 45–51):** "## Prior Insights from Literature"

**After:** Insert new section before "## Prior Insights from Literature":

```python
## Rich Input Bundles (paper + code + fibers/simulations)

If the input bundle includes observational data, reference code, and 
simulation inputs (fibers, mocks, baseline cosmology):

1. **Baseline cosmology decision:** Create an explicit decision for cosmological 
   parameters (e.g., fiducial_cosmology with options planck_2018, wmap9). Use 
   code value as default.

2. **Fiber/simulation inputs:** Declare fibers and baseline simulations as 
   explicit inputs (type: data or analysis). Document:
   - Which fibers are used
   - Baseline cosmology and parameters they assume
   - How they relate to observational inputs
   - Selection criteria (if some are excluded)

3. **Parameter overrides:** When code has default parameters (σ₈, Ωₘ, etc.) 
   that differ from paper-stated values:
   - Extract from code (implementation-notes will flag the discrepancy)
   - Create decision options if the paper mentions alternatives
   - Use code value in baseline universe

This is especially important for cosmological analyses where baseline 
parameters cascade through the entire pipeline.
```

---

### Priority: NICE-TO-HAVE

#### 5.5 Code exploration: baseline parameters (literature.py → summarize.py)

**File:** `/home/cdaley/code/LightconeResearch/Paper2ASTRA/src/paper2astra/prompts/summarize.py`

**Before (lines 118–137):** CODE_SYSTEM_PROMPT

**After:** Add sub-section before "What to produce":

```python
## Baseline parameters and simulation data

If the code defines baseline parameters (cosmological, statistical, or 
selection thresholds) or uses simulation inputs (fibers, mocks):

- **Baseline parameters:** Every default value that could be varied is a 
  decision point. Record: definition site (code line), current value, whether 
  alternatives are explored or mentioned in the paper.
- **Simulation inputs:** If the code uses fibers or baseline simulations, 
  document their cosmology, selection criteria, and relationship to 
  observations.

These feed into decision points in the specify phase; don't omit them as 
"just implementation details."
```

---

#### 5.6 Stale CLI command: verify-quotes (literature.py)

**File:** `/home/cdaley/code/LightconeResearch/Paper2ASTRA/src/paper2astra/prompts/literature.py`

**Before (lines 47–62):** Batch Verification Loop

**Status:** BLOCKED until ASTRA CLI tooling is clarified. This requires coordination with astra-spec maintainers.

**Action:** Verify whether `astra paper verify-quotes` is implemented in the current ASTRA CLI. If not:
1. Remove the batch verification loop instruction
2. Provide alternative guidance (manual verification, Python script, or wait for ASTRA CLI update)
3. Update the schema output format to reflect whatever verification method is used

---

## 6. Summary of Findings

| Category | Issue | Severity | Status |
|---|---|---|---|
| **Schema drift** | Reserved entity names not enforced at schema level | High | Prompt must check; schema should add pattern constraint |
| **Schema drift** | Conditional narrative requirement not explicit in schema | Medium | Prompts should reference; schema documentation clear |
| **Schema drift** | Evidence quote prefix/suffix optional in schema but required in prompts | Low | Prompts correct; schema should enforce |
| **Narrative skill** | SPECIFY unaware of narrative skill, five-key structure, paper-reproduction mode | Critical | Add sections to SPECIFY system prompt |
| **Narrative skill** | Anchor grammar unknown to spec builders | High | Add sub-analysis naming guidance to SPECIFY |
| **Narrative skill** | Paper-reproduction fidelity discipline not mentioned | Critical | Add fidelity constraints to SPECIFY |
| **Stale references** | CLI command `astra paper verify-quotes` unverified | High | Verify in ASTRA CLI; update or remove |
| **Coverage gaps** | No guidance on reconciling paper + code + fibers when all diverge | High | Add sections to SUMMARIZE and SPECIFY |
| **Coverage gaps** | Baseline parameters / cosmology not explicitly extracted | High | Add guidance to CODE_SYSTEM_PROMPT |
| **Coverage gaps** | Snakemake DAG awareness absent | Low | Nice-to-have; not blocking |

---

## 7. Audit Metadata

- **Prompt files audited:** 5 (summarize.py, extract_targets.py, literature.py, specify.py, review.py)
- **Schema files audited:** 3 (analysis.yaml, insight.yaml, universe.yaml)
- **Narrative skill references:** 2 (SKILL.md, paper-reproduction.md)
- **Gaps identified:** 12 (3 critical, 6 high, 2 medium, 1 low)
- **Must-fix edits:** 4
- **Nice-to-have edits:** 2
- **Blocked issues:** 1 (awaiting ASTRA CLI verification)

---

**End of audit.**
