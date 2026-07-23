# CLAUDE.md

## Project: UNIONS Paper II — pure E/B-mode validation (ASTRA reproduction)

ASTRA reproduction of the UNIONS B-modes paper (Daley et al. 2026, "Paper II"
of the UNIONS release series). The paper establishes that the fiducial
v1.4.6.3 leak-corrected catalog passes the configuration-space pure
E/B-mode null tests and characterises the residual COSEBIs / harmonic-Cℓ
behaviour. Cosmological inference is out of scope here — that's Papers III
(config-space) and IV (harmonic-space). Paper I (catalog construction) and
Paper V (simulations) are both Hervas-Peters et al.

This is the FIRST UNIONS reproduction in `LightconeResearch/Reproductions`,
and the first ASTRA spec where the analysis itself was AI-mediated. The
target audience is anyone (human or agent) who wants to interrogate the
paper's choices, alternatives, and quantitative claims.

### Scope

In: catalog preparation → 2pt measurement → pure E/B + COSEBIs + Cℓ_BB
decompositions → semi-analytical covariance + GLASS mock validation → PTE
matrices and combined null tests.

Out: cosmological-parameter inference (Papers III / IV). Their posteriors
are not findings of Paper II and do not appear in `findings:`. Paper II's
findings are the B-mode null-test PTEs and the version-comparison
conclusions.

### Sub-analysis decomposition (mirrors the snakemake DAG, not the paper section order)

```
unions_bmodes (root)
├── catalog              v1.4.6.3 leak-corr load, masks, n(z)
├── two_point            TreeCorr ξ± measurement (1k bins; 10k for COSEBIs)
├── pure_eb              pure E/B decomposition of ξ±  (config space)
├── cosebis              COSEBIs decomposition of ξ±   (config space)
├── cl_bb                pseudo-Cℓ BB measurement      (harmonic space)
├── covariance           CosmoCov + MC propagation, GLASS mock validation
└── null_tests           PTE matrices + combined statistics
```

Sub-analyses inherit shared decisions (fiducial cosmology, fiducial catalog
version, scale cuts) from the root via `from: ../<decision_id>`.

### Source of Truth

| Artifact | Location | Role |
|----------|----------|------|
| Paper TeX | `paper/unions_bmodes.tex` (symlink → `/automnt/n17data/cdaley/unions/pure_eb/docs/unions_release/unions_bmodes/`) | Authoritative for numerical claims, scope, prose voice |
| Paper macros | `paper/claims_macros.tex` | Auto-generated quantitative macros; the cleanest "what number does the paper report" |
| Source code | `/automnt/n17data/cdaley/unions/pure_eb/` | The pipeline that produced the paper |
| Snakemake DAG | `/automnt/n17data/cdaley/unions/pure_eb/workflow/` | What each output IS (rules, configs); rule names map to spec recipes 1-to-1 |
| Project context | `/automnt/n17data/cdaley/unions/pure_eb/CLAUDE.md` | Container/HPC details, code style, quick-reference commands |
| Memory (dense) | `/home/cdaley/.claude/projects/-automnt-n17data-cdaley-unions-pure-eb/memory/MEMORY.md` | Active work state, decisions made, gotchas |
| Fiber tree | `~/loom/.felt/` (project surface in `/automnt/n17data/cdaley/unions/pure_eb/.felt/`) | Lab-notebook record; primary source for prior_insights and decision rationale |
| Literature PDFs | `/automnt/n17data/cdaley/unions/pure_eb/docs/arxiv/<id>/` | Full sources for citation extraction |

### Input bundle (Paper2ASTRA-style staging)

| File | Purpose |
|------|---------|
| `work/notes/methodology.md` | Sub-analysis-by-sub-analysis methodology, decisions, recipes, dependencies |
| `work/notes/literature.yaml` | Prior insights with DOI + exact quote evidence |
| `targets/targets.md` | Replication targets: headline claims, figures, tables, quantities |
| `paper/` | Symlink to the paper source (read-only) |
| `implementation-notes.md` | Tricky-bits notes (placeholder; recipes are not executed at this stage) |

### Local skill sandbox

`.claude/skills/narrative/` is a project-local copy of the lightcone-ui
narrative skill (`add/narrative-skill` branch as of 2026-04-25). Edits made
during the ralph loop's reflective tail land here first; promotion upstream
is a manual cherry-pick once changes stabilise.

### Stage of execution

This project is at the **SPECIFY** stage (in Paper2ASTRA's nine-phase
lexicon). We produce `astra.yaml`, `universes/baseline.yaml`, the seven
sub-analysis specs, and their narratives. We do NOT run anything yet —
recipes are included on outputs as documentation of what the snakemake
rules do, but `prism run` / `astra validate --execute` is not invoked.

The IMPLEMENT phase (re-running the analysis end-to-end inside Prism) is a
follow-up effort once the spec stabilises.

### Stylistic conventions

- **Fidelity to paper voice.** Don't sharpen or soften claims relative to
  the paper TeX. Preserve hedges. The paper is the authoritative voice for
  scope, certainty, and numerical claims.
- **Use the macros for numbers.** Anything that has a `\<name>*` macro in
  `paper/claims_macros.tex` should appear in the spec with the same value;
  do not re-derive.
- **Snakemake rule names map 1-to-1 to recipe commands.** When in doubt
  about what produces an output, find the rule by output path in
  `workflow/rules/*.smk`.
- **Decisions are the spine.** Any choice where a different defensible
  value would shift a numerical result becomes a decision. Aim for 8–20
  root-level decisions plus per-sub-analysis decisions where they are
  locally scoped (mask version → catalog; binning → two_point; etc.).
- **Excluded options are first-class.** When a paper considers and rejects
  an option, declare it with `excluded: true` and `excluded_reason:`.
  Better than silently absent.

### Reading order for any agent starting fresh

1. This file.
2. The constitution at `felt show paper-ii-astra-conversion/constitution`.
3. The parent fiber `felt show paper-ii-astra-conversion` (current state and
   open work).
4. `work/notes/methodology.md`.
5. `targets/targets.md`.
6. `astra.yaml` (current state) and any `analyses/<sub>/astra.yaml`.
7. `paper/unions_bmodes.tex` — selected sections relevant to the current
   iteration.
8. The relevant fibers (search the loom by topic).
9. The narrative skill at `.claude/skills/narrative/SKILL.md` plus
   `references/paper-reproduction.md`.

### Spec-narrative invariant

Every declared decision, finding, output, and sub-analysis must be cited in
the narrative tree. Run `astra validate astra.yaml` after any spec change.
