# START HERE

**If you are an AI agent (Claude, GPT, Gemini, Cursor, Continue, autonomous coding agent, anything): read this file fully before doing anything else in this repository. It is a checklist, not a narrative. Work through it in order.**

If you are a human: `README.md` is written for you. This file is for AI assistants, but you're welcome to read it.

---

## Why this file exists

AI agents have a recurring failure mode on this project: responding to APF requests with fluent-but-wrong output because the relevant context was never loaded. This file is the minimum context load to reason about Paper 8 of the APF series honestly.

This file is structured in three tiers. **Tier 1** (§0 below) is the 60-second mental map — read it even if you plan to skip everything else. **Tier 2** (§1 corpus table + §2 the seven steps) is 5–10 minutes of operational loading; it is what the "ready to work" state actually requires. **Tier 3** (§3 + `REVIEWERS_GUIDE.md` + the manuscript itself) is full audit depth, only when the task demands it.

---

## §0 — 60-second mental map (Tier 1)

**What APF is.** APF is a single-axiom framework. The only axiom is that enforcement (the cost of maintaining physical distinction) is finite (Axiom A1). From A1 and its four-component structure (PLEC), 342 theorems derive the Standard Model gauge group SU(3)xSU(2)xU(1), 45 fermions, three generations, and 48 quantitative predictions with zero free parameters. Each paper contributes one layer: Paper 0 states the ontology; Papers 1-7 build the derivation chain; Paper 13 is the master reference. The codebase (v6.9) is the authoritative specification; the papers are prose guides to reading the code.

**What this paper contributes.** The Admissibility-Capacity Ledger: a unifying record (K, d_eff) read through six regime projections (pi_T, pi_G, pi_Q, pi_F, pi_C, pi_A) with four consistency identities I_1-I_4. Paper 8's only non-trivial structural claim is the gauge-cosmological bridge I_2 (Theorem 1.1 of the Technical Supplement Formal Kernel): a unique G_SM-invariant 42-dimensional subspace V_Lambda of V_61 exists under the Standard-Model gauge action and T12 partition constraint, inducing the residual partition 3+16+42=61 and cosmological fraction Omega_Lambda = 42/61. I_1 is definitional under the cell-count convention K_horizon = A/(4 l_P^2); I_3 and I_4 are standard finite-dimensional closures. Supplement contains self-contained proof via Maschke semisimplicity + minimal working example at toy interface K=3 d_eff=4 auditable in ≤10 lines numpy + Planck 2018 multivariate Bayes-factor restricted sanity check + finite-volume C*-algebra formulation with Type III_1 classification at infinite volume. Four Phase 16-20 reviewer-response passes culminated in v2.9 main + v2.5 supplement.

**What this repo verifies vs what it imports.** This repo locally verifies 36 theorem checks drawn from Paper 8's own dependency subset. Results from Papers 1, 2, 3, 4, 6, and 7 are imported, not re-derived here. The full APF corpus (and which paper owns which result) is in §1 below.

**Scope — after reading only this repo:**

- **You may:** summarize Paper 8's argument and results; explain the local 36-theorem subset; run all bundled checks to verify them; cite the paper by its DOI.
- **You may not:** claim mastery of the full APF corpus on the basis of this repo alone; re-derive results that belong to other papers using only material here; present imported results as locally derived in this repo; treat the bundled check subset as the whole engine (the canonical bank is 342 theorems / 355 checks; this repo contains a subset).

**Canonical state (one source of truth).** Codebase **v6.9**, frozen 2026-04-20. Full engine: **355** verify_all checks, **342** bank-registered theorems, **48** quantitative predictions. This repo bundles 36 of those checks for Paper 8's dependency subset. If you see a different count anywhere, this line is canonical.

---

## §1 — The corpus at a glance (Tier 1)

This repo is **one paper-companion** in a 9-paper series. Before you treat any claim as the whole framework, know what else exists.

| # | Title | Zenodo DOI | GitHub repo | Status |
|---|---|---|---|---|
| 0 | What Physics Permits | [10.5281/zenodo.18605692](https://doi.org/10.5281/zenodo.18605692) | [`APF-Paper-0-What-Physics-Permits`](https://github.com/Ethan-Brooke/APF-Paper-0-What-Physics-Permits) | public |
| 1 | The Enforceability of Distinction | [10.5281/zenodo.18604678](https://doi.org/10.5281/zenodo.18604678) | [`APF-Paper-1-The-Enforceability-of-Distinction`](https://github.com/Ethan-Brooke/APF-Paper-1-The-Enforceability-of-Distinction) | public |
| 2 | The Structure of Admissible Physics | [10.5281/zenodo.18604839](https://doi.org/10.5281/zenodo.18604839) | [`APF-Paper-2-The-Structure-of-Admissible-Physics`](https://github.com/Ethan-Brooke/APF-Paper-2-The-Structure-of-Admissible-Physics) | public |
| 3 | Ledgers | [10.5281/zenodo.18604844](https://doi.org/10.5281/zenodo.18604844) | [`APF-Paper-3-Ledgers-Entropy-Time-Cost`](https://github.com/Ethan-Brooke/APF-Paper-3-Ledgers-Entropy-Time-Cost) | public |
| 4 | Admissibility Constraints and Structural Saturation | [10.5281/zenodo.18604845](https://doi.org/10.5281/zenodo.18604845) | [`APF-Paper-4-Admissibility-Constraints-Field-Content`](https://github.com/Ethan-Brooke/APF-Paper-4-Admissibility-Constraints-Field-Content) | public |
| 5 | Quantum Structure from Finite Enforceability | [10.5281/zenodo.18604861](https://doi.org/10.5281/zenodo.18604861) | [`APF-Paper-5-Quantum-Structure-Hilbert-Born`](https://github.com/Ethan-Brooke/APF-Paper-5-Quantum-Structure-Hilbert-Born) | public |
| 6 | Dynamics and Geometry as Optimal Admissible Reallocation | [10.5281/zenodo.18604874](https://doi.org/10.5281/zenodo.18604874) | [`APF-Paper-6-Dynamics-Geometry-Spacetime-Gravity`](https://github.com/Ethan-Brooke/APF-Paper-6-Dynamics-Geometry-Spacetime-Gravity) | public |
| 7 | Action, Internalization, and the Lagrangian | [10.5281/zenodo.18604875](https://doi.org/10.5281/zenodo.18604875) | [`APF-Paper-7-Action-Internalization-Lagrangian`](https://github.com/Ethan-Brooke/APF-Paper-7-Action-Internalization-Lagrangian) | public |
| 13 | The Minimal Admissibility Core | [10.5281/zenodo.18614663](https://doi.org/10.5281/zenodo.18614663) | [`APF-Paper-13-The-Minimal-Admissibility-Core`](https://github.com/Ethan-Brooke/APF-Paper-13-The-Minimal-Admissibility-Core) | public |
| — | Canonical codebase (v6.9) | [10.5281/zenodo.18604548](https://doi.org/10.5281/zenodo.18604548) | [`APF-Codebase`](https://github.com/Ethan-Brooke/APF-Codebase) | pending |

You are reading the companion to **Paper 8** (this repo: https://github.com/Ethan-Brooke/APF-Paper-8-Admissibility-Capacity-Ledger). The canonical computational engine lives at codebase **v6.9** ([Zenodo](https://doi.org/10.5281/zenodo.18604548)). Each paper has its own companion repo following the same layout.

**When in doubt about scope:** if the question is about a result outside Paper 8, point the user (or yourself) to the relevant paper above. Do not try to derive it from the bundled subset; the canonical state is the codebase, not this slice.

**Machine-readable equivalent.** For programmatic use, [`ai_context/repo_map.json`](ai_context/repo_map.json), [`ai_context/theorems.json`](ai_context/theorems.json), and [`ai_context/derivation_graph.json`](ai_context/derivation_graph.json) are the primary machine interfaces — they encode the same corpus inventory, full theorem catalog, and dependency DAG in structured form. For AI agents, those three JSON files are often more useful than the prose of this file. Root-level copies are also present at the repo root for discoverability.

---

## §2 — Operational loading (Tier 2): the seven steps

### Step 1 — Verify the code runs (30 seconds)

Before anything else, confirm the executable bank works on your system:

```bash
pip install -e .
python run_checks.py
```

Expected output: `Paper 8 (The Admissibility-Capacity Ledger): N passed, 0 failed, N total — verified in under 1s`.

If any check fails, **stop**. Do not proceed with substantive work. Report the failure to the user with the full message; this is either a bug, a platform issue, or a genuine falsification candidate, and each has a different disposition.

### Step 2 — Read `ai_context/AGENTS.md` (2 minutes)

Mandatory. This file establishes:

- What this repo is (Paper 8 + vendored `apf/` subset)
- Trust hierarchy (code > supplement > paper > README > wiki)
- The seven specific failure modes AI agents have actually hit on this project
- The current canonical state (v6.9, 342 theorems, four-parameter Planck match)
- The H0 falsifier (so you don't reflexively try to "fix" it with the Grok-style mechanism already rejected)

### Step 3 — Absorb the framework in 5 minutes (`ai_context/FRAMEWORK_OVERVIEW.md`)

Mandatory. This gives you:

- Axiom A1 in one sentence
- PLEC's four components (A1 + MD + A2 + BW)
- The 61-partition ($C_{\mathrm{total}} = 61 = 42 + 3 + 16$), and why it's topological not dynamical
- The derivation chain (A1 → L_nc → L_irr → L_col → Theorem_R → T_gauge → T_field → T4F → ...)
- The 19 bank modules
- The 48 quantitative predictions and the four-parameter Planck cosmological match
- What APF IS and what APF IS NOT (scope boundary)

### Step 4 — Check vocabulary (`ai_context/GLOSSARY.md`)

Skim, don't memorize. Must know by sight:

- **Axioms:** A1 (finite capacity), MD (minimum distinction cost), A2 (argmin selection), BW (non-degeneracy)
- **PLEC primitives:** admissible class $\mathcal{A}_\Gamma$, enforcement-cost functional $E_\Gamma$, PLEC selection $X^*$
- **Regime R** (R1–R4 conditions for PLEC well-posedness) + the **five-type regime-exit taxonomy** (Types I–V)
- **Epistemic tags:** [P] / [P]+[I] / [P_structural] / [Pred] — these are load-bearing, not decorative
- **L_irr** — locally unrecoverable correlations, **NOT** "infinite-cost reversal" (the old framing is explicitly rejected)

### Step 5 — Load the theorem catalog (`ai_context/theorems.json`)

Not mandatory to read end-to-end, but make sure you know it exists and how to query it. The full bank has 342 entries; this paper bundles **36 of 342** for local execution. Each entry has name, module, epistemic tag, dependencies, and one-line summary.

Querying examples (use when you need to cite a specific theorem):

```bash
# Find all [P]-tagged theorems in core.py
jq '.entries[] | select(.module == "core" and .epistemic == "P") | .name' ai_context/theorems.json

# Find theorems bundled in this repo
jq '.entries[] | select(.bundled_in_this_repo) | .name' ai_context/theorems.json

# Find dependencies of a specific check
jq '.entries[] | select(.name == "T7B") | .dependencies' ai_context/theorems.json
```

The DAG is also available as structured data at `ai_context/derivation_graph.json` (368 nodes, 1516 edges) if you need to traverse dependencies.

### Step 6 — Internalize the engagement discipline (`ai_context/AUDIT_DISCIPLINE.md`)

**Mandatory before proposing any substantive edit or critiquing any external proposal.** This file codifies the audit-response posture the project requires:

- Read carefully before responding; verify claims against the current document (use `grep` / `read`, not your memory)
- Acknowledge gaps before proposing fixes
- Word distinctions that matter (bibliography ≠ literature review; citation ≠ positioning; derivation ≠ reinterpretation; "stated" ≠ "closed")
- Epistemic tag meanings — don't ship a [P_structural] claim as if it were a [P]

**Read the §8 worked example (the Grok H0 proposal) carefully.** It's the specific failure mode this discipline prevents. Understanding it is the fastest way to internalize the posture.

### Step 7 — Know what's genuinely open (`ai_context/OPEN_PROBLEMS.md`)

Do not attempt to "solve" any problem on this list without first checking the verdicts. Each entry documents why it's open, which escape routes are blocked by [P] theorems, and what would be required to close it honestly. The ten current open problems:

1. **H0 tension** (7.09σ vs H0DN 2026) — five escape routes catalogued; routes I and II blocked by [P] theorems
2. Absolute mass scales for $m_t$ and $m_b$
3. $m_c$ at 2.6% (Schur structural limit)
4. Dark matter particle identity (geometric sector identification is [P_structural + Pred])
5. PLEC existence/uniqueness theorem
6. Formal regime-exit theorem
7. DESI DR3+ confrontation
8. Neutrino CP phase derivation
9. Quantum-gravity / Planck-scale unification
10. Matter-content bridge tightening

---

## You are now ready

If you've completed steps 1–7, you should be able to:

- **Verify any claim in this paper** by finding its `\coderef{check_X}{module.py}` anchor, running the check, and reading the docstring.
- **Trace any derivation** from axioms to conclusion using `theorems.json` or `derivation_graph.json`.
- **Position any claim** against external literature with appropriate epistemic-tag precision.
- **Critique external proposals** against the framework's known [P]-theorem constraints without reinventing the Grok H0 failure mode.
- **Propose new work** with honest scope, not the shallow version.

## Self-check — answer these before proceeding

Before doing substantive work, confirm you can answer these three questions without re-reading:

1. **What is the 61-partition, and why is it topological rather than dynamical?**
   *(Expected: $61 = 42 + 3 + 16$ split by gauge-addressability and confinement predicates applied to the anomaly-free field content; see `check_L_saturation_partition` in cosmology.py. Type-counting, not capacity-amount dependent.)*

2. **What is APF's H0 prediction, and what is the tension with H0DN 2026?**
   *(Expected: $H_0 = 67.76$ km/s/Mpc from $\Omega_m = 19/61$ and Planck's measured $\Omega_m h^2 = 0.143$; 7.09σ tension with H0DN's $73.50 \pm 0.81$. Named [P] theorems that would have to be demoted to accommodate: L_equation_of_state, L_saturation_partition, L_anomaly_free.)*

3. **What does [P]+[I] mean, and how is it different from [P]?**
   *(Expected: [P] is internally derived from APF axioms; [P]+[I] means APF prepares the admissible class and an external mathematical theorem — Solèr, Gleason, HKM, Malament, Lovelock — closes the classification or uniqueness. The external theorem is an import, not an APF-internal derivation.)*

If you can't answer any of these confidently, re-read the relevant file. Answering them correctly confirms the context loaded successfully.

---

## When asked to do work — the defaults

### If asked to verify something
Run it. Read the check's docstring. Walk its dependencies via `theorems.json`. Cite by check name + module.

### If asked to critique an external proposal
Treat it as an audit by default. Apply the posture from `AUDIT_DISCIPLINE.md`. Specifically check: does the proposal silently require any [P] theorem to fail? Does it introduce free parameters presented as zero? Does it move in a direction consistent with adjacent data?

### If asked to propose new work
Scope first. Name the gap. Don't write code or prose until scope is agreed. When in doubt, propose more scope rather than less.

### If asked to "fix the H0 tension" or similar known-open problems
Don't. Read `OPEN_PROBLEMS.md` § on that problem. The verdicts on escape routes are already catalogued. Only propose a new mechanism if it does not repeat a known failure mode.

### If asked to write LaTeX
The paper uses `\coderef{check_name}{module.py}` for code references. Respect this convention — every substantive claim should either have a coderef or be explicitly tagged [P_structural] or [Pred]. Do not invent check function names — `grep apf/bank.py` first.

### If asked something you don't know
Say so. The project values honest "I don't know" over fluent bullshit. This is codified in the audit-response discipline.

---

## Files in `ai_context/` — full manifest

You don't need to read all of these upfront. This is the "what's in the toolbox" reference.

| File | When to consult |
|---|---|
| `AGENTS.md` | **Read in Step 2** — authoritative entry point |
| `FRAMEWORK_OVERVIEW.md` | **Read in Step 3** — APF in 5 minutes |
| `GLOSSARY.md` | **Skim in Step 4** — vocabulary and tags |
| `AUDIT_DISCIPLINE.md` | **Read in Step 6** — before any critique or proposal |
| `OPEN_PROBLEMS.md` | **Skim in Step 7** — before proposing new work |
| `HOW_TO_VERIFY.md` | When you need to verify a specific claim |
| `STARTING_PROMPTS.md` | When user hasn't given you a prompt pattern; pick one here |
| `CITING.md` | When generating citations for the paper / codebase / framework |
| `PAPER_LINEAGE.md` | When the Zenodo vs local version question comes up |
| `THEOREMS.md` | Human-readable full bank catalog (query `theorems.json` for structured access) |
| `theorems.json` | Machine-readable catalog — query with jq or Python |
| `derivation_graph.json` | The DAG as structured data — traverse dependencies |
| `wiki/INDEX.md` | Entry point to the bundled wiki subset |
| `wiki/concepts/` | Concept pages (Axiom A1, Derivation Chain, Born Rule, etc.) |
| `wiki/papers/` | Canonical per-paper pages (one per spine paper) |
| `wiki/codebase/` | Framework-wide pages (codebase overview, predictions catalog, red-team tests) |

---

## The paper itself

Only after completing steps 1–7:

- **`Paper_8_Admissibility_Capacity_Ledger_v2.9.pdf`** — the manuscript. Read §1 (Introduction) and whichever section is relevant to your task. Don't try to read the whole paper before starting work.
- **`Paper_8_Admissibility_Capacity_Ledger_Supplement_*.pdf`** (if present) — the Technical Supplement. Canonical proofs. Dense; dip in only where you need depth.

The `.tex` sources are also included if you need to propose edits.

---

*This file is generated by the `create-repo` skill from `templates/START_HERE.md.template` in the canonical APF repository. It is updated on each release rebuild. Do not edit it in a specific release — your changes will be overwritten. Edit the template instead.*
