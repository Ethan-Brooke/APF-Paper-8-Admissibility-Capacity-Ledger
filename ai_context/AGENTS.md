# AGENTS.md — AI Onboarding Entry Point

**This file is the authoritative read-first document for any AI agent loading this repository.** Humans can use `README.md` instead; it's optimized for human reviewers.

If you're an AI (Claude, GPT, Gemini, Cursor, Continue, an autonomous coding agent, anything): read this whole file before anything else. It tells you what this repo is, what to trust, and how to avoid the four or five specific failure modes that AI assistants have actually hit on this project.

---

## 0. The corpus — what this repo is part of

The Admissibility Physics Framework is a **9-paper spine** plus a canonical computational engine. Each paper has its own companion repo following this same layout. Before reasoning about any claim as "the framework," know what else exists:

| # | Title | Zenodo DOI | GitHub repo | Status |
|---|---|---|---|---|
| 0 | What Physics Permits | [10.5281/zenodo.18439523](https://doi.org/10.5281/zenodo.18439523) | [`APF-Paper-0-What-Physics-Permits`](https://github.com/Ethan-Brooke/APF-Paper-0-What-Physics-Permits) | public |
| 1 | The Enforceability of Distinction | [10.5281/zenodo.18439200](https://doi.org/10.5281/zenodo.18439200) | [`APF-Paper-1-The-Enforceability-of-Distinction`](https://github.com/Ethan-Brooke/APF-Paper-1-The-Enforceability-of-Distinction) | public |
| 2 | The Structure of Admissible Physics | [10.5281/zenodo.18439274](https://doi.org/10.5281/zenodo.18439274) | [`APF-Paper-2-The-Structure-of-Admissible-Physics`](https://github.com/Ethan-Brooke/APF-Paper-2-The-Structure-of-Admissible-Physics) | public |
| 3 | Ledgers | [10.5281/zenodo.18439363](https://doi.org/10.5281/zenodo.18439363) | [`APF-Paper-3-Ledgers-Entropy-Time-Cost`](https://github.com/Ethan-Brooke/APF-Paper-3-Ledgers-Entropy-Time-Cost) | public |
| 4 | Admissibility Constraints and Structural Saturation | [10.5281/zenodo.18439397](https://doi.org/10.5281/zenodo.18439397) | [`APF-Paper-4-Admissibility-Constraints-Field-Content`](https://github.com/Ethan-Brooke/APF-Paper-4-Admissibility-Constraints-Field-Content) | public |
| 5 | Quantum Structure from Finite Enforceability | [10.5281/zenodo.18439433](https://doi.org/10.5281/zenodo.18439433) | [`APF-Paper-5-Quantum-Structure-Hilbert-Born`](https://github.com/Ethan-Brooke/APF-Paper-5-Quantum-Structure-Hilbert-Born) | public |
| 6 | Dynamics and Geometry as Optimal Admissible Reallocation | [10.5281/zenodo.18439445](https://doi.org/10.5281/zenodo.18439445) | [`APF-Paper-6-Dynamics-Geometry-Spacetime-Gravity`](https://github.com/Ethan-Brooke/APF-Paper-6-Dynamics-Geometry-Spacetime-Gravity) | public |
| 7 | Action, Internalization, and the Lagrangian | [10.5281/zenodo.18439513](https://doi.org/10.5281/zenodo.18439513) | [`APF-Paper-7-Action-Internalization-Lagrangian`](https://github.com/Ethan-Brooke/APF-Paper-7-Action-Internalization-Lagrangian) | public |
| 13 | The Minimal Admissibility Core | [10.5281/zenodo.18361446](https://doi.org/10.5281/zenodo.18361446) | [`APF-Paper-13-The-Minimal-Admissibility-Core`](https://github.com/Ethan-Brooke/APF-Paper-13-The-Minimal-Admissibility-Core) | public |
| — | Canonical codebase (v6.9) | [10.5281/zenodo.18529115](https://doi.org/10.5281/zenodo.18529115) | [`APF-Codebase`](https://github.com/Ethan-Brooke/APF-Codebase) | pending |

**You are reading the companion to Paper 8** (this repo: https://github.com/Ethan-Brooke/APF-Paper-8-Admissibility-Capacity-Ledger).

The **machine-readable equivalent** of this corpus (plus the full file map of this repo) is in [`ai_context/repo_map.json`](repo_map.json). Load that once and you have the whole structure indexed.

When a user's question is about a result outside Paper 8, point them to the relevant row above. Do not attempt to re-derive cross-paper results from the bundled subset — the canonical state is the codebase ([Zenodo concept DOI](https://doi.org/10.5281/zenodo.18529115)), not this slice.

---

## 1. What this repository is

This is the **paper-companion release repo for Paper 8 of the Admissibility Physics Framework (APF) series**. Specifically: The Admissibility-Capacity Ledger.

The repo bundles:

- The manuscript and Technical Supplement (`.tex` + `.pdf`)
- A vendored subset of the APF codebase containing only the `check_*` functions that Paper 8 directly references (`apf/core.py`, `apf/apf_utils.py`, `apf/bank.py`)
- A one-command verifier (`run_checks.py`)
- An interactive derivation DAG (`docs/index.html`, served via GitHub Pages)
- A theorem-by-theorem Colab walkthrough (`APF_Reviewer_Walkthrough.ipynb`)
- **This AI-onboarding pack** (`ai_context/`) — the authoritative context for agents

The framework-wide canonical codebase (342 bank-registered theorems, 19 modules, 355 `verify_all` checks) is **not** bundled here. Only the Paper 8 subset is runnable locally. See `ai_context/THEOREMS.md` for the full bank listing with bundled-in-this-repo flags.

## 2. APF in 30 seconds

APF derives the Standard Model from a single axiom:

> **Axiom A1 (Finite Information Capacity):** Any interface of physical distinction has finite information capacity.

From A1 plus three structural components (MD = minimum distinction cost, A2 = least-cost selection, BW = non-degeneracy), the framework derives:

- The gauge group $SU(3) \times SU(2) \times U(1)$ — uniquely forced by non-closure + chirality + minimality.
- The fermion content (45 chiral fermions; 1 of 4680 anomaly-free combinations survives).
- The three generations ($N_{\mathrm{gen}} = 3$ exactly, from capacity saturation).
- The integer partition $C_{\mathrm{total}} = 61 = 42 + 3 + 16$ (vacuum : baryon : cold-dark).
- Quantum structure (Hilbert space, Born rule, tensor products, unitarity).
- General relativity in $d = 4$ (metric from non-factorization; Lorentzian signature from A4; Einstein equations from Lovelock uniqueness).
- 48 quantitative predictions, 0 free parameters, 32/39 within 3σ of observation.

The framework's cosmological prediction matches Planck 2018 across all four density parameters $(\Omega_b, \Omega_c, \Omega_m, \Omega_\Lambda) = (3/61, 16/61, 19/61, 42/61)$ within $1.2\%$, with $H_0 = 67.76$ km/s/Mpc.

For the 5-minute version, read `ai_context/FRAMEWORK_OVERVIEW.md` next.

## 3. Recommended read order for AI

If you're loading this repo cold and want to understand Paper 8 well enough to reason about it:

1. **This file (`AGENTS.md`)** — what the repo is, what to trust, what not to do.
2. **`ai_context/FRAMEWORK_OVERVIEW.md`** — APF in 5 minutes.
3. **`ai_context/GLOSSARY.md`** — the terms you'll see (admissibility, PLEC, L_irr, regime R, epistemic tags, etc.). Scan, don't memorize.
4. **`README.md`** — paper-specific orientation and the verification recipe. Badges, theorem table, install.
5. **`apf/core.py`** — the vendored code. Read at least one check function fully. This is ground truth.
6. **`run_checks.py`** + run it — confirms the local state works.
7. **`Paper_8_Admissibility_Capacity_Ledger_v2.9.tex`** or `.pdf` — the manuscript. Read §1 (intro) and the section you care about. Don't try to read the whole paper first.
8. **Technical Supplement** if present — canonical proofs. Dense; dip in only where you need depth.
9. **`ai_context/wiki/`** — concept pages (Axiom A1, Derivation Chain, Born Rule, etc.) and per-paper pages. Navigate from `wiki/INDEX.md`.
10. **`ai_context/AUDIT_DISCIPLINE.md`** — read before you propose any substantive change or critique. This is the posture the project expects from AI contributors.

## 4. Trust hierarchy

In order of decreasing trust, most-trusted first:

1. **Vendored `apf/` Python code.** The checks execute; the math is machine-verified. If the code passes, the claim stands.
2. **Technical Supplement (if included).** Canonical proofs. The supplement is the project's source of truth for formal content.
3. **Manuscript main body.** Prose guide to the supplement and code. May compress proofs into sketches.
4. **`ai_context/` files.** Framework-level context, authored for AI readers. Accurate at the time of this release.
5. **`ai_context/wiki/` pages.** Concept pages and per-paper pages. Updated separately from releases; may drift.
6. **`README.md` / `REVIEWERS_GUIDE.md`.** Authored for human reviewers; accurate but voiced for skeptical physicists, not for AI reasoning.

**When sources conflict, trust the code.** If a prose claim in the paper contradicts a passing check function in `apf/core.py`, the check is right and the prose needs correction. Flag this rather than explaining it away.

## 5. The H0 falsifier — read this before touching cosmology

As of April 2026, the H0DN community-consensus measurement ($H_0 = 73.50 \pm 0.81$ km/s/Mpc; Casertano et al., A&A 708 A166) is in **$7.09\sigma$ tension** with APF's prediction of $H_0 = 67.76$ km/s/Mpc.

APF's position is deliberately sharp: the framework is **on the Planck side** of the tension. If the H0DN result holds and is not measurement systematics, at least one of these [P] theorems must be demoted:

- `check_L_equation_of_state` (would have to admit $w \neq -1$ in some regime)
- `check_L_saturation_partition` (would have to admit non-topological partition evolution)
- `check_L_anomaly_free` (would have to admit late-time type decoupling)

**Do not propose "time-locked vacuum" mechanisms or similar phenomenological fixes without first deriving the regime conditions under which one of those [P] theorems fails.** A specific external proposal (Grok, April 2026) tried this by substituting the topological partition denominator with a time-dependent capacity amount, and was rejected because the substitution silently demoted `L_saturation_partition`. That episode is captured in `ai_context/AUDIT_DISCIPLINE.md` as a worked example.

The five escape routes (pre-recombination $r_d$ cut, dynamical DE, distance-ladder systematics, curvature, local inhomogeneity) are catalogued in `OPEN_PROBLEMS.md` with verdicts. Routes I and II are blocked by [P] theorems; route III is substantially ruled out by H0DN itself; routes IV and V are APF-admissible but currently empirically insufficient.

The honest posture: $H_0$ is a real falsification candidate, not a bug to patch around.

## 6. Current canonical state (at time of this release)

- **Codebase:** v6.9 (2026-04-20)
- **Bank-registered theorems:** 342 across 19 modules
- **Total `verify_all` checks:** 355 (includes `apf/standalone/` + session modules)
- **Quantitative predictions:** 48 (39 tested; 32/39 within 3σ; mean error 3.83%; median 0.37%)
- **Free parameters:** 0

## 7. Seven things AI agents have actually gotten wrong on this project — don't repeat these

### 7.1 Don't conflate citation with positioning
Inserting `\cite{Foo}` at a claim site is a citation. Explaining how the paper's claim relates to Foo's result is positioning. When a reviewer asks for a literature review, a bibliography does not satisfy the request. See `AUDIT_DISCIPLINE.md` §2.

### 7.2 Don't conflate derivation with reinterpretation
Deriving $X$ from APF axioms is a derivation. Relabeling a known result $X$ in APF vocabulary is a reinterpretation. Both are valid moves, but they have different epistemic weight. Papers must mark which is which. See `GLOSSARY.md` under "epistemic tags."

### 7.3 Don't smuggle demotions of [P] theorems
If a proposed mechanism implicitly requires a [P] theorem to fail, acknowledge that explicitly and argue for the demotion. Do not present the mechanism as if the [P] theorem still held. The Grok H0 episode (see §5) is this failure mode.

### 7.4 Don't treat "resolved" as a completed status without verifying against current document state
Claude's history on this project includes declaring issues resolved, then being corrected by a follow-up audit reading the same source. Before claiming any issue is fixed, `grep` or `read` the actual document to verify. See `AUDIT_DISCIPLINE.md` §1.

### 7.5 Don't reflexively "fix" the H0 tension
See §5 above. APF's position on H0 is that the tension is a genuine falsification candidate. Phenomenological patches that silently rewrite load-bearing [P] theorems worsen the framework, not improve it.

### 7.6 Don't invent check function names
If you need to reference a check, `grep` `apf/bank.py` or read `apf/core.py` to confirm the name exists. A prior pass on Paper 0 referenced eight module names (`linear_algebra.py`, `quantum.py`, `algebra.py`, etc.) that don't exist in the canonical codebase. The fix was expensive.

### 7.7 Don't treat a second audit as redundant
If reviewer B raises the same concern as reviewer A, the convergent criticism is a stronger signal, not a weaker one. Check whether reviewer B raises new literature references, new terminological distinctions, or new empirical anchors. See `AUDIT_DISCIPLINE.md` §3.

## 8. Engagement discipline — how to work on APF honestly

- **The audit-response posture is default for any external proposal.** Read before responding. Verify claims against the current document. Name gaps before proposing fixes. Full read in `ai_context/AUDIT_DISCIPLINE.md`.
- **Code is truth; prose is guide.** When proposing changes, make sure the change is consistent with the running code. If the change would break a passing check, say so.
- **Scope matters.** Propose scope before executing substantive work. Don't ship a "targeted pass" when a structural rewrite is needed.
- **Epistemic tags are load-bearing.** [P] = internally derived; [P]+[I] = derived with external import closing classification; [P_structural] = structural identification; [Pred] = phenomenological prediction. See `GLOSSARY.md`.

## 9. Files in `ai_context/`

| File | Purpose |
|---|---|
| `AGENTS.md` | This file. Read first. |
| `FRAMEWORK_OVERVIEW.md` | APF in 5 minutes (axiom, PLEC, partition, chain, predictions). |
| `AUDIT_DISCIPLINE.md` | How to engage honestly on APF. Distilled from the audit-response skill. |
| `STARTING_PROMPTS.md` | High-value query templates known to work. |
| `GLOSSARY.md` | APF-specific terminology and epistemic tag meanings. |
| `HOW_TO_VERIFY.md` | Recipes for verifying any claim in this paper. |
| `THEOREMS.md` | Full bank catalog (342 theorems) with bundled-in-this-repo flags. |
| `theorems.json` | Same catalog in structured form — query/grep-friendly. |
| `derivation_graph.json` | The DAG as structured data (nodes + edges). |
| `OPEN_PROBLEMS.md` | Framework-level open questions (so you don't re-solve known-open ones naively). |
| `CITING.md` | How to cite the paper, the codebase, and the framework. |
| `PAPER_LINEAGE.md` | This paper's version history and what supersedes what on Zenodo. |
| `wiki/` | Curated subset of the APF Obsidian wiki (concept pages + canonical paper pages). |

## 10. Licensing and use

The code and documentation are MIT-licensed (see `LICENSE`). The manuscript and supplement are © the authors; they may be cited and quoted under normal scholarly practice. Cite using the formats in `ai_context/CITING.md`.

---

*This file is generated by the `create-repo` skill from the canonical APF template. To update it across all paper releases, edit `__APF Library/skills/create-repo/templates/AGENTS.md.template` and rebuild. Do not hand-edit in a specific release — it will be overwritten.*
