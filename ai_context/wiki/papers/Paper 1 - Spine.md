---
type: paper
domain: apf
layer: 1-spine
created: 2026-04-14
updated: 2026-04-17 (Paper 0 audit edit pass)
sources: []
---

# Paper 1 - Spine
## "The Enforceability of Distinction"

**Status:** v4.0-PLEC (April 2026), **44pp** (audit-edit pass applied 2026-04-17; wayfinder convention added 2026-04-18). PLEC rollout complete — Enforceability of Distinction (ontological foundation) + PLEC (variational formulation) with four structurally necessary components A1, MD, A2, BW. Five interpretive bridge boxes added grounding the derivation chain in concrete physics (Pauli exclusion, uncertainty relation, Planck quantization, Stern-Gerlach, Aharonov-Bohm, coupled superconducting resonators). New "Operational primitives" wayfinder box at start of §2 mapping distinction / maintained / enforcement cost to their motivation, formalization, and measurement-protocol sites.

**Core claim:** Physical meaning is grounded in *enforceable distinction*; the [[Principle of Least Enforcement Cost]] (PLEC) is its canonical variational formulation in admissible regimes. PLEC has four structurally necessary components: [[Axiom A1]] (Finite Enforcement Capacity — PLEC upper bound), **MD** (Minimum Distinction — PLEC positive floor), **A2** (Minimum Cost — PLEC argmin selection), **BW** (Bound Wideness — PLEC non-degeneracy). Together they force the existence of the canonical object — a complex Hilbert space with the Born rule, finite-dimensional irreps of $SU(3) \times SU(2) \times U(1)$, and deterministic mass/mixing ratios. Zero free parameters.

**Key theorems (22 boxed results, PLEC framing):**
- §1.3 tcolorbox "The Principle of Least Enforcement Cost: four structurally necessary components" (A1 upper bound, MD positive floor, A2 argmin, BW non-degeneracy)
- §2.1 "Enforceability and its variational expression" with existing A1/A2/MD/BW tcolorboxes preserved; load-bearing line: "Enforceability of Distinction... is the ontological foundation; PLEC is its variational formulation; none of the four components is derivable from the others (see Technical Supplement essentiality proofs)"
- T_sep (Separation Biconditional), **L_ε* (Minimum Cost Floor — renamed to avoid collision with A2)**, L_cost (Cost Uniqueness), L_loc (Locality), NT (Non-Degeneracy), L_nc (Non-Closure)
- L_Δ (Superadditivity — the watershed), L_irr (Irreversibility)
- T1 (Order-Dependence — the bridge theorem), T_adj (Self-Adjointness), L_blk + T_alg (Noncommutativity)
- T2 (Hilbert Space), T_Born (Born Rule), T_Tsirelson (Tsirelson Bound)
- T1-NB (Non-Boolean second route), T_canonical, Conditional theorems, L_M/L_NT derived

**Content:**
- [[Axiom A1]] statement and physical grounding (§2)
- Definitions and structural lemmas — arena chain (§3)
- The formal bridge: from finite capacity to Hilbert space (§4)
- Consequences, applications, and physical content (§5)
- Scope, assumption inventory, identifications (§6)
- Three-act narrative arc with plain-English intros/conclusions for all 22 theorem boxes

**Implementation:** `paper1.py` (82 checks) + `supplement.py` (44 checks, strict subset)

**Theorem register:** 82 entries in 6 sections, with `†` markers for supplement.py checks.

**Changelog:**
- 2026-04-17: v4.0-PLEC audit-edit pass — 12 surgical edits driven by the Paper 0 ontology audit. Four copy-edits: §1.1 "The question is then immediate" → "With enforcement named, the framework's central question can be stated..."; §1.3 "The APF answers a prior question" → "one step upstream of these programs"; §2 "was quietly taking" → "had left implicit"; §3.2 "enforcement-theoretic origin of irreversibility" → "how irreversibility emerges from finite capacity --- not as a contingent statistical fact, but as a structural feature of the admissibility budget itself". Two figure fixes: Fig 3 Arena DAG C<∞ edge label repositioned (pos 0.3 → 0.72); Fig 5 Bridge DAG red→black on all right-column "classical model fails" labels + arrows, labels shifted x=6.0→7.0 to clear BW node, T_adj-targeted arrow rerouted with `bend right=25` (previously landed at phantom `Tadj.east+(3.2,0)` cutting through T1). Mead 2000 "identically" claim investigated against local book copy (§§4.5–4.8; equations 4.17, 4.18, 4.20, 4.22, 4.27, 4.31, 4.32 verified) — dynamical-identity claim retracted as overclaim (Mead's equations are dynamical, L_Δ is static superadditivity); restored as a physical-anchor tcolorbox with verified equation citations and structural-mapping discipline. Five interpretive bridge tcolorboxes added: (i) L_nc→L_Δ "The gap where classical models end" with Pauli-exclusion and uncertainty-relation examples; (ii) L_ε*→L_cost "From a floor to a ladder" with Planck's step; (iii) the T1/T_adj/L_blk/T_alg chain grounded in Stern-Gerlach three-apparatus sequence; (iv) T_2c "Why exactly ℂ" with Aharonov-Bohm anchor; (v) coupled superconducting resonators as L_Δ instance with Mead equations. Clean compile: 40pp → 44pp. Archive: `Old/paper1_v4.0-PLEC_pre-P0audit.tex`/`.pdf`.
- 2026-04-16: v4.0-PLEC — PLEC rollout (presentation-layer only, no theorems/proofs changed). Four edit sites on main body: (1) abstract expanded with Enforceability of Distinction ontological foundation + PLEC as variational formulation; (2) §1 intro subordination paragraph + four-component preview; (3) §1.3 tcolorbox retitled "The Principle of Least Enforcement Cost: four structurally necessary components," body rewritten to cover A1/MD/A2/BW; (4) §2.1 renamed "Enforceability and its variational expression," §2.0 preamble recast as "arena of PLEC," load-bearing line rewritten with subordination framing and forward-reference to Supplement essentiality proofs. Supplement companion: 6 edits including new `\label{sec:essentiality}` anchor and added A2 row in the "What is assumed vs. what is derived" table (previously missing). 40pp main, 81pp supplement. Clean compile both.
- 2026-04-16: v3.2-A2 — Dual-axiom rollout. A1 box title → "Finite Enforcement Capacity"; new §2.0 admissibility-structure preamble; new A2 box (formal statement, what A2 says/does not say, tie-breaking, non-derivability from A1); load-bearing line updated to "A1 and A2 together"; L_ε* renamed "Minimum Cost" → "Minimum Cost Floor" at three sites to avoid collision with A2. 40pp.
- 2026-04-15: v3.1 — Tables full-width, reader's guide removed, theorem register 41→82, DAGs updated (5 new nodes), narrative intros/conclusions for all 22 boxes, code files saved to codebase.
- 2026-04-14: v3.0 — Congruity pass, content updates, architectural split (37pp). Supplement 78pp.

## See also
- [[Axiom A1]]
- [[Enforcement Budget]]
- [[Paper 2 - Structure]]
- [[Derivation Chain]]
