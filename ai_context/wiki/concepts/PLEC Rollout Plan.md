# PLEC Rollout Plan — Papers 1 and 2 + Both Supplements

**Created:** 2026-04-16
**Status:** Approved, Phase 0 in progress
**Scope:** Papers 1, 2, Paper 1 Supplement, Paper 2 Supplement
**Deadline context:** May 15, 2026 ToE contest (29 days out at plan creation)

---

## Goal

Promote the **Principle of Least Enforcement Cost (PLEC)** to the foundational statement of APF across all four target documents. A1, MD, A2, BW become four named features of one variational principle, not four competing axioms.

**No mathematical content changes.** All of today's A2 edits stay. All supplement independence proofs stay — reframed as essentiality-of-feature proofs for PLEC's four features. Codebase and theorem counts untouched.

## The Principle

$$
\mathcal{R}(\Gamma) \;=\; \arg\min_{G \,\in\, \mathcal{A}(\Gamma)} \; \sum_{d \in G} \varepsilon(d)
$$

$$
\mathcal{A}(\Gamma) \;=\; \Big\{ G \;:\; \sum_{d \in G} \varepsilon(d) \leq C(\Gamma), \;\; \varepsilon(d) \geq \mu^* > 0 \Big\}
$$

**Four features:**

1. **A1 (upper bound):** Σε(d) ≤ C(Γ). Finite enforcement capacity.
2. **MD (lower bound):** ε(d) ≥ μ* > 0. Enforcement isotropy / positive floor.
3. **A2 (argmin):** realization at the minimum over 𝒜(Γ). Minimum cost selection.
4. **BW (non-degeneracy):** cost spectrum isn't trivially flat. Keeps the argmin interesting.

**Canonical lay statement:**
> *Reality is the minimum-cost expression of distinction compatible with finite capacity.*

**Name convention:**
- Primary: *The Principle of Least Enforcement Cost (PLEC)*
- Formal gloss: *Principle of Admissible Minimum Cost* (for contexts where "least" feels imprecise)
- Recommended: go with PLEC as primary name throughout.

---

## Phase 0 — Foundation (½ session)

Draft PLEC's canonical formal statement in plain prose for review: one-sentence lay form (today's epigraph), one-equation formal form, four named features. Apply to `Reference - APF Axiom Inventory v2.md` as the new primary entry. Archive previous Inventory.

**Deliverable:** Axiom Inventory with PLEC at the top. Canonical name and formal statement frozen.

---

## Phase 1 — Paper 1 + Paper 1 Supplement (1-2 sessions)

**Subordination discipline (applies throughout Phase 1).** Per the amended Axiom Inventory §0.0 (ontological hierarchy), Paper 1 retains Enforceability of Distinction as its opening claim. PLEC enters as the variational formulation of enforceability in admissible regimes, introduced *after* the enforceability framing is in place — not as a replacement for it. The locked claim *physical meaning is grounded in enforceable distinction* appears prominently in the abstract and/or §1 intro, with PLEC then introduced as its canonical variational expression.

### Paper 1 main body

- **Abstract.** Current opener ("This paper presents a quantum reconstruction program that begins from a single physical fact: the universe is finite") stays. Add one sentence after it explicitly naming Enforceability of Distinction as the ontological foundation. Add a subsequent sentence introducing PLEC as the canonical variational formulation of that foundation in admissible regimes. Epigraph stays.
- **§1 intro.** Current finite-universe / enforceability framing stays (unchanged spine). PLEC enters as a named variational expression after the enforceability case is made. Include the subordination sentence: *PLEC is the canonical variational formulation of the finite-enforceability framework; it does not replace the locked claim that physical meaning is grounded in enforceable distinction.*
- **§2 existing A1 + A2 material.** Promoted to a section titled *Enforceability and its variational expression* (not "The Principle of Least Enforcement Cost" alone — the title must signal the hierarchy). Structure: subsections for A1 (enforceability, with finiteness) and MD (enforceability floor) as the enforceability side; A2 (argmin) and BW (non-degeneracy) as the variational side; then PLEC as the unified variational formulation of both sides. A1/MD tcolorboxes and A2/BW tcolorboxes become components of this structure.
- **Load-bearing line.** "A1 and A2 together are load-bearing" → "Enforceability of Distinction (with its finite-capacity and positive-floor features) is the ontological foundation; PLEC is its variational formulation. Both are load-bearing."
- **L_ε*** stays named "Minimum Cost Floor" (rename completed earlier).
- **Body shorthand** "A1/A2/MD/BW" stays — each now denotes a structurally necessary component of PLEC, carrying enforceability content independently of PLEC's variational packaging (per §0.3 as amended).

### Paper 1 Supplement

Supplement is the natural home for deeper PLEC reframing — readers of the supplement are already committed to the spine.

- Section title "Inputs: Two Axioms, One Structural Regularity Condition, One Richness Premise" → "Enforceability of Distinction and the Principle of Least Enforcement Cost: Foundation and Variational Formulation."
- A2 tcolorbox countermodel (G₁/G₂ at C/2 and C) recast as **essentiality-of-component-A2** proof: without A2, admissible set contains both; with A2, only the cheaper one is realized.
- Same reframing for MD (essentiality-of-component-F2) and BW (essentiality-of-component-F4).
- Independence proofs preserved verbatim; relabeled as essentiality proofs for PLEC's structurally necessary components.
- Four-layer hierarchy from Inventory §0.0 transliterated into the Supplement intro.
- "Four inputs" → "one principle, four structurally necessary components, all grounded in layer-1 enforceability."

**Per pair protocol:** draft → approve → archive (`_pre-PLEC` suffix in Old/) → apply → compile (pdflatex × 2) → verify clean.

**Deliverable:** Paper 1 v4.0-PLEC (≈40pp), Paper 1 Supplement v2.1-PLEC (≈80pp), both compile-clean.

---

## Phase 2 — Paper 2 + Paper 2 Supplement (1-2 sessions)

### Paper 2 main body

- Abstract adds one sentence placing Paper 2 inside PLEC.
- §1 Paper 1 recap updated to name PLEC as foundational.
- Glossary / registry tables gain a PLEC row at top; A1/A2/MD rows relabeled as features.
- Today's 15 A2 attribution edits stay verbatim — "A2 selects" is correct because A2 is the argmin feature of PLEC.
- θ_QCD argument (Edit #13) stays; it's already a PLEC argmin argument.

### Paper 2 Supplement

- EC_CM_from_A1 lemma block (today's Diff 5): EC clause stays "A1 alone"; CM clause now reads "A2 applied to the admissible set — the argmin feature of PLEC."
- Premises table adds PLEC row; A1/A2/MD rows relabeled as features.
- "Correct characterisation" paragraph updated accordingly.
- Pre-existing undefined refs (`sec:vac`, `sec:equipartition`) remain deferred.

**Per pair protocol:** same draft→approve→archive→apply→compile. Paper 2 Supplement compiles via `/tmp` (filename has spaces).

**Deliverable:** Paper 2 v6.0-PLEC (≈47pp), Paper 2 Supplement v2.1-PLEC (≈96pp), both compile-clean.

---

## Phase 3 — Propagation (½ session)

- New wiki page: `wiki/Principle of Least Enforcement Cost.md` as canonical landing.
- Update: `wiki/Axiom A1.md` (recast A1 as upper-bound feature of PLEC).
- Update: `wiki/Paper 1 - Spine.md`, `wiki/Paper 2 - Structure.md` (v4.0-PLEC / v6.0-PLEC, PLEC changelog).
- Update: `wiki/Paper Update Work Plan.md` (new Phase 2.7 entry, marked complete).
- Update: `wiki/Index.md` (version tag updates).
- Update: `APF Reference Docs/Reference - APF Paper Update Work Plan v2.md` (Phase 2.7 detail block).
- Update: `APF_Paper_Index.md` (version bumps).
- New: `APF_Rollout_Timeline_v1.9.md` (successor to v1.8).
- Update: both `CLAUDE.md` files (rollout timeline ref v1.8 → v1.9).

---

## Phase 4 — Verification and signoff

- Final compile of all four PDFs.
- Side-by-side read of first 2 pages of each — verify consistent PLEC framing.
- Full signoff procedure.

---

## Sequencing rationale

Paper 1 first (spine, sets pattern). Main before supplement in each pair (supplements cross-reference their paper's framing). Phase 0 before Phase 1 (Inventory is the authoritative reference paper edits point to). Phases 1 and 2 can parallelize across days but Phase 2 reads smoother once Phase 1's patterns are set.

## Out of scope

- Papers 0, 3-7, 13, 14 (PLEC propagation deferred to post-contest Phase 2.8).
- Website / video / contest submission (Tracks A, D — separate).
- Codebase (PLEC is presentation only; theorems unchanged).
- Paper 2 Supplement pre-existing undefined refs (optional cleanup, defer).

## Risk register

| Risk | Mitigation |
|---|---|
| Feature-essentiality proof weaker than today's independence proof | Essentiality proofs must preserve or strengthen countermodel content; diff against today's text before commit. |
| Stanford-reviewer concern about merging MD into A1 | PLEC does NOT merge features — it names them together under one principle but keeps them structurally separate. Called out explicitly in Paper 1 Supplement. |
| Contest-deadline squeeze (May 15) | If Phase 2 runs long, ship Paper 1 + Supplement alone with PLEC; Paper 2 + Supplement ship post-contest. Phase 1 is a publishable standalone. |
| Cross-reference breakage | EC_CM_from_A1 label preserved for back-compat (7+ call sites). tcolorbox IDs preserved. |
| Paper 2 Supplement compile fragility | Compile via /tmp (filename has spaces). Same workflow as today. |

## Working rules that apply throughout

- Working Rule 7 (draft review before applying) at every "Draft … in plain prose" step.
- Working Rule 8 (compile after every edit) at every "Apply edits" step.
- Archive convention: `Old/` folder, `_pre-PLEC` suffix.
- Paper 2 Supplement: compile via `/tmp`.
- Independence/essentiality proofs: preserved and strengthened, never deleted.

## Progress tracking

- [x] Phase 0: Plan committed to wiki (2026-04-16)
- [x] Phase 0: PLEC canonical statement drafted and approved (2026-04-16)
- [x] Phase 0: Axiom Inventory updated to v2.1 with §0 PLEC consolidation + §10 rollout/placement tracking (2026-04-16)
- [x] Phase 1: Paper 1 main body (2026-04-16, v4.0-PLEC)
- [x] Phase 1: Paper 1 Supplement (2026-04-16, v2.1-PLEC)
- [x] Phase 2: Paper 2 main body (2026-04-16, v5.3-PLEC)
- [x] Phase 2: Paper 2 Supplement (2026-04-16, v2.1-PLEC)
- [x] Phase 2.8a: Paper 3 v1.1-PLEC + Paper 13 v8.3-PLEC (2026-04-17)
- [ ] Phase 2.8b: Papers 0, 4, 5, 6, 7, Enforcement Crystal — source-pending (awaiting editable .tex)
- [ ] Phase 3: Propagation (wiki concept pages, reconstruction-program docs)
- [ ] Phase 4: Verification and signoff

**Status:** Phase 2.8a complete (Paper 3 + Paper 13). Phase 2.8b blocked on source availability for remaining targets.
