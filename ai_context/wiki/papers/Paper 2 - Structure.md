---
type: paper
domain: apf
layer: 2-structure
created: 2026-04-14
updated: 2026-04-21 (Phase 12.5 I2 at L_count seam: v5.3-PLEC/v2 → v5.4/v2.1)
sources:
  - Paper_2_Structure_of_Admissible_Physics_v5.4.tex (48pp)
  - Paper_2_Structure_of_Admissible_Physics_Supplement_v2.1.tex (97pp, ~8,250 lines)
  - gauge.py + unification.py (codebase v6.9)
---

# Paper 2 - Structure
## "The Structure of Admissible Physics"

**Status:** **v5.4 main (48pp) + Technical Supplement v2.1 (97pp)** (Phase 12.5 I2 ingest at the L_count seam, 2026-04-21). Paper 2 establishes $C_{\mathrm{total}} = 61$ via $L_{\mathrm{count}}$; Phase 12.5 registers that this integer is the same integer T_ACC identity $I2$ pins (shared between $\pi_F$ gauge projection and $\pi_C$ cosmological denominator). Earlier: v5.3-PLEC/v2 (April 2026 PLEC rollout).

**Changelog:**
- **v5.4 main + Supplement v2.1 (2026-04-21):** Phase 12.5 I2 at L_count seam. Supplement: new §sec:I2_at_L_count at the end of §sec:capacity_count (L_count) titled "I2 (gauge-cosmological identity): the K=61 shared integer" with Proposition stating the integer identity between π_F reading at SM interface and π_C denominator (items (i)–(iv) linking L_count's 61 = 45 + 4 + 12 to the Paper 6 cosmological-partition denominator); prose paragraph on why the 1.2% Planck match is evidence for (not content of) the identity; subspace-refinement-via-bridge paragraph naming Paper 6 Supplement §sec:Bridge_compressed as strictly finer than I2; interactive figure reference to `Fig_61_Channels_Bridge_to_102_v6.html` in Paper 0 folder. Code-theorem concordance table gains I2 row. Main Summary gains paragraph placing I2 at the structural-enumeration / cosmological-partition seam with coderef to `check_I2_gauge_cosmological`. Incidental fix: added `\DeclareMathOperator*{\argmin}{arg\,min}` to main preamble (v5.3-PLEC relied on an auto-loaded macro absent from sandbox TeX Live — pre-existing issue inherited). Archive: `Old/Paper_2_*_pre-Phase12.5.{tex,pdf}`. Upstream reference: Three-Layer Ontology doc §§1 and 8.
- **v5.3-PLEC (April 16, 2026)** — PLEC rollout (presentation-layer only, no theorems/proofs changed). Main body: 4 edits — (1) epigraph block with canonical sentence + abstract reframed around "four structurally necessary components of PLEC" plus subordination statement; (2) §1 intro opener expanded with PLEC naming sentence + four-component load-bearing summary; (3) §1 "cost minimization under structural constraints" line replaced with full PLEC framing naming A1/MD/A2/BW roles; (4) end of §5 "single axiom about finite enforcement capacity" replaced with "Principle of Least Enforcement Cost---expressed through its four components." Technical Supplement: 4 edits. Compile: main 47 pages stable, supplement 96 pages stable. Archives: `Old/paper2_pre-PLEC.tex`/`.pdf`, `Old/Paper2_Formal_Supplement_v2_pre-PLEC.tex`/`.pdf`.

**Core claim:** Physical meaning is grounded in *enforceable distinction*; under the [[Principle of Least Enforcement Cost]] (PLEC) — expressed through its four structurally necessary components A1 (finite capacity / upper bound), MD (positive floor), A2 (argmin), BW (non-degeneracy) — the admissible set $\mathcal{A}(\rho, R)$ is well-defined and the realised gauge structure is the PLEC argmin over that admissible set. The Standard Model gauge group $SU(3) \times SU(2) \times U(1)$ and field content (45 fermions) drop out as the unique min-cost admissible configuration. Zero free parameters.

**Key theorems:**
- [[Non-Closure Theorem]] (`check_L_NZ`) — Admissible state spaces cannot be closed under all automorphisms
- L_gauge_template_uniqueness (`check_gauge_uniqueness`) — Only one consistent gauge structure avoids closure
- T_gauge (`check_T_gauge`) — Forcing $N_c = 3$ (color), $N_w = 2$ (isospin), $U(1)_Y$ hypercharge
- [[Gauge Uniqueness]] — Standard model group as the *only* admissible choice
- T4E/T4F — Flavor-capacity saturation: $E(3)/C_{EW} = 3/4$, mass hierarchy from capacity ordering
- T_theta_QCD (`check_T_theta_QCD`) — Strong CP dissolved: $\theta = 0$ admissible at zero maintenance cost; $\theta \ne 0$ admissible at $\varepsilon^* > 0$; **A2 selects the minimum** as the realised value
- **EC and CM (Paper 2 Supplement, lemma `EC_CM_from_A1`)** — EC follows from A1 alone via the derived structural separation; CM follows from A2 applied to the admissible set defined by A1 + MD
- T_vac,stab (`check_T_vacuum_stability`) — Vacuum absolute stability: unique enforcement well
- T_proton (`check_T_proton`) — Proton absolute stability: exact B conservation at Bekenstein saturation

**Content:**
- Proof of non-closure from enforcement cost growth
- Gauge template mechanics: irreps, constraints, forbidden closures
- Why $SU(3) \times SU(2) \times U(1)$ is unique
- Field content counting: 45 fermions (15 per family × 3 families)
- Baryon and lepton number emergence from $U(1)_Y$
- Generation argument strengthened with flavor-capacity saturation (T4E/T4F)
- Three structural predictions (§11): $\theta_{QCD} = 0$, vacuum stability, proton stability
- Comparison to other gauge unification proposals

**Dual-document architecture:**
- Main paper: `paper2.tex` (~1,580 lines, argument-first, 47 pages)
- Technical Supplement: `Paper2_Formal_Supplement_v2.tex` (~7,900 lines, proof-dense, 96 pages)
- Supplement §14 contains formal proofs for all three structural predictions (added April 2026)

**Implementation:** `gauge.py` (29 checks), `generations.py` (T4E/T4F)

**v5.3-PLEC changelog (April 16, 2026):**
- PLEC rollout (presentation-layer only, no theorems/proofs changed). Main body: 4 edits — (1) epigraph block with canonical sentence + abstract reframed around "four structurally necessary components of PLEC" plus subordination statement; (2) §1 intro opener expanded with PLEC naming sentence + four-component load-bearing summary; (3) §1 "cost minimization under structural constraints" line replaced with full PLEC framing naming A1/MD/A2/BW roles; (4) end of §5 "single axiom about finite enforcement capacity" replaced with "Principle of Least Enforcement Cost---expressed through its four components."
- Technical Supplement: 4 edits — (1) abstract reframed to introduce PLEC as variational formulation of enforceability-of-distinction foundation; (2) "All upstream results" inventory updated to name "the four PLEC components A1, MD, A2, BW; constitutive semantics SP, K3, FD1--FD6"; (3) "What Is Assumed, What Is Proved, What Is Deferred" tcolorbox — PLEC essentiality intro added, A2 inserted as item 3 (previously missing), PLEC role annotations on all four components, prose reference to Paper 1 Supplement §1 essentiality proofs; (4) "named assumption count" comparison paragraph includes A2 and uses PLEC component framing.
- Compile: main 47 pages stable (after resolving ToC-length oscillation), supplement 96 pages stable. Two pre-existing undefined refs in supplement (`sec:vac`, `sec:equipartition`) verified not introduced by this rollout.
- Archives: `Old/paper2_pre-PLEC.tex`/`.pdf`, `Old/Paper2_Formal_Supplement_v2_pre-PLEC.tex`/`.pdf`.

**v5.2-A2 changelog (April 16, 2026):**
- Dual-axiom rollout. Main body: 15 attribution edits A1 → A2 across §3 (Stage 3 cost closure), §6 (one doublet not realised by A2), §7.4 (θ_QCD now framed as A2 selecting min-cost admissible value), Appendix A registry (L_ε* row), Appendix B symbol index (new A2 row).
- Technical Supplement: 5 linked edits to `EC_CM_from_A1` lemma block — lemma statement reattributes EC to A1 alone and CM to A2; proof of CM rewritten as explicit A2 argmin step ($G' = G \setminus \{X\}$ admissible and strictly cheaper by $\varepsilon^*$); premises table gains A2 row; "correct characterisation" paragraph and "CM is not a selection principle" paragraph reattribute CM to A2; §1812 Theorem B footnote updated.
- All 7+ call sites of `lem:EC_CM_from_A1` elsewhere in the supplement still resolve correctly (label preserved).

**v5.1 changelog (April 2026):**
- Terminology pass: "Formal Supplement" → "Technical Supplement" throughout
- Code reference audit: 10 stale refs fixed, 9 missing refs added
- New §11: three structural predictions with full argument sketches
- Generation strengthening: T4E/T4F flavor-capacity saturation paragraph
- Appendix updates: 8 new theorem index rows, 3 symbol index entries, 4 summary table rows
- Supplement: 3 new formal theorem/proof sections (§14.1–14.3), concordance table updated
- Consistency audit: all paper↔supplement cross-references verified

## See also
- [[Non-Closure Theorem]]
- [[Gauge Uniqueness]]
- [[Paper 4 - Constraints]]
- [[Paper 1 - Spine]]
- [[T_ACC Unification]] — Paper 2's $C_{\mathrm{total}} = 61$ is the integer $\pi_F$ reads; $I2$ is the gauge-cosmological identity at this seam
- [[Interface-Sector Bridge]] — the subspace-level refinement of $I2$ (Paper 6 supplement)
