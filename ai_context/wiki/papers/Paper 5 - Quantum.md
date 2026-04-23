---
type: paper
domain: apf
layer: 5-quantum
created: 2026-04-14
updated: 2026-04-21 (Phase 12.3 I3 + quantum-horizon ingest: v2.0-PLEC/v1.0 → v2.1/v1.1)
sources: []
---

# Paper 5 - Quantum
## "Hilbert Space, the Born Rule, and Tensor Products"

**Status:** **v2.1 main (13pp) + Technical Supplement v1.1 (33pp)** (Phase 12.3 I3 + quantum-horizon coupling ingest, 2026-04-21). Paper 5's Hilbert space is now formalized as the quantum-regime reading $\pi_Q(\mathrm{acc}_{\mathrm{quantum}}(d))$ of the T_ACC ledger; measurement is reframed as a Sector A ↔ Sector B crossing via the interface-sector bridge theorem. Earlier: v2.0-PLEC/v1.0 (2026-04-19 evening).

**Changelog:**
- **v2.1 main + Supplement v1.1 (2026-04-21):** Phase 12.3 I3 + quantum-horizon coupling ingest. Supplement: new §sec:hilbert_acc_quantum Proposition tying the `acc_quantum(d)` factory to the Paper-5 Hilbert-space construction ($\mathcal{H} = \pi_Q(\mathrm{acc}_{\mathrm{quantum}}(d))$, $\dim\mathcal{H} = d = \pi_Q$, $\mathrm{ACC} = \ln d$, Born rule as $\pi_Q$-consistent probability assignment); new §sec:I3_thermo_quantum Theorem stating $I3$ on the Paper-5 Hilbert space ($S_{\mathrm{vN}}(\rho_{\mathrm{max-mixed}}) = \ln\dim\mathcal{H} = \mathrm{ACC}$) with two remarks (thermo-quantum content; entropy ≤ ACC bound reading); new §sec:sectors_and_measurement Corollary "Sector A ↔ Sector B crossings as measurement" framing unitary evolution as Sector-A-preserving target re-pointing and measurement as Sector A → Sector B crossing, with the record register identified at the interface-capacity level with $V_{\mathrm{global}}$ (epistemic tag [P + framing] — not a new proof, a new reading of existing $T_{\mathrm{decoherence}} + T_{\mathrm{interface\,sector\,bridge}}$ content); theorem index extended with 4 new rows. Main: new §sec:summary_TACC "Paper 5 inside the T_ACC ledger" paragraph in Summary; coderefs to $I3$ + bridge theorem + lemma. Archive: `Old/Paper_5_*_pre-Phase12.3.{tex,pdf}`. Upstream reference: Three-Layer Ontology doc §6.

**Core claim:** [[Axiom A1]] forces quantum mechanics from first principles. The structure of admissible distinction-enforcement necessarily gives rise to complex Hilbert space, the [[Born Rule]], and entanglement via tensor products. Not assumed; derived.

**Key theorems:**
- T2 (`check_T2`) — Hilbert space forced by operability constraints
- T3 (`check_T3`) — Noncommutativity of operators is necessary
- T_Born (`check_T_Born`) — [[Born Rule]] probability formula
- T_Tsirelson (`check_T_Tsirelson`) — Tsirelson bound on correlations
- T_CPTP (`check_T_CPTP`) — Completely positive, trace-preserving maps
- [[Tensor Product Structure]] — Multi-particle Hilbert spaces

**Content:**
- Operability bound → Hilbert space structure
- Why complex numbers, not real or quaternionic
- Proof of the Born rule from enforcement cost
- Superposition and entanglement as cost-efficient encodings
- Quantum non-locality without action-at-a-distance
- CPTP maps and the arrow from pure to mixed states
- Connection to quantum information and von Neumann algebras

**Implementation:** `core.py` (48 checks)

**Status notes:**
- Fully implemented in v6.9
- v2.0-PLEC main paper stable; Technical Supplement v1.0 built fresh 2026-04-19 evening as 31pp self-contained formal companion.
- Supplement contents: PLEC + imported lemma restatements; linearity proof from coherence closure; complex Hilbert via Frobenius trichotomy (R ruled out by Regime-R R4, H ruled out by L_loc commutativity); T_tensor; T_Born via Gleason; T_Hermitian + T_CPTP with Stinespring-dilation interface-mediator interpretation; T_decoherence via record-locking channels (Type III regime exit); L_spectral_action_internal with full Boltzmann-cutoff derivation and SM Lagrangian coefficients (a_0=61, a_2≈21.985, a_4≈87.201) from Seeley–DeWitt; L_QG_P1_closure UV completion as Type IV exit; six countermodels; red-team on H1/H2/H3 + Gleason meta-attack; theorem index + dependency diagram.
- **Phase 10 status:** supplement present, pending Stanford agentic review cycle.

## See also
- [[Born Rule]]
- [[Axiom A1]]
- [[Paper 1 - Spine]]
- [[T_ACC Unification]] — Paper 5 Hilbert space = $\pi_Q(\mathrm{acc}_{\mathrm{quantum}})$; $I3$ thermo-quantum realized here
- [[Interface-Sector Bridge]] — measurement as Sector A ↔ Sector B crossing
- [[Paper 3 - Ledgers]] — ACC formalization home (Supplement v1.2 §7)
