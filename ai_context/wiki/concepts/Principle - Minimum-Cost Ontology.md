# APF Framing Principle: What Is, Is Minimum-Cost Admissible

**Status.** Foundational reframe — applies to foundations, abstracts, and public-facing materials across all APF papers. Derivation-chapter prose retains active-voice scaffolding for pedagogy; see §Implementation.

**PLEC cross-reference (2026-04-17).** The canonical one-line statement of this principle — *"reality is the minimum-cost expression of distinction compatible with finite capacity"* — is now the **lay statement of the Principle of Least Enforcement Cost (PLEC)**, formalized in Papers 1 v4.0-PLEC and 2 v5.3-PLEC with the **four-component decomposition**: A1 (finite capacity; upper bound), MD (positive cost floor), A2 (argmin selection), BW (budget-window / non-degeneracy). Each component is structurally necessary, with explicit countermodels in the Paper 1 Technical Supplement §1 showing no three of the four imply the fourth. This page remains the canonical prose presentation of the ontological reading; PLEC is its canonical variational formulation. The "A1 admits / A2 selects" partition in §2 below is now more finely resolved by the four PLEC components, but the descriptive-framing discipline of this page is preserved in all PLEC-era documents.

**Relationship to the formalism.** This is prose alignment with existing formal content, not new physics. $L_{\mathrm{col}}$ is already an $\arg\min$; A1 is already an inequality; $L_{\mathrm{cost}}$ is already a proportionality. None of these are processes. The change is that the English now matches.

---

## 1. The Principle

**What exists is the minimum-cost admissible configuration. Nothing more is being said.**

There is no filter, no selection, no process, no mechanism. A1 is not an axiom about what happens; it is a description of what is. Enforcement cost is a property of configurations — like mass, or charge, or dimension. Some configurations have higher cost, some have lower cost, and some exceed the capacity bound entirely. The realized configuration is the one with the lowest cost inside the bound. The rest are conceptual comparisons — configurations whose cost we can compute, but not configurations that were considered and rejected, because no one is doing any considering.

## 2. Precise Statement

Let $\mathcal{A}$ be the set of configurations $G$ satisfying $\sum \varepsilon(d) \leq C$. The realized configuration is

$$G_* \;=\; \arg\min_{G \in \mathcal{A}} \sum \varepsilon(d).$$

"Enforcement" names the cost structure $\varepsilon$, not an act of enforcing. Non-realized configurations have higher cost than $G_*$ or exceed $C$; they are therefore not what is. They are not rejected. Rejection implies a rejector.

## 3. Two Descriptive Conditions on What Exists

**Condition 1 — Admissibility.** What exists has $\sum \varepsilon(d) \leq C$. Configurations violating this are physically inaccessible under A1; they are not part of the actual world.

*(Earlier drafts of this document compared inadmissibility to "a triangle with four sides." That analogy equated logical and nomological impossibility, which is a stronger claim than the framework needs. A1 is a physical postulate: inadmissible configurations are ruled out by A1, not by the laws of thought. "Physically inaccessible under A1" is the honest phrasing.)*

**Condition 2 — Minimality.** Among admissible configurations, what exists is the one with lowest total cost. Higher-cost admissible configurations are coherent to describe and to compute about, but they are not what is. No competition is occurring; no contest is being run. The minimum-cost configuration is simply the referent of "what exists."

Both conditions are descriptive facts about $G_*$. Neither involves an action, a choice, or a process.

## 4. The Stationary-Action Precedent (Promoted)

This is not a new move in physics. Classical mechanics underwent exactly this de-teleologization between Maupertuis and Lagrange. Maupertuis treated the principle of least action as evidence of divine design: God chose the path of minimum action, and the particle followed. Lagrange and Hamilton retired the teleology. What remained was a variational principle: the physical trajectory is the one satisfying $\delta S = 0$. Alternative paths — through Mars, wiggly, indirect — are not paths the particle considered and rejected. They are mathematical objects the analyst introduces to characterize the real trajectory by contrast. No one now says the particle evaluates all paths and picks one.

APF is the same move applied one level higher. Instead of "which trajectory does a particle follow," we ask "which gauge configuration / fermion spectrum does the universe instantiate." The answer has the same structure: *the one at the minimum of a well-defined cost functional, subject to constraints*. Alternative configurations — $\mathrm{SU}(5)$, $\mathrm{SO}(10)$, fermion templates with 48 or 63 Weyl DOF — are mathematical ghosts. They are configurations whose cost we can compute, whose admissibility we can check, whose properties we can contrast against what exists. They are not configurations the universe considered and declined. Nothing considers.

This is the strongest single defense against the "but what does the enforcing?" objection. The answer is the same as "what picks the classical trajectory?" — nothing. There is a cost structure; there is a minimum; the world is at the minimum. Two-centuries-old physics, already universally accepted.

## 5. Water Analogy (Supporting)

Water is at the lowest point available to it. Nothing pushes water downhill — gravity provides the potential, but "push" smuggles agency into a relation. Nothing selects the lowest point from alternatives. Nothing filters out higher points. Water is where it is, and where it is happens to be low. The alternative framings — "water chooses the lowest point," "gravity selects the minimum potential," "a filter rejects higher configurations" — are verbally active but the physics underneath is purely descriptive. There is a potential; there is a minimum; the water is at the minimum.

APF sits at the minimum of a cost structure the way water sits at the bottom of a valley.

## 6. The Chain, Fully Descriptive

| Step | Active framing | Descriptive framing |
|------|----------------|---------------------|
| Meaning | The universe maintains a distinction | A distinction is present |
| Robustness | It resists erasure | No admissible configuration lacks it |
| Enforcement | Resources are committed to it | It has nonzero enforcement cost |
| Finiteness | Resources are bounded | The capacity bound is finite |
| What is realized | The universe selects minimum structure | What exists has minimum cost |

The logic is unchanged. No verbs remain that imply agency, action, or process.

## 7. Where This Lives in the Formalism

The descriptive reading is already consistent with the existing formal content:

1. **$L_{\mathrm{col}}$** states $G_{\mathrm{realized}} = \arg\min_{G \in \mathrm{viable}} \dim(G)$. The $\arg\min$ operator returns an element, not a process.
2. **$L_{\mathrm{cost}}$** proves $\dim(G)$ is the unique cost measure (Cauchy uniqueness). A measure is a descriptive property of configurations.
3. **Lemma $EC_{CM}$ from $A1$** (Paper 2 Supplement) proves non-minimal configurations waste capacity, making them inadmissible. This is a descriptive fact about their cost, not an action taken against them.

Every use of active language in the papers is pedagogical convenience built on top of these flat descriptive theorems. The theorems contain no processes.

## 8. Worked Example

The Standard Model gauge group has $\dim = 12$. $\mathrm{SU}(5)$ has $\dim = 24$. $\mathrm{SO}(10)$ has $\dim = 45$. The total capacity is $C_{\mathrm{total}} = 61$. All three cost less than 61, so all three are admissible.

The realized gauge group has the lowest cost of the three. The realized gauge group is $\mathrm{SU}(3) \times \mathrm{SU}(2) \times \mathrm{U}(1)$. $\mathrm{SU}(5)$ and $\mathrm{SO}(10)$ are not what is. No further story is available, and none is needed. The question "why is the universe at the lower-cost option?" is not a further question — it is the definition of what the cost structure describes. Asking it is like asking why water is at the lowest point when the lowest point is where water is.

## 9. Cosmological Evolution, Descriptively

A cosmological trajectory is a one-parameter family of minimum-cost admissible configurations, parameterized by the *committed-capacity fraction* $s = k/C$ (the saturation level of the enforcement ledger). As $s$ takes different values, the minimum-cost admissible configuration takes different values. What we call "forces" and "evolution" is the functional dependence $G_*(s)$.

*The honest complication.* This reframing handles equilibrium beautifully but strains on the approach to equilibrium. If $s$ itself evolves, we need an equation of motion for $s$, and that equation of motion looks like a dynamical process — which is the language the descriptive framing was meant to retire. Two responses are available:

(i) **Slice-reading.** Treat $s$ as a purely mathematical parameter indexing a synchronic family of minima. The framework predicts the configuration at each $s$; it does not predict a temporal ordering among them. "Evolution" on this reading is a chart of the family, not a process traversing it. This is the cleanest descriptive reading but may be narrower than the framework's empirical ambitions require.

(ii) **Dynamical-emergent-reading.** Treat $s$ as a time-like parameter whose equation of motion is itself derived from the descriptive structure. In classical mechanics this works: $\dot q = \partial H / \partial p$ is a mathematical consequence of the Hamiltonian structure, not an additional law. If APF's committed-capacity dynamics emerge from the cost structure alone (e.g., via monotone non-decrease under $L_{\mathrm{irr}}$), "evolution" is the framework reading off its own time-parameter; no external process is required.

Paper 6 (cosmological evolution) is where this is decided. The descriptive framing commits the framework to showing that $\dot s$ is a consequence of $\varepsilon$ and A1, not a separate input.

The matching transition at $s_{\mathrm{crit}} = 1/102$ is the cleanest case: $G_*(s)$ is qualitatively different above and below this value. What looks like "the 42/19 partition crystallizing" is the minimum of $\varepsilon$ jumping to a different configuration as the constraint landscape crosses a threshold. No crystallization process is occurring; the function $G_*(s)$ has a discontinuity.

## 10. Why This Matters

1. **It removes every "who does it?" question.** Under active framing, "what enforces?", "what selects?", "what filters?" are natural questions with no satisfying answers. Under descriptive framing, these questions don't arise. They presuppose a process that isn't there.
2. **It aligns the prose with the formalism.** A1 is an inequality. $L_{\mathrm{col}}$ is an $\arg\min$. $L_{\mathrm{cost}}$ is a proportionality. None of these are processes. The formal content is descriptive; the prose should match.
3. **It preempts the anthropomorphism objection entirely.** There is nothing to anthropomorphize. The universe is not choosing, selecting, filtering, or enforcing. It is at the minimum of a cost structure, the way a classical particle is at the stationary point of the action.
4. **It sharpens the SM-vs-GUT explanation.** "Why $\mathrm{SU}(3) \times \mathrm{SU}(2) \times \mathrm{U}(1)$ and not $\mathrm{SU}(5)$?" has a one-line answer: *the first has lower cost and both are admissible.*
5. **It removes the temptation to reify the correlation space.** If there is no filter, there is no space of outputs the filter ranges over. There is one minimum, one configuration, one set of facts. The "space of admissible configurations" is a mathematical staging ground for the $\arg\min$, not an ontological layer.

## 11. Caveats and Limits

This framing makes commitments that deserve explicit acknowledgment:

- **Modal commitment.** Even "what is, is minimum-cost admissible" ranges implicitly over possibilities. The framework treats non-actualized configurations as mathematical objects (ghosts) whose costs are computable — a mild Platonism about possible-but-not-actual configurations. The alternative (strict actualism: non-actualized configurations don't exist in any sense) is also consistent but forecloses the arg-min language. The framework chooses the Platonist path for pragmatic reasons.
- **Admissibility is physical, not logical.** Inadmissible configurations are ruled out by A1, a physical postulate. They are not ruled out by the laws of thought. Comparison to logical impossibilities ("four-sided triangle") overstates the case.
- **Dynamics may re-import process language.** See §9. If Paper 6 ends up with an equation of motion for $s$ that cannot be read as a consequence of $\varepsilon$ alone, the descriptive framing breaks down there. The framing is strongest at equilibrium / at the late-time attractor.
- **This is presentational, not substantive.** The formalism already supports both readings. The descriptive framing is chosen because it aligns prose with formalism and preempts objections, not because it makes new predictions.

## 12. Implementation

- **The word "enforcement" stays.** It names the cost structure $\varepsilon$, not an action. The notation is intact; only the surrounding prose changes.
- **Tier the rollout.** Foundations, abstracts, Paper 0, Paper 13 §2, and public-facing materials go fully descriptive. Derivation chapters (Papers 1–7 main-text proofs) keep active-voice pedagogical scaffolding because the algorithmic narrative ("apply this filter; apply the next filter; the survivor is unique") is genuinely clearer in active voice, and the formal content is already descriptive regardless of the prose.
- **Active verbs to replace in foundations/abstracts:** "The universe enforces" → "the configuration has cost." "Selects" → "is" or "has minimum cost among." "Filters out" → "has higher cost than." "Rejects" → never appears.
- **Paper 13 §2 replacement (proposed).** *"A distinction present in a region has a positive enforcement cost. The total cost of all distinctions present in a region is bounded by the region's finite capacity. What is present in a region is the minimum-cost configuration consistent with that bound."* (Earlier drafts used "maintained distinction"; "present" is more consistent with the principle, since "maintained" smuggles back the activity being retired.)
- **"Selection framing" retires in favor of "descriptive framing" or "minimum-cost ontology."** The word "selection" still suggests a process.
- **Correlation space as ontological layer is dropped.** What A1 implies is a cost structure and a minimum, not a space of candidate configurations. Paper 8 is reframed accordingly.
- **Forces as descriptive.** Forces are the classical-limit appearance of the minimum shifting as conditions shift, not causes of change.

## 13. One-Liner for Public Materials

*The universe isn't being enforced. It's at the minimum of a cost structure — the same way a classical particle is at the stationary point of the action, or water is at the bottom of a valley — and everything else follows.*

---

## Appendix: Changelog from the original draft

- **§3 Condition 1.** Softened "not possible, the way a triangle with four sides is not possible" to "physically inaccessible under A1." The four-sided-triangle analogy conflated logical and nomological impossibility.
- **§4 Stationary-action analogy promoted to its own section.** In the original this was a single paragraph; it is the strongest single defense of the descriptive framing and deserves the weight.
- **§9 Cosmological evolution expanded.** The original punted on how $s$ evolves. The revised section names the slice-reading vs. dynamical-emergent-reading split and pins the resolution to Paper 6.
- **§11 Caveats added.** The modal commitment, the physical-vs-logical distinction on admissibility, and the equilibrium-vs-dynamics strain are now acknowledged rather than glossed.
- **§12 Implementation tiered.** Foundations go fully descriptive; derivation chapters keep active-voice scaffolding. Trying to purge active voice from proof narratives hurts clarity more than it helps rigor.
- **§12 Paper 13 §2 replacement.** "Maintained distinction" changed to "present distinction" — "maintained" is itself an active verb.
- **§1 framing.** The lede now says "prose alignment with existing formal content, not new physics." Honesty about what the document is.
