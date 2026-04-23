---
type: concept
domain: apf
layer: 5-quantum
created: 2026-04-14
updated: 2026-04-14
sources: []
---

# Tensor Product Structure
## "How Multi-Particle Hilbert Spaces Arise from Admissibility"

**Core claim:** Multi-particle quantum systems naturally decompose into tensor products $\mathcal{H} = \mathcal{H}_1 \otimes \mathcal{H}_2 \otimes \cdots$ because this is the unique cost-minimizing way to extend single-particle admissibility to composite systems.

**Mathematical structure:**

For $n$ distinguishable particles, the joint Hilbert space is:
$$\mathcal{H}_n = \mathcal{H}_{\text{particle 1}} \otimes \mathcal{H}_{\text{particle 2}} \otimes \cdots \otimes \mathcal{H}_{\text{particle n}}$$

**Why is this forced?**

From [[Axiom A1]], maintaining distinctions in a composite system should not cost exponentially more than in single-particle systems. The tensor product is the unique structure satisfying this:
- Cost to maintain $n$ distinguishable states in $\mathcal{H}$ scales as $O(\log d)$ where $d = \dim \mathcal{H}$ is finite (L_epsilon_star)
- For composite systems, this forces $d_n = d_1 \cdot d_2 \cdots d_n$ (product structure)
- Any other composition rule would either break finiteness or allow infinite cost growth

**Entanglement emerges naturally:**

Separable states $|\psi_1\rangle \otimes |\psi_2\rangle$ are easy to maintain (cost ∝ cost of parts).

Entangled states $\frac{1}{\sqrt{2}}(|00\rangle + |11\rangle)$ are slightly cheaper (fewer classical bits needed for measurement outcome correlation).

This is why entanglement persists: it's cost-efficient.

**Bosons vs. Fermions:**

The tensor product structure must be symmetrized for identical particles:
- Bosons: symmetric states $|\psi_A \psi_B\rangle = |\psi_B \psi_A\rangle$ (cheaper to maintain identical properties)
- Fermions: antisymmetric states (cost penalty for identical particles in same state, enforcing Pauli exclusion)

This distinction emerges from cost minimization in distinguishing identical vs. distinct particles.

**Bell states and non-locality:**

The tensor product structure naturally accommodates Einstein-Podolsky-Rosen (EPR) correlations:
- Bell state $|\Phi^+\rangle = \frac{1}{\sqrt{2}}(|00\rangle + |11\rangle) \in \mathcal{H}_A \otimes \mathcal{H}_B$
- Correlations exist even though subsystems are spatially separated
- T_Tsirelson (`check_T_Tsirelson`) proves quantum correlations cannot violate $S \leq 2\sqrt{2}$ (Tsirelson bound)

**Open questions:**

1. **Does locality emerge from tensor product structure?**
   - Current: Locality assumed; tensor product structure enforced by cost
   - Open: Can we derive locality (non-signaling) from admissibility alone?

2. **Can we reconstruct the metric from entanglement?**
   - Ryu-Takayanagi conjecture connects entanglement entropy to spacetime geometry
   - Does this emerge naturally from admissibility cost?
   - Relates to RP1 (Riemannian geometry derivation)

3. **What is the fundamental reason for the factor of 2 in $d_1 \cdot d_2$ scaling?**
   - In classical information, product states have multiplicative dimensions
   - Quantum tensor products preserve this; is there a cost reason?

**Computational note:**

Computing expectation values in composite systems requires careful tensor contraction. The codebase (particularly `toy_model.py`) uses efficient contraction algorithms to avoid exponential blowup in cost.

**Analogy to classical probability:**

Joint probability distributions obey $P(A,B) = P(A)P(B)$ for independent systems. Quantum tensor products are analogous, but preserve superposition structure. This is more efficient than classical composition.

## See also
- [[Born Rule]]
- [[Paper 5 - Quantum]]
- [[Enforcement Budget]]
