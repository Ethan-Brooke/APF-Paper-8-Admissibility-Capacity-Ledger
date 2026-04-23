---
type: concept
domain: apf
layer: 1-spine
created: 2026-04-14
updated: 2026-04-17
sources: []
---

# Enforcement Budget
## "The Cost of Maintaining Physical Distinctions"

**Definition:** The total operational cost required to maintain all physical distinctions in a system over a given timeframe. This is the core resource whose finiteness constrains all of physics (per [[Axiom A1]], the upper-bound component of [[PLEC Rollout Plan|PLEC]]).

**PLEC note (2026-04-17):** The enforcement budget $C$ is the quantity bounded above by A1, bounded below per distinction by MD (the positive cost floor $\varepsilon^* > 0$), minimized by A2 (argmin selection picks the cheapest admissible configuration), and kept non-degenerate by BW (budget-window ensures distinct-cost alternatives exist). The scaling laws below continue to hold without modification; the PLEC framing simply names which component of the variational principle each constraint represents.

**Formal structure:**

Let $\mathcal{D}$ be a set of distinguishable states. The enforcement cost $C(\mathcal{D})$ is the work required to:
- Prepare the ensemble of $\mathcal{D}$
- Maintain the distinctions against decoherence/perturbation
- Measure and verify the distinctions

**Key scaling laws (from `core.py`):**

- **L_epsilon_star:** $C \geq \epsilon_* \cdot \log N$ where $N = |\mathcal{D}|$ (minimum cost grows logarithmically)
- **L_cost:** For mixed ensembles, $C_{\text{mixed}} \geq \sum_i p_i C_i$ (convexity)
- **L_irr:** Reverse process cost $C_{\text{reverse}} > C_{\text{forward}}$ for thermodynamic processes (asymmetry)

**Physical consequences:**

1. **Hilbert space dimensionality** — Can't maintain infinite superpositions cheaply; forces finite-dimensional spaces
2. **Quantum indeterminism** — The [[Born Rule]] is the cost-optimal measurement rule
3. **Mass hierarchies** — Heavier particles cost more to maintain; drives fermionic mass ratios
4. **Dark energy** — Vacuum state enforcement cost sets the cosmological constant
5. **Irreversibility** — Entropy increase follows from cost asymmetry between forward/reverse processes

**Analogy:** Like a thermal energy budget in thermodynamics, but more fundamental—applies even at zero temperature, to pure quantum states.

**Measurement perspective:**
- Getting a definite measurement outcome requires "paying" in enforcement cost
- Superposition is cost-efficient (many states, low cost)
- Collapse is a cost spike (one state, high enforcement cost)

**Status:** Fully formalized in `core.py` and `supplements.py` (121 checks total)

## See also
- [[Axiom A1]]
- [[Paper 1 - Spine]]
- [[Non-Closure Theorem]]
