---
type: concept
domain: apf
layer: 2-structure
created: 2026-04-14
updated: 2026-04-14
sources: []
---

# Non-Closure Theorem
## "Admissible State Spaces Cannot Close Under All Symmetries"

**Theorem (L_NZ):** If a quantum state space $\mathcal{H}$ is admissible (satisfies [[Axiom A1]] with finite enforcement cost), then there cannot exist a continuous group $G$ such that every automorphism of $G$ maps $\mathcal{H}$ onto itself while preserving inner products. In other words, the state space does not form a complete representation of any classical gauge group.

**Statement (informal):** You cannot close the state space. Any attempt to add symmetries classically will violate the finiteness constraint.

**Proof sketch (`check_L_NZ` in gauge.py):**

1. Suppose $\mathcal{H}$ is closed under group $G$: all $g \in G$ yield automorphisms $U_g : \mathcal{H} \to \mathcal{H}$.

2. The cost to maintain all $g$-invariant states as a distinguishable ensemble grows with $|G|$ (or dimension of $G$ if infinite).

3. For a continuous group, this cost becomes infinite (uncountably many group elements).

4. Therefore, maintaining the closure structure costs $\infty$, violating finiteness (A1).

5. Conclusion: only local (spacetime-dependent) gauge structures avoid closure.

**Physical consequence:** [[Gauge Uniqueness]]

This forces the transition from classical global symmetries to quantum local gauge bosons. It explains why:
- We cannot have a GUT with a larger classical group
- The standard model gauge group is unique and necessary
- Additional fermions or exotics are forbidden (they would require closure)

**Comparison:**
- **Classic mistake:** Gauge theories are assumed; APF derives why they are forced
- **GUTs:** Assume a larger group then break it; APF forbids closure in the first place
- **Supersymmetry:** Would require closing under fermi-bose exchange; forbidden by non-closure

**Historical note:** The intuition that physics forbids infinite resources is old (from thermodynamics). APF makes this fully rigorous and derives gauge structure as a necessary consequence.

**Status:** Proven and formalized; central to [[Paper 2 - Structure]]

## See also
- [[Gauge Uniqueness]]
- [[Axiom A1]]
- [[Paper 2 - Structure]]
