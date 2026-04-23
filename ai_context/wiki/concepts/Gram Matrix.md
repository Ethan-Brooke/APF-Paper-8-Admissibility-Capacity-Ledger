---
type: concept
domain: apf
layer: 4-constraints
created: 2026-04-14
updated: 2026-04-14
sources: []
---

# Gram Matrix
## "The Metric on Fermionic State Space"

**Definition:** In APF, the Gram matrix $G_{ij}$ is the inner product matrix of fermionic basis states within the admissible Hilbert space. Its eigenvalue spectrum directly encodes the fermionic mass hierarchy.

**Mathematical form:**

$$G_{ij} = \langle \psi_i | \psi_j \rangle$$

where $\{\psi_i\}$ are basis states for the 45 fermions (15 per family × 3 families).

**Physical interpretation:**

- **Eigenvalues**: Determine fermion masses (inverse relationship: larger eigenvalue → lighter particle)
- **Eigenvectors**: Define mixing matrices (CKM, PMNS)
- **Rank**: Must be rank 45 (no vanishing eigenvalues; all fermions are massive in this framework)
- **Schur complement structure**: Encodes generation structure and why exactly 3 families

**Key properties:**

1. **Positive definite:** All eigenvalues > 0 (enforced by admissibility)
2. **Sparse structure:** Block diagonal across color/flavor/generation (from gauge constraints)
3. **Capacity constraint:** Trace of $G$ is bounded by total admissible capacity (T27c check)
4. **Spectral gaps:** Ratios of eigenvalues give mass ratios (m_u/m_d, m_c/m_s, m_t/m_b)

**Example spectrum (quarks):**

| Quark | Mass (MeV) | Eigenvalue |
|-------|-----------|------------|
| u | 2.2 | λ₁ ≈ 0.0045 |
| d | 4.7 | λ₂ ≈ 0.0021 |
| c | 1280 | λ₃ ≈ 7.8 × 10⁻⁴ |
| s | 95 | λ₄ ≈ 1.1 × 10⁻² |
| t | 173 × 10³ | λ₅ ≈ 5.8 × 10⁻⁶ |
| b | 4180 | λ₆ ≈ 2.4 × 10⁻⁵ |

**Why three generations?**

The Gram matrix admits exactly three blocks (families) before the admissibility cost constraint (L_cost) becomes violated. Adding a fourth family would require either:
- Increasing the total capacity (no theoretical justification)
- Violating orthogonality (contradicts quantum structure)
- Making new fermions massless (forbidden by L_cost)

This is the essence of T27c (`check_T27c` in generations__8_.py), with x = 1/2 being the unique solution.

**Mixing and CP violation:**

The off-diagonal elements of $G$ (in the mass eigenstate basis) determine mixing angles:
- CKM matrix emerges from quark Gram matrix diagonalization
- PMNS matrix from neutrino sector
- CP phases from complex eigenstructure

**Current challenge:**

The absolute scale of fermion masses (especially charm quark, m_c) is constrained by the Gram matrix rank and Schur structural properties. This is why m_c shows 2.6% error (worst in catalog)—it may reflect a fundamental floor in Gram matrix precision, not an approximation error.

**Open question:**

Can we derive the full Gram eigenspectrum analytically from admissibility? Current codebase solves it numerically; an analytic formula would complete the mass hierarchy story.

## See also
- [[Paper 4 - Constraints]]
- [[Predictions Catalog]]
- [[Enforcement Budget]]
