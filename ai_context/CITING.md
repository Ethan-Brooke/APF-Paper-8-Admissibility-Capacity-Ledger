# CITING.md — How to Cite This Paper, the Codebase, and the Framework

This release is one paper of a multi-paper series. There are three citation targets, depending on what you're citing.

---

## 1. Cite this specific paper

**APA-style:**

> Brooke, E. S. (2026). *The Admissibility-Capacity Ledger*. Version v2.9. Zenodo. https://doi.org/10.5281/zenodo.19721384

**BibTeX:**

```bibtex
@misc{brooke_paper8_2026,
  author       = {Brooke, Ethan S.},
  title        = {The Admissibility-Capacity Ledger},
  year         = {2026},
  version      = {v2.9},
  publisher    = {Zenodo},
  doi          = {10.5281/zenodo.19721384},
  url          = {https://doi.org/10.5281/zenodo.19721384}
}
```

The DOI `10.5281/zenodo.19721384` is the concept DOI; it resolves to the latest version. Specific versions have their own version-DOIs (see the Zenodo record).

If this release includes a Technical Supplement, cite the supplement inline where the formal proof is referenced:

> See Technical Supplement (Brooke, 2026, §S.N) for the full proof.

## 2. Cite the codebase (the executable theorem bank)

The codebase has its own Zenodo deposit — cite this if you're citing the machine-verification aspect or a specific check function.

**APA-style:**

> Brooke, E. S. (2026). *Admissibility Physics: Unified Theorem Bank & Verification Engine*. Version 6.9. Zenodo. https://doi.org/10.5281/zenodo.18529115

**BibTeX:**

```bibtex
@misc{brooke_apf_engine_2026,
  author       = {Brooke, Ethan S.},
  title        = {Admissibility Physics: Unified Theorem Bank and Verification Engine},
  year         = {2026},
  version      = {6.9},
  publisher    = {Zenodo},
  doi          = {10.5281/zenodo.18529115},
  url          = {https://doi.org/10.5281/zenodo.18529115}
}
```

Cite a specific check function as: `check_T7B` (from `apf/gravity.py`, APF codebase v6.9, Brooke 2026).

## 3. Cite the framework (the overall APF program)

For a general reference to the Admissibility Physics Framework as a whole, cite Paper 0 (the ontology paper) and the engine paper together:

- Brooke, E. S. (2026). *What Physics Permits: A Constraint-First Framework for Physics*. Zenodo. https://doi.org/10.5281/zenodo.18439523
- Brooke, E. S. (2026). *Admissibility Physics: Unified Theorem Bank & Verification Engine*. Zenodo. https://doi.org/10.5281/zenodo.18529115

For specific framework results, cite the relevant paper from the series:

| Paper | Topic | DOI |
|---|---|---|
| 0 | Ontology (*What Physics Permits*) | 10.5281/zenodo.18439523 |
| 1 | Enforceability of distinction (spine) | 10.5281/zenodo.18439200 |
| 2 | Structure and gauge uniqueness | 10.5281/zenodo.18439274 |
| 3 | Ledgers: entropy, time, accumulated cost | 10.5281/zenodo.18439363 |
| 4 | Constraints: field content, fermion spectrum | 10.5281/zenodo.18439397 |
| 5 | Quantum structure (Hilbert, Born) | 10.5281/zenodo.18439433 |
| 6 | Dynamics and geometry (spacetime, gravity) | 10.5281/zenodo.18439445 |
| 7 | Action, internalization, the Lagrangian | 10.5281/zenodo.18439513 |
| 13 | The minimal admissibility core | 10.5281/zenodo.18361446 |

## 4. Standalone derivation papers

Some specific quantitative results are published as standalone derivation papers, not part of the main spine:

- Brooke, E. S. (2026). *The Weak Mixing Angle as a Capacity Equilibrium*. Zenodo. https://doi.org/10.5281/zenodo.18603208

Cite these directly if your work specifically leans on their content (e.g., $\sin^2\theta_W = 3/13$ derivation).

## 5. When to cite which

**Citing a specific derivation or numerical result** (e.g., "$\Omega_\Lambda = 42/61$ as derived by APF"): cite the paper that derives it.

**Citing APF as a framework** (e.g., "We adopt the admissibility framework [ref] for ..."): cite Paper 0 + the engine paper together.

**Citing machine-verifiability** (e.g., "We verified this against the machine-checked codebase [ref]"): cite the codebase Zenodo DOI.

**Citing a specific check function** (e.g., "verified by `check_T7B`"): cite inline via the codebase DOI, with the check-function name.

**Citing the overall program** (e.g., "The APF program [ref]"): cite the `orcid` and the engine paper.

## 6. Author and ORCID

- **Author:** Ethan S. Brooke
- **ORCID:** https://orcid.org/0009-0001-2261-4682
- **Affiliation:** Independent researcher (Admissible Technologies)
- **Contact:** brooke.ethan@gmail.com
- **GitHub:** https://github.com/Ethan-Brooke

## 7. Zenodo version handling

Each paper's Zenodo deposit has a **concept DOI** (always points at the latest version) and one **version DOI** per release. When citing for strict reproducibility (peer-reviewed work, experimental confrontation), use the version DOI from the specific Zenodo record you downloaded. For general reference, the concept DOI is adequate.

The `zenodo.json` file in this release embeds the concept DOI as a `related_identifier` of type `isVersionOf`. When a new GitHub release + Zenodo deposit is minted, the new version-DOI is added under the concept.

---

*Generated by `create-repo`. To edit citation conventions, modify `templates/CITING.md.template` in the skill folder.*
