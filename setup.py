from setuptools import setup, find_packages

setup(
    name="apf",
    version="6.9.0",
    description="Admissibility Physics Framework — machine-verifiable theorem bank (v6.9 PLEC formalization)",
    author="E.S. Brooke / Admissible Technologies",
    packages=find_packages(),
    python_requires=">=3.8",
    install_requires=[
        "numpy>=1.20",
        "scipy>=1.7",
    ],
    # v6.9 (2026-04-20 interface-sector bridge pass + 2026-04-20 T_ACC
    #       unification pass + 2026-04-18 PLEC formalization):
    #   - 20 bank-registered modules (apf/plec.py, apf/unification.py both new);
    #     24 verify_all modules (includes apf/standalone/ sub-package);
    #   - 347 bank-registered theorems; 362 total verify_all checks;
    #   - 48 quantitative predictions, 0 free parameters;
    #   - 7 PLEC-infrastructure checks (2026-04-18): Regime_R,
    #     Regime_exit_Type_I..V, A9_closure (unified Lovelock-prerequisite chain);
    #   - 5 T_ACC-unification checks (2026-04-20 afternoon) in apf/unification.py:
    #     I1_holographic, I2_gauge_cosmological, I3_thermo_quantum,
    #     I4_action_thermo, T_ACC_unification;
    #   - 2 interface-sector bridge checks (2026-04-20 late evening) in
    #     apf/gravity.py: L_global_interface_is_horizon (tier 3 aux lemma),
    #     T_interface_sector_bridge (tier 4 theorem identifying T12's V_global
    #     stratum with T_horizon_reciprocity's Sector B — promotes the two
    #     "42"s in T11 and L_self_exclusion to a witnessed subspace identity,
    #     strictly finer than T_ACC identity I2);
    #   - depends on numpy + scipy.
)
