# Minimal Working Example

A 10-line numpy reproduction of Paper 8's core ledger identities, at the smallest possible structural scale.

## Run it

```bash
python3 toy_interface_numpy.py
```

Compare the output line-by-line against `toy_interface_expected_output.txt`. If the files match, the ledger machinery is working on your machine.

## What it proves

At $K = 3$, $d_{\rm eff} = 4$ (the smallest non-trivial interface):

- $\pi_F = 3$, $\pi_Q = 64 = d_{\rm eff}^K$
- $\pi_C = (1/3, 1/3, 1/3)$ — three equal residual sectors
- $Z(0) = 64 = N$ — partition-function limit at infinite temperature
- $\langle P_{\rm vac} \rangle = 16/64 = 1/4 = C_{\rm vac}/d_{\rm eff}$ — vacuum expectation on max-mixed state
- $\pi_F$ = denominator of $\pi_C$ — the bridge identity $I_2$ at the toy scale
- Slot / entropy / cosmological fractions all collapse to $K_X/K$ — Fractional Reading Equivalence
- Toy $\rho_\Lambda = 1/256$

All verified by hand. No `scipy`. No APF codebase. Just numpy and the `fractions` module.

## Why this is in the repo

The toy interface is the paper's answer to "can I audit this argument without learning APF?". It teaches the *shape* of the ledger machinery using only standard mathematical objects — tensor products of finite-dimensional complex vector spaces, max-mixed density matrices, projector traces. The SM-scale application (K = 61, $d_{\rm eff} = 102$) uses the same shape; you don't need to trust the SM-scale identifications to verify that the shape is well-defined.

See `WHY_THIS_MATTERS.md` for the pedagogical argument.

## What it does not prove

Running this script does not prove:

- That the SM interface has $K = 61$, $d_{\rm eff} = 102$.
- That the Standard Model gauge group is $SU(3) \times SU(2) \times U(1)$.
- That $\Omega_\Lambda = 42/61$ is the right cosmological fraction.

Those are structural claims that must be defended separately. The MWE only proves that the mathematical machinery Paper 8 uses is well-defined at small scale — a necessary precondition, not the content.
