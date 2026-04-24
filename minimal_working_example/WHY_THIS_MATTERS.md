# Why the Minimal Working Example Matters

A common failure mode when a reviewer (human or AI) encounters an abstract framework: they cannot tell whether the framework is *well-defined machinery* or *clever-sounding notation that hides the hand-waving*. The minimal working example is Paper 8's answer to that failure mode.

The toy interface at $K = 3$, $d_{\rm eff} = 4$ instantiates every piece of the ACC ledger:

- The ACC record $(K, \mathrm{ACC})$ as a literal 2-tuple.
- The six regime projections $\pi_T, \pi_G, \pi_Q, \pi_F, \pi_C, \pi_A$ as literal numbers (or "undefined" where applicable).
- The microstate space $H_{\rm micro} = (\mathbb{C}^4)^{\otimes 3} = \mathbb{C}^{64}$ as a literal Hilbert space.
- The bridge identity $I_2$ at the integer level: $\pi_F = 3 = $ denominator of $\pi_C$.
- The vacuum-expectation identity: $\langle P_{\rm vac} \rangle = 1/4 = C_{\rm vac}/d_{\rm eff}$.
- The partition-function limit: $Z(0) = N$.

If you understand the toy, you understand the shape. The SM-scale argument (Paper 8 Supplement §1 Theorem 1.1, plus everything in §A through §L) keeps the shape and fills in the specific $V_{61}$, $G_{\rm SM}$, and $T_{12}$ content. Whether those specific structural identifications are forced — that is the real argument, defended in Theorem 1.1 itself.

But if a reviewer rejects the shape, no specific identification can help. And a reviewer who sees only the SM-scale paper has no way to separate "this machinery is well-defined but applied to a contestable object" from "this machinery is wishful thinking".

The toy resolves that ambiguity by showing the machinery stands up on an object no one can contest — $\mathbb{C}^3$ with trivial gauge, uniform residual partition. Run it. It works.

## Onboarding sequence

If you are onboarding this repo for the first time:

1. **Run the MWE** (30 seconds). Confirm all identities pass.
2. **Read `ai_context/ARGUMENT_FLOW.md`** (2 minutes). The one-page argument spine.
3. **Read `ai_context/LOCAL_VS_IMPORTED.md`** (3 minutes). Know the boundary.
4. **Read Supplement §1** (Theorem 1.1, the formal kernel). 5 pages.
5. **Read `ai_context/CLAIMS_LEDGER.md`** (5 minutes). Row-by-row attack surface.
6. **Read `ai_context/DO_NOT_CLAIM.md`** before writing any summary.
7. **Only then** tackle the full supplement.

Total time to an honest working understanding: roughly 30 minutes.

## How this changes the conversation

A reviewer who has run the MWE can no longer say "this framework is not well-defined at any scale". The machinery demonstrably works at K = 3. The conversation shifts to the right question: "is the K = 61 identification forced?". That is the question Paper 8 Supplement §1 is designed to answer.

A reviewer who has NOT run the MWE may conflate these two distinct concerns. The purpose of this directory is to prevent that conflation.
