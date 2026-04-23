# STARTING_PROMPTS.md — High-Value Query Templates

These prompts are known to produce useful AI work on APF material. Use them as starting points; adapt to your question.

Grouped by intent:

---

## Verification — "show me this passes"

> Verify `ACC` in this repository. Run `python run_checks.py`, then read the function in `apf/core.py`, then explain in 200 words what it proves, what its dependencies are, and what its epistemic tag means.

> Walk through the vendored `apf/core.py`. For each `check_*` function, extract: the name, the claim (summary), the epistemic tag ([P] / [P]+[I] / [P_structural] / [Pred]), and the dependencies. Output as a markdown table.

> Run `python run_checks.py` and interpret the result. If any check fails, identify which and explain why — don't paper over it.

---

## Walkthrough — "explain the derivation chain to me"

> Using `ai_context/FRAMEWORK_OVERVIEW.md` and `ai_context/theorems.json`, trace the chain from Axiom A1 to the Born rule step by step. For each step, cite the check function that establishes it.

> What's the shortest derivation chain from A1 to `ACC`? List the intermediate theorems with their epistemic tags.

> Read `Paper_8_Admissibility_Capacity_Ledger_v2.9.tex` and produce a 500-word summary. Preserve all `\coderef{}{}` anchors. Distinguish internally-derived claims ([P]) from imported closures ([P]+[I]) and from structural identifications ([P_structural]).

---

## Audit / critique — "does this hold up?"

> I'm going to propose a mechanism for resolving [topic]. Before I describe it, review `ai_context/AUDIT_DISCIPLINE.md` and the §5 worked example (Grok H0 proposal). Then when I paste my proposal, check specifically: (a) does it contradict any [P] theorem in `ai_context/THEOREMS.md`? (b) does it introduce free parameters presented as zero-parameter? (c) does it move in a direction consistent with adjacent data?

> Here is an external review of Paper 8: [paste]. Read `ai_context/AUDIT_DISCIPLINE.md`. Verify each criticism against the current `.tex` source (use Grep / Read). Distinguish criticisms that are correct from those based on an older version. Propose full-scope remediation, not a targeted pass.

> Audit my proposal below against the framework's known constraints. Specifically check: [1] does it require demoting any check currently tagged [P]; [2] does it contradict `ai_context/OPEN_PROBLEMS.md` for known-open problems; [3] does it rely on a scaling law that is asserted rather than derived.

---

## Cross-referencing — "how does this relate to external work?"

> Position `ACC` against the mainstream comparator(s). Distinguish APF-internal derivation from reinterpretation. If the claim is an [I] import, name the imported theorem and its canonical reference.

> The Technical Supplement uses [Solèr / Gleason / HKM / Malament / Lovelock / other]. For each import, explain what APF-internal preparation makes the import applicable, and what the import closes. Use the language of "prepared class" vs "closure theorem."

> Compare APF's $\Omega_m = 19/61$ prediction to Planck 2018's $\Omega_m h^2 = 0.143$. Derive $H_0^{\mathrm{APF}}$. Compute the tension with the H0DN 2026 measurement. (Expected answer: $H_0^{\mathrm{APF}} = 67.76$ km/s/Mpc, 7.09σ tension.)

---

## Cosmology / open problems — don't re-solve the known-open

> Before proposing any dark-energy mechanism, read `ai_context/OPEN_PROBLEMS.md` §H0-tension. The five escape routes and their verdicts are catalogued there. If you propose a new mechanism, it must (a) not match the verdict on any already-ruled-out route, or (b) explicitly argue why a ruled-out route should be reopened.

> Read `ai_context/AUDIT_DISCIPLINE.md` §8 (the Grok H0 worked example). Then critique the following proposal: [paste]. Check specifically whether it repeats any of the six failure modes catalogued there.

> The APF four-parameter Planck cosmological match ($\Omega_b = 3/61$, $\Omega_c = 16/61$, $\Omega_m = 19/61$, $\Omega_\Lambda = 42/61$) is one of the framework's strongest empirical results. Explain why, at 0 free parameters, a four-parameter match within 1.2% is nontrivial. What is it evidence *for* and what is it evidence *against*?

---

## Writing / editing — "help me improve this paper"

> Read the paper's §N and identify every claim. For each: is it tagged correctly ([P] / [P]+[I] / [P_structural] / [Pred])? If the tag is wrong or missing, flag it. Do not change tags without a reason.

> Paper 8 cites [external work]. Check `ai_context/CITING.md` for the project's citation conventions. If the cited external work is a [I]-import (closure theorem), ensure the paper describes it as closing an APF-prepared admissible class rather than as a stand-alone result.

> Propose a Related Work section for this paper. Do NOT use a "reference list only" structure. Each comparator gets a paragraph positioning the paper's claim against the comparator's result. Distinguish APF-internal derivations from reinterpretations.

---

## Code extension — "add a new check"

> I want to add a new check to this paper's subset. Read `apf/core.py` for style. The new check should: (a) be a top-level function named `check_<name>`, (b) use `check()` and `CheckFailure` from `apf.apf_utils`, (c) return `_result(name=..., tier=..., epistemic=..., summary=..., key_result=..., dependencies=[...], cross_refs=[...], artifacts={...})`, (d) have an executable witness (not just assertions), (e) be added to `_CHECKS` at the bottom of the module and picked up by `register(registry)`.

> Is my proposed check consistent with the APF style guide? Specific things to verify: exact arithmetic where possible (use `fractions.Fraction`), no hidden free parameters, epistemic tag matches the actual derivation structure, dependencies list is complete.

---

## Meta — "what is this project, and what should I know?"

> I am an AI agent that just loaded this repository. Summarize in 10 bullets: what APF is, what this paper contributes, what's still open, what I should not try to re-solve, where the code lives, where the canonical truth lives, how to verify anything, what the epistemic tag system means, how to engage with external proposals honestly, and where to look first.

> Before I start substantive work on this paper, give me a 2-minute framework orientation. Use `ai_context/FRAMEWORK_OVERVIEW.md` as your source. End with the current open problems I should be aware of (from `ai_context/OPEN_PROBLEMS.md`).

---

## Anti-patterns — prompts that don't work

These look useful but produce low-quality output. Avoid.

- "Fix the H0 tension." — There is no internal fix; see `OPEN_PROBLEMS.md`. The right prompt is to state a proposed mechanism and have it audited.
- "Summarize APF in 50 words." — Compression below 200 words loses the [P]/[I] distinction that is load-bearing. Use `FRAMEWORK_OVERVIEW.md` §1 through §3 as a ~500-word summary instead.
- "Improve this paper." — Too vague. Constrain by section, claim, or audit concern.
- "Is APF correct?" — Not well-posed. Ask instead: is check X's witness correct; is the chain from A1 to Y gapless; is claim Z consistent with [data source].

---

*This file is generated by the `create-repo` skill. Editing it in a specific release won't persist across rebuilds; edit `templates/STARTING_PROMPTS.md.template` in the skill folder instead.*
