# AUDIT_DISCIPLINE.md — How to Engage Honestly on APF

This file exists because AI assistants have a recurring failure mode on this project: declaring issues resolved that aren't, accepting phenomenological proposals that silently rewrite load-bearing theorems, and conflating terms that the framework treats as distinct.

Read this before:

- Proposing any substantive edit to a paper, supplement, or code module.
- Responding to an external audit, reviewer comment, or lit review.
- Authoring content that claims to close, resolve, or discharge an issue.
- Evaluating an external AI-authored proposal (another model suggesting an APF mechanism, a derivation, a fix).

---

## The core discipline — three items, in order

### 1. Read carefully before responding

Read the whole audit/proposal/review. Pay attention to specific word choices. Distinct terms are usually distinct claims. If an external document distinguishes "bibliography" from "literature review," or "derivation" from "reinterpretation," treat them as different requests. See §2 below for the specific word-distinctions APF treats as load-bearing.

### 2. Verify against the current document, not your memory

For each criticism, `grep` or `read` the actual `.tex`, `.md`, or `.py` to check whether the concern is real. If the criticism is that something is absent, confirm it is actually absent in the current state. If Claude's earlier session claimed something was fixed, verify the fix actually landed — there is a history on this project of "resolved" claims not matching document state.

### 3. Name the gap before proposing the fix

If the criticism is correct, say so plainly. Do not segue directly into "here's the plan." First acknowledge what was missed, specifically; then propose scope. If a prior pass shipped a shallow version of the fix, name that pattern explicitly.

---

## Word distinctions that matter on APF

These pairs have been conflated in the past. Each pair is two different things.

| Looks similar | But means | And requires |
|---|---|---|
| **Bibliography** | List of references at the back | `\bibliography{}` entries |
| **Literature review** | Substantive section positioning the paper in the field | A dedicated `Related Work` section |
| **Citation** (`\cite{Foo}`) | Acknowledging prior work exists | A reference in `\bibliography` |
| **Positioning** | Explaining how this paper relates to that prior work | Prose comparison at the claim site |
| **Derivation** | Obtaining $X$ from APF axioms | A proof chain ending at $[P]$ |
| **Reinterpretation** | Relabeling a known $X$ in APF vocabulary | A `[P_structural]` or `[I]` tag |
| **Reference list** | Minimum for citing | Entries in `.bib` file |
| **Related work section** | Full positioning of the paper | Dedicated section, typically after intro |
| **Acknowledgment of deferral** | "This is deferred to Paper N" | A line of prose |
| **Delivery of substance** | Actually providing a sketch, bound, or partial result | Derivation content, with code |
| **Closed** (as in theorem status) | All prerequisites derived and unified | Formal check function, [P] tag |
| **Stated** | Written as a claim | No code anchor required |

## Epistemic tags — what each one means

All four tags appear in APF papers. Confusing them is a specific failure mode.

- **[P]** — internally derived. The claim follows from APF axioms (A1 + MD + BW + PLEC) using only APF-internal reasoning. No external classification or uniqueness theorem is invoked.

- **[P] + [I]** — derived with import. APF prepares an admissible class; an external mathematical theorem (Solèr, Gleason, HKM, Malament, Lovelock, etc.) closes the classification or uniqueness. The framework is responsible for preparing the class; the import closes it.

- **[P_structural]** — structural identification. The framework identifies an already-derived structure with a familiar physical representation or interpretive role. Weaker than a full derivation; often a definition-level result once the prerequisites are in place (e.g., T7B identifying the metric tensor once quadratic non-factorizing response is established).

- **[Pred]** — phenomenological prediction. An empirical claim whose status is downstream of the prior layers and must be evaluated against data. Tag for scaling laws, sector identifications (e.g., dark matter as $C_{\mathrm{ext}}$), and exclusion claims.

**Rule:** when a paper mentions a claim, check the tag. Do not treat a [P_structural] result as equivalent to a [P] result in argumentative weight. Do not cite an [I] import as if APF proved it. Do not ship a [Pred] claim as if it were a closure theorem.

## Response posture — two valid options, and how to choose

### Accept

The audit surfaces real gaps. Say so. Name the gaps specifically. Propose honest scope. If in doubt, propose MORE scope than LESS — shipping the shallow version of a fix creates the next audit.

### Push back

The audit makes a claim you can demonstrate is already addressed. Pushback is allowed, but requires:

- **Specific section citation** (§X.Y, tcolorbox label, line range, code module + function name).
- **Specific language from the current document** that addresses the criticism.
- **Acknowledgment** that the audit's concern was reasonable and the reviewer may have missed the addressing language.

Pushback without specific citations is not pushback. It is brushing off. The project has a history of this specific failure; don't repeat it.

## Scope discipline

When proposing remediation, the default is **full scope**, not minimum scope.

A "targeted citation pass" is almost never the right answer to a lit-review critique. A proper lit-review response typically requires:

- A dedicated `Related Work` or `Relation to prior work` section.
- Substantive positioning of each major claim against its mainstream comparators (not just a citation).
- Explicit acknowledgment of reinterpretation vs. new mathematics.
- Reference list expansion to include the audit-named gaps.
- Mirror positioning into any Technical Supplement.

If a lighter scope is genuinely right, make the case explicitly: name which audit items are duplicative of existing content, cite the existing content by section and line, and note the items where substantive work is still needed.

## Second audits — not redundant

If a second review raises concerns similar to the first, the default instinct is "this confirms the fix was right." Resist this. Look specifically for:

- **New literature references** the first reviewer did not cite.
- **Terminological distinctions** the first reviewer did not make.
- **Stronger derivational demands** on claims both reviewers flag (e.g., if both flag the same numerical result, the convergence is signal).
- **Empirical anchors** the first reviewer did not name.

Convergence across reviewers is a stronger signal than a single reviewer. Treat a second audit as potentially surfacing new content, not as confirming the first fix.

## Worked example — the Grok H0 proposal (April 2026)

This is the specific failure mode the discipline exists to prevent. Understanding it is the quickest way to internalize the posture.

**The proposal:** an external AI proposed a "time-locked vacuum" mechanism resolving the H0DN 7.09σ tension. The mechanism substituted the partition denominator ($C_{\mathrm{total}} = 61$, a type count) with $C_{\mathrm{active}}(a) = C_{\mathrm{total}} - L(a)$, a time-dependent capacity amount driven by L_irr locking. The resulting effective $\Omega_\Lambda(a)$ drifts upward, giving phantom $w(z) \approx -1.5$ near $z = 0$, yielding local $H_0 = 73.5$.

**Surface plausibility:** the proposal was well-structured, used APF vocabulary correctly (L_irr, capacity, admissibility), arrived at clean closed-form expressions, and presented itself as "zero free parameters, zero new axioms."

**What the discipline caught:**

1. **False zero-parameter claim.** The proposal has two free parameters (β and $L_{\max}$), both varied in its own tables. $L_{\max} = 9.71$ has no APF axiomatic origin shown. Picking β = 4 and $L_{\max} = 9.71$ is a fit, not a derivation.

2. **Wrong sign vs DESI.** DESI DR2 prefers $w_0 \approx -0.75$ (quintessence-like). The proposal predicts $w(0) \approx -1.5$ (phantom-like). These point in opposite directions from $w = -1$. The proposal would "resolve" H0 by worsening the DESI fit.

3. **Implicit demotion of `L_saturation_partition` [P].** The load-bearing move is substituting the type-count denominator with a capacity amount. `L_saturation_partition` explicitly proves the partition is type-count topological, not capacity-dependent: the same 61 types persist regardless of how much capacity is locked in records. The proposal requires this [P] theorem to fail, but doesn't say so.

4. **Asserted not derived scaling.** The proposal writes "$\rho_{\mathrm{DE}} \propto \Omega_\Lambda^2$ (derived from continuity equation)" but does not show the derivation. The entire analytic result depends on this scaling.

5. **Ad hoc functional form.** $L(a) = L_{\max}(1 - a^{-\beta})$ is chosen for analytic tractability, not derived from the dynamics of irreversible record formation.

6. **No end-to-end numerical demonstration.** The proposal asserts β = 4, $L_{\max} = 9.71$ produces $H_0 = 73.5$ but doesn't solve Friedmann with $\Omega_\Lambda(a)$ as defined and recover the number.

**Disposition:** reject the specific construction. Preserve as an artifact (the question it points at — H0 tension at 7.09σ — is real). Formalize the honest APF response: APF predicts 67.76 km/s/Mpc with zero free parameters and four-parameter Planck match; the tension is a falsification candidate, not a bug to patch around.

**The pattern to learn:** the proposal's fluency in APF vocabulary did not make its conclusion correct. The discipline requires checking whether the load-bearing step is a derivation or a fit, whether it contradicts any [P] theorem, and whether it's consistent with adjacent data (DESI). Fluency without this discipline ships harm.

## Important

- Do NOT claim an audit concern is "fixed" without verifying with `grep` or `read` against the current document state.
- Do NOT treat a second audit as redundant until it has been read in full and its specific novel claims enumerated.
- Do NOT segue from criticism straight to "here's the plan." Acknowledge the gap first, then propose the plan.
- Do NOT propose a "quick" or "light" scope if the audit is asking for structural work.
- Do NOT accept a phenomenological proposal that requires a [P] theorem to fail without acknowledging the demotion and making the case for it.
- DO name the pattern when repeating it.
- DO distinguish reviewer-version issues from current-version issues.
- DO propose more scope rather than less when in doubt.
- DO ask for scope confirmation before executing on any remediation larger than a few targeted edits.
- DO consult `OPEN_PROBLEMS.md` before proposing new work — many "open questions" have known-open research programs; re-solving them naively wastes effort.

---

*Canonical source: `__APF Library/skills/audit-response/SKILL.md`. This file is a repo-facing distillation; the full skill has more process detail for in-session work.*
