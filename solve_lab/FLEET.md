# FLEET — multi-agent solver campaign: roster, stop policy, restart procedure

This file is written and owned by the **coordinator session**, not by any solver agent.
Solver agents must not modify it. It exists so the campaign can be stopped at a session
limit and restarted afterwards without losing work or re-deciding the plan.

## Roster (10 independent solvers)

Every agent: reads `PROMPT.txt` as its mission, re-verifies the 39,026 partial itself,
treats prior lab conclusions as hypotheses, writes only into its own directory, runs no
git commands, and verifies every claimed improvement with `solve_lab/checker.py`.

| Agent | Work dir | Angle |
|-------|----------|-------|
| A | `agentA_work/` | Exact integer linear algebra: HNF / Smith normal form, LLL / BKZ on the lift lattice; retest "39,026 is optimal" under a larger movable-variable set |
| B | `agentB_work/` | Independent re-parse of the raw file; own model, ignoring the lab's atom/gate decomposition |
| C | `agentC_work/` | SAT / SMT / CP / MIP and algebraic solver encodings of the reduced core |
| D | `agentD_work/` | Large-scale stochastic search; fast exact incremental evaluator, parallel tempering / LNS across multiple restart partials |
| E | `agentE_work/` | Structural attack on the 7 failing equations: exact dependency neighborhood, computed free deformation space, simultaneous repair |
| F | `agentF_work/` | p-adic and multi-modular lifting: solve mod many primes, Hensel-lift, CRT-reconstruct |
| G | `agentG_work/` | Relaxation and rounding: homotopy, SDP / moment, MIP relaxations as a compass, exact integer snapping |
| H | `agentH_work/` | Decomposition: dependency graph, separators, treewidth, block DP over a small separator |
| I | `agentI_work/` | Build from scratch rather than repair: complete CDCL-style search, hard variables decided at the top of the tree |
| J | `agentJ_work/` | Attack the reduced parameterization (the claimed thirteen 296-bit numbers); verify it first, then solve directly in that space |

**Forbidden for all agents:** reverse-engineering how the instance file was generated — no
PRNG or coefficient-template forensics, no emission-order or index-ordering forensics aimed
at recovering a generator. Agent E originally held that angle and was redirected; artifacts
it produced under it are withdrawn, not built upon. The equations are to be analyzed as
mathematical objects only.

## Binding limits (measured, not assumed)

| Resource | Measured | Status |
|----------|----------|--------|
| CPU | 4 cores, 10 agents | **BINDING** — ~2.5x oversubscribed |
| Memory | 15 GB total, ~12 GB available | ample |
| Disk | ~30 GB available; `solve_lab/` ~350 MB | ample |
| Coordinator context | grows with each agent report | the real session limit |

CPU is the near constraint, so agents are compute-starved rather than blocked. That is
acceptable for search-style work but means wall-clock results arrive slower than a 10x
fan-out suggests.

## Stop policy

A stop is *appropriate* only when the agent's work is already durable. Order matters:

1. Send each agent a flush request (`SendMessage`): write its findings to
   `agent<X>_work/LOG.md` and `agent<X>_work/RESUME_<X>.md`, save any verified assignment
   JSON, and report back — within one tool round.
2. Wait for the reports, or for a short grace period.
3. Only then `TaskStop` anything still running.
4. Coordinator commits and pushes all agent directories.
5. Update the "Campaign state" section below.

Never `TaskStop` first: agents hold their best results in context, and an un-flushed stop
discards exactly the work the campaign is for.

Under CPU pressure specifically, prefer **thinning the fleet** over stopping it — stop the
agents whose angle has plateaued and let the rest have the cores.

## Restart procedure (after a limit resets)

1. Read this file, then `RESUME.md` and `STATE.json`.
2. Re-verify the deliverable: `python3 solve_lab/checker.py solve_lab/best/new_instance_partial_39026.json`.
3. For each agent to resume: read its `agent<X>_work/RESUME_<X>.md`, then relaunch a fresh
   agent with the same angle from the roster above, pointing it at that file so it continues
   rather than restarts. Agent work dirs are committed, so a fresh container recovers them
   from the branch.
4. Do not relaunch an angle that reported a definitive dead end — record it under "Retired
   angles" below and spend the core elsewhere.

## Campaign state

- Best verified partial: **39,026 / 39,033**, `best/new_instance_partial_39026.json`,
  failing `[12231, 12270, 12350, 14584, 18673, 22044, 29125]`. **No agent has beaten it.**
- Retired angles: generator inversion (dropped by user instruction, agent E redirected).

### Check-in 1 — fleet thinned 10 → 5

CPU was 5x oversubscribed (load 20 on 4 cores), so five agents were stopped flush-first.
Memory and disk stayed ample. Stopped, each having delivered its angle's result:

| Agent | Why stopped | What it delivered |
|-------|-------------|-------------------|
| B | angle complete | independent parser + model, 0/39033 mismatches; every eq is `scalar * L^k = 0` |
| E | below baseline (39,015) | clean parse + acyclic DAG, 8,365 free inputs, after losing its first angle |
| G | superseded | exact symbolic forward eval over F_p: 112 symbols, 57 non-constant checks |
| H | decomposition delivered | acyclic frame over all 38,748 vars, 8,747 free inputs — cleaner than prior sessions' |
| J | verdict delivered | reduced parameterization independently re-derived and CONFIRMED real |

### WITHDRAWN BY USER INSTRUCTION — do not resume, do not re-derive

The section below records a reading of the instance that the user has directed the fleet to
drop entirely. All agents have been told to delete its artifacts, purge its vocabulary from
their notes, and not re-derive it. **It is retained here only so a future session recognizes
the direction as closed rather than rediscovering it.** Do not act on it, cite it, or use it
to argue that the instance is unsolvable.

The reason for the change is sound independently of whether the reading was correct: every
agent had begun interpreting every result through a single lens, and that lens was being used
to conclude the instance cannot be solved — which narrowed what anyone was willing to search.
Work the equations as integer polynomials.

Note what survives translation, because it is most of the durable content: the exact affine
region models, the optimality lemmas, the effective-support reduction, the exact-gauge kernel
result, the "not affine" measurement of the boolean-to-residual map, the flat 5–6 cost profile
over discrete knobs, and every checker-verified score. None of those depend on the discarded
framing.

<details>
<summary>Withdrawn (historical record only)</summary>

### The finding that reframed the campaign

Agent I derives, from its own parse with each step reproducible, that the instance
**reduces exactly to a 256-bit ECDLP**: residual symbols trace to the short-Weierstrass
addition law after removing a K offset; the 512 conditional-pin constants lie on
`y^2 = x^3 + b`, a sextic twist of secp256k1 with prime 256-bit group order N; 185 of 219
table points have their double in the table, so the 256 selectors index a doubling ladder;
the instance asserts `k*G = P_target`. If step 6 of its chain holds — that there is no
mod-p freedom anywhere except the selector bits — then 39,026 is a **coding optimum, not a
near-miss**, and no search closes the remaining 7 equations.

Independent corroboration from agents that never shared a model: D's residual at 39,017 is
two mod-p pins on `(x3, y3)`; C reduces the system to three conditions; F confirms p is the
secp256k1 field prime with handles entering as `p*h`; A states the obstruction as two mod-p
congruences over 12 rank-7 rows.

### Surviving five, all retasked to falsify that claim from their own models

Independence is the point: none may import agent I's artifacts.

- **I** — harden the chain into a PASS/FAIL certificate; be adversarial about step 6, which
  is the only step carrying "unavoidable" rather than merely "present."
- **A** — are the two mod-p congruences the same object as the ECDLP? Does any integer
  combination escape them (HNF/enlargement)? An escape refutes step 6.
- **C** — decisive and cheap: in branch (1,1), set free inputs `x_22162, x_30213` to K2, K1
  and report the checker score. Closes the system, or breaks exactly as I predicts.
- **D** — did the search ever move `(x3, y3)` to a *different* valid curve point? Answerable
  from logs already held. If it moved freely, step 6 is wrong.
- **F** — independently verify N's primality (computed once, by one agent, via Cornacchia),
  and check for smooth/small-order structure. Also test the sharp prediction that lift
  obstructions concentrate at p alone.

A clean refutation would be the campaign's most valuable result; a clean confirmation lets
the situation be stated precisely instead of gestured at.

</details>

### Current tasking — after the withdrawal

All framing-dependent searches were stopped. The live work is algebraic and combinatorial,
and most of it needed no reframing because it was never framing-dependent:

| Agent | Task |
|-------|------|
| A | ISD / coset-leader search on the mod9118_0 basin: 89 affine rows, 65 knobs, rank 65, Q-consistent with a unique non-integral solution, so every integer point's violated set must contain a code support. Report the minimum support weight observed. |
| C | Fix the settable classifier so it reproduces the deliverable's true cost of 7, then re-run globalscan over all 3,349 settable handle-definition atoms with corrected pricing |
| E | Pin-feasibility scan per bit, then subsets: do the pin conditions constrain subsets independently, or only through their sum? |
| F | Multi-modular: solve mod many primes and prime powers, locate where lift obstructions actually concentrate, then Hensel-lift and CRT-reconstruct |
| G | Minimum-weight coset decoding in the exact equation-level model (6,613 linear + 161 nonlinear equations, 4,652 unknowns) — the only model in the fleet posed at equation level rather than atom level |
| H | Integer relations among the 512 load-pin constants and the two residual congruences, using close2.py's constructive cascade closer rather than search |
| I | Maximize cancellation: min number of nonzero equations over atom vectors in the image of the atom map, seeded from a 2-atom cut |
| J | Second off-manifold coding: 2 nonzero atoms vs the deliverable's 7, choosing handle lifts into as many equation-row kernels as possible |

I and J attack the same formulation from different models and are barred from importing each
other's work; A, C and G attack adjacent versions of the placement question the same way.

### Check-in 3 — the central question is answered

The open question above — that "all atoms zero" is sufficient but not necessary, so every
optimality argument lived inside one branch — is **closed**, from two directions.

**F computed it directly.** For the 39,033 × 39,033 equation–atom incidence matrix M
(525,982 nonzeros): **rank(M) = 39,033, dim ker(M) = 0.** Two independent methods — a
characteristic-free peeling certificate re-verifiable from M on disk (all pivots ±1 or ±2,
none divisible by any odd prime, so it holds over ℤ and every field of char ≠ 2), and
Wiedemann over a word prime, validated on controls of known rank before use. Therefore *any*
assignment satisfying all 39,033 equations must make all 39,033 atoms exactly zero:
all-atoms-zero is an **equivalence, not a restriction**, and no frame-optimality result in
this lab still depends on it.

**A's equation-level bound survives, and is tight.** Asked whether it was exposed to the
vacuity J found in the raw relaxation, A checked its own formulation rather than defending
it: its code is {N·u} over knob directions, so realizability is in the construction. No
single-equation atom appears in any window, and none of the 3,235 in the instance carries a
private variable. The sharpened statement: **neither condition alone gives 7** — condition (b)
admits 582 weight-6 supports at L=6 and 62 at L=16, the mod-p filter kills every one, and 7
is exactly the lightest weight the filter admits. The bound is tight, not merely unbeaten.

**Convergent closure.** eq8680 — the single row costing the deliverable its seventh equation —
was reached independently four times: I from the equation side, J from row structure, H from
the free-input side, G from its equation-level model (where it is a binary quadratic form of
discriminant exactly 0, i.e. a perfect square pinning one linear relation for one equation).
I then closed it by complete enumeration over its 18 atoms; H proved no direction anywhere
moves a22231 without moving a37887, making it a 1-for-1 trade.

**J's synthesis of why everything failed:** the instance keeps its free content
**equation-disjoint** from its residual cluster. Zero of the 83 confined knobs, and zero of the
900 free-knob atoms, share even one equation with the deliverable's cluster. Six dead ends,
one property.

### Live tasking after check-in 3

| Agent | Task |
|-------|------|
| A | Write THEOREM.md — the full statement, both conditions, tightness, and what is *not* proved (the L=344 affine ceiling) |
| C | Residue-aware pricer (structural pricing proved a dead instrument — no function of incidence structure alone can return 7) |
| E | The irreducible triple that every ≥2-bit subset reduces to |
| F | The two remaining conditional links — the **only** things between this lab and an unconditional result |
| G | Carry the span argument outward (enumeration stops being exhaustive past ~3e6 candidates) |
| H | Price all ~30 cascade pin atoms as alternative defect carriers |
| I | Finish the 66 pairs, then stop |
| J | The a23328 bridge — the one object not equation-disjoint from the residual |

**Limits at check-in 3:** 4 cores, load 17.3, 10 GB memory free, 29 GB disk, solve_lab 511 MB.
No limit near. Framing purge audited clean across all ten agent directories.

---

## THE INSTANCE, DECODED — read this before anything else

Agent F (`agentF_work/tree96.json`, `stage_law2.py`): **`EQUATIONS.txt` is a 96-stage binary
combination tree of depth 6 over 256 leaf constants.** A gate fires only when *both* its inputs
are live; with one live input the value passes straight through. Of the 72 stages with a full
six-tuple of free inputs, **all 72 obey the same degree-3 law with the same universal offset
constant — zero exceptions**, verified across every role partition and coordinate ordering with
two independent random draws required to agree. The remaining 24 are leaf-adjacent (one input
is a literal). One stage was decoded explicitly and its demanded output matched digit for digit.

**The task is: choose a subset of the 256 leaves whose fold through the fixed tree hits the
target at the root.** The space is exponential, not quadratic.

### The reversal that produced this

The campaign had converged on a floor of 7 from six directions with an infeasibility result
forming. F then **withdrew its own infeasibility argument**: its rigidity engine derived no
contradiction on same-tree selector pairs: instead an adder stopped being forced. That is an
accumulator, not a conflict. Had F defended the claim, this lab would have closed on a false
result and nobody would have looked for the tree.

**The standing caveat (agent E, written into `agentE_work/RESUME_E.md`)** — three independent
sightings in three frames: **the mod-p content of a row is not a property of the row; it is a
property of which selectors are on.** Coefficients, coprimality *and* the target residues all
move with configuration. Therefore **every mod-p rigidity / pinning / "coefficient divisible by
p" argument in this lab is conditional on a selector configuration that was not stated when the
argument was made**, and must be re-quoted with its configuration or it is not a claim about
the instance. This applies to results throughout `NOTEBOOK.md` and to several in this file.

**Broadened after E's fourth retraction — state your KNOB SET.** Four barriers were reported
and all four were retracted, every one caused by computing a property over a *filtered* knob
set and reporting it as a property of the instance: the triple (booleans excluded), §14
(affine-only), §15 (fixed residue multiset), §20 (dual-reaching knobs only — a knob moving just
one row is precisely an independent direction). The knob set, not the configuration, is the
variable that has actually been wrong every time. Any claim of the form "nothing can move X"
must state exactly which knobs were searched.

### What still stands unconditionally

- **ker(M) = 0** — rank 39,033 on the 39,033 × 39,033 incidence matrix, by three independent
  computations (peeling certificate over ℤ, re-verifiable from M on disk; Wiedemann at two
  distinct word primes). So all-atoms-zero is an *equivalence*: any full solution must make
  every atom vanish. Note agents' matrices differ in atom count (F 39,033, I 40,885, A 42,267)
  — compare decompositions before comparing kernel dimensions.
- **The deliverable: 39,026 / 39,033**, `best/new_instance_partial_39026.json`.
- **No infeasibility claim stands anywhere in this lab.**

### Verification rule

Some states carry values above Python's 4,300-digit cap, at which point **`checker.py` cannot
parse them and raises ValueError** — a bare "checker.py says" would be false. Use
`agentE_work/verifyE.py`, which raises only the digit cap and calls checker's own
`load_equations` / `load_assignment` / `evaluate_all` unmodified. Verified by the coordinator.

### Live tasking

| Agent | Task |
|-------|------|
| F | Fold evaluator over `tree96.json`; **test whether the stage law inverts** — that single fact decides whether meet-in-the-middle makes the tree attackable |
| E | Find a second independent residue class reaching (a28647, a20215) — the 41 dual-reaching knobs form exactly *one* class, and one class cannot satisfy two independent congruences |
| H | Enumerate detach sets (1-, 2-, 3-subsets) — the last unvaried axis; carrier class, selector count, region shape and knob budget are all priced and closed |
| G | Carry the span argument outward |
| I | Finish the 66-pair sweep |

---

## FINAL POSITION (coordinator context exhausted — handoff)

**Deliverable: 39,026 / 39,033**, `best/new_instance_partial_39026.json`, failing
`[12231, 12270, 12350, 14584, 18673, 22044, 29125]`. No agent beat it.

### Established
- **The instance is a 96-stage binary tree of depth 6** over 256 leaf constants; one degree-3
  law with one universal constant across all 72 fully-determined stages, zero exceptions; the
  law is **invertible in closed form**. Reachable space = 2²⁵⁶ − 1 non-empty leaf subsets.
- **ker(M) = 0** (three independent computations) ⇒ all-atoms-zero is an *equivalence*.
- **Two independent derivations of the same problem.** F from the circuit: root slots carry
  leaf supports **178 | 78**. E from the residual congruences: channels **178 | 41 + 21 + 16**.
  Same split, no contact between them. E's caveat stands — that the 16 inert booleans become
  live at some configuration is a *prediction*, untested, so the match is suggestive not
  established.
- **Contributions saturate**: a channel contributes at most once however many of its bits are
  on — the residual-side form of "a gate passes its input through when only one input is live."

### Refuted (all by the agents that reported them)
- The **infeasibility argument** — withdrawn by its author. **No infeasibility claim stands.**
- **Five barriers**, four of them E's, every one from computing a property over a *filtered
  knob set* and reporting it as a property of the instance.
- The **"rank > deficit" criterion** — supplied by H, propagated by the coordinator, then
  refuted by H: it fires on 3,781 of 3,889 detach sets, all of which zero one row against the
  witness's five. **Integer reachability of row targets is the binding quantity, not rank.**

### Highest-value next experiment
Finish F's 56 undecoded slot pairs and the 24 leaf-adjacent literals; build the fold evaluator;
**validate on ON-set {24601, 2081}** (must predict the fold of those two leaves, *not* the
target); then invert the target down the 78-side chain through the 88- and 50-support stages to
a node of leaf support ≤ 24, enumerate forward, match under 2²⁴. The window is populated — 66 of
96 stages sit in 2..24 (`agentF_work/stage_profile.json`).

### Read these, in this order
`agentF_work/RESUME_F.md` (136 lines, self-contained) → `agentE_work/RESUME_E.md` (92 lines) →
`agentA_work/THEOREM.md` → each other `agent*_work/RESUME_*.md`.

### Two standing rules
1. **Verification**: states above 4,300 decimal digits cannot be parsed by `checker.py` — use
   `agentE_work/verifyE.py` (raises only the digit cap, calls checker's own loader/evaluator
   unmodified). A bare "checker.py says" is false for those states.
2. **Any claim of the form "nothing can move X" must state its knob set AND its selector
   configuration.** Both have been wrong, repeatedly, and both change the answer.

---

## Check-in 6 — the channel/tree match is ESTABLISHED (agent M)

Coordinator context was reset; this entry is written from agent M's flushed report and its
committed directory, not from conversation. Deliverable re-verified at the top of this
check-in: `python3 solve_lab/checker.py solve_lab/best/new_instance_partial_39026.json` →
**satisfied 39026/39033 (7 failing) `[12231, 12270, 12350, 14584, 18673, 22044, 29125]`**.
No agent has beaten it. Nine agents (K, L, N, O, P, Q, R, S, T) are still running.

### What was open

The FINAL POSITION section above records the 178|78 split reached twice — F from the circuit
(root slot leaf-supports) and E from the residual congruences (channels 178 | 41 + 21 + 16) —
with E's own caveat that the match was **suggestive, not established**, because it rested on an
untested prediction: that E's 16 *inert* booleans become live at some configuration.

### What M measured

- **E's premise had a bug.** The "16 booleans that move nothing" are **14 inert + 2 already ON**
  (`x_1530`, `x_1603`); `channels.py` skips ON bits, so those two fell out of every class. They
  sit exactly where saturation predicts — `x_1603` in stage 19538, `x_1530` in stage 10649.
  Set both to 0 and the 16 split into **node 19538 exactly** and **node 10649 exactly**.
  The untested prediction is now tested and it holds.
- **339 configurations** (42-config sweep + 297 refinement: all-off, each of the 256 leaves
  alone, 40 random pairs). Maximal common refinement of the channel partition =
  **8 blocks: 178 | 41 | 21 | 6 | 3 | 3 | 3 | 1** (the last seven summing to 78).
  Of tree96's 88 non-root stages, **exactly 2 are ever split** (19538 → 6|3, 10649 → 3|1),
  always identically, and both splits fall on that stage's own two input slots.
  **86 stages are never cut; zero arbitrary crossings.**
- **Independent cross-check from F's `mux_wiring.json` alone** — selector vars → defining-atom
  cone → intersect with leaves, using none of M's own measurements: root 15298 inA == the
  178-block set-for-set and inB == the union of the other seven; 19538 == blocks 6 and 3
  set-for-set; 10649 == blocks 3 and 1; 21408's live class of 11 == its inB exactly.
- **The law:** *the channel partition at configuration C is the leaf-support partition induced
  by the tree's slot structure, cut at the deepest gate C saturates.* E's residual congruences
  resolve per subtree; the all-off signatures nest strictly by row, one row per level.

**Status change:** the two derivations are the same object. The channel model **is** the tree,
established rather than suggested. A corollary worth using: the channel measurement (0.6 s for
all 256 leaves at once) is a **direct oracle for the tree's binary slot structure**, including
slot pairs F has not decoded — subject to the resolution limit below.

### Negative results M established (these bound the tooling, not the instance)

- **E's engine cannot represent the deliverable.** `forward` from the deliverable's own free
  inputs returns **39,008, not 39,026** (23 derived vars differ); E's orientation zeroes the 8
  atoms the deliverable deliberately leaves nonzero and relocates the defect. `simsolve` from
  the deliverable's configuration also gives 39,008. Any enumeration run inside E's engine is
  therefore searching a space that does not contain the current best point.
- **Resolution limit of the residual oracle is 8 blocks.** It cannot see inside the 178, inside
  the 41, or inside 21408 past one slot split. Widening the signature to full bad-atom delta
  support splits all 256 (private pin atoms). The stage-wire oracle fails outright in E's engine
  (at all-off, 219/256 leaves change no stage wire) — E's forward is propagation-with-defect and
  does not realize intermediate stage values. **The residual side cannot replace F's decode.**
- Score attempts, both below baseline: single-leaf-through-`simsolve` 38,842; deliverable-
  neighborhood pair scan 38,944 (handles are tuned to the leaf pair, so swapping leaves without
  re-solving handles is a dead instrument).

### The lever this exposes

E's channel enumeration — "the empty set wins at 39,005", monotone in live channels — was
anchored at cfg0, **whose ON-set `{1530, 1603}` lies entirely under the root's inB slot**. Both
19538 and 10649 are 78-side, so **the root gate never fires anywhere E enumerated from**. The
deliverable's ON-set `{24601, 2081}` has one leaf per root slot and scores 21 higher. E's
monotonicity result does not price the configurations that matter, and the one configuration
class known to score better is exactly the class it never visited.

Knob set for every claim above (per the standing rule): the 256 boolean leaves of
`mcore.bools()` — booleans of the cone of ROWS `[7389, 10187, 20212, 20215, 28647]` — probed
singly 0→1, configurations named per row. M's `orient.pkl` is byte-identical to E's.

Artifacts: `agentM_work/LOG_M.md` (full narrative), `RESUME_M.md`, `blocks8.json`, `mcore.py`,
`xcompare.py`, `refine2.py`, `sweep1.py`, `fullsig.py`, `ancestors.py`.

### M re-tasked

M is continuing on its own highest-value experiment: re-run the channel/representative
enumeration from a base with the **root gate firing** (the deliverable's class), after fixing
the orientation so the engine can represent an 8-nonzero-atom residual of the deliverable's
shape — otherwise the enumeration prices a space the best known point is not in.

---

## Check-in 7 — the independent audit lands (agent P)

P was one of the two agents in the second fleet that are adversarial by design: audit the
tree claim from a parse of `EQUATIONS.txt` alone, importing nothing from any agent directory.
It read FLEET.md only for the statement of the claim under audit. Its own parser reproduces
the deliverable's failing set exactly — `[12231, 12270, 12350, 14584, 18673, 22044, 29125]` —
and it did not beat 39,026.

### Confirmed, independently, with zero exceptions

- **256 leaves / 512 constants.** Pure regex over the raw file: 512 triples
  `(selector, coord, K)`, `selector → #coords = {2:256}`, `coord → #K = {1:512}`.
- **One uniform law with one universal constant.** All 383 blocks are the same shape,
  spaced *exactly 43 apart* in the straight-line program; **382/383 match the template
  byte-for-byte with identical signs**, the 383rd being the root, which obeys the same law
  with its operands swapped (the law is commutative, checked 200/200). From six inputs:
  `A = i1−i2`, `B = i4−i3`, `E = i1+i2+i5+Q`, `N1 = E·A²−B²`, `N2 = A(i3+i6)−B(i2−i5)`,
  with three congruences per block whose 3×2 integer matrix is different in every block
  (382 distinct) but always rank 2, hence equivalent to `N1 ≡ N2 ≡ 0`.
- **Invertible in closed form, 300/300** random triples.
- **The pass-through mux** `out = (1−a)b·X + a(1−b)·Y + ab·Z`: 381/381. Liveness is fully
  determined by the selectors, so the configuration space is exactly **2²⁵⁶**.
- **Root split 178 | 78** — exact. This is now confirmed from three mutually independent
  parses (F's circuit decode, E/M's residual channels, P's audit).
- **The undecoded slot pairs hide nothing.** All 27 in P's decomposition resolve: they are
  the outputs of the 27 dead blocks, provably ≡ 0 mod P (empty leaf support).

### Refuted — F's stage counts are wrong

| | F / earlier FLEET sections | P, measured |
|---|---|---|
| stages | 96 | **383 law-blocks** = 255 merges + 101 pass-throughs + 27 dead |
| depth | 6 | **9** (178-side 8, 78-side 7, plus the root) |

255 merges over 256 leaves is a proper binary tree, so the *picture* is right and every
result that leaned on the picture survives; the two numbers do not. **F ran in an earlier
session and cannot answer, so this is recorded as unadjudicated**: P's numbers are the
measured ones and F's are superseded, with provenance attached to both. Any claim phrased
against tree96's node numbering is a claim about a coarsening of the real tree — including
check-in 6's "86 stages are never cut," which stands as measured but is scoped to the coarse
object. M has been told; nothing in check-in 6's measurements changes.

### New — the law is vacuous on the diagonal, and the deliverable is exploiting it

If a merge sees two **equal** live inputs then `A = B = 0`, so `N1 = N2 = 0` identically and
that block's output is **unconstrained mod P**, after which the root can be driven to the
target by inverting the law. The satisfying set is therefore

    { S : Σ_{i∈S} L_i = T }  ∪  { S : some merge sees two equal live inputs }

and **no optimality argument in this lab has ever accounted for the second family.** P
measured that the deliverable's fold matches its own at 376 of 382 blocks, the six mismatches
being the chain `277 → 330 → 357 → 370 → 377 → 380`, which forces `stage380 == stage381` —
which is precisely what the deliverable's 4 corrupted variables and 7 failing equations buy.
(The converse is a hard exclusion, not a freedom: any S whose intermediate sum has equal
first coordinate but unequal second gives `N1 = −B² ≠ 0`.)

**This is the first structural account of what 39,026 actually is**, and it reframes the
open question from "can the subset-sum be solved" to "is 7 the price of the degeneracy, or
the price of planting it at that particular block."

### Stated limits (P's own, kept verbatim in spirit)

Everything above is mod P. Each congruence carries a small multiplier `c` on the handle side,
so the integer condition is `c·P | R`. **P did not build the lift**, so "solve the subset-sum
⇒ full solution" remains a **conjecture**, evidenced only by the deliverable being such a lift
for its own configuration. Knob set for every determined/free statement: the 256 selectors,
liveness derived.

### A note on framing, for the record

P's derivation lands on a group law with an explicit addition formula — the reading the user
directed this lab to drop. It was dropped because it was being used to argue the instance is
unsolvable, and that narrowed what anyone would search. P's measurements are kept because
they are polynomial identities it verified, and they point the **opposite** way: P found a
*second family of solutions*, not a barrier. **No infeasibility claim stands anywhere in this
lab, and none follows from anything in this check-in.** Agents are to state these as the
identities they are and draw no unsolvability conclusion from them.

### Atom indices are NOT comparable across agents

Three agents name the deliverable's residual atoms three different ways — STATE.json
`[22229, 22230, 35758…35762]`, M `{23616, 23617, 36659…36664}`, P SLP positions
`36291…36297`. These are the same object under different decompositions (the fleet has
known differing atom counts: F 39,033, I 40,885, A 42,267, P 39,277). **Never compare atom
indices across agent directories without translating first.**

### P re-tasked

1. **Build the SLP-5497 carrier and check it.** P's sliding-window screen minimises
   `|equations touching [p, p+w−1]| − w` at **4** for every `w` from 1 to 12, at SLP position
   **5497** — against the deliverable's pocket, which this lab has never left. But this is a
   structural pricer, and check-in 3 already recorded that **no function of incidence
   structure alone returns the true cost of 7** (the deliverable's pocket touches 12 equations
   of which 5 cancel, and cancellation is residue content, not incidence). So the deficit of 4
   is a candidate, not a predicted 39,029. What settles it is an assignment, verified with the
   checker; a below-baseline result is to be reported as plainly as a win, because it would
   retire the last structural instrument.
2. **Price the degeneracy by placement** — enumerate where else a forced equality can be
   planted and what each placement costs. If a cheaper block exists, that beats 39,026
   *without* solving the subset-sum.

---

## Check-in 8 — three reports converge on what 39,026 is (agents S, K, M)

Deliverable unchanged: **39,026 / 39,033**, failing `[12231, 12270, 12350, 14584, 18673,
22044, 29125]`. S, K and M each re-verified it themselves; none beat it. S's own best
construction was 39,019; M's best over every base it priced was 39,008 before its engine fix.

### The convergence

Three agents working in three models arrived at the same account of the deliverable:

- **P (check-in 7)**: a merge seeing two equal live inputs has `A = B = 0`, so both its
  congruences vanish identically and its output is unconstrained. The deliverable's fold
  differs from P's at six blocks, forcing two blocks equal.
- **K, independently**: the deliverable's 7 failures are **one gate's off-pins**, and they
  work by forcing the **root's two inputs equal**, which makes the root check vanish and
  leaves the root output free — then it is set to the target. Verified that
  `2·leaf(24601) ≠ T` and `leaf(24601) ≠ T` while every root atom is zero.
- **M, in equation space** (a frame nobody had used): greedy keeps 79/86 rows and returns
  **exactly 39,026** — the deliverable already *is* that solve's optimum — stopped by a
  **divisibility obstruction on equation 29125**.

So 39,026 is not a search failure. It is the price of a degeneracy that cannot be had honestly.

### K — the degeneracy is unreachable by configuration (new closed negative)

**No stage anywhere can be made degenerate.** Interior stages: `|x−y| < 2ⁿ < N` forces
`x = y`, impossible on disjoint bit sets. Root: `x−y = ±N`, and because each bit belongs to
exactly one half the carry walk is deterministic — both directions end with a nonzero carry.
**Knob set as stated by K:** all 256 leaf selectors, all 2²⁵⁶ configurations, every other
variable free; also checked against all 256 single-exponent reassignments and **0 of 2000**
random 178/78 partitions. The adjacent hole (a half folding to the identity) is closed too.
This is consistent with P rather than contradicting it: the degeneracy is reachable only by
breaking pins, which is what the deliverable pays for.

**K's site-cost table** — the first thing in this lab that prices 7 against alternatives
instead of asserting it: deliverable's site **7** equations, cheapest leaf-pin pair **10**,
breaking the target pins **16**, and **the 12 atoms with footprint < 7 are all decoy
idempotency atoms**. That last clause bears directly on P's SLP-5497 window screen, which
selects precisely for low-footprint atoms; P has been asked to test the window against it
before building anything.

K also finished the decode (exactly 516 literals > 10²⁰: 512 leaf pins, 2 target constants,
`p`, `K` — nothing else), confirmed the 178|78 split independently, and **refuted its own
first two fold validations mid-session**: a plain cascade runs the target pin *backwards*
into the tree and returns "A = target X". Anyone building an evaluator must forbid those two
atoms. Open and not claimed: multi-leaf **B-half** folds are not yet reproduced by K's
closure (A-half is exact), so its uniqueness sentence is marked *strongly supported, not
proved*; and `k25_class.py` misses two of three idempotency atom spellings (logged as a bug).

### M — E's engine was structurally incapable of representing the deliverable, and is fixed

M's quadratic-branch hypothesis was **wrong**; the defect is structural. `harness._bootstrap`
assigns each derived variable a **definer atom** and `forward` solves that atom **to zero**,
so every definer atom is identically satisfied in every state E's engine can reach — and
**five of the deliverable's eight nonzero atoms are definers** (23616→x_7068, 23617→x_28730,
36659→x_29854, 36663→x_31864, 36664→x_642; atom 36663 *is* the expression `x_31864`, so the
engine defines that variable as zero while the deliverable sets it nonzero). The other 18
differing variables were downstream contamination.

**Fix (`engine2.py`): demote those 5 atoms, promote 5 variables to free. Zero variables now
differ from the deliverable; 39,026/39,033; same 7 failing equations; all 8 nonzero atoms
intact.** Guarded against overfitting: engine2 predicted 39,000 and 38,961 at two
non-deliverable points and `checker.py` confirmed both exactly. Those 5 variables drive 7 of
the 8 residual atoms **affinely** and were invisible to every solver before.

**Fleet-level consequence:** any result computed inside the unfixed engine was computed in a
space that **cannot contain the deliverable**. Check results that came out of it before
building on them.

**Monotonicity is an ARTIFACT.** M corrected its own check-in 6 reasoning: E's 178-block
representatives are A-side, so its rows *did* fire the root gate — but every root-firing
configuration E priced also carried cfg0's two B-side leaves, so **E never priced a clean
2-leaf root-firing configuration**. From a neutral base (all 2⁸ block-subsets): 38776, 38791,
**38804**, 38787, 38773, 38761, 38746, 38718, 38655 by live-block count — **unimodal, peaking
at 2, not monotone**. Root-firing beats 78-side-only 38,804 vs 38,791 overall and 38,804 vs
38,774 at equal count. The lattice optimum `(47, 490)` reproduces **the deliverable's slot
pattern**, recovered independently.

M's secondary for F is closed and negative: the 3 stages absent from `mux_wiring` are trivial
1|1 splits, so **the oracle does not advance the inversion attack**; `blocks8.json` stands.
M also fixed a bug in its own crossing test and downgraded its 32-block refinement from
"spurious" to "unvalidated" once P's audit removed the yardstick's authority.

### S — the filter that caused an earlier false barrier was still in place

S found it and removed it: E's channel measurement probed **boolean knobs only**. Over all
333 cluster-cone knobs there are **pure single-atom handles** — 1,815 free variables that move
exactly one atom, affinely, by exactly **±p** — so an atom is satisfiable iff its residual is
≡ 0 mod its handle step. a28647, a30787, a26958, a40306, a726 have **no** handle.

`lat3.py`: over the complete affine knob set (54 knobs, every single-row knob, every atom's
handle), eliminating the two target rows and solving the other 47 exactly over ℤ (feasible,
kernel dim 7), the reachable lattice on `(a20215, a28647)` is **exactly p·ℤ²** — so the
endgame condition is precisely `a20215 ≡ 0 (mod p)` **and** `a28647 ≡ 0 (mod p)`. A BFS over
configurations **terminated by exhaustion** at 48 mod-p tuples: a7389, a10187, a20212, a28647
all reach 0; **a20215 takes only 2 values and never 0**. The endgame therefore reduces to the
**single** condition `a20215 ≡ 0 (mod p)`.

S also confirms saturation is stronger than reported — it holds **across** classes, not just
within, and both non-boolean "knobs" are switches (1 through 10²⁰ give identical deltas) —
and refutes "the residual is a subset-sum" on the grounds that saturation makes it one-hot
selection, so **LLL on the bit vector is the wrong tool**. Whether that contradicts P's
subset-sum-plus-degenerate statement or merely re-describes it under a different
decomposition is **unresolved and flagged to both**.

### Standing cautions reinforced this check-in

- **Structural pricing is still a dead instrument.** S's §6 (769 footprint-1 atoms, so ≤6
  would beat the deliverable) and P's SLP-5497 window are the same kind of count, and K's
  table says the sub-7 footprints are decoys. Nothing here counts until it is an assignment
  the checker scores.
- **Three independent derivations have now landed on the same group law** (P, K, and earlier
  work). It is the reading the user withdrew, because it was being used to argue the instance
  is unsolvable. The measurements are kept as the verified polynomial identities they are;
  **no infeasibility claim stands anywhere in this lab, and none follows from this check-in.**
  K's own report says it plainly: "Not claimed anywhere: infeasibility. The instance is
  satisfiable." K has additionally been told not to call the final recovery step "mechanical"
  — every part of it is mechanical except the part that is the whole problem.

### Live tasking after check-in 8

| Agent | Task |
|-------|------|
| K | Close the B-half validation gap (instrument `CascadeP.close`, walk back from `x14853`); fix the idempotency-spelling bug |
| M | Characterise the divisibility obstruction on equation 29125 exactly; state the knob set, including the 5 newly-freed definer variables |
| P | Test the SLP-5497 window against K's decoy-idempotency finding **first**; then the carrier, only if it survives |
| S | Finish `lat5.py` (stopped 29/48, 0 feasible); reconcile its exhaustion with K's unreachability claim |

---

## Check-in 9 — the reduction is complete, and a factual conflict is open (agent L)

Deliverable unchanged: **39,026 / 39,033**, re-verified by L with plain `checker.py`. L did not
beat it; its own best is 39,018 (`agentL_work/assign_L1.json`, checker-verified), and L is
explicit that this is the canonical output of its constructor rather than a search result.

### The complete decoded reduction

All **383 nodes calibrated** — coordinate alignment from slot links, chord orientation solved
from each node's own three stage checks against the universal constant — **383/383, 0
failures**; all 256 leaf pin pairs extracted numerically with **0 conflicts**. For every ON-set
tried, the entire 39,033-atom system collapses to the **same 2 atoms**:

> The system is satisfiable iff some non-empty subset of the 256 leaves folds to
> `TARGET = (44859544763832475231923253825569092119321525945631045653619508440821028887,
> 36200939269128454586076546451607958467047992891178506183612554289882454126226)`
> (root coordinate order).

Everything else closes exactly over ℤ. The lift is **free**: every one of the 3,681 handle
variables appears in **exactly one atom** (knob set: all 8,747 free variables), so mod-p
suffices. This closes the conjecture P left open in check-in 7 — "solve the subset-sum ⇒ full
solution" is no longer conjectural in L's model.

L's node count of **383** independently corroborates P's recount against F's 96. L identifies
the unique root as **x9274**, with 384 leaves = **256 free booleans + 128 literal `:= 0`**, and
confirms F's configuration count of **2²⁵⁶ − 1**.

### The assigned case: no sum

L's task was the same-OR-group double-leaf case. Measured by exact re-propagation: `a` only →
the slot holds A exactly; `b` only → B exactly; **both → chord(A,B)** — not `A+B`, not A, not B.
All 10 local atoms vanish mod p in all three cases. **Mechanism:** F was right that both pins
fire — the leaf value-wire off-guard is `(1−leaf_bit)·w = p·h`, guarded by the *leaf bit*, not
the selector — but the mux coefficients are the mutually exclusive quadrants, so both firing
pins are multiplied by zero. Generalised to ON-sets of size 1, 2, 5, 73, 200, 256 (containing
same-group multi-leaf cases at every depth): all close. **F's same-OR-group sum caveat is
refuted; the fold law is unchanged.**

L also **retracted two of its own claims** mid-session: its 2^178 configuration count (F's
2²⁵⁶ − 1 stands) and its 178/78 measurement, which was a *sub*-forest under x8599/x21839 —
the lab's recurring 178|78 is the live-leaf split **at the root**.

Separated cleanly: **E's saturation is the pass-through branch; the double-leaf case is the
chord branch.** Different branches of the same mux, not the same phenomenon.

### OPEN — a factual conflict about the deliverable's ON-set

**L reports the deliverable has exactly one leaf ON (x24601)**, the trivial pass-through.
**M measured `{24601, 2081}`. K independently reports that at ON = {2081, 24601} the circuit
gives A = leaf(24601) and B = leaf(2081) exactly.** Three models, two answers.

This is load-bearing, not bookkeeping: K's and P's account of *why* the deliverable scores
39,026 — one gate's off-pins forcing the root's two inputs equal, killing the root check —
requires **two** live inputs at the root to be the mechanism at all. If L's reading is right
that account is wrong; if K/M are right, L's extractor is dropping a leaf. L has been asked to
settle it from the deliverable's own assignment and report which way it falls.

### Converged, from four models: beating 7 means searching cancellation, not support

L, knob set = the 378 parent-slot-wire pairs, incidence-only and therefore
configuration-independent: **the cheapest 2-atom cut anywhere costs 11**, while the deliverable
reaches 7 with *four* nonzero atoms — i.e. by **value cancellation inside shared equations**.
L's own inversion-based cuts priced 39–47. This is the same conclusion as K's site-cost table
(sub-7 footprints are decoys), M's equation-space result (8 atoms → only 7 failures because
they cancel inside equations), S's §5 (5 of 12 equations cancel), and check-in 3's finding that
no function of incidence structure alone returns the true cost of 7.

**Four independent models, one instrument nobody has built.** L has been re-tasked to build it,
since it has the only tooling that can: 383/383 calibration, exact `invchord`, and a
checker-verified constructor.

### Searches run and negative

|S| = 1 (256), |S| = 2 (32,640), |S| = 3 (2,763,520): **no hit**. |S| = 4 (174M) was left
running (~6.8 h). **Zero degenerate folds observed** — consistent with K's closed negative that
no configuration can make a stage degenerate. L has been told not to block on |S| = 4: 174M of
a 2²⁵⁶ space is not where the value is, and three empty levels are weak evidence about the
fourth and none about the rest. Root meet-in-the-middle is 2⁷⁸ and out of reach.

### Position of the campaign after this check-in

Four agents in four models now agree on the shape: the instance reduces to selecting a subset
of 256 leaf constants that folds through a uniform chord law to a fixed target, the lift off
that is free, and the deliverable's 39,026 is the price of a degeneracy at the root that
cannot be obtained by configuration. **No infeasibility claim stands anywhere in this lab**,
and none of the above is one — but the honest statement of the remaining gap is that every
route to 39,027+ now runs through either the cancellation search (unbuilt, L tasked) or M's
divisibility obstruction on equation 29125, and that the subset-sum itself has resisted every
enumeration attempted. Enumeration is not the lever; the fleet should not spend further cores
on level-by-level sweeps.

---

## Check-in 10 — the footprint screen is dead, and 7 is a placement price (agent P)

**P withdrew its own primary prediction**, in the first paragraph of its report, with the
measurement that killed it. Deliverable unchanged at **39,026 / 39,033**; P produced no new
assignment and says so plainly.

### The SLP-5497 carrier does not exist

SLP position 5497 is `x29741 − x17440 + x27926` — the `A = i1 − i2` wire **inside law-block
68** (leaf-support 2, both inputs leaves), footprint 5. Perturbing it changes that block's
`N1`/`N2`, which **forces the block's three law-congruence atoms at SLP 34872 / 34874 / 34876**
(footprints 13 / 12 / 11) to become nonzero. The realizable defect is therefore **4 atoms
touching 22 equations, not 1 atom touching 5**. The only escape is a configuration where block
68's gate is off — and then the perturbation never reaches the root, so it buys nothing while
still costing its 5 equations.

**The footprint screen is dead.** `min_p(|eqs touching [p, p+w−1]| − w) = 4` remains true as
incidence and is useless as a cost. This is the **third independent demonstration** that
incidence-only pricing fails (C at check-in 3, L's cheapest-2-atom-cut-costs-11, now P), and
the first to identify the error term exactly: **an atom's true price includes the congruence
atoms its perturbation forces downstream**, which is residue content, not structure.

### The decoy test, run first as instructed — and the window is *not* a decoy

K's decoy explanation is **corroborated at much larger scale**: in P's parse, **1,158 atoms
have footprint < 7, of which 1,152 are idempotency atoms** (`x − x²`; 1,145 at footprint 1).
K's "12" versus P's "1,158" is a decomposition difference, not a disagreement.

But **every atom in 5490–5511 is genuine law-block arithmetic** — products, squares, three-term
linears, `isidempotency = False` throughout. The window sits in the 6-atom non-decoy remainder.
So the carrier died of a **third failure mode** (downstream coupling), and **both explanations
are needed; neither subsumes the other.**

### 7 is the price of the PLACEMENT, not of the degeneracy

P answered this from material already built, without the long enumeration: cheapest **live**
merges touch 10–11 equations (block 279 → 10; blocks 2, 151, 193, 311, 330 → 11); the
deliverable's site touches 12 with 7 atoms; worst is 39. But **corrupting a live merge to hit
the target directly breaks three congruences, not two** — with `Z` forced to the target, `N1`
and `N2` are determined and no rank-2 combination `c_k1·N1 + c_k2·N2` vanishes. That is
precisely why the deliverable takes the vacuous route instead. **P found no placement below 7
and stopped there as directed.**

### A specific hole in K's unreachability negative — routed to K

P cannot adjudicate K's argument (it is stated in a vocabulary P's parse does not produce), but
identified one precise gap: both `|x−y| < 2ⁿ < N` and "at the root `x−y = ±N`" require the
modulus actually governing coordinate-pair equality to **exceed the largest signed subset
difference**, so that only one wrap (`k = ±1`) needs excluding. If that modulus is below 2²⁵⁶,
the enumeration runs over `x − y = kN` for several `k` and a carry walk covering only `k = ±1`
is **incomplete**. K has been asked to state which modulus bounds the walk and confirm it
exceeds 2²⁵⁶. Until then **K's negative is not settled**, and agents have been told not to
treat it as such.

P's second caution on it: K's argument rests on a premise about how the instance was built
(256 leaves ↔ distinct exponents in an order-N arithmetic), which sits inside the withdrawn
reading and inherits its caution. K has been asked whether the carry argument can be restated
using only the verified identities.

### P vs S adjudicated — subset selection, not one-hot

Decisive and cheap: **zero atoms in the file touch two or more distinct selector variables.**
Each of the 256 selectors appears in only 5–6 atoms, all local to its own coordinate load and
liveness fan-out. No cardinality atom, no one-hot tie, no cross-selector coupling. P's residual
is therefore **free independent subset selection over 256 booleans (2²⁵⁶)**. The one-line test
— *count atoms containing two or more selector variables* — has been sent to S to run against
its own parse. If S's parse also finds none, the one-hot reading is describing a different
object and the difference lies in what S's saturation measurement ranges over.

### The ON-set conflict — now three models to one

P states the deliverable's configuration as **leaves {21, 167} in its own numbering — two live
leaves**, joining M (`{24601, 2081}`) and K (A = leaf(24601), B = leaf(2081) reproduced
independently) against L's single-leaf reading. L has been told, and told the numbering caveat
cuts both ways: **compare by the pinned constants those leaves load, not by index.**

### P re-tasked — the caveat it has carried in every report

P has flagged in three consecutive reports that everything is mod P, that each congruence
carries a small multiplier `c` on the handle side so the integer condition is `c·P | R`, and
that "solve the residual ⇒ full solution" is therefore a **conjecture**. L now claims that
conjecture closed in its model, on the strength of **every one of the 3,681 handle variables
appearing in exactly one atom**. P — the agent that kept insisting it was unproved — has been
asked to test that property against its own parse and report the count either way. If it holds,
this lab's central reduction becomes unconditional; if it fails, L's reduction is over-claimed
and everything resting on "the lift is free" needs revisiting. P was told to test the property,
not to read L's artifacts.

### Vocabulary

P restated its §3 neutrally and unprompted: each block imposes the identities `E·A² − B² ≡ 0`
and `A(i3+i6) − B(i2−i5) ≡ 0 (mod P)`. What it verified are properties **of those identities** —
closed-form solution for `(i5, i6)` 300/300, operand-pair symmetry 200/200, order-independent
iteration on the 256 constant pairs 300/300. **No solvability conclusion is drawn, and P notes
its results point the other way.** This is the register the whole fleet is now writing in.
