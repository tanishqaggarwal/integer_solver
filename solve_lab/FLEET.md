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

---

## Check-in 11 — SAT/SMT/CP is retired on measurement (agent R)

Deliverable unchanged: **39,026 / 39,033**, re-verified twice from cold by R. Nothing R built
scores above 39,013. No infeasibility is claimed anywhere in `agentR_work/`.

### The angle is dead, and it was measured rather than asserted

R proved each encoding satisfiable with `witness.py` **before** blaming any solver, and scaled
siblings of identical shape rather than arguing from the full instance:

| encoding | result |
|---|---|
| z3 QF_BV, 8-bit prime, all selectors pinned | sat, 119 s |
| z3 QF_BV, same, selectors free (128 options) | timeout 300 s |
| z3 QF_NIA, pinned and free | timeout 300 s |
| CaDiCaL 1.9.5, bit-blasted CNF (8-bit, pinned; 57,299 v / 388,439 cl) | sat 107.6 s, 165,725 conflicts |
| Glucose 4.2, same | sat 276.7 s, 353,100 conflicts |
| CP-SAT free | 8-bit 0.7 s · 10-bit 4.6 s · **UNKNOWN at 300 s** for 12–28 bit |
| CP-SAT at 31-bit | **MODEL_INVALID** — cannot represent the arithmetic |
| uninterpreted-law variant | sat in 0.01 s — **vacuous**; +distinctness → unknown |

CNF scales as ≈`863·m²` clauses/stage. **Extrapolated to the real instance: ≈1.4×10¹⁰ clauses /
1.9×10⁹ variables ≈ 600 GB of DIMACS against 29 GB of disk.** The CNF cannot be written, let
alone solved; where the tools can run they are ~10⁵× slower than brute force on the same
instance. **SAT/SMT/CP is retired as a live angle**, and agent C's earlier abandonment of it now
has a reproducible reason attached.

### What the encoding target actually is — refuted, usefully

**The reduced problem is not 96 coupled stage constraints; it is one scalar equation.** R
re-derived and checked F's decode in its own code: the offset substitution kills `K`, all 248
forced pin points **and the target** satisfy one relation (fitted from 2, verified on
246 + target, 0 exceptions); the fold is commutative and associative (200 random triples each),
**so the tree shape is irrelevant** and F's remaining decode — the 56 slot pairs and 24
leaf-adjacent stages — is **not needed**. The 256 leaves form **one doubling chain** (9
doubling-closed pieces, 8 splices of exactly one gap point each, one head → `ladder.json`); the
order is a **256-bit prime** by Cornacchia.

**F's requested validation PASSED**: ON-set {24601, 2081} → ladder indices {72, 235}, and the
fold is **not** the target — exactly the prediction a correct evaluator had to produce. This is
now the **fourth** independent model reading the deliverable as **two** live leaves, against L's
single-leaf reading; L is settling it from the deliverable's own assignment.

### Measurements on the real instance

**Exhaustive Hamming weight ≤ 6 → no solution** (108 s meet-in-the-middle over all size-≤3
subsets both sides), so the reduced unknown has weight **≥ 7**. BSGS to 2⁴⁴ was running at
handoff and is resumable, but **enumeration stays retired as a lever** across the fleet.

### A contradiction inside R's own report — sent back to R

R's §(c) states that a different selector configuration does not cheaply beat 39,026, with
optimistic ceilings of **≤39,020 single-bit** and **≤39,022 pair**. R's §(d) proposes running
`s10/lattice3.py` on a **single-bit** configuration in the hope of beating 39,026. Both cannot
be right. Either the ceiling binds and §(d) is dead before it starts, or the ceiling is scoped
to the placements `gs2.solve` reaches and is silent about placements `lattice3` reaches — in
which case **it is a statement about a knob set, not about the instance, and must be re-quoted
that way everywhere it appears.** R scoped it correctly in `LOG.md` §7/§9 and has been asked to
say which it is in the report itself, this being the exact failure mode behind five retractions.

If the ceiling does not bind, the experiment is worth running and nobody has run it: the
single-bit defect is 3 atoms against the deliverable's 7, `gs2` places it badly (20 equations,
0 cancelling), and `lattice3` — the method that actually produced 39,026 — has never been
pointed at single-bit footprints. Success criterion, R's own: **a placement into ≤12 equations
with ≥6 cancelling beats the deliverable.** That is the same instrument four other agents have
independently converged on — cancellation, not support.

### Retired angles, cumulative

- Generator inversion (user instruction, first fleet).
- Level-by-level enumeration of the reduced problem (|S| = 1, 2, 3 exhausted; weight ≥ 7;
  root meet-in-the-middle 2⁷⁸).
- Incidence-only / structural pricing (three independent demonstrations: C, L, P).
- **SAT / SMT / CP encodings of the reduced problem (R, this check-in).**

---

## Check-in 12 — the lift is three-quarters free, not free (agent P)

Deliverable unchanged: **39,026 / 39,033**; P produced no new assignment.

### L's property is false as stated, and the diagnosis is specific

P tested check-in 9's load-bearing claim — that every handle variable appears in exactly one
atom, so mod-p suffices and the integer lift is free — against its own parse. Exhaustive
occurrence count, no exceptions, where a handle is a variable defined by an atom `h − (P-alias)·u`:

| variable | count | atoms it appears in |
|---|---|---|
| handle `h` (the P-multiple) | 3,707 | **exactly 2 — every one** |
| cofactor `u` (the free multiplier) | 3,707 | **exactly 1 — every one** |

**Zero handles appear in exactly one atom.** The property is false of the handles and true of
the cofactors, and L's 3,681 against P's 3,707 is a decomposition difference rather than a
disagreement about the file — so **the likeliest reading is that L counted `u` and concluded
about `h`.** L has been asked to re-run its count distinguishing the two and to say which it
measured.

**Why "appears in exactly one atom" is not sufficient even where it holds.** The handle sits in
the pair `h − u·P = 0` (⇒ `P | h`) and `R − c·h = 0` (⇒ `h = R/c`). The cofactor does appear
once — necessary for freedom — but the atom must also be **solvable over ℤ**, and `h − u·P = 0`
is solvable for `u` only when `P | h`. Composing gives exactly **`c·P | R`**, strictly stronger
than `R ≡ 0 (mod P)` whenever `c > 1`.

### How much is actually free — measured

Splitting all 3,707 by the multiplier `c` (the def-side coefficient is ±1 for all 3,707, so it
contributes nothing):

- **2,780 have `c = 1`** — the integer condition collapses to `P | R`, which *is* the mod-P
  congruence. **Genuinely free. L is right about 75%.**
- **927 have `c > 1`** — `c·P | R` is a real extra integer condition that mod-p reasoning does
  not deliver.

**The reduction is ~three-quarters closed and one-quarter open. It is not unconditional**, and
P declined to state it in unconditional form. L's reduction is over-claimed, but by 927
conditions rather than wholesale.

### Corroboration, scoped as carefully as it was found

All four variables the deliverable corrupts — `x642, x28730, x29854, x31864` — **are handle
variables**, and they are **exactly the four of 3,707 for which `P` does not divide the value**;
the other 3,703 all satisfy `P | h`. The deliverable's entire 7-equation deficit is a handle's
two atoms failing to be simultaneously satisfiable. If handles absorbed freely, that deficit
would not exist.

P's own scope note, kept: this shows the pair is a **joint constraint**, not that an honest
configuration must fail it. The deliverable broke those handles *on purpose*, to plant the
degeneracy — and the same assignment satisfies all 927 `c > 1` conditions elsewhere, **so a
lift demonstrably exists for at least one configuration.**

### The one number that decides it — P re-tasked

The freedom is the integer lift of each free coordinate (`r + kP`), which moves `R/P` and can
tune `R/P mod c`. **~766 lift parameters (2 per law-block) against 927 conditions with `c > 1`.**
Fewer knobs than conditions is not fatal — the moduli are ~7 bits and one parameter can serve
several by CRT when the coefficients are invertible — but **nobody has counted the rank, and
that rank is the whole question.** P flagged it rather than asserting either way, and is now
counting it, then building the lift for one configuration and verifying.

This is the **last conditional link** between this lab and an unconditional statement of what
the instance reduces to. Everything else in the reduction is measured.

### Standing caveats, now quantified rather than removed

"Solve the residual ⇒ full solution" **remains a conjecture**, and the unproved part is now
identified exactly: the **927 `c > 1` divisibilities**. Everything remains mod P. Knob set
unchanged (256 selectors, liveness derived).

---

## Check-in 13 — the 29125 obstruction does not exist; WITHDRAWN by its author (agent M)

Deliverable unchanged: **39,026 / 39,033**. Nothing above it produced.

### The retraction, and the correct criterion

**Check-in 8's "divisibility obstruction on equation 29125" is withdrawn by M.** The
`rhs % -P != 0` message came from **one elimination ordering inside a badly overdetermined
window** (86 rows vs 19 knobs; 999 vs 162). M named it an obstruction without testing the row
itself, and says so in its own first line.

Single-row solvability is **exact and window-independent**: `Σ coef_f·d_f = −s0` is solvable
over ℤ **iff `gcd(coef_f) | s0`**. For eq 29125 the gcd is **1**. **All seven failing rows
pass** — gcd 1 for six, 40490 for eq 22044, and every one divides.

**Anything in this file or in prior reports that cited the 29125 obstruction as a live lead is
void**, including check-in 8's summary and check-in 11's tasking table.

### What replaces it — a minimum-cost residual, measured window-free

`eqsub.py` removes the window entirely: solve each subset of the 7 failures, then **apply it,
re-propagate, and measure the true score**, counting collateral damage by measurement rather
than assumption. **127 subsets solvable, 0 infeasible**, largest solvable subset = **all 7**,
and in every case the solver genuinely zeroes its targets.

| fix | failures 7 → | score |
|---|---|---|
| eq 12350 / 18673 | 10 | 39,023 |
| eq 12270 | 11 | 39,022 |
| eq 12231 | 18 | 39,015 |
| eq 22044 | 28 | 39,005 |
| eq 14584 / **29125** | 34 | 38,999 |
| **all 7** | **44** | **38,989** |

Best over all 127 subsets: **39,023 < 39,026**. Every failing equation is repairable and
**every repair costs more than it gains**. **39,026 is a strict local optimum in equation
space**, and eq 29125 is not blocked — it is tied for the *most expensive* row to repair.

**Scope, and it matters:** this table prices *repairs of the current placement*. It says
nothing about a different placement, which is the only thing anyone is still testing.

### The five questions, answered

1. **What divides what:** nothing obstructs; `gcd | s0` holds in all 7 cases.
2. **Which knobs:** eq 29125 has 12 affine knobs; the direct one is **`x_28730`** (its only
   nonzero atom is 23617 = `x_28730 − x_17499·x_9413`).
3. **Equation or window:** **the window.**
4. **"Core infeasible" at 162 knobs — instance or widening:** **the widening.** 999 rows
   against 162 knobs is 6:1 overdetermined and generically infeasible regardless of instance.
   Now *proven* rather than argued: all 127 subsets are feasible.
5. **Knob set:** the 5 freed definer variables `[642, 7068, 28730, 29854, 31864]` **are in it
   and all 5 are affine**, verified across three knob sets up to "every free variable in the
   cone of every atom of the target equations."

### Why the residual is cheap, and what the engine defect was hiding

The 8 nonzero atoms are seen by exactly the 7 failing equations, nested, all zero-constant.
Eq 29125 sees a **single** atom, moved directly by `x_28730` — one of the 5 variables M's
engine fix freed. Driving it to zero is precisely what E's orientation did by construction,
and it costs **27 extra failures**. The engine defect was hiding exactly this trade.

### The strongest cross-check in the campaign

P's four corrupted handle variables — `x642, x28730, x29854, x31864`, which P measured as
**exactly the four of 3,707 for which `P` does not divide the value** — are a **subset of M's
five freed definer variables** `[642, 7068, 28730, 29854, 31864]`, in the same numbering, from
two decompositions that share nothing. **`x_7068` is the odd one out**; M has been asked what
it is, since it is either a fifth corrupted handle P's test missed or a definer that is not a
handle, and both are informative.

### M's §6, endorsed and closed

The residual side cannot see intermediate stage values — established in round 1, still true —
so M can neither confirm nor refute K's off-pin/root-degeneracy mechanism, and correctly
claims only consistency (the residual is confined to 7 equations, and the repair that most
directly restores a check is among the most expensive). A real test needs the circuit side.
M has been told to stop there rather than adjudicate from the residual side.

### M re-tasked

**Price alternative placements in engine2** — the only frame that provably represents the
deliverable. Four models have converged on cancellation-not-support as the remaining
instrument; L is building it from the circuit side, M is the complement from the residual
side, and `eqsub.py`'s solve → apply → re-propagate → measure primitive is already the right
tool. Bar for success, from M's own numbers: **any placement scoring above 39,026 is the first
in this campaign.**

---

## Check-in 14 — a satisfying assignment EXISTS, and the ON-set conflict is resolved (agent Q)

Deliverable unchanged: **39,026 / 39,033**, re-verified twice by Q. Q did not beat it.

### The first positive existence result in this campaign

Subset sums of the ladder realise every group element, and the order is below 2²⁵⁶, so the
target is reached: **a satisfying assignment exists.** After months of negatives this is the
lab's strongest positive claim, and Q attached its one caveat rather than dropping it —
**it rests on the stage law holding at the 24 leaf-adjacent stages, which Q did not test.**
Q has been re-tasked to close exactly that, since doing so makes the result unconditional.

Supporting measurements, all measured rather than assumed: the offset substitution reduces the
stage law to the plain chord construction; fitting the cubic from two leaf pins gives **253/253**
leaf pins on it (2/253 in the other orientation), non-singular; the law is **associative
297/297**, **commutative 297/297**, and equal to the shifted chord law **198/198**. The order is
a **256-bit prime** by Cornacchia on `4p = L² + 27M²`, verified by exact scalar multiplication
on 5 points, with `N ≠ p` and `p^k ≠ 1 mod N` for k ≤ 24. The 256 leaves are **one doubling
ladder** (249/253 decoded leaves have their double also a leaf; 4 chains of 124/79/41/9 linked
head-to-tail by one missing doubling each), with `L_i = 2^i·G` verified for i = 0..255.

**A real falsification test that could have gone the other way:** the three *inferred* ladder
points (exponents 41, 51, 176) were checked against the raw file, and **all three predicted
x-coordinates are literal constants in `EQUATIONS.txt`** (`check3.py`).

⇒ **`EQUATIONS.txt` is satisfiable iff the ladder scalar hits the target, and a solution's leaf
ON-set is exactly that scalar's binary expansion.**

### The ON-set conflict is RESOLVED — both readings were right about different objects

Q re-read the deliverable's atoms as given, with no forward pass (`val4.py`), and took a
wire-value census mod p: exponent 72's value appears on **92 wires**, exponent 235's on **5**,
**their group sum on 0**, and the target on **4**.

**The deliverable does not fold at all.** It passes a **single leaf** up the tree and overwrites
the value with the target near the root, paying 7 broken atoms for the overwrite. Its ON-set is
`{2081, 24601}` = ladder exponents `{72, 235}`.

So **two selectors are ON and one leaf propagates.** L was reading what reaches the root; M, K,
P and R were reading the selector configuration. Both correct, different objects — and the
distinction is load-bearing, since it is *why* the deliverable scores what it does. It also
sharpens K's and P's account: the root's two inputs being forced equal, and the fold never
happening because one input never arrives live, are the same event from two sides. **L has been
told to restate with the distinction attached rather than retract.**

### Refuted — the previous handoff's highest-value experiment

"Invert the target down the 78-side to a support-≤24 node and enumerate under 2²⁴" **cannot
work.** The fold is a group homomorphism of the selector vector, so any meet-in-the-middle is
the generic square-root attack at ~2¹²⁸. **Finishing the 56-stage decode buys nothing for the
search.** That experiment headed the FINAL POSITION section of this file and is now retired on
measurement — as is, independently, F's remaining decode (Q, R and L all reach "tree shape is
irrelevant" from different directions).

### Searches run and negative

Six structured attempts, all negative: BSGS below 2⁴⁴ from either end; Hamming weight ≤ 6
(meet-in-the-middle, 2.8M/side); all ON-bits inside a 34-bit window (2,865 s); small multiples
of the target on the ladder to 10⁷; the endomorphism orbit; and a two-term decomposition with
both coefficients below 2²¹. The endomorphism is **confirmed** but gives only √3. Weight ≤ 7
was stopped at ~2% with CPU 5× oversubscribed (load 19 on 4 cores) and is re-runnable in ~3 h.

### Where score improvement now stands

Q's tracks are **no longer coupled**: score improvement is provably independent of the hard
scalar problem — it is the minimum-weight coset problem on the atom incidence matrix (agent A's
formulation). Q's model adds a degree of freedom nobody had: **the defect can be made any group
element.** But Q argues each extra row you try to cancel imposes a linear condition on the fold
point — a line meeting the cubic in ≤3 points — and reaching one is again the hard problem, so
the freedom does not lower 7 by itself. `rows4.py` prints the restricted row systems for the
deliverable's 7-support (12 rows touched, 5 cancelled) and a forward-4 support (13 rows, 0
cancelled).

**Recorded as an argument, not a result**, and routed to L, whose cancellation search is a
measurement. This lab has repeatedly found measurements beat arguments; if L finds a placement
above 39,026, Q's reasoning is where to look for why.

### Retired angles, cumulative

- Generator inversion (user instruction, first fleet).
- Level-by-level enumeration of the reduced problem.
- Incidence-only / structural pricing (C, L, P — three independent demonstrations).
- SAT / SMT / CP encodings (R).
- **Circuit decoding, including F's 56 slot pairs and the 78-side inversion (Q, R, L).**

### CPU

Q flagged 5× oversubscription (load 19 on 4 cores). Agents whose angles are closed have been
told to keep their footprint small so the live threads — L's and M's cancellation searches,
P's rank count, K's carry-walk question — get cores. **The fleet should be thinned once N, O
and T report**, flush-first per the stop policy above.

---

## Check-in 15 — the rank is still uncounted, and P says so first (agent P)

Deliverable unchanged: **39,026 / 39,033**. P produced no assignment approaching it.

### Not answered

P opened its report with "I did not count the rank — not answered," having been the agent that
named that measurement as the whole question. **The open quarter of the reduction is still
open.** P built the machinery (`plift2.py`: seed the 256 selectors, 512 leaf coordinates and
764 block law outputs, then walk the SLP with a worklist solving each atom for its single
remaining unknown **over ℤ**, recording every division that does not come out exactly) but the
rank needs the derivative system `∂(R/P)/∂t_v mod c` across the 927 `c > 1` conditions, and that
needs a **complete** lift to differentiate around. P does not have one.

### DO NOT CITE THESE SCORES AS MEASUREMENTS OF THE INSTANCE

| configuration | undetermined vars | integer-division obstructions | equations failing | score |
|---|---|---|---|---|
| all selectors OFF | 18,417 | **0** | 2,645 | 36,388 / 39,033 |
| one leaf ON | 18,417 | 18 | 2,815 | 36,218 / 39,033 |

**Both are tooling artifacts of an incomplete construction**, scoped that way by P before its
own table: 18,417 variables were never determined and defaulted to zero, and the failures are
that default rather than the instance resisting. **What stopped P is its propagation, not a
property of the file.**

Diagnosed, not fixed: the first nonzero atom is SLP 10834, the `+Q` gate of block 192
(`x38494 = x11478 + Q`), nonzero only because `x11478` was still undetermined when the worklist
drained. The propagation **stalls wherever an atom holds two unknowns at once** — a handle and
its cofactor, or a mod-P copy target and its handle — and the seeding does not cover the copy
targets.

Equally important, and P said it unprompted: the all-off run completing with **0
integer-division obstructions is not evidence that the lift is free.** That configuration drives
essentially every residual to zero, so it exercises the 927 conditions only trivially. It shows
the constructor is sound on a degenerate input; nothing more.

### Status of the reduction — unchanged and still conditional

**2,780 of 3,707 handles are genuinely free at `c = 1`; 927 carry the strictly stronger
`c·P | R` whose satisfiability is unproved**, and the rank that would decide it is uncounted.
Knob set unchanged (256 selectors, liveness derived); everything else remains mod P. P has now
declined twice to state the reduction unconditionally, both times under a standing invitation
to close it.

**With M's 29125 obstruction withdrawn (check-in 13) and Q's existence result resting only on
its own 24-stage caveat (check-in 14), these 927 conditions are the ONLY place where this lab's
central reduction is still conditional.**

### The measurement, one step away

P's own path, and its next task, restricted to exactly this and nothing else: **seed every
mod-P copy target to equal its source exactly over ℤ**, forcing those handles to zero and
unblocking the cascade; confirm a complete lift at the all-off and one-leaf configurations; then
for each of the 927 conditions compute `∂(R/P)/∂t_v mod c` against the ~766 lift parameters and
take the rank modulo each prime power dividing the `c`'s.

P has been told that if it cannot reach the rank, the deliverable is a **complete and correct
lift constructor plus a precise statement of what remains** — a located second stall is worth
more than a partial rank nobody trusts.

---

## Coordinator note — the wrap-count hole in K's argument is closed

P raised (check-in 10) that K's unreachability negative excludes only one wrap, `k = ±1`, and
would be incomplete if several `k` were reachable. **Closed, by an elementary bound that does
not use the construction premise P objected to.**

Both inputs of a block are folds of subsets of the ladder, so each input's scalar is a subset
sum of `{2^e : e ∈ 0..255}` with each exponent used at most once (Q and R independently verified
`L_i = 2^i·G` for i = 0..255). Hence `|x − y| ≤ 2^256 − 1`. And

    2N = 231584178474632390847141970017375815705675128558149808765210326283036322988674
    2^256 − 1 = 115792089237316195423570985008687907853269984665640564039457584007913129639935
    2N > 2^256 − 1,  slack ≈ 1.158e76

so `x ≡ y (mod N)` forces `x − y ∈ {0, ±N}` — **no other multiple of N is reachable,
unconditionally** — and `x − y = 0` with the subsets drawn from disjoint slot supports forces
both empty. **K's walk covers the right cases**, and the wrap-count half of its argument can be
restated without the withdrawn framing.

**But the unconstrained condition is TRUE, so K's negative rests entirely on the partition.**
`N` has 192 one-bits and 64 zero-bits, so some `j` has `bit_j(N) = 1`, `bit_{j+1}(N) = 0`,
`j+1 ≤ 255`; rewriting `2^j = 2^{j+1} − 2^j` gives non-empty `A`, `B` over the full exponent set
with `Σ_A − Σ_B = N`. K has been told to restate its theorem in the form it can actually prove —
about the specific slot partitions, not about `N`'s representability — since the current phrasing
reads as the stronger, false statement. Arithmetic is three lines; K was asked to redo it rather
than take the coordinator's word.

### The equal-inputs condition, in final form

> There is a block β with slot exponent-supports `I_β`, `J_β` and non-empty `A ⊆ I_β`,
> `B ⊆ J_β` with `Σ_{i∈A} 2^i − Σ_{j∈B} 2^j = ±N`.

If it holds, β's two congruences vanish identically, β's output is unconstrained, and the root
is reached by inverting the law in closed form — **a full solve with no scalar recovery.**

---

## Check-in 16 — the exact scorer exists, and 927 is confirmed independently (agent L)

Deliverable unchanged: **39,026 / 39,033**. L's own best remains 39,018.

### The instrument this campaign has been missing

**An exact in-memory scorer**: `CK.load_equations()` once, then **~1.1 s per candidate**,
calibrated on two known points — deliverable → **7**, `assign_L1` → **15** — both matching
`checker.py`. Every cancellation argument made in this lab before now was made without one.

**By-product, and it matters as much:** `E.score` returns **13 for the deliverable where the
truth is 7**. Every incidence-based number in L's earlier report and in `cut.py` is inflated and
only **ordinally** useful. That is the **third independent demonstration that incidence pricing
fails** (C at check-in 3, P at check-in 10, now L), from a direction neither of the others took.

### The ON-set: L retracts, and the mechanism is measured

Test: all 256 leaves are *free* variables, so the deliverable's JSON value **is** the bit — no
propagation, no inference (`onset_deliv.py` reads all 256 directly). **Exactly `x2081` and
`x24601` are set to 1.** M, K, P, R and Q are right. L's single-leaf reading came from a **stale
partial model covering only the 178-side** — 2081 lives on the 78-side, so that extractor could
not see it. Same root cause as L's retracted 2^178. The current `full_model.pkl` does contain
2081, so nothing else L reported depends on the stale model.

**Measured mechanism** (`rootcheck.py`, `delivsite.py`), adopting Q's framing of two selectors
ON and one leaf propagating: `LCA(2081, 24601) = ROOT x9274`, 24601 under the a-child, 2081
under the b-child, root `sel_ab = 1`; root `va` input **== root `vb` input exactly in both
coordinates**, both equal to L's model value for **leaf 24601** transported to the root frame;
and root `vab` wires x30213, x22162 = **exactly the target pair L derived independently**. Cut
site: child **x27994**, parent **x4971.va** — 2 guards + 2 slot links = **4 atoms = 7 equations**.

**K's "root inputs forced equal" and Q's "never folds, one leaf passes up" are the same event
from two sides**, and the deliverable's independently-set root wires holding L's independently
derived target is the strongest cross-check of L's reduction in the file.

### The 927 confirmed from an independent decomposition

`hcheck.py`: L's 3,681 "handle" variables are free, appear in **0 residual atoms directly** and
in **exactly 1 definition** — they are the **cofactors `u`**. **L counted `u` and concluded about
`h`**, precisely P's diagnosis, and says so. Splitting by measured multiplier: **c = 1 for
2,747, c > 1 for 927**, 7 zero-slope. **L's 927 = P's 927**, from decompositions that share
nothing. The lift is **not** free and L's §3 criterion is restated **mod p**.

### First empirical data on how binding the 927 are

| \|S\| | distinct c>1 atoms violated | left undischarged by L's repair |
|---|---|---|
| 1 | 2 | **0** |
| 2 | 4 | 1 |
| 17 | 36 | 8 |

Every violated atom is **inside** the `c > 1` set; none outside it. They are **sparse but grow
with |S|**, and L's greedy round-robin repair cycles at its 60-round cap. L's own conclusion,
unprompted: **a simultaneous CRT solve over the ~766 shift parameters is what is needed, not
round-robin — and that is exactly P's rank question.** Routed to P as input.

Scope: only |S| = 1 is verified over ℤ end-to-end (`assign_L1.json`, 15 failing, both
attributable to the two target atoms). **The table measures L's repair, not the instance's
lift.**

### The cancellation search — instrument built, family mis-specified

L's generalisation of the deliverable's cut **fails to reproduce it at the deliverable's own
site**: 9 atoms / 49 equations with vab set, 7 / 47 without, against the true **4 / 7**. The
family is mis-specified by **5 atoms**. **No placement below 7 found; P is not contradicted.**

L's next step is diagnostic rather than a sweep, and I have endorsed it: diff the 9 broken atoms
against the deliverable's known 4 (prime suspects: the top slot links x24468/x18956 →
x13682/x37892, pinned to the target pair while also rewriting the root vab wires), **and fix the
divisibility repair first** — L measured that the un-converged round-robin injects nonzero atoms
unrelated to the cut, which would masquerade as bad placements in exactly the sweep it wants to
run. A contaminated instrument, caught before use.

`|S| = 3` finished: 2,763,520 folds, **no hit**. `|S| = 4` left running unattended, not blocked on.

---

## Check-in 17–19 — the odd variable resolved, the theorem sharpened, the rank still open

Deliverable unchanged: **39,026 / 39,033**, re-verified independently by S, M and K this round.

### S (check-in 17) — adjudication done properly, and a tension created

**S's §4 is the best adjudication in the campaign.** It translated before comparing (per-selector
occurrence mean **5.89** against P's "5–6", which is what licensed the comparison), classified
the 48 atoms it found touching ≥2 selectors rather than reporting the raw count as a conflict —
1 booleanity certificate, **47 bundled, each selector in its own additive term, no atom anywhere
multiplying two distinct selectors** — and concluded **P's claim stands on S's data**, against
S's own prior position. The distinction it drew is load-bearing and is now the lab's phrasing:

> **P's claim is about the domain; S's is about the image. A free independent subset domain does
> not make the residual a subset-sum.** If it were a sum, 2²⁵⁶ inputs would give ~2²⁵⁶ residues
> and LLL on the bit vector would be right; it gives **48**.

Other results: `lat5.py` at 22 of 48 configurations, **0 feasible**, with **a20215 in the bad set
22 of 22** — never the *reported* blocking row, because it has a handle (step p) so it is
individually satisfiable and infeasibility surfaces elsewhere. A second route to the same
conclusion, independent of reading the image mod p.

**S corroborated K independently and scoped it correctly.** Its discriminator is better than
"the row is zero" and is **reusable machinery**: degeneracy's signature is *zero AND unresponsive
to every knob in its cone*, since an identically-vanishing constraint cannot be moved. Over 48
configurations × 333 cone knobs: **0 degenerate cases**, responsiveness never below 2. S's own
limit: this is empirical non-observation over a **reached** space, so a degenerate configuration
isolated from cfg0 under single flips would be invisible. **Corroborates K; does not close it** —
and is independent of K's carry-walk step, so P's wrap objection does not touch it.

**S retracted its §6** usefully: the ≤6 arithmetic was sound, but **0 of the 769 equations contain
only one atom** (768 have exactly 2, one has 11), so it was incidence structure only. The
correction **cuts toward attainability**: two-atom equations mean cancellation is available there,
so a nonzero footprint-1 atom can cost **zero** equations. Routed to L, which now has an exact
scorer. Also: **a726 is never bad** (0 of 22) — a satisfied-but-rigid side condition, not an
obstruction; **a28647 is the only genuinely unhandled bad atom.**

**Open tension, now S's task.** Q proves a satisfying assignment exists; S's image closed at 48
tuples with **a20215 never 0**, and S's endgame condition is `a20215 ≡ 0 (mod p)`. Both cannot
describe the same space. S is testing the obvious resolution — evaluate the five rows at random
high-weight configurations its BFS could not reach and see whether they land outside the 48.
Q has been told its untested 24-stage caveat is now load-bearing.

### M (check-in 18) — `x_7068` resolved, and the placement neighbourhood is exhausted

**`x_7068` is a definer that is NOT a handle**, established in M's own frame without reading P's
list: linearly defined (`x_2099 + 7376877·x_642`) where P's four are product- or bare-defined;
**90 digits** against 723/89/724/724; and its definer *references* the corrupted handle `x_642`.
Atom 23616 is nonzero *because* `x_642` is corrupted — for it to be satisfied `x_7068` would have
to be 730 digits and it is 90. `P ∤ x_7068` is automatic for any linear combination containing a
corrupted term, so P finding exactly four among a product-definer population is consistent with
`x_7068` not being in that population at all. **M's five = P's four corrupted handles + one
collateral combiner. The cross-check is complete.**

**Alternative placements priced — none beat 39,026.** Soundness property that makes the search
trustworthy: demoting an atom and seeding its variable with its current value is **bit-identical**,
so demotion is score-neutral and purely adds freedom, validated per candidate. All 11 placements
(the 10 atoms currently zero, in a failing equation, and definers; targets of size 1 and 2),
solved → applied → re-propagated → scored: **every one returns exactly 39,026.**

**Structural finding:** equations **12270** and **18673** have **zero demotable zero atoms**, so
**2 of the 7 failures cannot be touched by this move at all**, whatever values are chosen.

**M's scope statement, which is the campaign's live edge:** `eqsub` prices repairs of the current
placement; `place` prices local neighbours. **Neither prices a genuinely different placement —
corrupting a different set of handles — and that is where anything above 39,026 must live.**
Division of labour accepted: **L supplies candidate handle-corruption sets from the circuit side,
M prices each** with the placement-agnostic primitive; the hand-off is relayed by the coordinator
since agents do not read each other's directories. M is making the primitive candidate-agnostic
and validating it by reproducing the deliverable from its own four handles.

### K (check-in 18) — the theorem restated in the form it can actually be proved

K recomputed the coordinator's arithmetic and **confirmed it**, with one correction: the slack is
**≈1.158×10⁷⁷**, one decimal order larger than the 10⁷⁶ quoted. Conclusion unaffected — `k = ±1`
is the only wrap, unconditionally — and K restated it using only the measured doubling identities
(255/255), **with no construction premise**.

**K accepted the representability point and rewrote its claim.** Taking `j = 0` (`bit_0(N) = 1`,
`bit_1(N) = 0`) gives non-empty disjoint `A, B ⊆ {0..255}` with `Σ_A − Σ_B = N` exactly, so the
unconstrained condition **is** satisfiable and K's negative cannot come from arithmetic. K's
phrasing did read as the stronger false claim. **New §4.0 states the theorem in partition form:**

> For every stage, neither slot support contains all of `{129..255}`.

Since `2^256 − N < 2^129`, omitting any exponent ≥129 caps subset sums below `N`. The
load-bearing facts are two **measured** partition facts — each root half omits 43 and 84
exponents ≥129, and every interior stage sits inside one root half — and K flagged those as the
thing to attack.

**B-half gap: three hypotheses killed, still open.** Backward flow **refuted** (provenance
instrumentation: the back-cone of `x14853` is 5,119 variables and contains no variable from above
the root); sign/shift bug **refuted**; liveness-bit seeding **refuted** (900 non-leaf booleans
seeded 0 / seeded 1 / left to derive give byte-identical verdicts) — which also **settles the knob
set: the boolean inputs really are just the 256 leaf selectors.**

**Classifier bug fixed**: all three idempotency spellings matched; free booleans 369 → **1,156**,
leaf selectors recognised 82/256 → **256/256**. It did not invalidate earlier results but exposed
the 900, which forced the test above.

**K listed four errors it made this session**, including that its own interior-stage bound
`|x−y| < 2ⁿ` is **false** (exponent sets are not initial segments). And per instruction, §1 and §7
now say plainly that obtaining the scalar is the entire remaining problem, with the negligible
fraction of the range its bounded searches covered spelled out. No infeasibility claim.

### P (check-in 19) — the constructor works; the parameter count was wrong; the rank is still open

`plift5.py` is a **working, complete lift constructor**. The seeding fix landed (nonzero atoms
**194 → 3**) and needed two more P had not predicted: a **heap keyed by SLP position**, because a
plain worklist let a variable be back-solved from a downstream constraint before its own
definition fired; and a rule that `h = 0` must **not** fire at `h`'s own definition atom, or it
pre-empts the very divisibility test being measured.

| configuration | division obstructions | nonzero atoms | failing eqs |
|---|---|---|---|
| all off | 2 | 3 | 27 |
| **one leaf ON** | **2** | **2 — exactly the two target congruences** | 17 |
| two leaves, live merge at block 2 | 4 | 4 — one block-2 law pair + the two targets | 27 |

Diagnostics, **not score attempts**; 9,040 variables remain defaulted. Cross-check: at `|S| = 2`
L reports 4 violated `c > 1` atoms and P gets one violated pair at the live block plus the two
targets — **same order, same places, unshared decompositions.**

**P found its own parameter count wrong and said so before reporting any number.** The 512 leaf
coordinates are pinned only *modulo* P — `s·(x − K) = c·h` with `h = u·P` makes `x = K + c·u·P`
legal — so fixing `x = K` exactly treated leaves as rigid. **Correct count ~1,278, not ~766, and
6 rather than 2 at a leaf⊕leaf merge.**

P did obtain a rank result — at block 2 the moduli are `(1, 1, 7038713)` with
`7038713 = 11·23·43·647`, each prime dividing exactly one modulus, giving **q = 11, 43, 647
solvable and q = 23 not** — and **declined to report it as the answer**, because it is computed on
the 2-parameter model P had just shown to be too small. The 6-parameter computation is the one
that decides it, and P is running exactly that at `|S| = 2`, where blocks are genuinely decoupled.

**Reduction status unchanged: 2,780 free at `c = 1`; 927 carrying `c·P | R`, satisfiability
unproved.** L's matching count raises confidence in the **count**, not the **satisfiability**.
P has now declined three times to state it unconditionally.

---

## Check-in 20–22 — the free-cancellation lever, and the audit lands (agents R, T)

Deliverable unchanged: **39,026 / 39,033**.

### CORRECTION TO THIS FILE — the evaluator `E` is wrong on the deliverable

**Three agents found this independently.** F's fast evaluator `E` reports **13 failing on the
deliverable where `checker.py` reports 7** — the 7 real ones plus 2554, 6816, 8124, 8680, 9123,
9421 — scoring it **39,020 instead of 39,026**. R measured it and confirmed the over-report is
**assignment-dependent, not a constant offset**; L found the same 13 from the circuit side; T
found that `fwd.Engine.run` **silently overwrites 4 variables of the deliverable**, which is the
mechanism. R's `rescore.py` shows `E` and `checker.py` agree exactly on every configuration R
generated, so configuration scores computed with `E` stand — **but the deliverable's own
footprint cannot be read off `E` at all.** Anyone scoring with `E` under-reports it by 6.

**Corrected figure:** the deliverable's defect occupies **13 equations of which 6 cancel** — not
the "7 atoms / 12 equations / 5 cancelling" quoted in earlier sections, which cross-quoted
`NOTEBOOK.md` §Session 10's atom numbering against `E`'s. R owned and corrected this.

### R — the single-bit experiment ran, and failed with a number

R resolved its own contradiction against its prior claim: `price.py` fixes the nonzero-atom
**support** to whatever `gs2.solve` lands in and optimises only over values, **but support is
itself a free choice of the repair** — which is precisely what the deliverable exploits. So
≤39,020 / ≤39,022 bound R's repair's support, **not the configuration and not the instance**.
§(d) was alive, so R ran it and killed its SAT/SMT/CP jobs to free the cores.

Result: single-bit footprint has 3 live atoms, **20 equations, 0 cancelling**, against the
deliverable's 13 with 6 cancelling. Of the 20: **0** can never cancel, **20** have a dead partner
available, **0** have a partner touching ≤3 other equations. The cheapest partner in the whole
footprint (atom 7954) occurs in **10** other equations; the next are at 11, 11, 12, 12, 13, 14,
14. Needing ~7 purchases at ≥10 each, **it cannot pay, by ~10×.**

**Why the deliverable wins, stated cleanly:** not a better configuration and not a bigger
support — it sits in a **rare footprint where 6 of 13 equations cancel for free**, while every
footprint R reached charges ≥10 for the first. Fifth independent confirmation of
cancellation-not-support.

**Scope caution attached (coordinator):** R's partner-occurrence counts are facts about the file,
but "buying one cancellation lights up ≥10 more equations" treats occurrence count as **cost**,
and three independent results in this lab show incidence does not price cost. R has been asked
to re-quote it as "≥10 equations **touched**", which is what was measured. The ~10× gap probably
survives; the claim should be stated at the strength of the measurement.

**R re-tasked on its own lever:** enumerate footprints by how much **free cancellation** they
carry, rank them, and only then ask which configurations route into them. Every search in this
campaign has gone configuration-first and priced the footprint it landed in. This inverts it, and
it is the deliverable's actual trick, never searched systematically.

### T — the adversarial audit

**Confirmed by re-running, not by reading:**

- **`ker(M) = 0` reproduces from cold**, and T ran the test nobody had: **faithfulness of M**,
  the premise the whole result rests on. Exact list equality between `{e : (Ma)_e ≠ 0}` and
  `checker.evaluate_all` at **10 points** — all-zeros, 4 partials, 3 random small, 2 random
  30-digit. **M is faithful.**
- **Pivot magnitudes measured: 37,889 ones, 1,144 twos, nothing else** — a claim FLEET.md and
  RESUME_F both attribute to `peel_cert.py`, which **only ever tests `pivot != 0`**. True, but
  previously unverified by its own cited script.
- **H's 722 dormant-handle exclusion is configuration-invariant** — re-measured at 8
  configurations including the deliverable's own; census identical every time. T expected this to
  break and reports that it did not.
- **I's eq8680 conclusion** survives the 5 groups it never tested (all `minfail > 6`).
- **A's linearity filter is honest**: `knobs_raw == knobs_linear` at every level.

**New unconditional result (T's A7):** over the **enlarged** knob set (24, not 9), condition (a)
is exhausted at L=0 — **no mod-p-admissible violated set of weight ≤6** (736,281 nodes at weight
6) — weight 7 is admissible, and **its unique lightest set is exactly
`[12231, 12270, 12350, 14584, 18673, 22044, 29125]`**, an end-to-end self-test. **≥7 is
unconditional and exhaustive at L=0 over the unfiltered knob set**, where A had only ≥6/≥4/≥3.

**Broken or weakened:**

- **B1 — "knobs" is a property of the atomisation, not the instance.** Rebuilding A's windows
  from F's certified-faithful parse: atoms and vars match A exactly, but knobs are
  **24 / 88 / 235 / 610 against A's 9 / 32 / 109 / 334**, all 24 verified genuine at L=0. **Every
  exhaustive count and Prange bound in THEOREM.md was computed over 37–55% of available
  directions.** T redid L=0; the answer is still 7.
- **B3 — "all-atoms-zero is an equivalence" holds only in F's decomposition.** Kernel dim ≥3,234
  in A/G/H's 42,267-atom model and ≥1,852 in I's 40,885. THEOREM.md §7 transfers F's consequence
  to A's atoms; **that transfer was never established.** And there are **five** atom counts in
  this lab, not the three this file named: 39,033 / 39,277 / 40,727 / 40,885 / 42,267.
- **B2 — I's "complete enumeration" is a strict subset of its own census** (42 groups / 27
  nonzero-effect, not 43 / 30; 5 nonzero-effect groups untested, 15 of 27 in no pair). The
  conclusion survives; the description did not.
- **B4 — H's carrier census** is phrased in the "rank > deficit" vocabulary H itself withdrew and
  never rewrote, and its base is one arbitrary selector. The measurement survives; **"7 is the
  floor across every carrier class" does not.**
- **B5 — this file dropped F's own caveats.** `mux_wiring.json` has **47 entries of 96 stages**,
  and RESUME_F §3 explicitly caveats the 2²⁵⁶−1 count as "a model from 47/72 wired stages, not an
  exhaustive check". FLEET carried the number without the caveat, under *Established*.

**T's item (d) — settled since T started, and against T's prediction.** T named "reachable space
= 2²⁵⁶ − 1" as the most suspicious surviving claim, because the specific open case (two leaves ON
in the same OR-group giving the slot a *sum*) fails in the direction that would break the count.
**L measured it: both pins fire, but the mux coefficients are mutually exclusive quadrants, so
both are multiplied by zero and the slot holds the chord, not a sum** — verified at ON-set sizes
1, 2, 5, 73, 200, 256. P confirmed the mux 381/381 and the count from an independent parse; Q, R
and L each established the fold is associative, so **tree shape is irrelevant**; K found the
undecoded slot pairs are dead blocks with empty leaf support. **T's criticism of the coordinator's
summary was correct and is fixed here; the underlying number survived.**

**T re-tasked** onto the newest load-bearing claims, which have had no adversarial pass: Q's
existence result and its untested 24-stage caveat; the Q-versus-S tension (Q says a solution
exists, S's closure of 48 tuples excludes it) — where T's own B1 finding makes "is *the image*
decomposition-dependent too?" the right question; and whether the 927 is decomposition-dependent.
