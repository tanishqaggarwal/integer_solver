# RESUME_T — agent T, the campaign auditor.  Self-contained.

Mandate: audit the claims that CURRENTLY survive, by re-running them, looking for the fleet's
documented failure mode — *a property computed over a filtered knob set, or at an unstated
selector configuration, reported as a property of the instance.*

All my files are in `solve_lab/agentT_work/`.  I modified nothing outside it and ran no git.

---------------------------------------------------------------------------------------------
## A. CONFIRMED BY RE-RUNNING (these gained weight)

### A1. Deliverable 39,026 / 39,033 — CONFIRMED
`python3 solve_lab/checker.py solve_lab/best/new_instance_partial_39026.json`
-> `satisfied 39026/39033 (7 failing)`, failing `[12231,12270,12350,14584,18673,22044,29125]`.
No file anywhere in `solve_lab/` claims a checker-verified score above it.

### A2. ker(M) = 0 — CONFIRMED, and for the first time its *premise* was checked
`python3 agentF_work/peel_cert.py` reproduces from cold:
`certificate verified: True, atoms forced 39033 of 39033, rank(M)=39033, dim ker(M)=0`.
The peel is a genuine triangularity certificate (each step's row has one not-yet-zeroed column,
with a nonzero pivot), so `M a = 0 => a = 0` by induction.  Valid over Q and Z.

Two things nobody had checked:

* **Pivot magnitudes.**  RESUME_F and FLEET.md say "all pivots are ±1 or ±2, none divisible by
  any odd prime".  `peel_cert.py` **does not test this** — it only tests `pivot != 0`.
  I measured them (`t_pivots`, inline): **37,889 pivots = 1, 1,144 = 2, nothing else.**
  The claim is TRUE; it had simply never been verified by the script cited as verifying it.
  (Immaterial to ker(M)=0 over Q/Z — any nonzero pivot suffices — it only supports the
  "char != 2" adornment.)

* **Faithfulness of M.**  ker(M)=0 is worthless unless `eq_e = 0 <=> (M a)_e = 0` for the real
  instance.  Nobody had tested it.  `t_faithful.py` / `t_faithful2.py`: evaluate all 39,033
  atoms at a given assignment (**no forward re-derivation** — F's own `fwd.Engine.run` silently
  overwrites 4 variables of the deliverable and reports 13 failing instead of 7, so that route
  is not a test), then compare `{e : (M a)_e != 0}` with `checker.load_equations` /
  `checker.evaluate_all`.  **Exact list equality at 10 points:**
  all-zeros (11,684 failing), the 4 saved partials (7 / 9 / 12 / 20), 3 random small-integer
  points (39,005-39,008), 2 random 30-digit points (39,033).
  **M is faithful.  all-atoms-zero <=> full solve is CONFIRMED for F's decomposition.**
  Bonus: at the deliverable exactly **7 atoms** are nonzero and they yield exactly 7 nonzero rows.

### A3. Agent H's "dormant handle" exclusion — CONFIRMED (I expected this one to break)
`handles.py` drops 722 of 1,865 one-check-atom free inputs as "dormant" (zero effect), measuring
the effect at the **all-zero** state `St({})`; `hsweep.py` then prices the survivors at a
**different** state (the closed 1-selector state, `BITS['A'][0]` = x_47).  Textbook setup for the
fleet's failure mode.  `t_dormant.py` / `t_dormant2.py` re-measure at 8 configurations —
all-zero, the closed hsweep base (score 39,018), three other single selectors, the deliverable's
own free-input values (score 39,020 in H's frame), and three random 8-selector configurations.
**The census is identical at every one: `{dormant: 722, p: 1143, other: 0}`; zero handles change
class.**  The solo-handle class is configuration-invariant.  H's exclusion is safe.

### A4. Agent I's eq8680 hunt — its GAP is now closed, and I's conclusion survives
See B2 for the gap.  `t_hunt_gap.py` runs the 5 nonzero-effect knob groups I never tested,
through **I's own** `eq8680.build` / `minfail_bnb`:
```
X21279  |E|=47 knobs=10  266s  minfail > 6      X23754  |E|=15 knobs=10   4s  minfail > 6
X34600  |E|=50 knobs= 9   20s  minfail > 6      X3629   |E|=53 knobs= 9  50s  minfail > 6
X8976   |E|=47 knobs= 9   37s  minfail > 6
+ all five paired with the cheapest compensator X19964: 4 of 5 minfail > 6, X19964+X21279 TIMEOUT at 600 s
```
No hit.  I's headline stands; its coverage claim did not (B2).

### A5. Agent I's 66-pair sweep — CONFIRMED COMPLETE
FLEET's live tasking still lists "I — finish the 66 pairs".  They finished: `hunt.log` has
22 singles + 66 pairs = 88 rows and terminates with `done`; `grep BEATS hunt.log` is empty.

### A6. Agent A's linearity filter is a NO-OP on A's own parse — CONFIRMED
`regsolve2.pick_knobs` *is* literally a loop that discards knobs until every monomial has <=1
knob, which looks like it contradicts THEOREM.md §1.2 ("no knob is ever discarded to preserve
linearity").  It does not: running A's own `eqwin.py` gives `knobs_raw == knobs_linear` at every
level 0-6 (9, 15, 32, 58, 80, 107, 109) and `nonlinear_atoms = 0`.  Nothing is discarded.
**THEOREM.md §1.2 is accurate.**  The A/B split is also honest: Theorem A is genuinely the
exhaustive part, Theorem B genuinely carries the Prange miss probability, and §6 states the
scope limitation (80% of window variables excluded at every depth) rather than hiding it.

### A7. NEW unconditional result — A's bound of 7 is *tight and exhaustive* at L=0
(this is the audit's one positive contribution to the mathematics)
See B1: the same L=0 window carries **24** genuine knobs in F's decomposition, not A's 9.
`t_conda.py` redoes A's binding condition (a) — minimum-weight mod-p-admissible violated set —
over all 24, exhaustively:
```
L=0 rows=30 knobs=24 rank_p(N)=12 left-kernel dim w=18
 weight 0..6 : NONE admissible   (1 / 31 / 465 / 4,495 / 31,465 / 169,911 / 736,281 nodes)
 weight 7    : ADMISSIBLE, D = [12231,12270,12350,14584,18673,22044,29125]
```
The weight-7 hit is **exactly the deliverable's failing set** — an end-to-end self-test of my
implementation.  Consequence, with no probabilistic content and no filtered knob set:
> **Every integer assignment agreeing with the deliverable outside K_0 (all 24 knobs free)
> fails at least 7 equations, and 7 is attained.**
A's own Theorem A gives only >=6 / >=4 / >=3 (L=2/6/16) unconditionally; >=7 was Theorem B,
probabilistic.  At L=0 scope, >=7 is now unconditional and exhaustive over the *un*filtered
knob set.  Supporting structure (`t_struct.py`): at L=0/2/4 with the enlarged knob sets,
rank_Q(N) = |K| exactly, the system is Q-consistent and mod-p INCONSISTENT — so A's Lemma
hypotheses all survive the enlargement too.

---------------------------------------------------------------------------------------------
## B. BROKEN OR WEAKENED, with the exact exclusion

### B1. **A's knob set is decomposition-dependent, and it is the smaller choice.**  (the main find)
A's `K_L` = "variables all of whose atoms lie in `A_L`" is computed in A's **42,267-atom**
decomposition.  `t_window.py` rebuilds the identical windows from F's **39,033-atom**
decomposition — the one I certified faithful in A2:

| L | atoms (T / A) | vars (T / A) | **knobs (T / A)** |
|---|---|---|---|
| 0 | 24 / 24 | 56 / 56 | **24 / 9** |
| 2 | 88 / 88 | 202 / — | **88 / 32** |
| 6 | 235 / 235 | 537 / 537 | **235 / 109** |
| 16 | 611 / 611 | 1,312 / — | **610 / 334** |

Atoms and vars agree exactly; **knobs do not**.  The extra knobs are not an artefact:
`t_knobs.py` perturbs each of the 24 at L=0 against the *raw equation values* and finds
**0 escape R_L and 0 make any row non-affine** — all 24 are genuine zero-collateral, exactly
affine directions.  A's atom set is *finer*, so a variable meets more atoms, so more variables
get disqualified as knobs.  Nothing is wrong with A's arithmetic; the point is that
**"knobs" is not a property of the instance, it is a property of the atomisation** — the same
failure mode as the five retracted barriers, one level up.
Consequence: **every exhaustive count, Prange probability and "lightest weight the filter
admits" in THEOREM.md was computed over 37-55% of the available directions**, and none of them
covers the enlarged system.  §6.1's "excluded" column is numerically right but attributes the
whole exclusion to window membership.
Mitigation: at L=0 I redid it (A7) and the answer is still 7.  L>=2 is not redone — C(104,6)
is out of reach in Python.

### B2. **Agent I's "complete enumeration" for eq8680 is a strict subset of its own census.**
RESUME_I: "the complete candidate set is ... **43 groups**, of which **30** have nonzero net
effect ... That census is **an enumeration, not a sample**", then "**ALL 22** single knob-groups
give minfail > 6".  Parsing `eq8680.log` directly: **42 groups, 27 with nonzero net effect** —
neither number matches the prose.  `hunt.py`'s `CANDS` is a hardcoded list of 22, and its pairs
are `itertools.combinations(CANDS[:12], 2)`.  So:
* **5 nonzero-effect groups were never tested at all**: X21279, X23754, X34600, X3629, X8976
  (X21279's net effect on eq8680 is a ~1,000-digit number — the largest in the census);
* **15 of the 27** appear in no pair; the "66 pairs" are C(12,2) over the 12 cheapest.
The word "complete" describes the *census*; the *test* is a hand-picked subset, and the prose
does not distinguish them.  I ran the 5 (A4) — no hit — so the conclusion survives, the
description of it did not.  Also: I's "a knob is a variable all of whose atoms are in the
support" makes the candidate list complete only for moves that keep every atom outside the
support at exactly zero; that restriction is not stated where "complete" is claimed.

### B3. **"all-atoms-zero is an equivalence" is TRUE only in F's decomposition.**  FLEET states
it unqualified, under "What still stands **unconditionally**", with only a parenthetical
"agents' matrices differ in atom count ... compare decompositions before comparing kernel
dimensions".  That parenthetical is doing load-bearing work it cannot do:
* In **A's / G's / H's** 42,267-atom decomposition the incidence matrix is 39,033 x 42,267 and
  its kernel has dimension >= 3,234; in **I's** 40,885 it is >= 1,852.  There, all-atoms-zero
  is strictly *stronger* than solving the instance.  A's own §6.3 proves this from the other
  side: 3,235 atoms occur in exactly one equation, so a single equation can be satisfied with
  two of its atoms nonzero and cancelling.  I documents the same thing (926 pairs).
* THEOREM.md §7 nevertheless writes "**Consequence:** any assignment satisfying all 39,033
  equations must make **every atom** exactly zero" while working in the 42,267-atom model.
  That transfer is **not** established by F's computation.  (It happens not to damage A's
  theorem, which is about knob moves, not atom vectors — but the sentence is false as written.)
* FLEET names three atom counts.  There are **five** in the lab: 39,033 (F, K), 39,277 (P),
  40,727 (E), 40,885 (I), 42,267 (A, G, H).  Only I reconciles its own arithmetically
  (40,885 - 39,033 = 1,852 = 2 x 926 pairs); A explicitly declines to; E, G, H, P never mention it.

### B4. **H's carrier census is not "complete", and its stated reason is a withdrawn criterion.**
RESUME_H §9: "**Carrier census, complete** ... The witness is the only carrier **in the
instance** whose realizable knob-image rank exceeds its balance deficit.  7 is the floor across
every carrier class."  Three problems, in decreasing severity:
1. The criterion "rank > deficit" is **withdrawn by H itself** in §11 ("the rational rank of the
   knob image is not the binding quantity").  §9's headline sentence is stated in exactly that
   retracted vocabulary and was never rewritten.  What survives is the *measurement* (best score
   over 1,147 priced carriers = 39,017), not the inference "7 is the floor".
2. `hsweep.py`'s base is **one** state — `close_trace(St({}).set_free({BITS['A'][0]: 1}))`, a
   single arbitrary selector, scoring 39,018.  H's own §10 re-ran the *pin* sweep from 2- and
   3-selector bases but **not** the 1,143-handle sweep.  So "in the instance" means "at x_47 on".
3. The 722 excluded handles.  This one I could not break (A3) — it is the sole exclusion in the
   lab that I re-measured and found configuration-invariant.
Also: FLEET's live tasking says "price all **~30** cascade pin atoms"; RESUME_H reports **20**.

### B5. **FLEET's "Established" tree facts drop F's own caveats.**
FLEET FINAL POSITION: "The instance is a 96-stage binary tree of depth 6 ... **Reachable space =
2^256 - 1 non-empty leaf subsets**", listed under *Established*.  RESUME_F §3 states the caveat
in its own words and FLEET does not carry it: "*this is a model from 47/72 wired stages plus the
verified quadrant structure, not an exhaustive check.  What remains genuinely open is what
happens when two leaves in the SAME OR-group are both ON ... That case is not covered by the
count.*"  Measured on disk: `mux_wiring.json` = **47** entries, `stage_roles.json` = 72,
`stage_profile.json` = 96.  So the decode covers **47 of 96 stages (49%)**; 56 stages have an
undecoded slot pair and the 24 leaf-adjacent stages' literals are unresolved.  A section headed
"**THE INSTANCE, DECODED** — read this before anything else" is 49% decoded.
What IS solid: `stage_law2.log` on disk backs "72 of 72 stages, same universal K, zero
exceptions" verbatim, correctly scoped to the 72 (the 24 leaf-adjacent stages were *skipped*,
so LOG.md §24's "So all 96 stages run ONE law" is an extrapolation over 25% of the tree —
FLEET's phrasing, "all 72 fully-determined stages", is the accurate one).

---------------------------------------------------------------------------------------------
## C. NOT CHECKED, AND WHY
* **A's L=2 / L=6 / L=16 exhaustive counts over the enlarged knob set.**  C(104,6) = 1.4e9 rows
  at L=2 alone; my L=0 weight-7 pass already cost 789 s for C(30,7).  Needs C/Flint, not Python.
* **A's Prange miss probabilities.**  Randomised; re-running samples the same distribution and
  cannot falsify them.  The 1.6e-4 at L=16 is honestly reported as the weak one.
* **F's 72/72 stage law and the invertibility (200/200 triples).**  Re-derivation is a ~322 s
  run plus a decode I would have to trust; I verified the *logs* match the prose instead.
* **E's saturation law** ("2 or 3 bits of one class give exactly the delta of ONE bit").  The
  number of measurements behind it is never stated in RESUME_E; sample size unrecoverable
  without re-running `channels.py`.  Flagged, not broken.
* **The one X19964+X21279 pair that timed out** at 600 s in A4.  Re-run with a longer limit.

---------------------------------------------------------------------------------------------
## D. THE SINGLE MOST SUSPICIOUS SURVIVING CLAIM
**FLEET's "Reachable space = 2^256 - 1 non-empty leaf subsets", and the whole
"THE INSTANCE, DECODED" section, presented as Established.**  It is the claim the fleet's
designated next experiment is built on; it is derived from a 49%-complete decode; its author
attached an explicit caveat that the coordinator's summary dropped; and the specific open case
(two leaves ON in the same OR-group) is precisely the one where a slot would see a *sum* rather
than a single value — i.e. where the "at most 2^256 - 1 reachable folds" count would fail in the
direction that matters.  F retracted an infeasibility argument once already for the analogous
reason.  **Finish the 56 undecoded slot pairs and settle the same-group double-leaf case before
any meet-in-the-middle is built on the 2^256 - 1 figure.**

Runner-up: **B1**, because it is live — the fleet is actively told to read THEOREM.md, and its
knob set is 37-55% of the available one at every depth above L=0.

---------------------------------------------------------------------------------------------
## E. FILES
`t_faithful.py` `t_faithful2.py` (M faithfulness, 10 points) · `t_dormant.py` `t_dormant2.py`
(H's dormant exclusion, 8 configurations) · `t_hunt_gap.py` (I's 5 untested groups, uses I's own
solver) · `t_window.py` (A's windows from F's parse) · `t_knobs.py` -> `window_L0.json` (the 24
knobs verified genuine + affine) · `t_struct.py` (rank/consistency on the enlarged sets) ·
`t_conda.py` -> `t_conda_L0.log` (**A7**, the exhaustive weight-<=6 result) · `t_minfail.py`
(integer BnB; superseded by t_conda and stopped).
Reproduce the headline: `cd solve_lab/agentT_work && python3 t_faithful2.py && python3 t_conda.py 0 8`.

=============================================================================================
# SECOND PASS — audit of agent Q's existence chain and the ROUTING LAYER
(coordinator check-in 22/23; item 2 dropped — agent S resolved the Q-vs-S tension itself)

## F. Q's LADDER AND ARITHMETIC — CONFIRMED, AND STRONGER THAN Q CLAIMED  (`t_ladder.py`)
Q's step 5 ("subset sums realise every k; N < 2^256; therefore T is hit") is airtight only if the
leaf set is exactly {2^i G : i = 0..255}, every exponent present and distinct.  Drop one exponent
and subset sums cover 2^255 of 2^256 and T surviving is a coin flip, not a theorem.  Measured:
* **all 249 checkable consecutive doublings `L_{i+1} == 2*L_i` verify exactly, 0 failures**;
* the 253 decoded ladder points are **253 distinct** points, all on the cubic;
* `L_i == 2^i G` for **all 253**, no mismatch;
* the 3 "missing" exponents are **not inferred** — vars x18184, x22579, x33434 are decoded leaves
  and equal `2^51 G`, `2^176 G`, `2^41 G` exactly.  So all 256 exponents are backed by a decoded
  leaf.  (RESUME_Q §1(d)'s "4 chains linked by one missing doubling each" understates its own
  result; the coordinator reports Q has since re-derived this as a single 256-chain.  Agreed.)
* `N*G == O`; N prime-ish; `N-1 <= 2^256-1`, so [1,2^256-1] does cover every residue mod N.
**Steps 4 and 5 of the existence chain survive.**  Also: all 256 selectors are genuine booleans
(each carries an `x*(x-1)` atom) and **all 256 are free variables** — so an arbitrary ON-set is at
least *syntactically* assignable.  That premise survives too.

## G. THE ROUTING LAYER — 'liveness is determined by the selectors' is ASSUMED, and as stated it
##    is FALSE in agent E's parse  (`t_sel.py`, `t_live.py`)
Q verifies the law each gadget enforces *as a function of its four input coordinate wires*.  What
carries a leaf to those wires is the selector/mux layer.  Measured:
* The two wire vars agent Q associates with each leaf (`qleaf[sel][2:4]`) are **FREE variables —
  253/253 for w1 and 253/253 for w2**.  They are not computed from the selector.
* Flipping each of the 256 selectors, at **four different bases** (deliverable; deliverable with
  all selectors off; agent S's triple8_seed; all-zero): every selector is a live knob (**0 inert**,
  median ~50 wires move), but the number that **make their own leaf's coordinate appear anywhere
  in the circuit is 0 of 256, at every one of the four bases.**
* Directly: 12 sampled leaves, selector set to 1 alone, from two bases — **0/12 arrived**, on their
  own wires or anywhere else.
* The deliverable has exactly 2 live leaves (72, 235) and **still has exactly those 2 with all 256
  selectors forced to 0** — its leaf values are sitting on free variables it assigns, not routed.
So forward evaluation from the selectors does **not** realise an ON-set.  The routing is a
*constraint* (the pins), not a propagation.  This is measurement, at four bases, not one.

### G1. A hypothesis of mine that BROKE under its own test — recorded because it did  (`t_pins.py`)
The pin for exp 16 is `x_32872 * (x_34615 - X_leaf) - 4949965 * x_5923`, i.e. atom = 0 gives
`w1 = X_leaf + 4949965*x_5923`, not `w1 = X_leaf`.  I expected this to be an additive knob that
un-pins the leaves and makes Q's ladder only the slack=0 section of a larger reachable set.
Census over all 256 leaves: **588 pin atoms — 82 exact `sel*(w - CONST)`, 506 with a slack term —
and 0 of the 506 slack variables are free.**  They are all defined.  **The leaves are pinned; my
crack is not one.**  Q's ladder is not weakened by the pin shape.

### G2. WHAT Q's END-TO-END TEST WILL MISS (the coordinator's actual question)
Q is running "set selectors in the real DAG at random subsets, check the root equals the
independently computed fold".  Three things that test cannot see:
1. **It cannot be driven from the selectors alone, so it risks assuming its conclusion.** Because
   w1/w2 are free (G), the harness must *choose* leaf-wire values.  If it sets them to the leaf
   constants for ON leaves and 0 for OFF leaves, it has hard-coded "OFF leaf = identity" — which
   is exactly the proposition under test.  The non-circular version must **solve** w1/w2 from the
   pin atoms and check consistency, never assign them.
2. **Random subsets never enter the regime that matters.**  Weight is Binomial(256, 1/2): within
   [104,152] with probability > 99.7%.  So random subsets never test weight <= 6 or >= 250.  At
   weight ~128 essentially every gadget sees two live inputs and takes the **chord** branch; the
   **pass-through branch (exactly one input live)** is barely exercised at depth — and pass-through
   is precisely what makes an OFF leaf behave as the group identity, the linchpin of "fold = group
   sum".  Note the sting: weight <= 7 is exactly the regime of Q's own lottery-ticket sweeps
   (`lowwt.py`, `wt7.py`), which computed the fold **in Q's group model** and never checked the
   circuit agrees there.  If routing misbehaves at low weight, those sweeps' "clean miss" verdicts
   are not evidence.  **Test low weight explicitly — 1, 2, 3, 5, 7 — not just random.**
3. **The degenerate/doubling branch.**  The chord `l=(b_y-a_y)/(b_x-a_x)` is undefined when the two
   live inputs are equal; the group law needs the tangent there.  Q's associativity check used
   *random* leaf triples and never hits it; agent P's "degenerate family" is exactly this case and
   is known to occur.  Distinct leaves cannot collide (`2^i G = +-2^j G` is impossible for i != j,
   N prime), but **intermediate fold values can**, and a random-subset test hits that with
   negligible probability.  Enumerate collisions on the actual fold tree instead of sampling.
4. Related, and unresolved by me: with the selector OFF, `x_16886 = 1 - sel` and the atom
   `x_16886 * x_17479` forces that leaf's **y-wire to 0**.  The curve has prime (odd) order N, so
   it has **no 2-torsion** — `(w1, 0)` is not a curve point.  An OFF leaf is therefore *not* the
   identity as a curve point; identity behaviour has to come from the mux coefficients, not from
   feeding a point into the chord law.  I did not settle whether it does.  **This is the specific
   thing to check**, and it is the same object L measured (both pins fire, coefficients are
   mutually exclusive quadrants) in a different model.

## H. NOT DONE THIS PASS, STATED PLAINLY
* **Item 3 (is the 927 decomposition-dependent?) — not started.**  P and L agreeing from unshared
  decompositions is real evidence, but per B1 a count matching across models is not the same as a
  count being intrinsic; the test is to rebuild 927 in F's certified-faithful 39,033-atom parse.
* **Agent L's "cancellation is a value property, not a support property" — not audited.**  Still
  the basis on which M is pricing 378 candidates.  The obvious first check is whether the 12
  cofactor/handle variables L names are free, and whether setting them to the deliverable's values
  in L's own reconstruction moves 13 -> 7.
* I did **not** verify Q's Schwartz-Zippel gadget census (383/383, 89/89 leaf-adjacent) — reported
  by the coordinator, not present in agentQ_work when I read it.
* `t_bfs_audit.py` (item 2) reached only generation 1 before item 2 was dropped; its depth-1 finding
  (257 assignments -> 7 distinct 5-tuples, quotient valid on the 6 reps tested) is consistent with
  S's own later self-correction and should not be cited for anything more.

## I. NEW FILES
`t_ladder.py` (Q's ladder, doublings, order) · `t_sel.py` (256 selectors x 4 bases, liveness) ·
`t_live.py` (leaf wires free; leaves do not arrive) · `t_pins.py` (588 pin atoms; slack census —
the hypothesis that broke) · `t_fold.py` / `t_fold2.py` (fold-vs-group-sum probe + the control that
reproduces Q's §5b: C1 on 4 wires, group sum on 0) · `t_spaces.py` (Q and S are on the identical
256 booleans — the "different spaces" hypothesis is refuted) · `t_bfs_audit.py` (item 2, dropped).

=============================================================================================
# THIRD PASS — audit of agent L's cancellation result (coordinator check-ins 33 / 37)

## J. ACCEPTED CORRECTION TO MY OWN SECOND PASS (section G)
Q's symbolic solve of the mux layer explains my measurement: a leaf pin is `sel*(w - C) - z`,
so the coordinate lands on the wire only once `z` is separately forced to 0.  **Routing IS
determined — by a simultaneous system rather than by propagation.**  My 0/256 and 0/12 numbers
stand as measurements, but they measure the weakness of unit propagation, not an absence of
determination.  **"Liveness is not determined by the selectors" is NOT established and I withdraw
it**; the correct statement is the narrower one, "forward evaluation from the selectors does not
realise an ON-set".  My §G2 item 4 (OFF leaf's y-wire = 0, and `(w,0)` is not a curve point) is
also resolved by Q: the identity value is `(0,0)` and only ever passes through, never entering a
chord, because `cC = ab = 0` whenever a child is dead.  §G2 items 1-3 (the design objections) are
unaffected.

## K. L'S CANCELLATION RESULT — CONFIRMED, from an independent direction  (`t_cancel.py`)
Tested from the DELIVERABLE side rather than through L's constructor, so the result does not
inherit the un-converged divisibility repair L itself flags.  Zero L's 12 cofactor vars in the
deliverable and measure both the exact score and the nonzero-atom support (support read in F's
decomposition, certified faithful in T2):
```
deliverable as given            7 nonzero atoms   FAILING  7  [12231,12270,12350,14584,18673,22044,29125]
same, L's 12 cofactors zeroed   7 nonzero atoms   FAILING 12  [+2554,6816,8124,9123,9421]
support IDENTICAL: True
```
* **All 12 are free variables** (12 of 12, agent E's parse).
* **The support is byte-identical and the cost differs by 5.**
> **"Cancellation is a value property, not a support property" is ESTABLISHED.**  The search
> really is site x handle-values, and M's premise is sound.

### K1. Two corrections to the numbers M should be using
1. **The gap is 5, not 6, and the far side is 12, not 13.**  From the deliverable, zeroing the
   12 gives **12** failing, not 13.  L's 13 is its *own* build2's score; **one of L's 13 is not
   explained by the 12 cofactors** and is almost certainly the un-converged repair L flagged.
   Price against 7 -> 12.
2. **Eight of the twelve do nothing — the list of 12 is really a list of 4.**  Zeroing each alone:
```
   x1329  +3     x9413  +4     x10903 +3     x17325 +4        <- the real cancellation knobs
   x105 x3387 x5081 x5676 x11436 x14393 x14768 x22820  ->  +0 each
```
   The 8 no-ops are **already 0 in the deliverable**, so "the deliverable sets them to specific
   nonzero integers" is false for two thirds of the list.  Only x1329, x9413, x10903, x17325 are
   nonzero there.  The cancellation degree of freedom is **4-dimensional, not 12**.

### K2. **M IS PRICING THE WRONG VARIABLES** — flag this before the lattice solve  (`t_cofactor.py`)
The coordinator reports M pricing an exact lattice target on **x642 and x28730**, "two of the
twelve cofactor variables in L's claim".  Neither is in L's twelve, and neither is a cofactor:
```
   x642    free=False   occurs in 2 atoms      x17325 (its cofactor)  free=True  1 atom
   x28730  free=False   occurs in 2 atoms      x9413  (its cofactor)  free=True  1 atom
   x31864  free=False   occurs in 2 atoms      x10903 (its cofactor)  free=True  1 atom
   x29854  free=False   occurs in 2 atoms      x1329  (its cofactor)  free=True  1 atom
```
x642/x28730 are the **P-multiples h** — defined wires, each appearing in two atoms (its own
definition and the guard above it).  The free knobs are the **cofactors u**: x17325 and x9413.
L's own §6c table states the pairing correctly (`h=x28730 u=x9413`, `h=x642 u=x17325`); the
conflation is downstream of L.  **A lattice solved over x642/x28730 as if they were free is
solving over the wrong coordinates** — and note they are exactly the two whose cofactors carry
the largest single-variable effect (+4 each), so the error is in the most load-bearing place.

### K3. THE PREMISE UNDER THE 15-ATOM FILTER — CONFIRMED ACROSS ALL 3,681  (`t_cofactor.py`)
L's criterion `equation e contains atom a  <=>  u_a in vars(e)` rests on "every residual atom has
exactly one free cofactor u, occurring nowhere else".  Checked in F's certified-faithful parse
against `checker.load_equations()`'s own varsets:
```
   cofactors that are free variables            : 3681 / 3681
   cofactors occurring in exactly ONE atom      : 3681 / 3681   (violations: 0)
   cofactors with eqs(u) == eqs(atom_u) exactly : 3681 / 3681   (mismatches: 0)
```
**The criterion is sound, and so is the "of 3,681 atoms exactly 15 are incident" filter that M is
enumerating against.**  This is the one place in this pass where a premise I expected to be soft
held completely.

### K4. THIRD CALIBRATION POINT for L's exact in-memory scorer
L calibrated on two points (deliverable -> 7, assign_L1 -> 15).  A third, verified with the
`checker.py` CLI, not just in memory:
> **deliverable with x105 x1329 x3387 x5081 x5676 x9413 x10903 x11436 x14393 x14768 x17325 x22820
> removed -> `satisfied 39021/39033 (12 failing)`,
> failing `[2554, 6816, 8124, 9123, 9421, 12231, 12270, 12350, 14584, 18673, 22044, 29125]`.**
My in-memory scorer (`checker.load_equations` + `evaluate_all` in-process) agreed exactly.

## L. STILL NOT DONE
* **Item 3 (is the 927 decomposition-dependent).**  Not started; ranks below the above.
* **Whether tuning the 4 real cofactors can beat 7.**  That is M's search, not an audit; but note
  the space is 4-dimensional over free vars whose defining atoms the deliverable already breaks,
  so x642/x28730/x31864/x29854 are *also* effectively assignable in a partial assignment — the
  live space is larger than 4 and M should establish its true dimension before solving a lattice.
* Q's Schwartz-Zippel census artifacts were still absent when I last looked.

## M. NEW FILES (third pass)
`t_cancel.py` (L's 12, support + score, one-at-a-time) · `t_cofactor.py` (3,681-cofactor premise;
h-vs-u).  Reproduce: `cd solve_lab/agentT_work && python3 t_cancel.py && python3 t_cofactor.py`.

=============================================================================================
# FOURTH PASS — item 3 (the 927) and item 2 (K's withdrawal)   [coordinator check-in 47]

## N. THE 927 IS INTRINSIC — CONFIRMED, and robust to a 3.5x change of decomposition
`t_927.py`, `t_927b.py`.  Rebuilt in F's 39,033-atom parse (certified faithful in T2).  Structure:
a cofactor `u` sits in exactly one atom, of shape `(h - (P*u))` with `P` a p-valued wire; `h` then
appears in exactly one further atom, the **guard**; `c` is the literal multiplying `h` *inside the
guard* (bare `h` => c = 1).  `c > 1` is what makes `c*p | R` strictly stronger than `R == 0 mod p`.

**Run 1 — borrowing only L's cofactor list, re-deriving every multiplier from F's atom text:**
```
   c == 1 : 2,754        c > 1 : 927        total 3,681
```
**Run 2 — borrowing NOTHING**, deriving the family from F's parse alone by shape:
```
   family 13,092 atoms (9,626 cofactors -- 2.6x looser than L's 3,681)
   c > 1 : 927      and all 927 multipliers are DISTINCT (one per handle)
```
> **927 is reproduced in a third, independent, certified-faithful decomposition, and is
> unchanged when the handle family is delimited 2.6x more loosely.**  It is a property of the
> instance, not of the atomisation.  Contrast B1, where the knob count moved by 2.7x under
> exactly this kind of re-decomposition — the test discriminates, and here it comes out the
> other way.

**The 2,747 vs 2,754 gap reconciles exactly.**  L reports `c==1` for 2,747 "plus 7 with zero
slope"; F's parse puts all of them in `c==1`, giving 2,747 + 7 = **2,754**.  Independently, the
set difference between L's cofactor list and my F-only shape family is `L \ F-only = 7` — the
same seven.  Both directions agree.  **No discrepancy remains; L's and P's 927 stand.**

## O. K's WITHDRAWAL — the null does look like it measured the aliasing  (`t_alias.py`)
Third independent check, from L's calibrated model rather than from K's or Q's code.  For every
parent node, compare each slot wire against the corresponding child's output wire:
```
   parent-slot / child-output links examined : 764   (768 unresolvable in L's model)
   child output IS the parent slot wire      :   0
   child output is a DIFFERENT wire (aliased): 764      -> direct rate 0.0%
```
**Q's 0/383 is confirmed independently: no slot output ever feeds a slot directly.**  Of 400
sampled aliased links, **272 have a single atom containing both wires**, and the shapes are the
additive alias Q describes:
```
   ((x-x)+x)   71        ((x-x)-x)   65        ((x-x)-(1730409*x))  1   ...
```
i.e. `parent_slot - child_out -+ third_wire = 0`, matching `x_24468 = x_13682 + 12354891*x_34243`.
> **A search for a DIRECT composition would return 0 of 38,748 by construction.**  K's null is
> consistent with the aliasing layer and is not evidence that the composition is absent.  The
> withdrawal looks premature; K should re-run against the alias form, as it is already doing.
**Stated limit:** 128 of the 400 sampled aliased links have **no single atom** containing both
wires, so the alias is not always one hop — some are chains.  I did not chase those, so "every
link is a one-atom alias" is NOT established; "no link is direct" is.

## P. FINAL STATUS OF MY QUEUE
Everything the coordinator assigned me has now been run.  Nothing of consequence is left open on
my side except the two limits stated above (the 128 multi-hop aliases; Q's Schwartz-Zippel census,
which I still have not read now that the artifacts exist).

## Q. NEW FILES (fourth pass)
`t_927.py` + `t_927.json` (927, L's list borrowed) · `t_927b.py` (927, nothing borrowed) ·
`t_alias.py` (0/764 direct, alias shapes).

=============================================================================================
# FIFTH PASS — the last gate: are the six shared slack factors forced to zero?
(coordinator check-in 51; answered from F's certified-faithful parse, independent of L's model)

## R. ANSWER: **NO — and nothing should force them to zero.  They are forced to p.**
`t_slack.py`, `t_slack2.py`, `t_slack3.py`, `t_slack4.py`.
The six shared factors `x_4116, x_16153, x_1962, x_12682, x_19049, x_15616` are not slack that
ought to vanish.  **They are the modulus p, replicated as wires.**

```
   1. deliverable values:  all six are EXACTLY p = 2^256 - 2^32 - 977   (6 of 6)
   2. copy-equivalence class under atoms of shape (xA - xB):
         the six lie in ONE class, rooted at x26064, of 220 wires
         all 220 are set to exactly p in the deliverable (220 of 220)
   3. atoms anywhere in the whole instance containing the literal p:  EXACTLY ONE
         (x26064 - 115792089237316195423570985008687907853269984665640564039457584007908834671663)
```
**This is a proof, not a measurement.**  M is faithful (audit T2) and `ker(M) = 0` (F's peel
certificate, re-verified in T2), so **every atom is zero in any full solution**.  Then
`(x26064 - p) = 0` forces `x26064 = p`, and each copy atom `(xA - xB) = 0` propagates it across
the whole 220-wire class.  The six are forced — **to p, not to zero.**

### R1. WHAT THAT MEANS FOR THE REDUCTION — the hand-off is exact MOD P, not over Z
The "slack" wires are the products `(w - (P*u))` with P in the p-class — **3,707 such atoms** —
i.e. `w = p*u`.  So the slack term in the coordinate hand-off is **p times a free cofactor**:
* **mod p it vanishes identically** -> the alias is exact, the hand-off does follow the measured
  tree, and **the reduction closes mod p on measurement**, which is exactly what L's §3 already
  stated ("MOD P, the 39,033-equation system reduces to ... the target congruence");
* **over Z it does not vanish** -> the residue is precisely the extra integer conditions
  `c*p | R`, i.e. **the 927** confirmed intrinsic in my fourth pass.

> **Q's "the slack is not pinned" and L's/P's "927 extra integer conditions" are the same
> phenomenon seen from two sides.**  The slack *is* pinned — to `p*u`, not to 0.  Q was right to
> decline to close the existence result on the shape of the alias alone, and right that nothing
> forces those factors to zero; the resolution is that they are forced to p, so the closure is
> mod p and the 927 are exactly what is left over Z.  **Nothing downstream needs restating,
> provided every statement of the reduction says MOD P** — L's does; any that says "over Z"
> without discharging the 927 does not.

### R2. Honest discrepancy
3,707 atoms of shape `(w - (P*u))` with P in the p-class, against L's 3,681 residual atoms and my
own fourth-pass family of 3,681 — a gap of **26**.  I did not chase it.  It does not touch the
927 (which reproduced exactly, twice, including from a 13,092-atom family), but whoever uses 3,681
as a closed census should reconcile the 26 first.

## S. THE CHAINED-ALIAS LIMIT — resolved, and the answer is NOT the six  (`t_slack4.py`)
Folding in my fourth-pass caveat, over the full set rather than a 400 sample:
```
   aliased parent/child links : 764     one-atom aliases : 486     multi-hop : 278
   multi-hop links with a p-class wire in an incident atom :   0 of 278
```
**The 278 multi-hop aliases do not terminate in the six shared factors, or anywhere in the
p-class.**  They route through something else, which I did not identify.  So: "no parent/child
link is direct" is established (T18, 0 of 764); "every alias is one hop" is false (486 of 764);
and **what the other 278 route through is open** and is the one loose end I am leaving.

## T. NEW FILES (fifth pass)
`t_slack.py` (unary-pin census over the six) · `t_slack2.py` (every atom, definer, values) ·
`t_slack3.py` (p-class, the unique literal-p atom, the slack products) · `t_slack4.py` (forcing
chain + the 278 multi-hop aliases).
Reproduce the answer: `cd solve_lab/agentT_work && python3 t_slack3.py`.

=============================================================================================
# SIXTH PASS — the two loose ends, both closed   [coordinator check-in 57]

## U. LOOSE END 1 CLOSED — the 278 were MY artifact.  The hand-off covers all 764.
`t_278.py`, `t_final.py`, `t_cross.py`.
The 278 "multi-hop" aliases were an artifact of my own pairing, not a feature of the instance.
I matched `OUT[n][j]` against `OUT[child][j]` — index to index.  But L's `calib2` already measured
a per-node parent/child **coordinate alignment** (188 orient=1, 67 orient=0), so at a flipped node
the correct partner is `OUT[child][1-j]`.  Decoding a concrete path showed it immediately: for
`x23450.va` the path ended `((x34166-x28992)-x33628)`, i.e. the child output aliases to the node's
*other* coordinate wire.  Re-pairing to allow the cross:
```
   aliased via SAME coordinate index : 486
   aliased via CROSSED index         : 278
   still no one-atom alias           :   0
   TOTAL                             : 764
   slack wire is a p-handle (= p*u)  : 764 of 764
```
> **Every one of the 764 parent/child links is a one-atom affine alias whose slack is exactly
> `p*u`.**  The mod-p closure is **not** scoped to 486 — it covers all 764.  **Q's hand-off result
> is complete**, and the qualification I flagged as possible does not exist.  My T18 "0 of 764
> direct" stands; my "278 multi-hop" is **withdrawn** — I was measuring my own index convention.

## V. LOOSE END 2 CLOSED — the 26 reconciles exactly, and it costs M three atoms
`t_26.py`, `t_final.py`.  My p-handle family vs L's 3,681 census:
```
   mine \ L : 33     L \ mine : 7     33 - 7 = 26      3,707 - 33 + 7 = 3,681  exactly
```
* **The 7 L counts that I do not are not p-handles at all** — in each, *neither* operand of the
  product is in the p-class (`(x14163-(x13271*x11852))` etc.).  They are a different family that
  L's census sweeps in.
* **The 33 I count and L does not ARE genuine p-handles** (`h = p*u`, u free), whose guards are
  stage checks and leaf pins rather than slot links — e.g. `((x34600-x30108)+x23642)`,
  `((x22579*(x19965-9843406673...))-x3178)`.  L's census shape appears to be scoped to slot-link
  guards.

### V1. **CONSEQUENCE FOR M: the enumeration space is 18 atoms, not 15**
Three of the 33 L omits are **incident to the baseline-failing set**, and all three satisfy the
exact criterion L verified on its 3,681 (free cofactor, occurring in exactly one atom,
`eqs(u) == eqs(atom_u)`):
```
   u=x10422  (x23642-(x8173*x10422))    guard ((x34600-x30108)+x23642)          -> 12231 12350 14584 29125
   u=x15120  (x18253-(x4339*x15120))    guard ((x13502*x3629)-x18253)           -> 12231 12350 14584 29125
   u=x35531  (x37720-(x14466*x35531))   guard ((9994531*(x13502*x8976))-x37720) -> 12231 12350 14584 29125
```
> **L's "of 3,681 atoms exactly 15 are incident" is missing 3.  The true count is 18, and M's
> enumeration space is 2^18 = 262,144 candidates, not 2^15 = 32,768.**  The filter logic is sound
> — my K3 check confirmed the criterion on all 3,681 — but it was applied to a census that omits
> the stage-check/leaf-pin guarded handles.  Re-run the incidence filter over the full p-handle
> family (3,707, or 3,714 counting both operand orders) before enumerating.
Note this is the *opposite* failure mode from B1: not a count that changed under re-decomposition,
but a family delimited by guard *shape* when the defining property is `h = p*u`.

## W. NEW FILES (sixth pass)
`t_26.py` (the 33/7 split) · `t_278.py` (path search over the multi-hop set) ·
`t_final.py` (the 3 omitted incident atoms; concrete path decode) · `t_cross.py` (**the close**:
764/764 one-atom aliases, 764/764 slack = p*u).
Reproduce the close: `cd solve_lab/agentT_work && python3 t_cross.py`.

=============================================================================================
# SEVENTH PASS — agent L's |S| = 2 closure over Z   [coordinator check-in 65]

## X. THE PREMISE NOBODY TESTED — and it HOLDS.  First instance-level verification.
`t_S2.py`, `t_S2b.py`, `t_deg.py`.
Every number in L's |S|=2 result is computed inside **L's own engine** (`E.run`, a 9,032-residual
-atom model).  `solve927.py` **prints and dumps no assignment**, so the closure had never been put
in front of `checker.py`; `assign_L2.json` predates the run by two hours and is not its output.
L's |S|=1 result *was* checker-verified; the |S|=2 closure was not.

**Reproduced from cold** (driver in agentT_work, agentL_work read-only):
```
   greedy fixpoint, round 0: 1 stuck
   SOLVED c=6672769  deg=2  wire x24908  t=2990790     <- identical to L's log
   ALL c>1 CONDITIONS DISCHARGED after 1 outer round
   IN L's MODEL: 0 undischarged, 2 nonzero atoms of 9,032 = the two target congruences
```
**Dumped the assignment and checked it against the real instance:**
```
   checker.py  ->  satisfied 39018/39033  (15 failing)
   failing [4573,7123,7469,9648,11854,16622,17726,21382,25539,28653,29437,31061,32894,32916,34517]
```
**And in F's certified-faithful 39,033-atom parse** (independent of L's 9,032-atom engine):
```
   nonzero atoms: 2   -- exactly ((x24468-x13682)-(12354891*x34243)) and ((x18956-x37892)-x32237)
   equation footprint of those 2 atoms: 15
   footprint == checker's failing set:  EXACTLY
```
> **L's |S|=2 closure is confirmed, at instance level, for the first time.**  Every one of the
> 927 integer conditions really is discharged: the only atoms left nonzero anywhere in the
> instance are the two target congruences, and the 15 failing equations are precisely their
> footprint — nothing unexplained.
`assign_L1.json` (|S|=1) gives the **identical** 2 atoms and the **identical** 15 equations, so
the two ON-sets are indistinguishable at equation level; |S|=2 closing is a statement about the
integer lift, not an improvement in score.

## Y. THE THREE PREMISES THE COORDINATOR NAMED
1. **"0 undischarged" is NOT "all 927 were checked."**  The count is
   `stuck = [a for a in relift(vv) if r[E.residx[a]] % p == 0]` — bad-list entries whose residual
   is **not** 0 mod p are silently **excluded**.  Measured at |S|=2: **2 such exclusions**, and
   they are exactly the two target congruences, so benign here.  But the metric carrying the
   weight is the companion **"2 nonzero atoms of 9,032"**, which *is* a complete check, and which
   my F-parse pass now extends to the whole instance.  **Report the nonzero-atom count, not the
   stuck count** — at |S|=17 the two differ ("1 undischarged, 3 nonzero").
2. **Direct recomputation: genuine.**  `solve_one`'s guard is
   `val = probe(vv,i,[w],[t]); if val % (c*p) == 0`, and `probe` calls `E.run(vv)` with the shift
   actually applied — the fitted polynomial is used only to *propose* the root.  P's second guard
   is properly inherited; no sign-bug of P's kind can survive it.
   **Limitation, and it is the one that blocks |S|=17:** the guard verifies **only the target
   atom**.  Applying `vv[w] += p*t` can disturb others, which is exactly the oscillation on the
   shared wires x23238 / x10261.  Per-condition verification does not establish global
   correctness; only the final whole-model nonzero count does.  At |S|=2 that final check passes,
   so the result is safe — at |S|=17 it is the thing still to prove.
3. **The degree bound is real, and it is NOT load-bearing for soundness.**  `fit()` samples
   t = 0..4, exactly 5 points, so it can only *see* degree <= 4 — a true degree-5 polynomial would
   be silently aliased.  Re-fitted the |S|=2 condition at deg <= 6, 8 and 10 (up to 11 points) on
   all six influencing wires:
```
   x24908 -> 2 2 2 2     x16742 -> 2 2 2 2     x14853 -> 3 3 3 3
   x30213 -> 1 1 1 1     x22162 -> 1 1 1 1     x12186 -> 3 3 3 3
```
   Same top degree at every sample size: **the degrees really are <= 3, not an aliasing artifact,
   and P's bound is independently confirmed a third time.**  Moreover, even if the bound were
   wrong, the direct-recomputation guard would reject the resulting root — **a bad degree bound
   can only cause a missed solution, never a false verified one.**  It is load-bearing for
   completeness and cost, not for correctness.

## Z. SCOPE — what this does and does not establish
It establishes the integer lift closes for **one** ON-set of size 2, verified against the
instance.  It does **not** establish that it generalises: the |S|=17 joint solve still ends
**1 undischarged, 3 nonzero atoms**, and the obstruction is the shared-wire simultaneity in
premise 2 above, not the degree bound and not the fit.  **The decisive experiment remains
|S| = 17.**  I would also want |S| = 3, 5 and 8 checked *and dumped*: three more points would
distinguish "closes for small |S|" from "closes generally", and each is minutes of compute.
**Whatever is run next should dump the assignment and pass it through `checker.py`** — that step
cost nothing here and is what turned a model-internal claim into an instance-level fact.

## AA. NEW FILES (seventh pass)
`t_S2.py` (reproduce the closure + dump the assignment) · `t_S2_assign.json` (the artifact L's
run never produced; checker-verified 39,018/39,033) · `t_S2b.py` (F-parse atom check; footprint
== failing set) · `t_deg.py` (degree bound at 5/7/9/11 sample points).

=============================================================================================
# EIGHTH PASS — agent O's eq8680 Lemma   [coordinator check-in 69]

## AB. THE LEMMA SURVIVES.  Two numbers in its statement are wrong; the conclusion is not.
`t_eq8680.py`, `t_eq8680b.py`.  Checked against the raw `EQUATIONS.txt` line and against F's
certified-faithful parse as a third source.

**CONFIRMED — the substance:**
* The equation really does factor as a perfect power of a single **affine** form.
* That form `S` is **affine in ALL 43 of its variables — 0 non-affine**, tested by second
  differences on every one.
* `dS/dx_4432 = +1`, `dS/dx_28730 = -1`, `dS/dx_19964 = -1`, measured exactly.
* `a23618 = x_4432 - x_19964 - x_28730` enters at coefficient **exactly +1** — it is the first
  of the 18 terms.
* **F's certified parse independently gives the same decomposition**: 18 `(coef, atom)` entries
  with coefficients `[-27,-21,-14,-13,-5,-4,1,1,1,6,15,17,20,23,25,25,28,35]`, identical to what
  I flattened out of the raw text (the only set difference is F's `x4432` vs the file's `x_4432`).
  **Three sources agree on the decomposition** — this is the case where my five-atom-counts rule
  is satisfied rather than violated.

**CORRECTION 1 — the equation is `S^4`, not `S^2`.**  The nesting is two levels deep:
`LHS = T*T` with the two factors textually identical, and `T = S*S`, again identical.
Numerically, at 4 random points, `LHS == S^k` **only for k = 4** (and `LHS == T^k` only for k=2).
**CORRECTION 2 — the linear form has 18 atoms, not 20.**  Confirmed twice (raw flatten, F's parse).

**The two corrections expose an internal inconsistency in the statement as written.**
"`eq8680 = T^2`" and "`dT/dx_4432 = +1`" cannot both be about the same object: the thing that
squares to eq8680 is `T = S^2`, which is *not* linear, and whose derivative I measured as
`dT/dx_4432 = 38046996267 = 2S+1` at my test point.  The object with derivative `+1` is `S`, and
`eq8680 = S^4`.  O conflated one level of the nesting.

## AC. THE CONCLUSION IS UNAFFECTED — and the modulus worry does not arise
`checker.py` evaluates each equation's LHS as an **exact integer** and requires `== 0`.  So the
constraint is `S^4 == 0` **over Z**, and Z is an integral domain, so `S^4 == 0 <=> S == 0`.
* **No modulus is in play at equation level**, so nothing needs p to be prime or squarefree.
  The coordinator's flagged risk — "`T^2 == 0 (mod something)` rather than over Z" — **does not
  materialise here.**  O applied "a square has a single zero locus" to the right object in the
  right ring.
* And the conclusion is **robust to the exponent entirely**: `S^k = 0 <=> S = 0` for every k >= 1
  over Z, so O getting the power wrong could not have broken it.  Same shape as my |S|=2 finding
  that a bad degree bound cannot produce a false verified result.

> **Verdict: agent O's Lemma holds.  `S = 0` is forced in every satisfying assignment, with no
> knob set, no frame, no configuration and no divisibility condition.  Fix "T^2" to "S^4" and
> "20 atoms" to "18"; the seven-way trade, the death of delta_0, and M's incidence argument all
> stand.**  This is the first result in this lab that is both unconditional and audited.

## AD. CROSS-LINK — eq8680 corroborates my sixth-pass finding
All three cofactors I found in the sixth pass to be genuine p-handles that L's 3,681 census omits
**and** incident to the baseline-failing set appear as terms of `S`:
```
   +25  (x_18253) - ((x_4339)*(x_15120))       <- x15120
    +1  (x_37720) - ((x_14466)*(x_35531))      <- x35531
   +23  (x_23642) - ((x_8173)*(x_10422))       <- x10422
```
So O's own equation independently confirms they are incident.  **Caution against a numerical
coincidence:** `S` has **18 atom terms**, and M's enumeration exponent is also now **18** (my
15 -> 18 correction).  These are different 18s and must not be conflated.

## AE. NEW FILES (eighth pass)
`t_eq8680.py` (factorisation, affinity, derivatives, the modulus argument) ·
`t_eq8680b.py` (18-term flatten; cross-check against F's parse; the x10422/x15120/x35531 link).

=============================================================================================
# NINTH PASS — agent N's detach exhaustion   [coordinator check-in 73]
Auditing the largest **reported** row in my own ledger, per my own closing flag.

## AF. N'S REDUCTION IS CONFIRMED — AND IT IS A PROOF, NOT AN ENUMERATION  (`t_detach.py`)
N's claim: `make(D)` gives detached pool members their witness values; only 4 of the 65 pool
variables have witness != gate, namely `{642, 28730, 29854, 31864}`; so detaching the other 61 is
a no-op and the 2^65 lattice has exactly 16 states.  Checked in **F's certified-faithful parse**,
not N's, via an identity N does not state:

> a pool variable `v` is defined by an atom `(v - RHS)`, so
> **`witness(v) != gate(v)` at the deliverable  <=>  that atom is NONZERO at the deliverable.**

I already knew from audit T2 that F's parse has exactly **7** nonzero atoms there.  Measured:
```
   pool size 65; all 65 have an identifiable defining atom (xV - RHS)
   defining atom NONZERO at the deliverable:  x642, x28730, x29854, x31864   -- exactly 4
   equals N's witness set exactly?  TRUE
```
**The other 61 have their defining atom zero, so re-attaching them is a no-op at the witness state.**

### AF1. The gap N's argument leaves, and it closes
"No-op at the witness state" is *not* the same as "no-op at all 16 states": if a non-witness pool
variable's RHS depended on a witness variable, its gate value would change when that witness
variable is re-attached, and the lattice could exceed 16.  N does not address this.  Measured:
```
   of the 61, DIRECTLY referencing a witness variable                    : 0
   of the 61, reaching one transitively WITHIN the pool                  : 0
   of the 61, reaching one anywhere in the FULL definition DAG (30,001)  : 0
```
> **Zero, over the whole instance.**  The 61 gate values are independent of `D`, so `make(D)`
> depends only on `D & {642,28730,29854,31864}`.  **The 2^65 detach lattice has exactly 16 states
> by proof, and the 16 measured `(R,b)` signatures are complete by construction rather than by
> having happened to be reached.**  N's conclusion — OPT = 5 for all 16, best 39,026, axis closed —
> stands on a reduction that is now verified from an independent parse.
The row is promoted in LEDGER.md from *reported* to *T verified independently*.

## AG. THE 924/924 p-OBSTRUCTION IS **NOT** O'S LEMMA — and N got there first
The coordinator's hypothesis was that N's `924/924 obstruction denominator divisible by p` might
be measuring O's `S = 0`, making the two results one.  **It does not, and N had already
established why** (RESUME_N, "O's Lemma confirmed, and it is a SECOND obstruction, not mine"):
* `T = 0` already holds at the witness;
* **eq8680 is exactly the one equation that detaching `x_28730` buys** — so O's Lemma *is* the
  39,025 -> 39,026 step, which is a much sharper statement than "the same result twice";
* the witness region **excludes 8680**, so the 924/924 p-obstruction is independent of it;
* in the 13-row region, max rows zeroable subject to 8680 being zeroed is **0**, so the knobs
  cannot reach `T = 0` at all and detaching `x_28730` is the only route.
Nothing for me to correct.  Cross-link worth recording: `x_28730` is simultaneously one of N's 4
witness variables, one of the h-wires in L's cancellation set, and the variable entering O's `S`
at `dS/dx_28730 = -1`.  **Three threads are describing one wire.**

### AG1. ONE ACTIONABLE ITEM FOR N
N writes `eq_terms[8680] = (1, True, [(1, 37887)])` and "`optN.inner` returns the inner form,
never its square, so my model already carried `T`".  Per my eighth pass the nesting is **two
levels deep — the equation is `S^4`, not `S^2`** — so an `inner` that strips **one** level yields
`S^2`, which I measured to be **non-affine in all 43 of its variables**, not the linear form.
* If N only uses the zero locus, nothing changes: `S^2 = 0 <=> S = 0`, and N's conclusions hold.
* If N **linearises** that row anywhere, it is linearising a quadratic.  **Strip twice.**
This is the "wrong number, right result" pattern again — flagging it as a code note, not a defect
in N's conclusion.

## AH. NEW FILES (ninth pass)
`t_detach.py` (65-pool witness/gate check in F's parse; full-DAG independence of the 61).
Reproduce: `cd solve_lab/agentT_work && python3 t_detach.py`.

=============================================================================================
# TENTH PASS — agent M's enumeration engine   [coordinator check-in 79]
Question 2 (are the 30 optima supersets?) was closed by M's own `verifysup.py` before I started,
with a refinement that falsifies the neighbouring claim — I did not duplicate it.
**Question 1 was the whole audit: is M's scorer exact, checked OUTSIDE M's parse?**

## AI. M'S ENGINE IS EXACT — verified independently, and it reproduces the deliverable exactly
`t_meng.py`, `t_meng2.py`.  M's six calibration gates include "incremental == full engine3,
0 vars differing", but that is M checking M, and the entire value of a 4,096-subset exhaustive
verdict rests on the scorer.  Drove `ieng.tune` on H12's witness subset, **materialised the
assignment M's engine actually scores**, and put it in front of `checker.py` and F's parse:
```
   M engine on {642,28730,29854,31864}: base_score 39008 -> score 39026, 5 knobs, 5 vars changed
   CHECKER (independent)              : satisfied 39026/39033, 7 failing
   failing == the deliverable's exact [12231,12270,12350,14584,18673,22044,29125] :  TRUE
   M's reported score == checker's score                                          :  TRUE
   F's certified parse                : exactly the deliverable's 7 nonzero atoms
```
**And the assignment M's engine produces is byte-identical to the deliverable: 0 of 38,748
variables differ.**  M's engine does not merely reproduce the score — it reconstructs the file.

### AI1. The footprint/failing gap is the cancellation, not a defect
The 7 nonzero atoms touch **12** equations of which only **7** fail: **5 cancel to zero**, namely
`[2554, 6816, 8124, 9123, 9421]`, and **no failing equation lies outside the footprint** — nothing
unexplained.  Cross-link worth recording: those are **exactly** the 5 equations that appeared as
new failures in my third pass when I zeroed L's cofactors (`t_cancel.py`, 7 -> 12).  **The five
equations that cancel are precisely the five that break when the cancellation handles are
zeroed** — L's mechanism and M's engine agreeing exactly, from opposite directions.

### AI2. Exactness away from the calibration point — 9/9
A scorer can be exact where it was calibrated and wrong elsewhere, and the enumeration's value is
its verdict on the *other* 4,095.  Spot-checked 9 subsets spanning the range, each materialised
and scored by `checker.py`:
```
   (642,28730,29854,31864) 39026=39026   (28730,) 39009=39009   (642,29854) 39011=39011
   (642,1844,9629,28730,29854,31864) 39026=39026   (31864,35619) 39009=39009
   (23754,35619) 39008=39008   (642,1844,35619,37413) 39009=39009
   (18253,23642,31864,35619) 39008=39008   [+1 duplicate subset]
   agree 9 / disagree 0
```
> **M's engine is exact at every point tested, including the maximum, the base and four
> non-optimal interior points.  The row is promoted to *T verified independently*.**
**Scope, stated plainly: 9 of 4,096 subsets.**  "Nothing above 39,026" rests on the scorer being
exact *everywhere*, and I have verified it at 9 points, not 4,096.  What I can say is that the
engine is exact wherever I could check it and that it reconstructs the known answer bit-for-bit —
which is the strongest cheap evidence available, and considerably more than a self-calibration.

## AJ. NEW FILES (tenth pass)
`t_meng.py` (materialise + checker + F parse + identity with the deliverable) ·
`t_meng2.py` (9-subset exactness spot-check) · `t_meng_assign.json` (the engine's own output,
checker-verified 39,026).  Reproduce: `python3 t_meng.py && python3 t_meng2.py`.

=============================================================================================
# ELEVENTH PASS — L's closure sweep, run by T   [coordinator check-in 83]
Scripts handed over in `agentT_work/from_L/`; `closeS4.py`'s global guard (accept a shift only if
the TOTAL nonzero-atom count strictly decreases, verified by direct recomputation).
Driver `t_sweep.py`: data load with `agentL_work` as cwd, then **chdir back so the dumps land in
my directory**.  Launched detached, PID in `job.pid`, liveness by `kill -0`.  Artifacts named
`close_T2ctl / T3 / T5 / T8.json` — deliberately distinct from L's invalid `close_S3/S5/S8.json`.
`agentL_work` was not written to.

## AK. RESULT — **the integer lift does NOT close generally.  It closes at 2, 3, 5 and breaks at 8.**
```
   |S|  tag     nonzero atoms of 9,032   checker      wall     closes?
    2   T2ctl            2               39,018      160 s     YES   <- CONTROL, matches audit T24
    3   T3               2               39,018      171 s     YES
    5   T5               2               39,018      179 s     YES
    8   T8               3               39,002      289 s     ** NO **
```
**Control passed**: exactly 2 nonzero atoms and 39,018, reproducing the number I established in
audit T24 — so the global guard does **not** have the leak the scoped guard had (which gave 8).

At `|S| = 2, 3, 5` the only nonzero atoms are the two target congruences, and all three give the
**identical** failing set `[4573,7123,7469,9648,11854,16622,17726,21382,25539,28653,29437,31061,
32894,32916,34517]` — the 15-equation footprint of those two atoms, exactly as at `|S|=1,2`.

At `|S| = 8` a **third** atom survives:
```
   ((x21408*x10138)-(15333171*x658))        c = 15333171 = 3 * 7 * 19 * 83 * 463
```
and the score drops to **39,002 / 31 failing**.  The solver also worked visibly harder (289 s vs
~175 s) and still did not clear it.  Note the c factors into small primes, so this is **not** the
large-prime-factor cost case — root-finding was cheap and the obstruction is genuine.

> **Answer to the campaign's last open question: closure is a small-`|S|` phenomenon, not a
> general one.  The boundary lies between 5 and 8.**

### AK1. What this does and does not establish — read before quoting
* **One ON-set per size**, drawn by L's own `random.Random(7)` convention, not an exhaustive
  sweep.  `|S|=8` failing at *this* ON-set does not prove every 8-leaf ON-set fails; `|S|=3,5`
  closing at *these* ON-sets does not prove every one closes.  **The honest statement is: the
  first observed failure is at 8, and 3 and 5 were observed to close.**
* It is **"this solver did not close it"**, not "it cannot be closed".  `closeS4` stops when no
  single-wire shift strictly decreases the global nonzero count; a residue that needs two wires
  moved *together* would look identical.  Given that the |S|=17 obstruction was already diagnosed
  as shared-wire simultaneity, that is the live hypothesis for |S|=8 too.
* **Next, and it is cheap:** the nested ON-sets mean `T3 ⊂ T5 ⊂ T8`, so re-running `|S| = 6, 7`
  on the same nested prefix would localise the break to a single added leaf — about 6 minutes,
  and it would say whether the failure is a property of size or of one particular leaf.

## AL. PARKED — O's seven-way trade, partial, NOT a verdict  (`t_oktrade.py`, `t_oktrade2.py`)
Superseded in priority by the sweep; recording what I measured so it is not lost.  **This is not
an audit verdict on O and must not be quoted as one.**
* **`K` is frame-dependent.** Rebuilding O's `K` in the *default* orientation from F's parse gives
  **12** free S-carriers (O says 26), **11** free inputs reaching a nonzero region atom (O says
  15), union **23** (O says 34), overlap **0** (O's numbers imply 7).  O explicitly scopes its
  theorem to "frame B's orientation", which promotes defined variables to free — the most likely
  and innocent explanation.  **I did not reproduce frame B, so this is a flag, not a defect.**
* **The seven-way uniformity is NOT structurally forced.**  I tested the natural reduction — if
  every knob that moves a failing row also moved `S`, then "every purchase costs exactly eq8680"
  would be one fact seen seven times.  It is false: **7 knobs** (`x1329, x7068, x8731, x9118,
  x9413, x10903, x17325`) move a failing row with `dS = 0` exactly.
* **But that does not contradict O**, and the distinction matters: *moving* a row is not *buying*
  it.  Restricted to those S-preserving knobs, 6 of the 7 failing rows are individually solvable
  over ℤ with `S` held at 0 (eq29125 is not) — **but that check ignores collateral**, and O's
  "buyable" means the score does not drop.  So my result is a **necessary-condition** statement
  only.  **Conclusion: O's uniformity is a genuine search result resting on its collateral
  accounting over `K`, not a restatement of N's "eq8680 is what detaching x_28730 buys". Whether
  that accounting is right is UNAUDITED.**

## AM. NEW FILES (eleventh pass)
`t_sweep.py` + `t_sweep.log` + `job.pid` · `close_T2ctl.json`, `close_T3.json`, `close_T5.json`,
`close_T8.json` (all checker-verified above) · `t_oktrade.py`, `t_oktrade2.py` (parked O work).
