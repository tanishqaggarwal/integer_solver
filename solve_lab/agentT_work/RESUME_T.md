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
