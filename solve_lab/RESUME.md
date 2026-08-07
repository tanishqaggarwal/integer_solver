# RESUME — read me first

## STATUS: deliverable **39,026 / 39,033** (re-verified with `checker.py` this session)
Verify: `python3 checker.py best/new_instance_partial_39026.json` -> 39026/39033 (7 failing).
Second branch, independently verified: `s11/data/finish3_named.json` -> **39,018** (x5647 channel).

> **READ `S11_PART4_LOCAL.md` FIRST.**  It replaces the "rigidity"/"one factor of p" language of
> Parts VI-VIII (already retracted in Part IX) with exact algebra: what the obstruction *is*,
> that it is *reachable*, and what the new wall is.  Then `S11_SEMANTICS.md`, `S11_PART2.md`.

### The 90-second version of the newest work
1. **The checkpoint's defect is nine atoms and fifteen equations, and nothing else.**  Nine
   variables occur in no atom outside them.  Eliminating those nine knobs gives the reachable
   set exactly, as a coset of congruences (`s11/local1.py`).  Enumerating drop-sets proves
   **7 is exactly the local optimum** — a proof, replacing eleven sessions of "it looks rigid".
2. **The whole obstruction is three congruences**: `p | x9118`, `p | x8731`, `p | x28730`
   (the fourth, `7376877*p | x7068 - x2099`, follows from two gates).  `local2.py`: the 15x9
   coefficient matrix has rank 9, so all fifteen hold only when all nine atoms vanish.
3. **They are reachable.**  x9118 and x8731 are FREE inputs; x28730 is defined by the very gate
   that demands its divisibility.  `s11/fix2.py` clears all nine atoms — the first time the
   checkpoint's obstruction has actually been removed.  It relocates rather than vanishing.
4. **Second, independent route.**  `x4287 = 1` makes `x21279 = 1` hence `x7075 = 0`, and since
   `T = 5113045*x7075*x9118`, `U = x7075*x8731`, BOTH congruences evaporate.  The price is three
   new ones on the x21279 channel, which collapse to `x27177 = 0` and `x4306 = 0` (mod p) —
   affine in x9118 and x8731, hence solvable; `s11/sw6.py` achieves all five simultaneously.
5. **Both routes hit the same wall.**  Either way x7068 must change residue mod p, and x7068 is
   copied into free "mirror" inputs across the circuit.  Repairing the mirrors fans out and, in
   four independent searches from both routes, always terminates at
   `a19297`, `a19299`, `a30984`.  Each has a free p-handle one gate up (x30317, x5146, x2936),
   so they reduce to `p | x11150`, `p | x25739`, `p | x37758` — **the same three-congruence
   obstruction, recurring one level out on GATE-COMPUTED values instead of free inputs.**  Their
   collateral is only two single-equation checks, so meeting them would jump the score.
6. **That wall is channel-specific.**  In the x5647 channel (`s11/data/finish3_named.json`,
   score 39,018) `x15298 = 0` and all three wall checks are vacuous.  That branch's own defect
   is three atoms / fifteen equations needing `8640431*p | x12000`, `p | x12926`, `p | x21364`.

### NEWEST OF ALL: the obstruction is CONSERVED QUANTITIES of the message — `S11_PART6_BITS.md`
Each obstruction certificate y gives `INV_y = sum_a y_a * r_a`, constant under every continuous
move, computable from any message in 0.08 s (3.9 ms with the cone evaluator).  A full solve needs
all six to vanish, so the objective stops being "count failures" (plateaued at 4, no gradient) and
becomes "drive six explicit GF(p) numbers to zero".

* **inv5 depends on only 18 of the 256 bits.**  Enumerated EXHAUSTIVELY over all 2^18 = 262,144
  subsets (`s11/bits24.py`, 19 min): **232 distinct values, zero never attained.**
* At the checkpoint's own message the certificates are exactly conserved (all 14 knobs that touch
  their rows are annihilated, `s11/bits25.py`), so **that message provably cannot be completed** —
  the first actual reason for the plateau rather than another failed search.
* At a sibling message only 12 of 14 are annihilated, so across messages this is a strong screen,
  NOT a proof.  Channel A is very likely dead; channels B and C are where to go.
* Bit landscape mapped exactly: 256 real bits (900 of the 1,156 free booleans are provably inert),
  trees 88/90/37/41, each bit owning two private load pins.  Exhaustive weight-2 scan: minimum 4
  failing checks, reached by exactly three messages — {x24601,x2081} (the checkpoint),
  {x24601,x4287}, {x24601,x13195}.

---

### The barrier is FOUR numbers in GF(p) — read `S11_PART5_MODP.md`
Mod p every handle dies (they are all `free * wire`, wire = p), so the instance is a plain GF(p)
circuit that forward-evaluates **in 0.08 s** (`s11/gmp1.py`).  From the best inputs, every gate is
satisfied and exactly SIX checks fail; two of them (a35759, a35760) are cleared for free by x9118
and x8731 with zero collateral.  The barrier is the other four: a7930, a29539, a40826, a41512.

Established this session, all by exact computation over the WHOLE instance:
* no continuous move reaches them — 1,726 live free inputs (full response matrix cached in
  `s11/data/resp_modp.pkl`) plus 725 gate-purchased knobs, still inconsistent;
* the message bits are frozen against continuous motion — each bit is the *only* knob touching
  its own load-pin row, so any system containing the pins forces its coefficient to 0;
* bits are non-additive, no single flip improves the count, and exactly 256 of the 1,156 free
  bits are real (the other 900 are provably inert);
* 7 failing equations is optimal at this base over all continuous knobs + gate purchases + drops
  (`s11/gmp33.py`), so the checkpoint is the best the continuous structure permits.

**Next attempt: the 256 discrete bits.**  0.08 s per exact evaluation means ~40,000 patterns an
hour.  The three channels are U*V (checkpoint, 4 failing), (1-U)*V (39,018 state, 3 failing) and
U*(1-V) (turn x2081 off — reachable, never explored).

---

### Earlier step (see `S11_PART4_LOCAL.md` Parts XII-XIII for the full chain)
`localopt2.py` has already priced the x5647 branch: drop 0-6 all infeasible, and every
ripple-response certificate at 3 hops contains a26719, so that branch is ~11 — worse than 7.

The live lead is the **core-side chain**.  The three wall congruences of §17 ARE solvable
(two affine knobs x22162, x30213; rank 2 and consistent), which clears a19297 and a30984
outright and leaves the bill on the core trio a688/a1618/a40608 — the one part of this instance
with mature machinery (Part I: rank 2, reduces to a cubic mod p, solved exactly by
`s11/solveA.py`; `a40608 = (W-C)^2` is never independent).  Concretely:

  1. Re-derive the x22162 / x30213 correction as a point ON the core's cubic variety rather than
     as a free shift — single-knob solves fail the p-divisibility (`s11/core2.py`), but on the
     variety the two core quantities move together, which is the only motion p-quantisation allows.
  2. Close a19299 with its extra factor 6672769, by the `s11/sw7.py` recipe (shifts of x22162,
     x30213 by multiples of p preserve all the mod-p work and move x25739 mod 6672769).
  3. Close a21617 via its free variable x14623 (exact solve already verified).

If a688, a1618, a19299 and a21617 all close, the residue is four single-equation checks —
**4 failing, score 39,029**.

---

## (previous RESUME follows)


## STATUS (session 11): deliverable **39,026 / 39,033** (unchanged, re-verified)
Verify: `python3 checker.py best/new_instance_partial_39026.json` -> 39026/39033 (7 failing).
Session-11 best in a NEW branch, independently verified: `s11/data/finish3_named.json` -> **39,018**.

> **READ `S11_SEMANTICS.md` THEN `S11_PART2.md`.** Session 11 stopped doing atom algebra and
> decoded the circuit as a program; Part II then priced every branch and localised the
> obstruction. They supersede S9/S10 wherever they conflict.

### The 60-second version of session 11
1. **The checkpoint sits in channel U=V=1 with BOTH mirror gates off** (only 2 message bits on:
   x_2081, x_24601). `a23000 = (1-U)(1-V)` with `U = OR(x_8599,x_21839)`, `V = OR(x_7304,x_25956)`.
2. **`a40608 = (W - C)^2`** exactly — never an independent obstruction (earlier price lists
   double-counted it).
3. **Each "core" is rank 2 in two quantities**, and eliminating a variable turns group 1 into a
   **CUBIC mod p**. Cubic roots are invisible to Jacobian/Newton/beam — which is exactly why ten
   sessions called the core rigid. Cantor-Zassenhaus solves it instantly.
4. **Best branch: U=0, V=1 with bits (490, 91)** — first core dead, group-2 mirror dead, group-1
   mirror already satisfied. Drives `a688 = a1618 = a40608 = 0` EXACTLY.
5. **The mirror trio needs an extra divisibility by 8640431 = 53 x 163027.** `gamma(k,l)` under
   p-shifts of `x_31339`, `x_33708` has bidegree (2,3); exact interpolation + CRT gives gamma = 0.
6. **THE DEFICIT IS EXACTLY 2, AND IT IS TOPOLOGICAL.** Exhaustive scan of all 7,253 non-locked
   free inputs: `a14445` and `a34580` have exactly ONE non-bit control each and it is the SAME
   variable `x_33129`; `a27139` and `a33796` likewise share only `x_37088`. `x_33129` is the free
   variable of a14445 *and* feeds `x_15111 -> x_20541 -> x_10170`, the other side of a34580. No
   message choice removes this — it is circuit topology.
7. **Why 39,026 still wins.** The deficit is 2 in every channel. What differs is the price of the
   absorbing set: in the (490,91) branch the cheapest is the mirror trio at **15** equations; in
   the checkpoint's channel it is the `x_2099` ladder at **7**. Cheapness of the absorber decides.


### PART III (same session): the deficit is PROVED = 2 (bipartite matching)
- **Bits are genuinely boolean**: `a1430 = x_490^2 - x_490` and one such check for EVERY one of
  the 256 message bits (13-14 equations each). So they are not continuous controls -- this
  closes the one loophole that would have removed the deficit outright (`s11/boolform.py`).
- **Maximum matching over the constraint/control graph** (all 7,273 free inputs scanned, bits
  removed): **14 constraints, matching 12, DEFICIT = 2** (`s11/hall.py`). Unmatched: the two
  mirror residuals. Explicit **Hall violator**: 9 constraints
  `{a14445,a21050,a26839,a29539,a34580,a40065,a688,mirror3719,mirror25118}` whose whole control
  neighbourhood is only 8 variables `{14515,16441,18751,19750,28955,31339,33129,33708}`.
- **Absorbers priced** (`s11/pairprice.py`): cheapest pair `a688+a1618` = 15 equations, mirror
  trio = 15. So 15 is the branch floor -- achieved and verified (39,018).
- No single-atom equation-space compensation exists (`s11/compensate.py`: 0 proportional
  shadows for any absorber), and the constrained equation-space solve over the full
  173-equation region with 26 exact-linear handles returns NONE (`s11/realise3.py`).

> The obstruction is now a THEOREM about this instance, not an observation: a Hall violator of
> size 9 over 8 controls. The score is decided entirely by the price of the absorbing set --
> 15 in the (490,91) branch, 7 in the checkpoint's channel.

- **Equation space offers no escape** (Part IV). Of the six knob atoms of the 15-equation
  region only THREE have a private handle (`a26719<-x24175` step -8640431p, `a26721<-x4615`
  step -p, `a26723<-x13992` step -p); the other three are gates with no freedom. The honest
  system is 15 equations in 3 LATTICE knobs, and no assignment satisfies even |S|-8 of them,
  at any saved state (`s11/eqopt.py`, `s11/eqopt2.py`). Steps are multiples of p, so they can
  only repair equations whose atom-sum is already 0 mod p; the failing ones are not.

- **Gate-breaking priced (Part V).** Running `fw.forward` on the checkpoint DESTROYS it (37
  failing): its score depends on five broken GATE atoms. Breaking a gate frees its output as a
  new control; 817 gate atoms live in <=8 equations, cheapest `a41332` [1 eq] -> x_24453, then
  `a36244` [4 eq] -> x_3432, and 12 cheap gates move the mirror (which has no non-bit control).
  Breaking those two costs 5 equations and would give 39,028 IF they were independent — but the
  joint 6x6 Newton is singular at all 12 starts (`s11/joint6.py`): they buy ONE dimension, not
  two, so the total returns to 15. The deficit survives the one attack the rest of the session
  had structurally excluded.

### PART VI (same session): treated as an INTEGER PROGRAM — the obstruction is ONE NUMBER
Scripts `s11/ip1.py` .. `s11/ip11.py`; full write-up in `S11_PART3_IP.md`.
- **IP #1** min-cost defect placement (exact over all 2^14 subsets): optimum **15** equations for
  the (490,91) channel — certifying the 39,018 construction as optimal *for that channel*.
- **IP #2** global lower bound: absorbers must lack a private handle; cheapest pair = **2
  equations**, so no score above **39,031** is possible anywhere. Also found that the
  checkpoint's 7 atoms span 12 equations but only 7 fail — **cancellation is real**, so the
  objective must be over EQUATIONS, not atoms.
- **IP #3/#4** minimum-weight coset `min ||b + G k||_0` + integer Newton: independently returns
  **15** for the same channel. Two structurally different exact methods agree.
- **IP #5/#7/#8** the checkpoint as a RAW integer program with no circuit orientation (needed,
  because forward-evaluating destroys it): 19 exact-linear core vars, 123 collateral equations,
  50 compensators, system 130x69. **allow = 0, 1 and 2 are ALL infeasible.**
- **IP #9/#10** THE RESULT: the system is consistent over Q, its solution is supported on
  exactly the seven x_2099 ladder variables, and the least d making `M x = d*rhs`
  integer-solvable is **d = 2458959 * p** (= 3 * 819653 * p) — every proper divisor fails.
  > **The entire obstruction at 39,026 is a single divisibility by 2458959*p.**
  At the checkpoint 0 of 7 failing values are divisible by p (gcd = 1); same at 39,018.

> **NEXT TARGET, stated exactly:** reach a state where the failing equation values are
> = 0 (mod p). The obstruction then collapses to divisibility by **2458959**, a 7-digit
> modulus — attackable by the same CRT/quadratic-form method that cleared 8640431
> (`s11/quad3.py`).

- **IP #12/#13 — the p-factor is UNIVERSAL.** Computed the invariant at every saved state:
  7 of 7 are consistent over Q and 7 of 7 have invariant divisible by p; cofactors D/p are only
  {1, 2458959, 8640431} (both small ones are handle multipliers and both are CRT-clearable — Part
  II cleared 8640431). At `closehit2` the cofactor is **1**, i.e. the invariant is EXACTLY p.
  > **THE TRAPDOOR IN ONE SENTENCE: every reachable state leaves a residual integer program that
  > is solvable over Q and whose sole integrality obstruction is a single factor of
  > p = 2^256 - 2^32 - 977.**
  This is what earlier sessions called p-quantisation / the conserved obstruction / the deficit
  of 2 / "7 is an invariant". A full solve requires REMOVING the p from the invariant, i.e.
  reaching a state whose failing right-hand side is already p-divisible.

- **IP #14/#15/#18 — the p-factor cannot be removed locally (loop closed).** A full solve needs
  a p-divisible failing right-hand side. Asked with all variables and no constraint it looks
  solvable (`ip14`) — but applying it turns satisfied equations into nonzero multiples of p
  (`ip15`: 28 failing -> 6,097). Asked properly, INSIDE the integer kernel that keeps every
  satisfied equation exactly zero (`ip18`): kernel dim 2 at the checkpoint, 10 at 39,018, and
  Stage A is **NOT solvable** at either — you cannot make even ONE failing value p-divisible
  while preserving the rest.
  > Sole obstruction = one factor of p; absorbing it needs a p-divisible RHS; p-divisibility is
  > unreachable inside every preserving move. That is why local search cannot shift the score.

### CORRECTION (read before trusting Part VI/VII): the obstruction proofs were LOCAL
`s11/perm.py`, `s11/hensel.py`, `s11/relax.py`, `s11/closure.py`, `s11/newtonp.py`.
1. **Nothing is permanently unfixable.** Every failing equation has 6-30 variables with genuine
   mod-p leverage, at both the checkpoint and 39,018. The wall is COUPLING, not rigidity.
2. **The "exact-linear" filter deleted exactly that leverage.** Requiring
   f(u+2)-f(u) == 2*(f(u+1)-f(u)) keeps only variables entering LINEARLY and rejects every
   quadratic one. With the filter the system is unsolvable even mod p at every region size
   (130x69 .. 1400x1014); with the TRUE symbolic Jacobian the mod-p region system solves in 11s.
3. **Regions were never closed**: closing the checkpoint's failing region reaches 26,598 of
   39,033 equations and 28,232 of 38,748 variables. The problem DOES NOT LOCALISE.
> So "the sole obstruction is one factor of p" and "the p-factor is universal" are statements
> about restricted move sets, NOT about the instance. Withdrawn as global claims. The 39,018
> construction and the per-channel optima stand (correctly scoped).

### CORRECTED NEXT STEP (highest value)
Global mod-p solve: find an assignment with EVERY equation = 0 mod p, using Newton with the
true symbolic Jacobian (`s11/newtonp.py` has it) and the circuit's triangular gate structure to
make the step tractable. Then the residual is p*r and the p-quantised handles absorb it exactly
-- a clean two-stage route to a full solve.

### Do NOT redo
- The clean all-zero frame, the MUX/OR-tree decode, `a40608 = (W-C)^2`, the core rank-2
  reduction, the cubic, the 8640431 CRT step, the channel taxonomy, the control scans.
- Newton/beam/local search on the cores. They cannot find cubic roots mod p.
- Turning on all four of a,b,c,d — it lights BOTH mirror cores for nothing.

### START HERE NEXT SESSION
1. The only lever left is **the price of the absorbing set**. Enumerate, per channel, every
   2-deficit absorbing set and its equation cost; the checkpoint's 7 is the best known. Look for
   a channel where two 1-equation checks can absorb the deficit -> that would score ~39,031.
2. `s11/quick.py` (cone-restricted, ~170x faster than a full forward) makes such sweeps cheap.
3. If attacking the collision itself: find any second non-bit control for `a34580` or `a33796`
   (i.e. any free input reaching `x_10170` or `x_6858` that is not `x_33129`/`x_37088`).
   `s11/last4.py` is the tool; it currently returns none.

### Toolchain
`cd solve_lab/s9 && python3 atomize.py && python3 poly.py && python3 gates.py && python3 fwd.py`
(rebuilds caches; validates against the raw file -- 0 mismatches over 39,033 equations).
Session-11 scripts are in `solve_lab/s11/`; they import `s9/eff/lib.py`, which chdirs to `s9/`,
so write outputs with absolute paths (see `s11/uv01.py`).

### Git
Branch `claude/math-problem-solving-2y9sl3`.

---
# RESUME — read me first

## STATUS (session 10): best verified **39,026 / 39,033** — and PROVED optimal for this defect placement
Deliverable: `best/new_instance_partial_39026.json`
Verify: `python3 checker.py best/new_instance_partial_39026.json` -> `satisfied 39026/39033 (7 failing)`
Failing lines: `[12231, 12270, 12350, 14584, 18673, 22044, 29125]`.
**Read `S10_EXACT_RESIDUAL.md` first** — it supersedes `S9_STRUCTURE.md` where they conflict.

### The 60-second version (re-derived from the file this session, not inherited)
At the delivered witness there are **7 nonzero atoms** — `22229, 22230, 35758, 35759, 35760,
35761, 35762` — living in **exactly 12 equations** (5 satisfied, 7 failing). Every other one
of the 39,033 equations is already exact in Z. The 12x7 coefficient matrix has **rank 7**, so
all 12 hold iff all 7 atoms vanish. Writing `A = (a22229, a22230, a35758..a35762)`,
`D = x_7068 - x_2099`, `K2 = x_28730 mod p`, the knobs realise **exactly**

    (1)  A1 + 7376877*A7  ==  D   (mod p)      [mod 7376877*p before the free k*p shift, see below]
    (2)  A2              ==  K2   (mod p)
    (3)  A3, A4, A5, A6 free

Both residues verified 0 at the witness. **That is the entire remaining problem.**

### The accounting rule — this sets the endgame
`dim ker(M_S) = 7 - rank(M_S)`, and `c` independent mod-p congruences need `c` free parameters:

| binding congruences | max equations satisfiable | score |
|---|---|---|
| 2 | 5 of 12 | **39,026 (current)** |
| 1 | 6 of 12 | 39,027 |
| 0 | 12 of 12 (`A = 0`) | **39,033 — full solve** |

> **There is no partial credit between 39,027 and a complete solution.** One congruence is
> worth exactly one equation; killing both solves the instance.

Optimality of 39,026 is **proved, not searched**: `s10/lattice3.py` enumerates all 2^12
subsets and tests exact integer solvability (integer kernel by column HNF, then a 2-row
integer linear system for the two congruences). Sizes 12..6: 0 solvable. Size 5: 300+.

### NEW (session 10): the ripple's repair rule is too weak — re-audit the "forced chain" verdicts
`x_7068 += k*p` *appears* to break atoms 29539 and 40826. It does not: both close through
their handles (`x_29967`, then `x_30163`), and the nonzero-atom set returns to exactly the
same seven for k = 1, 2, 7, -3228258 (`s10/repairD.py`). Consequence: `D mod 7376877` is
free, so congruence (1) is only mod p.
> **`lib.ripple` repairs an atom only through its canonical output variable with exact
> division, so it silently misses handle-based repairs. Session 9's §14.4 "the chain is
> FORCED", `chase.py` and `solve_branch.py` all rest on that weaker rule. Re-audit them
> with `s10/tools.solve_lin` (effective-linear solve over ALL variables of an atom).**

### Where the two congruences come from, and the one weak link
* (2) binds **only** because `x_28730` is not free: moving it drags `x_4432`, which breaks
  atom **7930** = `9367949*(x_24548 - x_25442) - x_7927` (15 eqs) and 41512 (`s10/isolate.py`).
  Everything else — `x_642, x_17325, x_9413, x_1329, x_10903, x_29854, x_31864, x_9118,
  x_8731` — is fully free with zero collateral.
  **If atom 7930 can be closed while `x_28730` moves, the score is 39,027 immediately.**
* (1) is the core congruence on `D mod p`.

### Global handle census (session 10) — why both congruences are rigid
`s10/handles.py` over all 42,267 atoms: **1,249** free inputs occur in exactly one atom.
Of those handles **1,240 have granularity exactly `p`**, 9 are dormant/rigid, and
**zero** are unquantised or have any other granularity.

> Every solo handle in the instance is exactly p-quantised, so solo-handle moves shift any
> atom only by multiples of p and **every residue mod p is invariant under them**. None of
> the 12 residual equations contains an atom with a free (Z) handle.

That is the trapdoor's construction verified exhaustively, not inferred from failed searches.

### Beam search under the strong repair rule — and the criterion that matters
Strong repair (effective-linear solve over ALL variables of a broken atom), beam 200,
depth 10, seed forbidden as a repair choice (`s10/beam7930.py`, `s10/beamD.py`).
> **CRITERION: collateral empty AND the residue actually moved.** "Collateral closed" alone
> is not enough — `lib.ripple` recomputes gate outputs and will silently restore the seed.
> `x_2099 / x_37158 / x_22542 += 1` all report a clean close that is an **artefact**
> (ripple rebuilt `x_2099` from definer 29090; `D mod p` and the score are unchanged).

Real results: `x_28730 += k*p` closes (via `x_7927` then `x_11052`) but leaves `K2` fixed;
every `d` with `d % p != 0` on either congruence fails, the collateral walking a ladder
(11625 -> 11624 -> 11621 -> 30238 -> 24948 ... ; 27314 -> 29539/40826 -> 19482 -> 19480 ...).
So Session 9's verdict survives the stronger rule — only the *evidence* in §3 was wrong.

### PART II (same session): the GLOBAL attack — and the first crack in the design

**Forward-eval frame** (`s10/forward.py`): take the witness's FREE INPUT values and
forward-evaluate every gate. Result: **6 nonzero atoms, all CHECKS, zero broken gates**
(37 failing, 38,996). Each contains a FREE input:

    a7930   x_24548 == x_25442 (mod p)     a35759  x_9118 == 0 (mod p)
    a29539  x_14853 == x_1308  (mod p)     a35760  x_8731 == 0 (mod p)
    a40826, a41512   big checks, 1 equation each

This explains the 39,026 witness: it violates five GATE atoms on purpose so `x_1308`
and `x_25442` land on `x_14853`/`x_24548`.

**The point is RIGID** (`s10/ad.py`, `s10/closure.py`, `s10/rankdef.py`). Exact
reverse-mode AD mod p (validated against finite differences) + the closed constrained
system gives **rank 79 of 79 columns — full column rank, zero null space — with 6
independent inconsistencies**, and 0 degenerate rows (so not Session 9's square bug).
With all 256 message bits relaxed to GF(p) the closure is 2,352 x 710 and STILL
inconsistent. A full 1,156-way single-bit scan with genuine forward-eval finds nothing.
=> local and first-order methods are definitively dead.

### THE CRACK: the p-wire is not rigid (`s10/wirekernel.py`)
Every handle enters as `wire * handle` with `wire` one of **220 variables equal to p**,
held by ONE bare pin `a37694 = x_26064 - p` in only **12 equations**. Setting the wire
off p flips the whole census: **1,240 handles go from granularity p to granularity 1.**

Writing `w_u = p + d_u`, every wire-identity atom is linear and homogeneous in d, so the
219 equations containing them give `M d = 0` in Z^220:

>  **rank(M) = 217 of 220  ->  KERNEL DIMENSION 3.  The wire can deform for free.**
>  161 of 220 members have gcd 1 (can take ANY value), including handle multipliers
>  `x_11360, x_28599, x_17499, x_22665, x_28961`.  The root `x_26064` has gcd 0 — FIXED.

`wire = 1` reaches **39,020 with only TWO nonzero atoms** (`a37694`, `a39417`);
its 13 failing equations contain only wire-copy atoms and boolean pins, no free inputs.

**Why it does not pay yet** (`s10/deform2.py`): applying a kernel vector and re-solving
handles restores 3,346/3,349 product gates; of 235 broken atoms **215 are wire copy atoms
whose equations still cancel by construction**, so the genuine cost is **20 atoms** —
13 of them multi-wire monomials `w_i*w_j`, whose invariance is QUADRATIC in d
(`p(d_i+d_j) + d_i d_j = 0`) and so invisible to the linear kernel. Net ~38,981, ~39,018
after closing the checks.

**Deformation route CLOSED, and exactly why** (`s10/multiwire.py`): the multi-wire
monomials are almost all SELF-pairs `w_i*w_i`, coming from the degree-4 SQUARE checks.
Their condition `w_i*w_j = p^2` degenerates to `w_i^2 = p^2 => w_i = +-p`, and the
constraint graph is a SINGLE non-bipartite component over all 220 members with
self-loops — so every member is forced back to +-p, and -p keeps granularity p.
The crack is real but bounded: the LINEAR wire-identity system allows a 3-dim
deformation; the QUADRATIC square checks pin it. Each such check is in exactly ONE
equation (hence 38,981, ~39,018 after closing the six checks — still above 7 failing).

**SETTLED, and the route is comprehensively closed:** counting multiplicity over all
42,267 atoms, **220 of 220 wire members appear SQUARED somewhere — none is unsquared.**
So every member is pinned to +-p by a degree-4 square check. The 3-dimensional linear
kernel is real but is entirely absorbed by the quadratic constraints. The p-quantisation
of all 1,249 handles stands, and with it Part I's two congruences and 39,026's optimality.

### START HERE NEXT SESSION
1. The wire is pinned by ~20 degree-4 SQUARE checks, each living in exactly ONE equation.
   That is the cheapest guard found in any session. Attack those squares directly in
   equation space: a square check `E^2` in one equation costs 1 equation if violated, so
   breaking k of them to free the wire costs k. Compute, for each squared wire member,
   how many square checks must break to release it, and compare against the 7 equations
   the current branch pays.
2. Equation-space compensation for `a37694`'s 12 equations using deformed copy atoms
   (same linear algebra, pin's row moved to the RHS).
3. LLL-reduce the 3-dim kernel lattice (current basis has ~325-digit entries).

### PART III (same session): the BUDGET attack — the trapdoor priced exactly

**Reframe.** The current branch pays 7 failing equations, so any structural violation
costing <= 6 beats it. Price list (`s10/pricelist.py`, price = equations a check lives in):
a degree-4 square check = **1 equation**; `a40826`/`a41512` = **1 each**; a wire copy
atom = 12-14; the wire root pin `a37694` = 12; a boolean check = 13-15.
(Earlier greedy searches scored by NUMBER OF NONZERO ATOMS -- wrong objective. Fixed.)

**The six certificates** (`s10/certs.py`, via `[A | b | I]` elimination):

    rank(A) = 79/79        INCONSISTENCY CERTIFICATES: 6
    cert 0: cheapest (1, 36602) (1, 37887)      cert 3: cheapest (1, 41507) then 10
    cert 1: cheapest (10, 2423) (10, 21617) (10, 31670) (11, 19297)   <== NO cheap member
    cert 2: cheapest (1, 41400) then 10          cert 4: cheapest (1, 41827) then 10
    cert 5: cheapest (1, 11007) (1, 25676) (1, 39800) (1, 42245)

> **Five of six certificates cost 1 equation to hit. Certificate 1 cannot be hit for
> less than 10, and its members 2423/21617/31670/19297 are exactly the section-15.2
> trapdoor chain.** Min-cost hitting set = **15 equations vs a budget of 7** — the
> design carries a margin of exactly 8. Note 2423, 31670, 19297 each hit FOUR of the
> six certificates, so the optimum has the shape `cost(hub) + 1 + 1`.

**Hub compensation does NOT pay (a claim I made and then corrected).** `truecost.py`
first reported hub a31670 at true cost 1 (9 of 10 equations "compensable"). That
heuristic counted the EXISTENCE of a helper per equation, not that a helper's value is
one number shared across all of them. Exact computation (`hub31670.py`): the region is
41 equations x 16 atoms with FULL column rank 16; with a31670 != 0 forced there are ~2
free parameters against 10 equations so at most ONE is savable, and recruiting the other
helpers drags in 31 more equations. Hub cost ~9, and the 15-equation hitting set stands.

**Also priced and closed:** wire deformation re-tested with square-check repair
(`deform3.py`) -- zero admissible moves across all 6 kernel directions, so section 18
survives the stronger (necessary-condition) test. The forced OR gate is a non-lever:
`x_9274 = 1` holds automatically at all-bits-zero and that branch is 38,871; only two
boolean free inputs are set at the deliverable (x_2081, x_24601).

> **The productive question is now exactly one: can certificate 1 be hit for under 9
> equations?** Everything else in the instance is already cheap.

### PART IV: every route PRICED — the margin is uniformly 6

**Correction to Part II:** 161 wire members have kernel-gcd 1, but reachability != usable
magnitude. Hitting `d_u = 1-p` needs kernel coefficients ~10^250, blowing other coordinates
to ~10^575. Measured (`s10/deform_solve.py`): raw kernel directions give |w_3915| =
|w_11360| = **325 digits** -- handle granularity far WORSE than p -- and score 38,990.
Only a SHORT kernel vector would matter.

**Complete price table** (`s10/memberprice.py`; cost = identity eqs + square-check eqs):

    cheapest member overall   x_15413   13   (10 identity + 3 square)
    uniform wire shift                  13   (root pin 12 + a39417 1)   -> 39,020
    cheapest USEFUL multiplier x_3915   15   ( 9 identity + 6 square)
    x_11360 36   x_22665 48   x_14466 54   x_15616 56   x_28961 208   x_28599 292
    kernel deformation                 ~20   (and useless magnitudes)
    certificate hitting set             15
    ---------------------------------------------------------------
    current branch (give up)             7   <== the deliverable

> **No route into the design costs less than 13 against a give-up cost of 7. The margin
> is 6, and it is UNIFORM across every attack surface measured: wire, certificates, hubs,
> bits.**

### THE ONLY REMAINING QUESTION
Hit certificate 1 for under 9 equations (cheapest member: 10). Two open sub-questions:
1. Is there a SHORT vector in the 3-dim wire kernel (LLL over ~325-digit entries)? Support
   <= 5 would cut the deformation cost from ~20 to ~5 and take the wire route under budget.
2. Can certificate 1 be hit by a variable the closed 79-column system never reached?

### PART V: the INVARIANT — 7 is pinned

**No sparse kernel vector** (`s10/sparsekernel.py`). A kernel vector with support 2-3 would
cost ~3-6 and win. Viewing the kernel as a 220x3 matrix K (row = 3 basis values per member),
a vector supported on S exists iff all rows outside S are coplanar in Q^3. Measured:
kernel-zero rows 3; 215 distinct directions among 217 nonzero rows; largest rank-1 cluster
2; **largest coplanar set 4** => sparsest kernel vector has support **>= 213 of 220**.
Every free wire deformation moves >= 213 members and breaks essentially all their square
checks. The wire route's floor of 13 (uniform shift) stands.

**Seven is an INVARIANT, not a property of the placement** (`s10/eighth.py`,
`s10/invariant.py`; exact integer subset enumeration each time):

    extra atoms          region   satisfied   FAILING
    (none)                  12         5         7
    35756                   15         8         7
    35754                   17        10         7
    35756 + 35754           18        11         7

> **Every extra free parameter buys exactly as many equations as it drags in.** Only two
> adjustable atoms in the entire instance even touch the twelve equations (35756 overflow 3,
> 35754 overflow 5); both tested alone and together. This is the exact reproducible form of
> what earlier sessions called the "conserved obstruction".

**Consolidated:** give-up 7; invariant 7; uniform wire 13; cheapest member 13; certificate
hitting set 15; kernel deformation ~20. Three independent lines -- equation-space lattice,
GF(p) closure + certificates, wire geometry -- all return 7, every route costs >= 13.

### PART VI: a real gap found in Part I's model — and the invariant SURVIVED it

**Forensics on the setter's constants: no backdoor** (`s10/forensics.py`). gcd of all
2,817 large literals = 1; gcd of all pairwise differences = 1; no constant equals either
residue mod p; `D0/K2 mod p` is a full 253-bit number and `D0 != k*K2 (mod p)` for k < 60.
Both residues ARE quadratic residues -- the only structure found, and not exploitable.

**THE GAP.** Part I called `x_28730` "not free" because moving it drags `x_4432` and breaks
atom 7930. That was an artefact of HOW I moved it -- together with `x_4432`, to hold
`a22231 = 0`. **`a22231` need not be zero.** Moving `x_28730` ALONE:

    a22230 changes by +d ; a22231 changes by -d ; x_4432 UNTOUCHED, no collateral
    (verified for d = 1, 2, p)

and a22231's ten equations lie ENTIRELY inside the twelve -- zero overflow. Correct model
has EIGHT atoms and ONE paired congruence:  `A1 + 7376877*A7 == D0 (mod p)`,
`A2 + A8 == K (mod p)`, A3..A6 free. Exact optimum: **6 of 12**, not 5. Constructed and
verified end to end (`s10/build27.py`): all eight atom values realised EXACTLY, x_4432
untouched, six equations satisfied -- the first time the region passed 5.
**But `a37887 = R^2` then lights up and breaks eq 8680, restoring 7.**

Extending further (`s10/kill37887.py`): R is a linear atom combination
(`a22231 + 6*a22232 + 15*a22233 - 21*a22234 - 13*a22235 + ...`), and a22232..a22235 move
in pairs, so R = 0 is reachable (a22231 = 9d + 34e, gcd(9,34)=1). 12-atom model with R=0
imposed exactly: region 12 -> 16, max satisfied 5 -> 9. **Failing: 7.**

    model                          atoms  region  satisfied  FAILING
    Part I baseline                   7      12       5         7
    + 35756                           8      15       8         7
    + 35754                           8      17      10         7
    + both                            9      18      11         7
    + a22231  (the gap, fixed)        8      12       6         7   <- a37887 costs 1 outside
    + a22231,a22232..35, R=0         12      16       9         7

> **Six independent placements -- including one built specifically to exploit a genuine
> error in the earlier model -- all return exactly 7.** Every degree of freedom added is
> matched, to the equation, by the equations it drags in. The invariant is not an artefact
> of the defect set: I found and fixed a real gap in that choice and the number did not move.

### PART VII: number theory CLOSED; and the root pin costs 1, not 12

**secp256k1 hypothesis refuted** (`s10/curve.py`): (D0,K2) is not on y^2=x^3+7, neither is
a valid x-coordinate, n/G_x/G_y are absent (p itself IS present), 7870/15734 constants have
(c mod p) a valid x-coord vs random expectation 7867, and 507/7999 multipliers are prime vs
~470 expected. Exactly random on every axis -- p is a convenient modulus, not a curve.

**Rational reconstruction: no structure** (`s10/ratrec.py`). Every residue (D0, K2, D0/K2,
K2/D0, D0*K2, D0+-K2, 1/D0, 1/K2, HUGE mod p, C1 mod p) returns MAXIMAL 38-39 digit a and b
-- right at the sqrt(p/2) bound, i.e. no small rational. gcd(HUGE,C1)=1; HUGE//p and C1//p
are unremarkable 12-13 digit numbers; HUGE != k*C1 (mod p) for k<200. Of 2,815 constants
exceeding p, ZERO have residue < 2^80 and ALL 2,815 residues are distinct. The seven
residual equation values have gcd 1.
> **No arithmetic backdoor exists. The number-theoretic line is closed.**

**NEW: the wire root frees for ONE identity equation, not twelve** (`s10/rootfree.py`).
e_root lies in the identity row space; writing e_root = y0^T M, supp(y0) = {eq 37257} -- a
SINGLE equation. Equation 37257 is the unique identity equation whose wire content is the
root pin ALONE; in the other eleven, a37694 sits beside copy atoms that absorb it under a
non-uniform deformation. Constructed (`s10/freeroot.py`): dropping 37257 gives rank 216 and
a 4-DIMENSIONAL deformation space, all four directions moving the root, with all 218 other
identity equations holding. **Still does not pay:** entries are ~324 digits with support
217, so 17 non-copy atoms break -> 38,984. Identity cost 1 + square-check cost ~12 = 13,
the same floor from a third independent direction.

### PART VIII: every door opened and priced

**A counting error, corrected.** Part III said "6 independent inconsistencies". b is a
SINGLE column, so rank([A|b]) - rank(A) <= 1: there is ONE obstruction and the six were six
witnessing rows of it. That raised a real hope (one dropped row might suffice). Tested
exhaustively (`s10/singledrop.py`): **all 128 single-row drops fail; all pairs among the 30
cheapest fail.** The obstruction survives removal of any one row.

**The region, closed exhaustively** (`s10/regionknobs.py`). `eighth.py` defined "adjustable"
as carrying a SOLO handle -- exactly why it missed a22231. Correct definition: a variable is
a region knob if moving it changes no equation outside the twelve. Scanning every variable
in every atom of the twelve:

    variables with ENTIRE footprint inside the twelve : 9
      x_642, x_1329, x_8731, x_9118, x_9413, x_10903, x_17325, x_29854, x_31864
    they reach exactly [22229,22230,35758..35762] -- the Part I atoms
    next cheapest: x_28730 at 1 outside (eq 8680 via a37887); everything else >= 3

> Nine knobs, exactly the seven atoms. No hidden freedom. Region CLOSED.

**Boolean branches, exhaustive in the WITNESS frame** (`s10/bitwitness.py`). All 1,156
boolean free inputs flipped in the deliverable's frame with exact repair: best 20 all give
failing = 7 with identical region; x_4287 -> 34, x_24601 -> 83, x_2081 -> 106.
> No boolean flip improves on 7. Branch structure EXHAUSTED.

**NEW FREEDOM FOUND -- and it is inert** (`s10/cycles.py`, `s10/slide.py`). fwd.py covers
only 29,675 of 31,475 defined vars: 1,800 sit in gate CYCLES. A cyclic block is a system,
and a singular one has a solution FAMILY that forward-eval silently collapses to a point.
Measured: 40 non-trivial SCCs, all size 2, **local Jacobian rank 1 of 2 in every one =>
kernel dim 1 each = 40 free parameters invisible to any local method.** Sliding along all
40: no new nonzero atoms, failing stays 7, D0 and K2 UNCHANGED (8 are literally inert).
The freedom is real and ORTHOGONAL to the obstruction.

    door                                   status                cost
    region knobs beyond the 9              closed exhaustively    --
    boolean branches (1,156, witness)      closed exhaustively    >= 7
    single-row / cheap-pair sacrifice      closed exhaustively    none work
    cyclic freedom (40 params)             OPEN but inert         0 gain
    number theory / curve / ratrec         closed                 no structure
    wire (uniform/member/kernel/root)      closed                 >= 13
    certificate hitting set                closed                 15
    give up (the deliverable)              --                     7

### PART IX: the sacrifice question, answered exactly

**Reformulation.** Dropping rows S leaves A_{-S}x = b_{-S}; its left null vectors, extended
by zeros on S, are exactly the y in leftnull(A) with supp(y) disjoint from S. So
**consistent after dropping S  <=>  t in colspace(Y[:,S])**, with Y a basis of leftnull(A)
and t = Y.b. Each test becomes a 49 x |S| rank check instead of a 128 x 80 elimination
(`s10/budget6fast.py`).  closed system 128x79, rank 79, leftnull dim 49, t nonzero.

**The minimum sacrifice is exactly THREE rows** (sizes 1 and 2 impossible):
`{a3578, a26731, a35759}` = setter load pin (price 14) + mirror (16) + one currently-failing
check (7). Union = **37 equations**, i.e. score 38,996 -- exactly the forward-eval floor.
The cheapest *sized* solution is the most expensive kind.

**Budget <= 6 exhausted** over all 46 rows priced <= 6, cost-pruned:
size 1: 46 sets -> none; size 2: 1,081 -> none; size 3: 16,261 -> none;
size 4: 179,446 -> none; size 5: 1,550,200 -> none; size 6: DFS (see runs/).
> Too few rows is impossible; cheap enough is unreachable. Sacrifice route CLOSED.

**Every door now has a number:** give up 7; invariant 7 (6 placements); min sacrifice >=3
rows costing 37; wire routes 13; certificate hitting set 15; kernel deformation ~20; region
knobs beyond the nine: none; boolean flips (1,156, witness frame) >= 7; cyclic freedom (40
params) real but inert; number theory: no structure.

### PART X: the message space, EXHAUSTIVELY CLOSED

**Global rigidity, tested bluntly** (`s10/randomize.py`). Randomising the non-boolean free
inputs (1 / 10 / 100 / 1000 / ALL 6,117) gives 37 / 148 / 1084 / 5219 / 7355 failing; best
over every randomisation = 37 = the base. **The four core checks 7930, 29539, 35759, 35760
fail from EVERY starting point.** The residual is pinned against the non-boolean inputs
GLOBALLY, not merely locally.

**The "256-bit codeword" collapses to 5 dimensions** (`s10/bitgroups.py`). Exact AD gradient
of every failing check w.r.t. all 1,156 boolean free inputs:

    boolean inputs moving ANY failing check : 128 (not 256)
    distinct signature vectors              : 5
    multiplicities                          : 75, 50, 1, 1, 1
    => reachable message states = 76*51*2*2*2 = 31,008

Within a group bits are interchangeable, so only the COUNT matters. Enumerable in a second.

**Swept, with the model validated where it holds.** A first sweep claimed 2 of 6 zeroable;
CONSTRUCTING it refuted that (`msgverify.py`: 62 failing, nothing zeroed) -- its two bits
were x_2081 and x_4287, the structural MUX controls, where b*(X-HUGE) has X depending on b
so linearity fails. Re-run properly (`msgvalid.py`): linearity VALIDATED exactly on both
large groups (bits x_91, x_47); sweeping all 76*51 = 3,876 states gives histogram
**{0: 3876}** -- the 125 ordinary load bits cannot zero a SINGLE failing check. The only
bits with leverage are the 3 structural controls x_2081, x_4287, x_13195 = the branch flips
already measured (34 / 83 / 106 failing in the witness frame).

**Sacrifice route exhaustively closed** (`budget6fast.py`): 10,917,019 within-budget sets
tested, NONE restore consistency; minimum sacrifice is 3 rows {a3578, a26731, a35759} at
cost 37; sizes 1 and 2 impossible.

### Do NOT redo
- The MUX branch (`x_4287 = 1`). It **does** zero all seven residual atoms simultaneously
  (`s10/muxzero.py`) but leaves 8 collateral atoms -> 44 failing; best repair 38,991. Its own
  load pins have only p-quantised handles (`x_27676 = p*x_6504`, `x_7574 = p*x_26658`), so
  `x_31861`, `x_14865` stay pinned mod p: 4 mod-p conditions vs 2 free residues — the same
  deficit of 2, relocated.
- Greedy equation-space repair from the MUX state (it just switches the MUX back off -> 39,008).
- Treating `x_28730` as free: `s10/lattice.py` finds a 6-subset and `s10/construct.py`
  realises the target atom vector **exactly**, but it scores 39,011. Good model check, dead end.

### Next experiments, in priority order
1. Re-audit every session-9 forced-chain verdict with handle repair (see the warning above).
2. Atom 7930: enumerate everything that moves `x_24548` and `x_25442`, look for a joint move
   that lets `x_28730` float (`s10/atom7930.py`). Worth exactly one equation -> 39,027.
3. The two residues are the whole problem — attack them, not assignments:
   `D0 mod p = 61705020361863629770768910187978745858728889529652486596432934143473517757811`
   `K2      = 33310166114805471624282140578459083391052142224394967852279417483154815501175`

### Toolchain
`cd solve_lab/s9 && python3 atomize.py && python3 poly.py && python3 gates.py && python3 fwd.py`
(rebuilds all caches; `atomize.py` validates the decomposition against the raw file —
**0 mismatches over 39,033 equations**). Session-10 scripts are in `solve_lab/s10/`
and import `s9/eff/lib.py`; `s10/tools.py` has the stronger `solve_lin` repair.

### Git
Branch `claude/read-prompt-xpo2kf`.

---
# RESUME — read me first

## STATUS (session 9): best verified **39,026 / 39,033**
Deliverable: `best/new_instance_partial_39026.json`
Verify: `python3 checker.py best/new_instance_partial_39026.json` → `satisfied 39026/39033 (7 failing)`
Independent check: `python3 s9/verify_ast.py best/new_instance_partial_39024.json` (AST walk, no eval/regex).
Failing lines: `[12231, 12270, 12350, 14584, 18673, 22044, 29125]`.

**How 39,026 was reached:** use EFFECTIVE footprints, not syntactic ones — a variable multiplied by a
currently-zero variable has no effect. That admits `x_9118` and `x_8731` as knobs, and by CRT
(`gcd(5113045, p) = 1`; `x_8731` moves atom 35761 in steps of 1) two of the four binding congruences
dissolve: lattice invariant factors go `[1,1,P,1,P,1,P,7376877P]` → `[1,1,P,1,1,1,1,7376877P]`,
so max-simultaneously-zeroable rises 4 → 6 and failing drops 9 → 7. Exhaustive over all 2^13
subsets: no 7-subset is integer-solvable, so 7 is optimal for this defect placement.

> **RETRACTED:** an earlier version of this file claimed 39,022 was a proved local optimum.
> That was wrong — see `S9_STRUCTURE.md` section 6. The proof assumed every atom outside the
> defect set must vanish. **It must not.** An equation is zero iff its *linear combination of
> atoms* is zero, so any atom whose whole equation footprint already lies inside the failing set
> is a FREE knob. Atom 22230 (`x_28730 − x_17499·x_9413`) is exactly that — `x_9413` and
> `x_28730` appear in no other atom — which frees `x_28730` from the lattice `p·ℤ`. Extending
> along the `35754…35762` ladder gives 5 knobs over 13 equations, 4 of them simultaneously
> zeroable: 11 + 2 − 4 = 9 failing. **Work in EQUATION space, not atom space.**
**Read `S9_STRUCTURE.md` first** — it supersedes the older analyses below on every point of conflict.

### The 60-second version
Exactly **3 atoms** are nonzero at the partial (22229, 22231, and the square 37887 whose root
contains 22231). Since `x_28599 = x_17499 = p = 2^256−2^32−977` exactly, the
entire residual system is two congruences:

    x_7068 ≡ K1 (mod p)      and      x_4432 ≡ K2 (mod p)

with `x_17325`, `x_9413` as free quotient handles and `K1, K2` the constants that atoms 3576/3578
pin `x_6418`, `x_12553` to. Everything else in the 39,033 equations is already exact in ℤ.

### Three results from session 9 that change the picture
1. **METHODOLOGICAL BUG FOUND — re-check anything built on a Jacobian.** 783 check atoms are
   perfect squares `E²`; at the partial `E = 0`, so a finite-difference row for `E²` is *quadratic*
   (`c²δ²`), not linear. Feeding squares into a linear/Newton/null-space solve solves the wrong
   system. `s9/roots.py` extracts all 783 roots. After the fix the reported mod-p inconsistency
   moved off a *spurious* square (42245) and onto the *genuine* core atom (19297). Earlier
   sessions' "inconsistent Jacobian" verdicts were reading the artefact.
2. **The core has a SECOND branch, never previously recorded.** Core ⟺ `S ≡ T ≡ 0 mod p`, and
   eliminating `w` gives `u²·(A·c² − B²) ≡ 0`. So it is *not* only `u ≡ w ≡ 0`: the alternative
   `A·c² ≡ B² (mod p)` frees `u` completely. Blocked only by the mod-p pins on `x_22162`,
   `x_30213`, `x_16742`. **This is the one door opened and not closed — start here.**
3. **39,022 is a local optimum, proved (not just observed).** Handles absorb exactly the
   p-multiples, so the reachable defects are `A ∈ D1+pℤ`, `B ∈ D2+pℤ`. Each failing equation is
   `m·(c₁A + c₂B)`, so it can vanish only if `c₁D1 + c₂D2 ≡ 0 mod p` — checked for all nine
   (none), and eqs 8680/29125 need `B = 0` exactly (impossible). All six alternative defect
   placements cost 23–29 failing equations vs 11. **Improving the score requires cracking the core.**

### BIGGEST NEW LEAD (session 9): the core is NOT a wall — both cores can be zeroed
`x_8599 = 1` (88 of the 1,156 boolean free inputs do it, keeping `x_21839 = 1`) reroutes
`x_12186` from the computed `x_30454` to the **free input `x_5096`**. Set `x_5096 := K1`,
`x_14853 := x_12186`, let the gates close C1/C2 → `u = w = S = T = L1 = 0` exactly: **the core
(19297/19299/30984) is satisfied simultaneously with C1 and C2 for the first time.**
It lights a second core of identical shape (26733/28438/32342, gated by `x_38170`), and *that one
is zeroable too*: ~90 boolean free inputs shift its controls by exact complements mod p
(`u' + δ = p` and `w' + δ' = p` exactly). Two-bit constructions (e.g. `x_2527` + `x_1502`,
`s9/construct3.py`) leave **both cores clean**, with the entire residual being activated load pins.
Closing those pins moves variables in the first core's cone and re-lights it (11→16→13, stalls).
**So the real invariant is the pin/mirror cascade, not the core.** Attack that next:
enumerate each bit's pin set and search for a bit subset whose pins are jointly satisfiable.

The core's other branch (`A·c² ≡ B² mod p`) is **CLOSED**: it needs `A` to be a quadratic residue;
`A` is fixed by the pins and is a non-residue for both reachable `x_12186` residues (a ~50% sample
of shifts *would* be residues — the setter's pin of `x_14853` to `K1` is what closes it).

### PROBLEM SIZE / SOLVER SUITABILITY → read `REDUCED_PROBLEM.md`
Irreducible residual = **2 congruences mod p** (512 bits). MEASURED: only **2 of 1,156** boolean
free inputs move those residues (the two MUX controls) and their deltas are identical → **rank 1 of
2** over GF(p); only 4 of ~6,117 integer free inputs move them, all pinned mod p in a complete
solution. **The 256-bit message is decoupled from the verification — there is no subset-sum kernel.**
So the hardness is arithmetic, not combinatorial: quantum annealers are the wrong machine (needs
~10⁵–10⁶ logical binaries and ~512 bits of coupler dynamic range vs ~4–6 effective bits available).
Right tools: Legendre/QR conditions, LLL/CVP, or the §12 identity.

### THE LOCK, IN CLOSED FORM (read S9_STRUCTURE.md section 12)
`x_12186 ≡ x_22649 (mod p)` exactly 1:1 and `x_22649` is FREE — it is the circuit's output wire
(found with `s9/modtrace.py`). Moving `{x_22649, x_22152, x_14853, x_7068}` together by
`δ = K1 − x_12186` **kills chain 1 outright** (atom 22229 = 0) while holding the core's `u = 0`.
What stops it is load pin **31670**, gated by bit `x_24601`, which pins the output wire to exactly
the original `x_12186` — a constant differing from `K1` by precisely `D1`.

Flipping the MUX control **`x_4287 = 1`** makes `x_2099`/`x_19964` read the FREE inputs
`x_9118`/`x_8731` instead of the pinned `x_6418`/`x_12553`, so **atoms 22229 AND 22231 are both
zero simultaneously** — the first time. Its price is three new loads that collapse (3-from-2) to
`x_4306 ≡ x_27177 ≡ 0 (mod p)`, which IS solvable — but `∂x_27177/∂x_8731 = 0`, so `x_9118` mod p
is uniquely forced to `33371159155735472537534252650716501592825364489306217536352743247010353604716`,
while the mirror+core chain needs it to equal the pinned `x_12186`
(`82007976…230357` if `x_24601=1`, or `0` if `x_24601=0`). **Neither matches — that mismatch IS the
trapdoor.** Any further attack should target that single scalar identity, not search assignments.

### THE ALL-ZERO BRANCH IS CLOSED (rigorously) — `S9_STRUCTURE.md` §14
`x_2081 = 0, x_24601 = 0` is the only quadrant where the MUX output and the circuit's output wire
are consistently both 0. Driven through the staged construction (`zero.py` → `zero9.py`, `chase.py`)
it reaches **3 atoms / 17 failing** with `A = B = 0`, `x_15298 = 0` (core dead), `u = w = 0` and the
forced OR gate satisfied — every condition that blocks the main branch. It stalls because:
(1) only the 256 **gate-bits** can satisfy the forced OR `x_9274 = 1` (all 900 pin-free booleans
fail, exhaustively verified); (2) only `x_47` makes the drivers `x_16742 → x_18956` and
`x_14681 → x_24468` live, i.e. only `x_47` lets setter pins 688/1618 close; (3) the mirror chain
that closing them requires is **FORCED** — every target has exactly one 1:1 driver, no branching;
(4) it terminates exactly on `x_24221`/`x_25477`, **`x_47`'s own pinned pair**.
The only bit that unlocks the setter pins is the bit whose own pins the unlocking chain must violate.

### Do NOT redo
- Greedy ripple-repair from the canonical orientation (converges to the core in 3 rounds, 39,013).
- Single-bit flips: all 1,156 boolean free inputs scanned; only `x_2081`, `x_24601` deactivate the
  core (`x_15298 = 0`), both strictly worse locally (17/19 residuals vs 4).
- Repairing via `x_7068`/`x_4432`/`x_6418`/`x_12553` (strategies A/A1/A2/B1/B2 in `s9/drive.py`):
  39,002–39,013, all worse than 39,022.
- Mod-p linear solve over all 7,273 free inputs at fixed bits — inconsistent at S0 *and* S1 even
  with the corrected root model (270-row certificate).

### Next experiments, in priority order
1. **The pin/mirror cascade** (see the lead above) — with both cores now zeroable this is the only
   remaining invariant. Per activated bit, enumerate its load pins `bit·(x_B − HUGE) = s·x_C`, then
   search for a bit subset whose pins are jointly satisfiable (exact-cover over the pinned free
   inputs). `s9/pinclose.py` is the naive greedy version and it stalls; the set-cover view is new.
2. `s9/construct3.py` reached both-cores-clean; extend it rather than restarting.
3. Full two-bit flip scan (~667k pairs). Session 9 only ran 88x14 = 1,232 targeted pairs.
4. Quadrant re-solve at `x_2081 = 0` / `x_24601 = 0` with the corrected root-based residual model.

### Toolchain (`s9/`, self-contained, caches regenerable)
`cd solve_lab/s9 && python3 atomize.py && python3 poly.py && python3 gates.py && python3 fwd.py`
then e.g. `python3 state0.py`, `python3 drive.py A`, `python3 bitscan.py`, `python3 newton.py`.
`atomize.py` validates the whole decomposition against the raw file (0 mismatches over 39,033 eqs).

### Git
Branch `claude/read-prompt-324ju4`.

---
# RESUME — read me first

## ⚠️ CURRENT (re-randomized) INSTANCE — 39,013/39,033; core REDUCED but not cracked
The EQUATIONS.txt in the repo is a NEW re-randomized instance (39,033 eqs). Full analysis in
**`NEW_INSTANCE_STATUS.md`** and **`CORE_REDUCTION.md`** (read both). Best verified partial:
**39,013 / 39,033** (`best/new_instance_partial_39013.json`, quadrant (1,1)).

### BREAKTHROUGH (latest session): the 20-equation core is fully reduced.
All 20 remaining verifier squares = integer combos of three monsters M1,M2,M3, wiring-defined by
two base gates S=x_35389, T=x_6671. Core ⟺ M1=M2=M3=0 ⟺ **S≡0 and T≡0 mod p** (quadrant 1,1),
then set private quotient handles x_30317,x_2936,x_5146. See CORE_REDUCTION.md for the full chain
down to control differences x_29322=x_14853-x_12186, x_3558=x_24908-x_16742.

### WIRE ESCAPE (the actionable path — read CORE_REDUCTION.md's wire section)
Agent B PROVED the wire=p mod-p system is a rigid isolated point: rank(J_sat)=3035/3036 active
cols, null space dim 1, 19/20 core conditions directly contradict the wiring. The ~5547 "dead"
free inputs feed ONLY products against the p-wire (wire=p≡0 mod p → wire·handle≡0). This is WHY
the witness is unreachable on the wire=p branch.
THE ESCAPE: the 220-var identity wire (root 38100, forced to p by x_26064's single-var atom) is,
only "meant to vanish" — the witness lives on the wire≠p branch. Set the whole
wire = sign·1: then wire·handle = handle ≢ 0 mod p, ACTIVATING all ~5547 quotient handles → a huge
new null space. Core collapses: M1=L1+x_30317→x_30317=−L1 (trivial), M3→x_2936=537773·L3 (trivial),
M2→x_5146=L2/6672769 (needs 6672769|L2 — a 2^23 modulus; L2 mod 6672769 is message-controllable:
4239005 at 39013, 2032135 at 39018). Only ~13 "active unpackings" break (wire members as standalone
terms + (x_26064−p) checks + wire·x_31342 products): [8429,11166,11915,12594,23869,25313,26785,
31400,32300,36106,36767,37257]. NEXT: build the wire=1 global solve over the activated handles
(Dixon lift), heal the 13, set L2≡0 mod 6672769, verify. Agents B (wire=1 consistency via
tangent-linear) and E (wire=1 construction) are on it. Best partial: 39,018 (best_agentD_39018.json).

### OBSTRUCTION on the wire=p branch (superseded by the escape above): residues are pinned.
The sparse wiring solution is unique (only 30 slack inputs nonzero, rank 30). Sparse-witness
null-space solves up to the FULL closure (6,114 inputs, 7,119 constraints) are INCONSISTENT —
S,T residues are linearly pinned by the wiring; no local move reaches the core. Both quadrants
(1,1) and (0,0) have the SAME 20 hard squares via different branches; multi-role control
variables couple the core to the whole system. Remaining paths: a global nonlinear solver
(basin-hopping / large mod-p Gaussian) or the setter's witness. NEXT: try a global linear solve
from the all-zero (0,0) point (cleaner linearization: all products vanish → residue conditions
become linear), or re-attempt with a genuinely different activator/quadrant that unpins S,T.

Definitive findings (exhaustive):
- Gate DAG fully ACYCLIC; forward-eval from free inputs satisfies all wiring automatically.
- The ONE large identity wire (220 vars) is PINNED to **p = 2^256-2^32-977 = the field
  prime** (x_26064=p, appears 13x); it is the twist multiplier. No free wire exists (unlike the
  solved instance, whose wire was FREE — that was the whole solve).
- The remaining obstruction is a **256-bit boolean codeword message** (disjoint 178+78-bit cones
  of control bits x_7715,x_34554). x_9274=OR(controls)=1 is FORCED, so activation is mandatory.
- Each set bit triggers a huge additive load; the GF(p) load matrix over the 256 bits is FULL
  RANK (no bit self-cancellation). Data can absorb mod p but the ℤ-lift imposes p-divisibility
  carries -> every wrong-message data solve is ℤ-inconsistent (SNF pivot=p); iterative repair
  diverges (27->300+). Inhomogeneous GF(p) message solve inconsistent (rank 3). No small vinegar
  linearization (21,922 vars bilinear). => structured GF(p) codeword/MQ; needs the setter secret.
- Key tools: build_twist.py (activate+route MUX), newton2.py (simultaneous absorber solve),
  p_message.py (GF(p) message solve), scc.py, localize.py, scan_bits.py, p in huge_consts.json.
- Untried heavy attack: lattice/LLL on the mod-p codeword (standard, likely resisted by design).

---
# RESUME — read me first

## (historical) prior best partial
Best verified partial was **39,019 / 39,031** (`best/best_partial_39019.json`).

## SESSION 7 — SLACK-ACTIVE SOLVER BUILT; obstruction reduced to R=0 (read NOTEBOOK Session 7 tail)
The slack-active evaluator EXISTS now: `slack_active.py` (freeze x_24026:=x_18274-x_35186,
x_27116:=x_17728-x_1642 with x_12779=1 via a single 22-side bit e.g. 1858). It makes BOTH twist
halves hold by construction — the state plain forward-eval cannot represent. Activating the
slack ripples into ~18 verifier CHECK atoms; SA-with-square-roots (`slack_sa.py`, replaces the
deg-4 squares a40782/a39550 by their deg-2 roots Q=0 via `check_square.try_sqrt`) drives the
frustrated core 18 -> 6. Run the 4-way fleet: `python3 slack_sa.py <activator> <seed> <out.json>`
with activators in {1858,26947,27512,30104,5443,...}.
CRISP OBSTRUCTION: for verifier square a40782, satisfying it AND a1817 reduces to R=0 where
R = 28*x_10783 + (ripple terms), x_10783=x_16644*x_17301, all fixed by the RIGID 3183-slack
(a44271: x_3183=x_17728, so x_27116=x_17728-x_1642 is pinned). The continuous knobs x_24026 and
the FREE var x_31302 (df=None) CANNOT change R (Q40782 slope in x_24026 is 0 once a1817 held).
So the witness = a DISCRETE 233/22-bit choice whose rigid-slack ripple self-annihilates in every
verifier square (a knapsack). The div-wire escape (x_8821=x_17810*x_27292 in {-2,-1,0,1}) lets
x_18274/x_17728 leave their g2/h2 lattice but only onto (base/2)*Z, still coprime => degenerate.
NEXT: keep the slack-active SA fleet running; or attack R_i=0 across the 530 squares as a system
(linearize ripple monomials); or find the setter's 233-bit knapsack solution (LLL blocked by
numerator nonlinearity — 7/50 linear). NOTE: single 22-side bits give x_12779=1 (not 2).

## TRAPDOOR MECHANISM — fully reverse-engineered (Session 6, read NOTEBOOK Session 6 tail)
The obstruction (atoms 1817,30378,40782,44271) is the twist x_9770=x_18274 & x_3183=x_17728.
KEY: the confluent forward-eval QUANTIZES both sides to COPRIME units and ZEROS the slack
products, so it can NEVER represent the (feasible) witness — this is why every forward-eval
search (SA/mitm/greedy/pairs/enum/local) plateaued. Specifically:
- Under forward-eval: x_9770=m*g, x_3183=m'*h, x_18274=m2*g2, x_17728=m2'*h2 (g=119182..,
  g2=91416..; gcd(g,g2)=1, gcd(h,h2)=2). Rigid twist => degenerate 0 only. (codewords.py, quant_structure.py)
- BUT the wire DEFS carry product slacks: x_9770 = x_35186(=m*g) + x_3368, x_3368=x_12779*x_24026;
  x_3183 = x_1642(=m'*h) + x_10466, x_10466=x_12779*x_27116. Both gated by x_12779=x_23380*x_36336.
- forward-eval sets x_12779=0 (slacks off) -> quantization. The WITNESS activates x_12779 (22-side
  bit pairs give x_12779=2) AND x_24026/x_27116 (deeper, via x_38215) so
      x_9770 = m*g + x_12779*x_24026 = x_18274 = m2*g2   (bridges the coprime gap).
- So the TRUE solve = search WITH the slacks active. With slacks on, x_9770 is NOT limited to 27
  values and CAN equal x_18274; the decoupling (x_9770<-22 only) is a slacks-OFF artifact.

NEXT-STEP for a solver: build an evaluator/search that DRIVES x_12779, x_24026, x_27116 nonzero
(find their activating bit cascades: x_12779<-{1858,2795,5443,10652,19520,26947,27512,30104,...},
x_24026<-x_38215<-...), then solve m*g + x_12779*x_24026 = x_18274(B) (coupled product match).
Do NOT rely on the all-0 forward-eval regime — it structurally excludes the witness.

## Earlier (still true) reduction
- `A` = the 22 control bits `BITS22`; `B` = the other 233 bits.
- `x_18274 = x_6773/x_8821`, `x_17728 = x_17233/x_8821` (SHARED denominator x_8821).
- `x_8821` is **exactly linear** in the 233 bits; numerators are high-degree.
- best_partial_39019 sets ALL 255 control bits = 0.
- twist eqs: 1817 = 6033033*(x_9770-x_18274)+x_26977; 44271 = x_3183-x_17728;
  30378 = x_3183-x_9982-x_17728. (x_26977, x_9982 identically 0.)

## How to evaluate (the correct model)
`confluent_eval5.build5()` -> (A_atoms, kind, info, seq, bestval, ncyc). Build `seq`:
```python
order = json.load(open('eval_order.json'))['order']
defset = set(v for v in kind if kind[v] != 'const')
seq = [v for v in order if v in defset and v not in (9770,3183)]
seq += [v for v in (9770,3183) if v in defset]
seq += [v for v in defset if v not in set(order) and v not in (9770,3183)]
```
`make_forward(kind,info,seq,bestval)` -> Z solver `solve(list(bestval), setbits)`;
`make_forward(...,mod=P)` -> mod-P solver. forward_Z([]) violates exactly {1817,30378,40782,44271}.
The forward-eval satisfies every ORIENTED gate/load/div atom by construction for ANY bit set;
only the twist "check" atoms float — so it is a valid oracle for x_9770/x_3183/x_18274/x_17728.
NOTE: integer forward-eval is *lossy* (leaves a stale value when a division isn't exact) — use
the mod-P solver for any linearity/degree probing.

## Highest-EV next experiments
1. `runs/tab22_full.log` — full 2^22 (x_9770,x_3183) mod two 31-bit primes; saves
   tab22_9770_{p}.npy / tab22_3183_{p}.npy. When done: confirm B=0 fails; hash S and inspect
   structure (common factors, moduli). S then lets you INVERT the 22-side in O(1) (lookup).
2. Residue-pool identity: `extract_huge.py` -> huge_network.json (865 huge atoms; 512 simple
   loads bit*(x_B-HUGE)=s*x_C). Check whether x_9770(A) and x_18274(B) are combinations of the
   SAME HUGE residues => matching becomes combinatorial, not brute 2^233.
3. MITM/lattice via x_8821 (the linear coordinate on the 233 side) — see NOTEBOOK Session 6.

## Exhausted this session (do NOT redo)
- SAT/SMT (user directive: custom heuristics only; z3/cvc5 return unknown anyway).
- v4 evaluator / anything freezing x_18274 (fixed in v5).
- Linear algebra / lattice: `linalg255.py` (CORRECT, over all 255 bits, mod-P) has RANK 255/255
  and forces ALL bits = 0. The witness (!= all-0) is OUTSIDE the linear neighborhood of all-0, so
  linear/lattice/subset-sum provably cannot reach it. Supersedes the Session-5 "slaved-B" claim.
- B=0: ruled out (full 2^22 scan, 0 matches).
- Modulus (gcd residues=1), residue-lattice relation (none), slack vars (x_26977/x_9982 rigid).
- Local search / greedy / SA / pairs / triples from all-0 — all plateau (all-0 is the local min).

## The ONE remaining avenue (unimplemented)
A custom NONLINEAR solver / backward circuit-inversion: pin x_18274=N1, x_17728=N2 for a chosen
(N1,N2) from the 2^22 table S, and propagate/search backward through the 233-side acyclic circuit
(residue-load selects + product/sum gates) to determine the bits. Big build, uncertain (z3 failed
the analogous forward CSP). This is the only path not proven dead — everything else is exhausted.
The instance is a genuine obfuscated-circuit trapdoor; a full witness likely needs the setter's
secret or a cryptanalytic break of the specific 233-side residue circuit.

## Git
Branch `claude/read-prompt-5t2raw`. Commit+push after meaningful experiments.
