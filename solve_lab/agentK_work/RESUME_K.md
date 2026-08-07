# RESUME_K — agent K. Self-contained. Everything here was measured on this box, from the file.

No named-problem framing anywhere. Everything below is integer / modular arithmetic on
constants read out of `EQUATIONS.txt` and verified against the equations themselves.

--------------------------------------------------------------------------------------------------
## 0. SCORES

- Shared baseline **39,026 / 39,033**, re-verified by me at the start of this session:
  `python3 solve_lab/checker.py solve_lab/best/new_instance_partial_39026.json`
  -> `satisfied 39026/39033 (7 failing)`, failing `[12231,12270,12350,14584,18673,22044,29125]`.
- **I did not beat it.** Nothing in `agentK_work/` is a better partial.
- **No infeasibility is claimed.** The instance is *satisfiable*; see §5.

--------------------------------------------------------------------------------------------------
## 1. THE INSTANCE IS NOW FULLY DECODED — closed form, no undecoded slots left

Agent F's decode is completed and corrected here. The whole 39,033-equation system is
equivalent to a **single 256-bit statement**.

**Constants (all read out of the file, `k15_leaves.py` finds every literal > 10^20):**
there are exactly **516** large literals — 512 leaf-pin constants, 2 target constants,
`p`, and `K`. Nothing else.

```
p     = 115792089237316195423570985008687907853269984665640564039457584007908834671663   (x26064)
K     = 97553848499418123410591666447050222001188385549510401465815187079080512838891    (x24453)
shift = K/3 mod p
b     = 64019533680030876408443198762210829058751700634554282185987325820393598524794
N     = 115792089237316195423570985008687907852837564279074904382605163141518161494337
```

**Composition law.** Each stage's three checks are exactly
`out_x = l^2 - a_x - b_x - K`, `out_y = l(a_x - out_x) - a_y`, `l = (b_y-a_y)/(b_x-a_x)`.
Substituting `X = x + shift` removes `K` completely, giving `out_X = l^2 - a_X - b_X`.
**Under that substitution all 256 leaf constant pairs and the target satisfy the same cubic
`Y^2 = X^3 + b` (a = 0), 256/256, verified** (`k16_points.py`). The law is therefore the
classical chord composition on that cubic: commutative and associative, so **the fold of a
leaf subset does not depend on the tree shape at all** — the 96-stage tree is irrelevant to
the value, only to the wiring.

**The 256 leaves are one doubling chain.** 255 of the 256 have their double also in the set;
one leaf (index 28) is nobody's double, one (index 243) has no double in the set. Following
the chain labels every leaf with an exponent `e ∈ 0..255` and leaf(e) = 2^e · leaf(0).
(`k18_struct.py`, `k19_chain.py`, `chain.json`.)

**Group order.** `4p = L^2 + 27M^2` via Cornacchia gives six candidate traces; exactly one
`N = p + 1 - t` annihilates the chain base, checked by composing the chain itself
(`k21_order.py`). N is prime (Miller-Rabin, 40 rounds).

**Therefore, writing `G` for leaf exponent 0 and `T` for the target point:**

> the system is satisfiable exactly by the leaf ON-sets `S` with
> `sum_{e in S} 2^e ≡ dlog_G(T) (mod N)`,
> and since `0 < dlog < N < 2^256` and every exponent 0..255 is available exactly once,
> **the solution is unique: S = the binary support of `dlog_G(T)`.**

--------------------------------------------------------------------------------------------------
## 2. THE FOLD EVALUATOR — built and VALIDATED (F's priority 2, done)

`fold.py` — extracts the 256 points and the target from the equations (`build_points()`),
composes them (`add`, `fold`). `points.json` is its output.

Validation is against the real equations, not against a model:
`cascadep.py` closes all 39,033 atoms **mod p** from a chosen leaf ON-set in ~1–3 s
(`k26_drive.py: drive(on)`), and `rootpair(v)` reads the root's two input pairs
`(x12186,x16742)` and `(x14853,x24908)` back out of the circuit.

**CRITICAL — forward-only closure.** A plain cascade runs the *target pin backwards* down the
tree (it derives `x24468` from the target literal, then `x13682 -> x38045 -> x22162 -> the root
check -> the A-side input`), which silently produces "A = target X" and makes every fold check
look wrong. `drive()` therefore marks the two target-pin atoms as underivable
(`k26_drive.FORBID`). **Any future evaluator must do the same** — this cost me two wrong
readings before I caught it (A.x came back equal to the target's X, which is what exposed it).

| ON-set (by exponent) | prediction | measured from the equations |
|---|---|---|
| `{2081, 24601}` (the deliverable's) | A = leaf(24601), B = leaf(2081) | **both exact** |
| `{e163, 2081}` | A = leaf(e163), B = leaf(2081) | **both exact** |
| `{0,1,3}` | A = 3G, B = 8G | **both exact** |
| `{0,1,2,3,5}` | A = 7G | **A exact**; B not reproduced |
| `{4,6,7,10,12}` | A = 208G | **A exact**; B not reproduced |

So the composition law is confirmed as the group sum on the **A half for 1-, 2- and 3-leaf
folds**, and on the **B half for single leaves**. Multi-leaf B-half folds are *not yet
reproduced by the closure* — see `bside.log`. That is a closure/wiring gap on the 78-side, not
evidence against the law (the A half uses the identical stage law and matches exactly); but it
is **open**, and §1's uniqueness statement inherits that caveat until it is closed.

The validation the handoff asked for **passes in the required direction**: the evaluator
predicts `fold({24601,2081})`, and that value is **NOT** the target
(`fold.py` prints `equals target? False`).

Root leaf split, measured: **A side 178 exponents, B side 78** (`rootsupport.json` +
exponent 163 settled empirically in K26, not guessed). Independently equals F's 178|78.

--------------------------------------------------------------------------------------------------
## 3. WHAT THE 39,026 DELIVERABLE ACTUALLY IS (previously unexplained)

Its 7 failing equations come from exactly 7 nonzero atoms, all in one cluster:

```
def  x29854 - p*x1329        res  5113045*x7075*x9118 - x29854
def  x31864 - p*x10903       res  x7075*x8731 + x31864
def  x642   - p*x17325       res  (x7068 - x2099) - 7376877*x642
def  x28730 - p*x9413
```
`x7075 = 1 - x21279`, `x21279 = x4287*x2081` is a gate's liveness bit and `(x8731,x9118)`
is that gate's output pair. So the four conditions are the gate's **off-pins**:
`x7075*x9118 ≡ 0`, `x7075*x8731 ≡ 0 (mod p)` plus two slot congruences.

**Why it scores 39,026:** at the root the deliverable has `inA == inB` **exactly**
(both equal `leaf(24601)`), which makes `b_x-a_x = 0` and `b_y-a_y = 0`, so the root check
`x35389 = (out_x+a_x+b_x+K)(b_x-a_x)^2 - (b_y-a_y)^2` vanishes *identically* and the root
**output becomes completely free**. It is then set to the target. Confirmed:
`2*leaf(24601) != T` and `leaf(24601) != T`, yet every root atom is zero.
The deliverable pays for forcing the B-side chain to carry `leaf(24601)` instead of the
honest `leaf(2081)`; that lie surfaces only at the gate-21279 off-pins, cost 7 equations.

**So a degenerate stage has a free output. That is a real, exploitable hole — and it is closed.**

--------------------------------------------------------------------------------------------------
## 4. THE DEGENERACY ROUTE IS CLOSED (a complete, non-heuristic negative)

A stage is degenerate iff its two children fold to the *same* point (equal x AND equal y;
`a_x=b_x, a_y=-b_y` gives `x35389 = -(2a_y)^2 != 0`, so that case does not work).
For a stage whose children own disjoint exponent sets `J1, J2`:

```
x = sum_{S1} 2^e  (bits in J1),   y = sum_{S2} 2^f  (bits in J2),   need x ≡ y (mod N)
```
`2^255 < N < 2^256`.
* **interior stage** (`|J1|+|J2| = n < 256`): `|x-y| < 2^n <= 2^255 < N`, so `x = y`, and
  disjoint bit sets force `x = y = 0`. **Impossible.**
* **root** (`J1 = IA` 178 bits, `J2 = IB` 78 bits, together all of 0..255): `|x-y| < 2^256 < 2N`,
  so `x - y = ±N`. Every bit position belongs to exactly one side, so the schoolbook addition
  `x = y + N` has **no free choice at any bit** — it is a deterministic 256-step carry walk.
  `k22_dp.py` runs it in both directions: **both terminate with a nonzero final carry.**
  **Impossible.**

Knob set for this claim: **all 256 leaf selector bits, over all 2^256 configurations, with
every other variable free.** Configuration: any. It is a statement about the exponent
arithmetic, not about a particular assignment.

**Robustness (this matters, because my 178/78 split is measured, not proved):** re-running the
DP with any single exponent moved to the other half — all 256 variants — still fails, and over
**2,000 random 178/78 partitions the DP succeeds 0 times**. The obstruction is a property of
`N`'s bit pattern (bit 255 of `N` is 1 and the carry can only be cleared on an A-bit where
`N_i = 0`), not of the particular split. So the conclusion survives even if the split is
slightly wrong.

Consequence: 39,026's trick cannot be completed honestly, and there is no cheaper stage at
which to play it.

--------------------------------------------------------------------------------------------------
## 5. WHY 7 IS HARD TO BEAT (measurement, not a proof)

`k27_sites.py`: equations-per-atom over all 39,033 atoms.
Minimum footprints: `1,2,3,4,5,5,6,6,6,6,6,6` — **all twelve are idempotency atoms
`(x*(x-1))` of decoy booleans**, which are equation-disjoint from the fold path (J's finding,
re-measured). The lightest *load-bearing* atom sits in **7** equations, and the deliverable's
cluster is exactly one of those 7-equation atoms with its partners cancelling inside the same
7 rows. Costs I measured for alternative injection sites:

| site | equations broken |
|---|---|
| deliverable's gate-21279 off-pins | **7** |
| cheapest leaf-pin pair (`sel x5090`, `x22106`) | 10 |
| breaking both target pins (`x24468`, `x18956`) | 16 (10 shared, ≤2 cancellable ⇒ ≥14) |

So beating 7 requires *cancellation below a single atom's footprint*, not a cheaper site.
I did not find one. **This is a measurement over the sites I enumerated (leaf pins, slot pins,
target pins, all single residual atoms); it is not a proof that no coding does better.**

--------------------------------------------------------------------------------------------------
## 6. TOOLS BUILT (all run from cold in this directory)

| file | what it does | runtime |
|---|---|---|
| `cascade.py` | exact-**integer** cascade closure over all 39,033 atoms, no DAG orientation assumed; 3-point linearity test, solves only when the coefficient divides exactly | 0.5 s |
| `cascade2.py` | incremental version: propagate first, guess only when nothing is forced | 0.5 s |
| `cascadep.py` | the same closure **mod p** — this is the fast driver | 1–3 s |
| `k25_class.py` | role-based variable classification -> `varclass2.json` (369 bools / 3747 handles / 4631 wires; **all 256 leaf selectors are free**) | 30 s |

**Known bug in `k25_class.py`, fix before trusting `varclass2.json`'s bool column:** it detects
idempotency atoms only in the form `(xA*(xA-1))`. They also occur as `(xA*(1-xA))` and
`((xA*xA)-xA)` — e.g. `x4287` is a boolean but is classified as a "wire". This does not affect
the leaf-selector list (that comes from the leaf pins, independently) and so does not affect
any result above, but it is the first thing to fix when extending the driver.

| `fold.py` | leaf points + target extraction, group composition | 5 s |
| `k26_drive.py` | `drive(on_set)` -> full mod-p state; `rootpair(v)` -> the root's two input pairs | 3 s/run |
| `k19_chain.py` `k21_order.py` `k22_dp.py` | chain labelling, N, the degeneracy DP | fast |

Data: `points.json` (256 points + target), `chain.json` (exponent labels + selector map),
`rootsupport.json` (178/78 split), `order.json` (N), `varclass2.json`, `bigliterals.json`.

Key facts to re-derive cheaply if anything looks wrong:
- **handles absorb everything over Z**: `k9_handles.py` seeds every handle to 0 and the
  integer cascade closes with **0 conflicts**, leaving only the root-check atoms. So the
  binding content of the instance is mod p and nothing else.
- The 4 "big" free variables per stage are output/slot wires, not knobs.

--------------------------------------------------------------------------------------------------
## 7. NEXT, IN ORDER (for whoever picks this up)

1. **The only thing between here and a full solve is one integer**: `k` with `k·G = T` in the
   chord group on `Y^2 = X^3 + b (mod p)`. Given `k`, the solve is mechanical: set leaf
   selector `e` on iff bit `e` of `k` is 1, run `k26_drive.drive()`, then the integer cascade
   (`cascade2.py`) with handles free to lift to Z, then `checker.py`.
   Bounded searches already run and **negative**: `k` is not < 2^40 (BSGS), `T` is not a leaf,
   not a sum of 2 or 3 leaves, `2P_i != T` for all i, no pair-sum collisions among the leaves.
2. `k28_shots.py` extends the BSGS window; N is prime so there is no exponent split.
3. If instead the goal is to beat 39,026: the only opening left is **cancellation** — find an
   atom set whose image under the incidence matrix has weight ≤ 6 *and* is realizable. Use
   `k27_sites.py`'s footprint table as the starting point; the 12 sub-7 atoms are decoys, so
   any solution must cancel 7 rows of a load-bearing atom against something.
4. Do **not** redo: the tree decode, the slot wiring, the leaf/target extraction, the
   178/78 split, the stage-degeneracy question. All settled above.

--------------------------------------------------------------------------------------------------
## 8. VERIFICATION RULE (fleet standing rule, respected here)

States above 4,300 decimal digits cannot be parsed by `checker.py`; use
`solve_lab/agentE_work/verifyE.py`. Every score quoted in this file came from
`checker.py` on the 39,026 file, whose largest value is 909 digits, so it is in range.
