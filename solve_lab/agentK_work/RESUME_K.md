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

> writing `k` for the integer in `[1,N)` with `k·G = T` under the composition law above,
> the system is satisfiable exactly by the leaf ON-sets `S` with
> `sum_{e in S} 2^e ≡ k (mod N)`,
> and since `0 < k < N < 2^256` and every exponent 0..255 is available exactly once,
> **the solution is unique: S = the binary support of `k`.**

**Status of that statement.** The constants, the cubic, the chain, `N` and the 178/78 split are
all *measured*. The step "fold = composition sum" is validated end-to-end on the A half up to
3-leaf folds (§2) and rests on F's independent 72/72 stage-law check; the B half is validated
for single leaves only. Treat the uniqueness sentence as **strongly supported, not proved**,
until §2's B-half gap is closed.

**What the reduction does and does not do — read this before quoting it.** It converts 39,033
equations into *one* unknown integer `k`. That is a real simplification of the search space
(from 38,748 coupled variables to 256 bits with a closed-form consistency test), and it is what
makes §4's negative provable. It **does not finish the instance**: every step after `k` is known
is bookkeeping, but **obtaining `k` is the entire remaining problem**, and this reduction
relocates that problem rather than solving it. Concretely, the recovery is a 256-bit one, and
what I actually ruled out is small: `k < 2^40` (baby-step/giant-step), `k` of Hamming weight
≤ 3, `T` equal to a leaf or to `2P_i`, and any pair-sum coincidence among the 256 leaves. That
leaves essentially the whole 256-bit range untouched. I know of no shortcut and I did not find
one. **No claim is made here, in either direction, about whether `k` is obtainable** — the
instance is satisfiable, a satisfying assignment exists, and nothing below argues otherwise.

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

So the composition law is confirmed as the composition sum on the **A half (178 leaves, depth 6)
for 1-, 2- and 3-leaf folds**, and on the **B half for single leaves**. Multi-leaf B-half folds
are *not reproduced by the closure* and this is **still open**. What I established about it this
session, so nobody repeats it:

* **The "backward flow" hypothesis is REFUTED.** `k29_trace.py` records, per variable, which
  atom derived it. The back-cone of `x14853` is 5,119 variables, contains **no variable from
  above the root**, and `x14853` is derived by its own forward slot pin. The B half *is* being
  computed forward from its 78 selectors. The only seeded things in the cone are the 78
  selectors and 662 handle factors (all of shape `xY - (xA*xB)`), which is correct behaviour.
* **The value is not a near-miss.** `k31_whatisB.py`: for `ON={e3,e5}` the B value is not the
  composition, not the difference either way, not the negation, not a coordinate-wise sum, not
  `a·P1+b·P2` for any `|a|,|b| <= 4`, not any leaf, and **not on the cubic**. For `ON={e3,e10}`
  it *is* on the cubic but is still none of those. So it is not a sign or shift bug.
* **`k34_diverge.py` is not a usable diagnostic as written** — it compares every wire against a
  predicted point, but most wires are chord-arithmetic intermediates (the shallowest "failing"
  one, `x7340`, is defined as `x32101 + x12781`, a partial sum), so mismatches there mean
  nothing. It also uses the *inflated* supports. Rewrite it to compare only slot wires, using
  `k36_tight.py`'s corrected supports.
* **The liveness-bit hypothesis is REFUTED too.** `k35_otherbools.py` drives the same ON-sets
  with the 900 non-leaf free booleans seeded 0, seeded 1, and left for the closure to derive.
  **All three give byte-identical verdicts** (`{0,1,3}`: A and B both match in all three;
  `{0,3,10}` and `{3,10}`: B fails in all three). So the 900 are not the cause.

So three hypotheses are dead — backward flow, sign/shift, and liveness-bit seeding — and the gap
is still there. Whatever it is, it is specific to the 78-side wiring at ≥2 live leaves.

None of this is evidence against the law — the A half runs the identical stage law and matches
exactly, and F verified the law on 72/72 stages independently — but §1's uniqueness statement
keeps its caveat until it closes. Note that §4's negative does **not** depend on it: §4d needs
only `N > 2^256 - 2^129` and the root split, neither of which involves the B-half closure.

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
## 4. THE DEGENERACY ROUTE IS CLOSED — as a claim about the PARTITION, not about N
### (sub-sections run 4.0, 4a, 4b, 4d, 4c — 4d is the argument that stands; 4b/4c are the repair history)

A stage is degenerate iff its two children carry the *same* coordinate pair (equal x AND equal
y; `a_x=b_x, a_y=-b_y` gives `x35389 = -(2a_y)^2 != 0`, so that case does not work).
For a stage whose children own disjoint exponent sets `J1, J2`:

```
x = sum_{S1} 2^e  (bits in J1),   y = sum_{S2} 2^f  (bits in J2),   need x ≡ y (mod N)
```

### 4.0 THE EXACT FORM OF THE CLAIM — read this before the rest of §4

The negative is **a statement about the partition, not about `N`**. The arithmetic-only version
is **false**, and I verified that rather than assuming it:

> `N` has 192 one-bits and 64 zero-bits. Take `j` with `bit_j(N)=1, bit_{j+1}(N)=0` (`j = 0`
> works); rewriting `2^j = 2^{j+1} - 2^j` gives non-empty disjoint `A, B ⊆ {0..255}` with
> `Σ_A 2^e - Σ_B 2^e = N`. Checked numerically: holds exactly.

So *some* pair of disjoint exponent sets satisfies the equal-inputs condition. What is being
claimed, and all that is being claimed, is:

> **THEOREM (partition form).** For every stage of this instance, neither of its two slot
> supports contains all of `{129..255}`. Since `2^256 - N < 2^129`, a set omitting any exponent
> `>= 129` has subset sums strictly below `N`. Hence `|x - y| < N` at every stage, so `x = y`,
> so (disjoint supports) both sides are empty — which the liveness gate forbids.

The load-bearing facts are therefore **measured partition facts**, not properties of `N`: the
root's two halves each omit many exponents `>= 129` (43 and 84 witnesses), and every interior
stage's support is a subset of one root half. Attack those, not the arithmetic.

My earlier phrasing did not make this explicit and could be read as the stronger, false claim
that no disjoint pair represents `N`. That reading is wrong and I am correcting it here.

### 4a. Which modulus bounds the walk, and how many wraps must be excluded

**The governing modulus is `N`**, the order of the chain base: measured, not assumed —
`k21_order.py` composes the chain with itself `N` times and lands on the law's identity, and
`N` is prime, so the order is `1` or `N`, and the base is not the identity.

**`N` does NOT exceed 2^256.** It is just below it:
`2^256 - N = 432420386565659656852420866390673177327`. So if the argument needed
"modulus > largest signed subset difference", it would be **broken**, and P is right to press.

It needs something weaker, which does hold, by a factor of about two. **Stated without any
premise about how the instance was built** — it uses only the measured facts that each input's
scalar is a subset sum of `{2^e : e = 0..255}` with each exponent available at most once
(that is the doubling identities of §4c, checked 255/255) :

```
max possible |x - y|  =  2^256 - 1
2N                    =  231584178474632390847141970017375815705675128558149808765210326283036322988674
2N - (2^256 - 1)      =  115792089237316195423570985008687907852405143892509244725752742275123193348739   > 0
```

so the only multiples of `N` reachable are `k = 0, +1, -1`, **unconditionally**. `k = 0` needs
`x = y` on disjoint bit supports, i.e. `x = y = 0`, excluded because both halves must be live.
**`k = ±1` therefore exhausts the cases** and a two-direction walk is complete. The right
condition to quote is `2N > 2^256 - 1`, i.e. `N > 2^255`, **not** `N > 2^256`.

*(Recomputed independently after the coordinator supplied the same bound. Same conclusion; the
slack is `≈ 1.158 × 10^77`, one decimal order larger than the `10^76` quoted to me. The
conclusion is unaffected — the slack is about half of `2N` either way.)*

### 4b. My interior-stage argument was WRONG. Corrected here.

I had written: "interior stage with `n < 256` leaves has `|x-y| < 2^n <= 2^255 < N`". **That is
false.** A stage's exponent set is an arbitrary subset of `{0..255}`, not an initial segment —
a stage owning exponent 255 has `x >= 2^255` no matter how few leaves it has. Only the root case
was argued correctly. P's challenge is what surfaced this; it is exactly the "case the author
had not enumerated" pattern. The repair:

`|x-y| = N > 2^255` forces one side's set `J` to have `sum_{e in J} 2^e >= N`. But
`sum_{e<=254} 2^e = 2^255 - 1 < N`. **So a degenerate side must own exponent 255** — only stages
containing that one leaf are candidates, and the walk gains a third per-bit case: positions
outside `J1 u J2` are available to neither side, so the walk can now die mid-way, not only at
the final carry.

* `k32_allstages.py` — the 12 nested supports owning exponent 255, both directions each,
  the unassigned leaf tried on either side: **all 12 fail.** But that pairing comes from my
  tree recovery, and K32's "chain" contains two *distinct sets of the same size*, so it is not
  provably a nest and the pairing is not provably the real stage structure.
* `k33_allpairs.py` — drops the tree: tests *every disjoint pair* of recovered supports (every
  real stage pair is among them, so extra pairs only make the test stricter), with an exact
  prune (`x - y = N, y >= 0` forces the "+" side's mask value `>= N`; 161 of 6007 supports
  qualify). **RESULT: the test FAILS — 262 pairs tested, many admit `x - y = N`.**

### 4d. CLOSED — and the carry walk turned out to be unnecessary

`k36_tight.py` fixed the support recovery: identify liveness/boolean variables by fixpoint, and
in a gated term `(xA*xB)` take only the operand that is *not* one of them. Sanity gate passes —
root comes back **A = 177, B = 78, disjoint** (the 178th is exponent 163, measured separately in
K26). With supports no longer inflated, **0 of 5480 recovered supports can even serve as the
"+" side**, and the reason exposes a far simpler and more robust closure than the carry walk:

```
2^256 - N = 432420386565659656852420866394968145599  <  2^129 = 680564733841876926926749214863536422912
```

So `sum_{e in J} 2^e >= N` requires `J` to omit **nothing** of value `>= 2^129` — i.e. `J` must
contain **every** exponent in `129..255`. Then:

1. `|x - y| <= max(maskval(J1), maskval(J2))` for any stage's two sides.
2. A nonzero multiple of `N` needs `max(maskval) >= N`, hence needs one side `⊇ {129..255}`.
3. The root's two halves are disjoint, and **each omits many exponents `>= 129`** — measured:
   the B half owns 43 of them, the A half owns 84. So neither half contains `{129..255}`.
4. Every interior stage's side is a subset of one root half, so it omits them too.
5. Therefore `|x - y| < N` at **every** stage, so `x = y`; disjoint bit supports force
   `x = y = 0`; both halves would have to be dead. **Impossible.**

`maskval(IA) >= N : False`, `maskval(IB) >= N : False`, verified directly.

**Why this version is worth more than the walk.** It needs only two measured facts, both with
enormous margin: `N > 2^256 - 2^129`, and *each root half omits at least one exponent >= 129*
(43 and 84 witnesses respectively, where one would do). It does not depend on the tree shape, on
my 178/78 split being exactly right, on which side exponent 163 sits, or on any carry
bookkeeping. The earlier deterministic walk (`k22`, `k32`) and the 0-of-2000-random-partitions
measurement agree with it and are kept as independent confirmations.

**Audit trail, kept deliberately.** This negative was wrong three times before it was right:
the interior case (`|x-y| < 2^n`, false — §4b); the tree-free test in `k33`, which *failed to
close* because my supports were over-approximations (§4d); and the **phrasing**, which read as
the stronger claim that no disjoint pair of exponent sets represents `N` — that claim is false
and §4.0 exhibits the witness. All three were found by taking outside challenges seriously
rather than defending the claim.

**If it is challenged again, attack step 3** — the measured claim that each root half omits an
exponent `>= 129`. That single fact, plus "every interior stage sits inside one root half", is
the entire content. Everything else is arithmetic that holds unconditionally.

### 4c. Restated without any premise about how the instance was built

P notes, correctly, that "256 leaves at distinct exponents in an order-`N` arithmetic" would be
a construction premise. It is not assumed here — every ingredient is a **verified identity among
the file's own constants**:

| identity | check | result |
|---|---|---|
| the law is commutative | 3000 random pairs of file constants | 0 failures |
| the law is associative | 3000 random triples | 0 failures |
| `leaf(e+1) == law(leaf(e), leaf(e))` | all 255 steps | 0 failures |
| `N`-fold composition of the base == identity | direct | holds |
| every leaf and the target satisfy `Y^2 = X^3 + b` | 256 + 1 | all satisfy |

Associativity + commutativity is what licenses "fold(S) depends only on the *set* S", and the
doubling identities are what turn a subset into the integer `sum 2^e`. Nothing about the
instance's provenance is used.

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

**The adjacent hole is closed too.** If a half folded to the group identity the chord law would
have nothing to represent and its output might likewise go free. That needs
`sum_{S} 2^e ≡ 0 (mod N)` for `S` inside one half, i.e. `sum = N` exactly, i.e.
`supp(N) ⊆ IA` or `⊆ IB`. Measured: `popcount(N) = 192`, and 55 of those bits fall outside
`IA`; `IB` has only 78 bits. **Impossible**, and impossible at every interior stage for the same
size reason.

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

**Where this table is weakest, stated so it can be attacked.** Two gaps, either of which would
revise it rather than confirm it:
* I classified the twelve sub-7-footprint atoms as decoys **by inspection** — they are all
  idempotency atoms of booleans whose only other appearance is a single trivial definition
  (`x20650 - x157*10101`, etc.), which is why I judged them unable to carry a defect that
  changes the fold. I did **not** prove they are inert. If someone exhibits a *realizable*
  defect carried by a low-footprint atom, the "load-bearing atoms all cost ≥ 7" line is wrong
  and this table is the thing that gets corrected, not the counterexample.
* My footprint numbers are **unions, i.e. upper bounds with cancellation ignored**, except for
  the deliverable's cluster where I checked cancellation explicitly. A site I priced at 10 or 16
  could be cheaper if its atoms cancel. I did not search for that.
Anyone screening atoms by footprint (e.g. minimising equations-touched against window size)
will land on exactly those twelve decoys first; that agreement is not corroboration, it is the
same measurement. The open question is realizability, which footprint alone cannot see.

**Independent corroboration, and its limit.** An independent decomposition of the same file
(different directory, different parse) measures **1,158 atoms with footprint < 7, of which
1,152 are idempotency atoms, 1,145 of them at footprint 1**. My twelve and that 1,152 are the
same phenomenon seen through different atom decompositions — the counts differ because the
decompositions differ, not because the instances differ. So "the sub-7 footprints are
idempotency atoms" holds at ~100x the scale I measured it.

**But decoy-ness is not the only way a light carrier dies, and my table only knows about
decoy-ness.** That same independent work found a window where *every* atom is genuine law-block
arithmetic — no decoys at all — and the light carrier there still failed, for a different
reason: **downstream coupling** (the carrier's value is forced by what consumes it). My table
cannot see that failure mode at all; it prices sites by equation footprint only. Treat the two
explanations as complementary, neither subsuming the other: a candidate site must survive
*both* "is it load-bearing?" and "is its value free once its consumers are accounted for?".
If a realizable low-footprint carrier is ever exhibited, **this table is what gets revised**,
not the counterexample.

--------------------------------------------------------------------------------------------------
## 6. TOOLS BUILT (all run from cold in this directory)

| file | what it does | runtime |
|---|---|---|
| `cascade.py` | exact-**integer** cascade closure over all 39,033 atoms, no DAG orientation assumed; 3-point linearity test, solves only when the coefficient divides exactly | 0.5 s |
| `cascade2.py` | incremental version: propagate first, guess only when nothing is forced | 0.5 s |
| `cascadep.py` | the same closure **mod p** — this is the fast driver | 1–3 s |
| `k25_class.py` | role-based variable classification -> `varclass2.json` (369 bools / 3747 handles / 4631 wires; **all 256 leaf selectors are free**) | 30 s |

**Classifier bug — FOUND AND FIXED this session.** `k25_class.py` originally matched
idempotency atoms only as `(xA*(xA-1))`. The file also writes them `(xA*(1-xA))` and
`((xA*xA)-xA)`. All three spellings are now matched, and the counts move a long way:

| | before fix | after fix |
|---|---|---|
| free booleans | 369 | **1156** |
| leaf selectors recognised as boolean | 82 / 256 | **256 / 256** |
| free "wires" | 4631 | 3844 |

**This did not invalidate any result above** — `drive()` takes its 256 selectors from the leaf
*pins* (`points.json`), independently of the classifier, and seeds them ahead of everything
else, so §2's validations were driven correctly either way. But it is a clean example of the
failure mode worth naming: a knob set filtered by a regex that misses a spelling will be
reported as a property of the instance. Re-run `k25_class.py` before using `varclass2.json`.

**Consequence, now SETTLED.** The fix exposed **900 free booleans that are not leaf selectors**,
which would have meant §1's knob set was too small. `k35_otherbools.py` settles it behaviourally
rather than structurally: driving the same ON-sets with those 900 seeded **0**, seeded **1**, and
**left for the closure to derive** gives byte-identical results in all three modes, on every
ON-set tried. They do not move the fold. So **§1's knob set stands: the 256 leaf selectors.**
(`k30_decoys.py` is the structural version of the same test — forward cone under the most
generous, *undirected* notion of influence — if you want a second line of evidence.)

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

0. **Close the B-half validation gap** (§2) — it is the only thing still gating §1. The
   instrumentation exists (`CascadeP.close` now records `trace`/`deps`), and the two obvious
   hypotheses are already dead: it is **not** backward flow (K29) and **not** a sign/shift bug
   (K31). Start from `k35.log`, then rewrite `k34_diverge.py` to compare **slot wires only**,
   using `k36_tight.py`'s corrected supports rather than the inflated ones.
   (`k25_class.py`'s regex bug is fixed; §6 has the before/after.)
1. **What remains is one integer `k` with `k·G = T` under the chord composition on
   `Y^2 = X^3 + b (mod p)` — and finding it is the whole job, not a formality.** Do not read
   "one integer" as "nearly done": *given* `k` the rest is bookkeeping (set leaf selector `e` on
   iff bit `e` of `k` is 1, run `k26_drive.drive()`, lift to Z with `cascade2.py` leaving the
   handles free, check with `checker.py`), but obtaining `k` is a 256-bit recovery and every
   cheap avenue I tried came back empty: `k < 2^40` (baby-step/giant-step, `k28_shots.py`),
   Hamming weight ≤ 3, `T` a leaf, `T = 2P_i`, pair-sum coincidences among the 256 leaves —
   all negative. Those rule out a negligible fraction of the range. Budget accordingly, and
   do not spend a session assuming the reduction has almost closed the instance.
2. `k28_shots.py` extends the search window; `N` is prime, so the exponent does not split into
   smaller independent pieces, and that route is closed.
3. If instead the goal is to beat 39,026: the only opening left is **cancellation** — find an
   atom set whose image under the incidence matrix has weight ≤ 6 *and* is realizable. Use
   `k27_sites.py`'s footprint table as the starting point; the 12 sub-7 atoms are decoys, so
   any solution must cancel 7 rows of a load-bearing atom against something.
4. Do **not** redo: the tree decode, the leaf/target extraction, the 178/78 split, the constant
   inventory, or the stage-degeneracy question (§4d — closed, and closed by a size bound with
   43-and-84-witness margin, not by a fragile carry walk).
   Do **not** trust: `k34_diverge.py` as written, `k33_allpairs.py`'s verdict (it uses the
   inflated supports and its hits are artifacts — `k36_tight.py` supersedes it), or any support
   set that does not pass K36's root-must-be-178/78 sanity gate.

--------------------------------------------------------------------------------------------------
## 9. THINGS I GOT WRONG THIS SESSION, KEPT ON PURPOSE

Four, all caught here rather than downstream. The pattern is the same each time — a claim that
looked finished until someone asked which case it had not enumerated.

1. **Two fold "validations" that were reading the target back at me.** The cascade derived the
   root inputs *from* the target pin. Symptom: `A.x` came back exactly equal to the target's X.
   Fix: `k26_drive.FORBID`.
2. **The interior-stage degeneracy argument** (`|x-y| < 2^n`). False — exponent sets are not
   initial segments. Surfaced by a challenge from an independent parse (§4b).
3. **The tree-free replacement for it** (`k33`) *failed to close*, and the reason was my own
   support recovery inflating sets by unioning both operands of gated terms (§4d).
4. **The classifier** missed two of three idempotency spellings, undercounting free booleans
   369 vs 1156 and mis-typing 174 leaf selectors (§6).

Only (4) was harmless. (1) would have produced a fake model, (2) and (3) a fake barrier.

--------------------------------------------------------------------------------------------------
## 8. VERIFICATION RULE (fleet standing rule, respected here)

States above 4,300 decimal digits cannot be parsed by `checker.py`; use
`solve_lab/agentE_work/verifyE.py`. Every score quoted in this file came from
`checker.py` on the 39,026 file, whose largest value is 909 digits, so it is in range.
