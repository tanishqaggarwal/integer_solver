# The constraint catalogue on `k` — everything this campaign has measured about the unknown scalar

**Agent X.** Every row below was re-derived or re-executed here unless the *Verified by X* column
says otherwise; rows I could only audit rather than re-run are marked as such, and rows I could not
reproduce would be marked **UNVERIFIED** (there are none).

## 0. What `k` is

Four independent parses agree that `EQUATIONS.txt` is satisfiable **iff** some subset `S` of the
256 leaf selectors folds to a fixed target: leaf `i` carries the scalar `2^i` on a curve of prime
order `N`, and with `k = Σ_{i∈S} 2^i` the condition is `k·G = T`. **`|S|` is exactly the Hamming
weight of `k`**, and `k ↦ S` is the binary expansion, so a constraint on `k` is a constraint on the
ON-set and conversely.

**Re-verified here from agent Q's raw instance-derived files (`xsetup.py`), not taken on trust:**

| fact | measurement |
|---|---|
| `p = 2^256 − 2^32 − 977` | **True** — the secp256k1 prime |
| `N = 0xFFFFFFFF…BAAEDCE6AF48A03BBFD25E8CD0364141` | **True** — the secp256k1 group order |
| 256 leaf points on `y² = x³ + b` (`a = 0`) | **256 / 256** |
| `L_i = 2^i·G` by independent repeated doubling from `G = leaf(x2779)` | **256 / 256, 0 bad** |
| `N·G = O`, `N·T = O`, `T` on the curve | **True** |
| `G`, `T` identical to the values Q's six searches used | **True** |

So the group is **isomorphic to secp256k1** (`b = u⁶·7`). That is an independent confirmation of
"prime order, non-anomalous, no small embedding degree, no exploitable structure" — the reduction
lands on a curve whose hardness is the most-studied in the world.

---

## 1. THE CATALOGUE

| # | constraint on `k` | who | model | exhaustive? | verified by X | what would falsify it |
|---|---|---|---|---|---|---|
| 1 | **weight(`k`) ≥ 9** | X | group, mod `p` | **exhaustive** for weight ≤ 8: table = all 177,589,056 subsets of size ≤ 4, scan = all 177,588,800 of size ≤ 4, plus \|S\| ∈ {0,1} direct | **re-run from cold, mine** | a wrong ladder (`L_i ≠ 2^i·G`), a wrong `T`, or a truncation collision missed — excluded: 64-bit keys, expected false-positive count 0.085 over the whole sweep, and every planted target was found |
| 1b | **weight(`k`) ≥ 10** | X | group, mod `p` | scan of all 8,809,549,056 size-5 subsets against the same size-≤4 table — **in progress** | in progress | as above |
| 2 | `k > 2^44` and `N − k > 2^44` | Q | group, mod `p` | exhaustive, BSGS 2²² baby × 2²² giant, both signs | **code-audited**; superseded by row 3 | a bug in Q's `jadd`/`batch_affine_x` |
| 3 | **`k > 2^52` and `N − k > 2^52`** | **X** | group, mod `p` | exhaustive, BSGS 2²⁶ baby × 2²⁶ giant, both signs, C engine | **re-run from cold, mine**; planted `k₀ = 5·2²⁶+1234567` recovered at exactly `i = 5` | as row 1 |
| 4 | ON-bits are **not confined to any 34-bit window** | Q | group, mod `p` | exhaustive over `k = a·2^s`, `a < 2³⁴`, `s ≤ 222`; 2,865 s | **code-audited only, not re-executed** — the enumeration and the `sa == sb` join are correct on inspection | a bug in the shift bookkeeping; re-running costs ~1 min with X's C engine |
| 5 | `k ≠ 2^i / m mod N` for `m ≤ 10⁷` | Q | group, mod `p` | exhaustive in `m`, all 256 `i` | **code-audited only, not re-executed** | as row 4 |
| 6 | `T ∉ { ±λ^j·2^i·G }` (the 1,536-point endomorphism orbit) | Q | group, mod `p` | exhaustive | **re-verified independently** (`xverify_lam.py`): `φ(X,Y) = (βX, Y)` **is** multiplication by `λ`, `λ³ = 1 mod N`, checked on 5 ladder points; `T` **not** in the 1,536-point orbit | a wrong `β`/`λ` pairing — excluded by the 5-point check |
| 7 | no `k = a + bλ` with `\|a\|,\|b\| < 2²¹` | Q | group, mod `p` | exhaustive over that box | **code-audited only, not re-executed** | as row 4 |
| 8 | the endomorphism buys **√3 only** | Q | group | — | **re-verified**: `\|⟨λ⟩\| = 3`, so 2¹²⁸ → 2¹²⁷·². It does not move the figure | — |
| 9 | **the solution is essentially unique** | coordinator | arithmetic | exact | **re-computed**: `2²⁵⁶ − N = 432420386565659656852420866394968145599 ≈ 2^128.35`. `k·G = T` fixes `k mod N`; a second ON-set exists iff `k₀ < 2²⁵⁶ − N`, probability `2^-127.65` | nothing — this is arithmetic |
| 10 | **NEW — Q's slot-collision caveat is VACUOUS for every `\|S\| ≤ 42`** | **X** | instance | exact | **computed here**: two children of a slot coincide iff `Σ_{S1} 2^i − Σ_{S2} 2^i = ±N` for disjoint `S1, S2 ⊆ S`; that is a signed-binary representation of `N` with `\|S1\|+\|S2\|` nonzero digits, and by Reitwiesner's theorem the minimum such is the **NAF weight of `N`, which is 43** (`weight(N) = 192` in plain binary). So no two children can ever coincide when `\|S\| ≤ 42` | an error in the collision criterion itself, not in the arithmetic |

**Row 10 matters more than its size suggests.** Q listed the collision criterion as *open* in its §4
because it "needs the particular scalar". It does not: it is bounded uniformly by a property of `N`
alone. Everything the fleet is doing at low `|S|` — my sweep, T's integer lifts at `|S| = 2,3,5,6,7,8,17`,
V's handle-less census — sits far inside the safe region, so the degenerate-branch escape hatch
(where a stage gadget imposes **nothing** because `dx = dy = 0`) **never opens** there.

---

## 2. DO Q'S SEARCHES HAVE INSTANCE-LEVEL STANDING? — my ruling: **YES**, for the negative direction

Q withdrew its six searches (§15) on the ground that they computed the fold inside the group model
without checking the circuit agrees, restored them (§24) once the chain was measured mod `p`, and
its thread closed before it was asked to rule on the later evidence. Here is the ruling.

The searches are **negative** results, and a negative needs exactly one implication:

> a satisfying assignment with ON-set `S` **⟹** `k_S·G = T`.

Its contrapositive — no `k` in family `F` satisfies `k·G = T`, therefore no satisfying assignment has
ON-set in `F` — is what makes the search a statement about the instance. That implication needs
four links, and **all four are measured**:

1. **All 39,033 equations = 0 over ℤ ⟹ every gate atom = 0.** F's incidence matrix has
   `ker(M) = 0` by three independent computations; T re-verified the certificate from cold and
   separately established **faithfulness** of the decomposition by exact list equality against
   `checker.evaluate_all`. Note this link needs nothing about moduli: the assignment is integral and
   the equations vanish over ℤ, so `M·a = 0` over ℚ and `ker(M) = 0` gives `a = 0`.
2. **Atoms zero ⟹ the leaf, gadget, mux and tree identities hold.** Q: 256/256 leaves on the
   cubic and 256/256 satisfying `L_i = 2^i·G`; 383/383 chord gadgets Schwartz–Zippel-verified
   against the real sub-DAG with orientation `(+1,+1)`; 383/383 muxes implementing
   identity / pass-through / **sum** in both coordinates with identical coefficient wires;
   382 parent←child edges, one root, all 256 selectors beneath it.
3. **The coordinate hand-off.** Agent L: every slack wire is a constant multiple of `p` times a
   free variable, **3,681/3,681, zero exceptions**, so `parent_input ≡ mux_out (mod p)`
   unconditionally.
4. **The degenerate branch cannot be used to escape.** Row 10 above: for `|S| ≤ 42` no slot ever
   sees two equal live children, so no gadget is vacuous.

**Because links 2–4 are identities mod `p`, and ℤ-truth implies mod-`p` truth, the modulus is not a
gap in this direction.** `p` is precisely the modulus the negative needs. The 927 `c > 1`
divisibility conditions over ℤ sit on the **converse** (existence) and are irrelevant to a negative.

So I rule: **Q's six searches, and mine, are statements about the instance** — stated as
*"established mod `p`"*, which for a negative is the same thing.

**The one residual risk, stated plainly.** Link 1 rests on **F's** decomposition. This lab holds five
mutually inconsistent atom counts — 39,033 (F), 39,277, 40,727, 40,885 (I), 42,267 (A/G/H) — with
kernel dimensions 0, ≥1,852 and ≥3,234 respectively. `ker(M) = 0` is a fact about F's parse; T's
faithfulness check was **10 sampled points**, which is evidence, not a proof, that F's `M·a` equals
`evaluate_all` as functions. If F's decomposition splits or merges atoms the way Q's own 47,198-term
parse did, link 1 weakens. **That, and not the modulus, is where I would push if someone wanted to
attack these negatives.** Concretely: re-run T's faithfulness check at 10⁴ random points rather
than 10, and the last soft spot closes.

---

## 3. IS THERE *ANY* PER-BIT INFORMATION? — measured, and the answer is no

The coordinator's structural observation is **exactly right**, and here is the measurement behind it.

### 3.1 Quantitatively, from the searches: the posterior on any single bit moves by `< 2^-201`

Every constraint in §1 is of the form `k ∉ F` for an explicitly enumerated `F ⊆ [0, 2²⁵⁶)`. Their
sizes:

| family | size |
|---|---|
| weight ≤ 9 | 11,711,713,815,280,289 = 2^53.38 |
| BSGS both ends, 2·2^52 | 2^53.00 |
| `k = a + bλ`, `\|a\|,\|b\| < 2²¹` | 2^43.58 |
| 34-bit window | 2^41.80 |
| small multiple `m ≤ 10⁷` | 2^31.25 |
| endomorphism orbit | 1,536 |
| **union (upper bound)** | **< 2^53.9** |

A single-bit fiber `{k : bit_i(k) = b}` has `2²⁵⁵` elements. The searches remove at most `2^53.9`
of them, a fraction **`2^-201.1`**. Both fibers of every bit are hit essentially equally. **No
search performed by this campaign has moved the posterior on any individual bit of `k` by more than
`2^-201`.** That is the precise sense in which every constraint is global.

### 3.2 Structurally, from the group: per-bit information *is* the whole problem

The fold `S ↦ Σ_{i∈S} 2^i·G` is a group homomorphism of the selector vector, and `k ↦ k·G` is
random self-reducible: for any `s ≠ 0`, `r`, the point `sT + rG` has discrete log `sk + r`, uniform
over `ℤ_N`. So an algorithm that reads one bit of the log of an arbitrary target reads one bit of
`sk + r` for freely chosen `s, r` — and the standard hidden-number-problem reduction turns that into
the whole log. **Any method that extracts information about a single bit position of `k` from `T`
solves the ECDLP on secp256k1.** This is also why meet-in-the-middle is generic here: the only
handle is the *weight*, which is a symmetric function of the bits, and every family in §1 is either
symmetric under permuting positions (weight, small multiple, orbit) or defined by an interval
(window, BSGS) — none is a per-bit predicate.

### 3.3 From the instance: no selector is distinguished by its value

Measured directly on `EQUATIONS.txt` (`xperbit.py`), over the 256 leaf-selector wires:

* **0 / 256** appear in an equation that mentions no other wire — **no selector is pinned**.
* Q's unit propagation, at ON-weights 1, 2, 3, 5, 7 and 128, solved **0 / 256** selectors and found
  **0 contradictions** — the routing is determined by a simultaneous system, not by propagation.
* Q measured the mux at **383 / 383 slots** implementing identity, pass-through and sum. Every slot
  can pass, drop or add, so **every subset of the 256 leaves is realisable by the routing layer**;
  no ON-set is structurally forbidden.
* The selectors *are* structurally distinguishable — footprints run from **77 to 185 occurrences**
  across **30 to 51 equations**, and the gadget census splits the leaves **178 paired / 78
  pass-through**. But that is *where a leaf sits in the fold tree*, fixed by the circuit and
  identical for every candidate `S`. It carries no information about **whether that bit is ON**.

### 3.4 The one place where non-uniform information could still appear — and it is not per-bit

The integer lift is indexed by `|S|`, not by position. T has closed it at `|S| = 2, 3, 5, 6, 7, 8, 17`
and found a solver-coverage gap (handle-less atoms) at `|S| = 32`. If the lift genuinely failed at
some weight, that would be information about **the weight of `k`** — a global quantity, and
exactly the quantity my sweep bounds from below. It would still say nothing about any individual
bit. **So the structural claim survives even in the one direction that could have broken it:
there is no per-bit information to be had, and the fleet should stop looking for it.**

---

## 4. WHAT THIS LEAVES

`k` is a 256-bit scalar with weight ≥ 9 (≥ 10 pending), not near either end of `[0,N)` by 2⁵²,
not window-confined, not in the endomorphism orbit, not a small multiple — and **essentially
unique**, so there is no "many needles" advantage. The residual search space is `2²⁵⁶ · (1 − 2^-202)`.
Generic DLP on secp256k1 is ~2¹²⁸ with a √3 discount. Nothing in the circuit reduces it, because the
fold is a homomorphism of the selector vector.

**Be blunt about what the sweep is.** For a uniformly random 256-bit scalar, `P(weight ≤ 9) ≈ 2^-203`.
It is a lottery ticket. It is worth buying only because it is cheap and because a hit is a complete
solve. A miss is *"weight ≤ w exhausted, no solution"* — a real, citable bound — and **nothing more**.
