# LIFT_SUPPORT — agent AF. Is the integer lift `|S|`-dependent, analytically?

**Answer: NO. Neither condition family is a function of `|S|` in any way that can obstruct.**
Stated as Theorem 5 below, with hypotheses inside the statement.

Everything here was re-derived in this directory from `EQUATIONS.txt` with **my own parser**
(`af1_parse.py` … `af24_check.py`). No other agent's code, model or pickle was imported. Where I
reuse another agent's *result* it is attributed inline and, where I could, spot-checked.

**Baseline verified before starting:**
`python3 solve_lab/checker.py solve_lab/best/new_instance_partial_39026.json`
→ `satisfied 39026/39033 (7 failing) [12231, 12270, 12350, 14584, 18673, 22044, 29125]`. ✔

**I did not beat 39,026.** Nothing here is a better partial. No infeasibility is claimed.

---

## 0. TL;DR for the coordinator's four questions

| your question | answer |
|---|---|
| **(1)** extract the conditions symbolically | done, from my own parse; the `927 / 2780` split reproduces **exactly** (§1) |
| **(2)** is every condition small-support, so `|S|` is invisible? | **No — the short argument fails.** 310 of the 927 have support > 1, up to **256** (the root). §8 does **not** die by support (§3) |
| **(3)** does `c·P \| R` get harder as leaves switch on? Does no-wraparound kill it? | **It kills it**, and for a reason *stronger* than no-wraparound alone: the system contains **no bounded data variable at all**, so no growth argument can bite. No-wraparound then supplies the second half — it makes the per-block ledger *identical at every block for every `S`* (§4, §5) |
| **(4)** does the irregular shape (178\|78, 101 pass-through, 27 dead) make a condition see the *count*? | **No.** The live-block law `Σ L_v = \|S\|−1` is **shape-independent** (proved, §2.4). The 128 constant-gate blocks are `\|S\|`-blind by construction. Reconciles with Z (§6) |

**And a fact the fleet's map is out of date on:** `|S| = 128` is **not** a stall. Two artefacts
(`agentT_work/close_T128s59.json`, `close_T128s7fix.json`) are checker-verified at **39,018/39,033**
with the **byte-identical 15-equation footprint** of the `|S| = 8/17/32/64` closures, and I evaluated
all 3707 lift conditions on them with my own parse: **0 violated `c>1` conditions** (§7). Agent T
found this and recorded it in `RESUME_T.md` §BD/§BE; **AB's `UPPER_BOUND_MAP.md` §S5 still rests on
the stall and should be corrected.**

---

## 1. The conditions, extracted (re-derivation, not a re-reading)

Pipeline: recursive-descent parse → additive peel into atoms → definition DAG → alias closure →
constant propagation → locate `P` → the `W = M·u` slack family → their consumers.

* `P = x_26064 = 115792089237316195423570985008687907853269984665640564039457584007908834671663`
  `= 2^256 − 2^32 − 977` ✔, with a **220-member alias class**; independently validated — all 220
  are assigned exactly `P` in `new_instance_partial_39026.json`.
* **3707 slack wires** `R` with a unique definition `R = P·u`, `u` a **free** variable
  (0 definitions, 3707/3707).
* Each `R` has exactly one *consumer* atom `Expr − c·R`. Therefore each condition is exactly

  > **`c·P | Expr`**, `u = Expr/(c·P)` free.

* **`c = 1` for 2780; `c > 1` for 927.** All 927 multipliers distinct, range `193 … 16 718 259`,
  893 distinct primes, 202 non-squarefree. **This reproduces the coordinator's 927/2780 exactly
  from an independent parse.**

### 1.1 The five families (my classification of all 3707)

| family | shape of `Expr` | n | of which `c>1` | active when |
|---|---|---|---|---|
| **congruence** | `L · B`, `B = α·N1 + β·N2` | **1149** = 3·383 | **288** | `L = 1` |
| **off-pin** | `(1−L) · i_j`, `i_j` FREE | **766** = 2·383 | **192** | `L = 0` |
| **leaf pin** | `s_i · (w_i − C_i)` | **512** = 2·256 | **256** | `s_i = 1` |
| **difference** | `X − Y` (ungated) | **766** = 2·383 | **191** | always |
| **product** (ungated) | `A·B` | **512** | 0 | always |
| other | `V − K` | 2 | 0 | — |

`1149 + 766 = 1915` congruence+off-pin atoms — **identical to agent W's counts**, reached here
without reading W's code. `288 + 192 + 256 + 191 = 927` ✔.

### 1.2 A parser bug I made and caught — the same one agent U documented

My first `peel_sum` walked the whole left-deep `+/−` chain, so when a body's **leading summand was
itself a subtraction** `(x_a − x_b)` it was emitted as a spurious *copy* atom. That over-merged the
alias classes and made **254 of 383 liveness gates evaluate to the constant 0** — an absurdity loud
enough to catch it. Fix: a chain step is *only* `(const)*(ATOM)`; peeling stops at the leading
summand. Copies `5833 → 3738` (U's corrected parser: 3749).
**Blast radius, measured before reporting:** `P`-alias count (220), lift-atom count (3707) and the
`927/2780` split were **unchanged**; what changed was the gate algebra, which went from nonsense to
`128 const-0 + 255 AND`, and the deliverable's nonzero-consumer count `33 → 3`. This is agent U's
§7 failure mode reached independently; it is worth a permanent warning line in FLEET.md.

---

## 2. The gate algebra — the load-bearing structure

### 2.1 Every gate expands to a monotone read-once formula
Expanding all 383 gates through the definition DAG (`af11_gates.py`) terminates with **0 opaque**:

* **128 gates are identically `0`.** (= the 101 pass-through + 27 dead blocks.)
* **255 gates are `AND`.** Operator census inside the 255 formulas: `255 and`, `1604 or`,
  `2114 leaf`, **zero `not`**.
* Every one of the 255 is exactly `L_v = OR(I_v) ∧ OR(J_v)` with `I_v ∩ J_v = ∅`,
  both arms **pure OR-trees over distinct leaves** (255/255, 0 exceptions).

> **Gate law.** `L_v(s) = [S ∩ I_v ≠ ∅] · [S ∩ J_v ≠ ∅]`.

All 256 selectors are forced boolean: exact polynomial test finds an atom equal to `k·(s²−s)` for
**256/256** (my first shape-matcher found only 165 — the shape matcher was wrong, not the instance).

### 2.2 The 510 slot supports are a binary tree — independently reproduced
* 510 distinct slot supports, **0 laminarity violations**.
* Exactly **2 maximal sets, of sizes 178 and 78**, disjoint, covering all 256.
* `I_v ∪ J_v` is a node of the family (or the root) for **255/255** gates; the 255 unions are
  distinct and are *exactly* the 255 internal nodes. Node census **511 = 256 leaves + 255 internal**.

This reproduces agent U's `511 / laminar / 178|78` from a third independent parse (U's own, L's
model, K's — now mine).

### 2.3 Cross-check against every existing artefact
For all **20** `agentT_work/close_*.json` plus the 39,026 deliverable, the gate law predicts the
live-block count and it is right in **20/20** cases (`af23`, `af24`). E.g. the deliverable has
exactly 2 ON selectors → **1** live merge block, matching U §6 ("exactly 1 of the 383 chord stages
has non-zero inputs") and W §3 ("1 degeneracy block, gate on") by a completely different route.

### 2.4 THEOREM 1 (live-block count) — proved, not sampled
> Let `T` be **any** binary tree with leaf set `Λ`, and for an internal node `v` with children
> `a, b` put `L_v(S) = [S∩Λ_a ≠ ∅]·[S∩Λ_b ≠ ∅]`. Then for every `S ⊆ Λ`,
> **`Σ_v L_v(S) = max(|S| − 1, 0)`.**

*Proof.* Induct on the subtree. Let `g(v)` be the number of live nodes strictly inside `v` and
`n_v = |S ∩ Λ_v|`; claim `g(v) = max(n_v − 1, 0)`. Leaf: `g = 0 = max(n−1,0)` for `n ∈ {0,1}`.
Internal with children `a,b`: `g(v) = g(a)+g(b)+[n_a>0][n_b>0]`. If both `n_a,n_b > 0`,
`= (n_a−1)+(n_b−1)+1 = n_v−1`. If exactly one is `0`, `= (n_v−1)+0+0`. If both `0`, `= 0`. ∎

**The tree shape does not enter.** So the 178|78 root split, the 101 pass-throughs and the 27 dead
blocks are irrelevant to the *count*; they only change *which* leaves a given block separates.
Verified numerically over `|S| = 0 … 256`, 51 random sets: deviation from `|S|−1` is 0 everywhere.

**Consequence (the only channel by which `|S|` is visible at all):**
`#active congruence conditions = 3(|S|−1)`, `#active off-pin conditions = 766 − 2(|S|−1)`,
`#active leaf pins = 2|S|`. **`|S|` *is* recoverable from the gate pattern.** So §8 is not killed
by blindness; it has to be killed by the ledger.

---

## 3. THEOREM 2 (support) — the short argument you asked me to try FIRST, and it fails

Exact selector-support of all 927 conditions with `c > 1` (`af21_supp.py`):

| sub-family | n | support size |
|---|---|---|
| congruence @ constant-0-gate block | 99 | **0** (identically vacuous: `c·P \| 0`) |
| difference, `X ≡ Y` by an explicit copy atom | 180 | **0** (identically vacuous) |
| off-pin @ constant-0-gate block | 71 | **0** (free wire, no selector in it) |
| difference, free LHS | 11 | ≤ 4 |
| **leaf pin** | 256 | exactly **1** |
| **off-pin @ merge block** | 121 | 2 … **178**, median 3 |
| **congruence @ merge block** | 189 | 2 … **256**, median 3 |

Distribution of the 310 gated `c>1` conditions on merge blocks:
`>1: 310, >2: 191, >4: 97, >8: 56, >16: 30, >32: 15, >64: 6, >128: 3, =256: 1`.

> **THEOREM 2.** 361 of the 927 conditions have selector support `∅` and are identically vacuous or
> identically satisfiable; 256 have support exactly one selector. But **310 have support ≥ 2 and one
> has support 256.** Therefore **the locality argument does not close §8**: some conditions can, in
> principle, see `|S|`.

This is a negative result and I report it as one. §8 has to be attacked at the arithmetic.

---

## 4. What the conditions actually *are* — the per-block ledger

### 4.1 Off-pin side: measured, exhaustive
* **766/766** off-pin residuals are **free variables** (zero definitions).
* Each appears in only **3–6 short atoms**, all inside its own block's cluster, and
  **766/766** appear in exactly one mux atom `· − (L_v · i_j)`.
* So when `L_v = 0`: the block's 3 congruences are `c·P | 0` (vacuous) and its 2 off-pins are
  `c·P | i_j` with `i_j` free — **satisfied by `i_5 = i_6 = 0`**, which is also the mux's own value
  (`L_v·i_j = 0`). **Cost zero, at every `S`.**

### 4.2 Congruence side: the Jacobian, re-derived here
Expanding a merge block's two residual differences as multivariate polynomials in the free wires
(`af20_law.py`, exact sparse expansion, no truncation):

* `N1`: 13 monomials, total degree 3, **linear in `i5`** with a **3-monomial degree-2** coefficient,
  **degree 0 in `i6`**.
* `N2`: 6 monomials, total degree 2, **linear in both**, each with a **2-monomial degree-1**
  coefficient.

A 3-monomial degree-2 form in `{i1,i2}` is `(i1−i2)² = A²`; 2-monomial degree-1 forms are
`i4−i3 = B` and `i1−i2 = A`. Hence, **from my own parse**,

> `∂(N1, N2)/∂(i5, i6) = [[A², 0], [B, A]]`,  `det = A³`.

which is agent W's Jacobian, obtained here by counting monomials rather than by trusting the
`E·A²−B²` / `A(i3+i6)−B(i2−i5)` formulas. (I did **not** independently re-derive rank-2 of the 3×2
`(c_k1,c_k2)` matrix — my reconstruction could not pair the six per-block difference wires. That
fact is **attributed to W**: exhaustive over 383 blocks, all six 2×2 minors nonzero,
`max|minor| = 260 582 651 319 840 < 2^48 ≪ P`.)

### 4.3 The arithmetic content — measured over all 383 blocks
* **At most ONE congruence row per block has `c > 1`** (288 blocks have exactly 1, 95 have 0).
* **At most ONE off-pin row per block has `c > 1`** (192 have 1, 191 have 0).
* For all **288/288** `c>1` congruence rows: **`gcd(c, α) = gcd(c, β) = 1`.**
  (In particular `c` never divides out — the conditions are genuine, not decoration.)
* `gcd(a, P) = gcd(c, P) = 1` throughout (`P` prime, `c ≤ 2^24`).

### 4.4 THEOREM 3 (the ledger) — the same at every block, for every `S`
> Fix a merge block `v` with `L_v = 1`, and write `N1 = P·n1`, `N2 = P·n2` once the mod-`P` law
> holds. Then, by CRT (`gcd(c,P)=1`), block `v`'s five conditions are exactly
>
> * **two** conditions mod `P`: `N1 ≡ N2 ≡ 0 (mod P)`. The 3 rows collapse to these two because the
>   3×2 coefficient matrix has rank 2 mod `P` (W). They are solved **uniquely** by `(i5, i6) mod P`,
>   because `det ∂(N1,N2)/∂(i5,i6) = A³ ≢ 0 (mod P)` — **this is exactly what U's partition theorem
>   buys**: no proper slot support has masked value ≥ `N`, so no merge block ever sees two equal or
>   opposite inputs, so `A ≢ 0 (mod P)` at **every** merge block for **every** `S`.
> * **at most one** condition mod `c`: `α n1 + β n2 ≡ 0 (mod c)`.
>
> The residual freedom after the mod-`P` step is `(i5, i6) → (i5 + P t5, i6 + P t6)`, `t5,t6 ∈ ℤ`,
> which moves `α n1 + β n2` by `(α A² + β B) t5 + β A t6 (mod c)`.
>
> **`A` and `B` are freely adjustable mod `c`.** `B = i4 − i3` is a difference of the two input
> *y*-values; every leaf `y`-pin has `c = 1` (measured: all 256 leaves carry `m = 1` on `Y` and
> `m > 1` on `X`), so a leaf `y` is `≡ D (mod P)` and, `P` being invertible mod `c`, **free mod `c`**;
> an interior input is a child's free chord wire, likewise free mod `c`. Since `gcd(β,c) = 1`,
> `α A² + β B` can be made a **unit** mod `c`, and then `t5` alone discharges the single congruence.
>
> **Ledger per live block: 2 free knobs, 2 mod-`P` conditions + ≤1 mod-`c` condition, one knob
> spare.** Per dead block: 2 free knobs, 2 conditions, both discharged by `0`.

---

## 5. THEOREM 5 (main) — the lift is `|S|`-blind

> **Hypotheses.**
> **(H0)** Atom-level semantics (`every atom = 0`), which is the *sound* direction: an atom-level
> solution is an equation-level solution, so a **positive** result at atom level transfers.
> (A negative one would not — see §8.)
> **(H1)** U's partition theorem: no proper slot support has masked value ≥ `N`; equivalently
> `A ≢ 0 (mod P)` at every merge block, for every `S`. *(U, exhaustive; not re-derived by me.)*
> **(H2)** W's rank-2 of the 3×2 congruence matrix mod `P`, all 383 blocks. *(W, exhaustive.)*
> **(H3)** Measured here, exhaustively: 766/766 off-pin residuals free and block-confined;
> ≤1 `c>1` congruence and ≤1 `c>1` off-pin per block; `gcd(c,α)=gcd(c,β)=1` for all 288;
> all 256 leaf `Y`-pins have `c = 1`.
> **(H4)** Measured here: **no data variable in the instance is bounded.** The only variables
> carrying a booleanity atom are the 2340 boolean control bits; **0 of the 766 chord-output wires
> and 0 of the 512 leaf coordinate wires** carry one, and no atom anywhere imposes a range.
>
> **Conclusion.** Under (H0)–(H4), **no condition in the 927-family and none in the 766 off-pin
> family is a function of `|S|` in any way that can obstruct.** For every `s ∈ {0,1}^256`, the whole
> system of 3707 lift conditions is dischargeable, block by block bottom-up, with a ledger that is
> **identical at every block and independent of `|S|`**. The only surviving `S`-dependence in the
> instance is the mod-`P` layer, whose content is precisely `Σ_{i∈S} 2^i·G = T` — a statement about
> **which** `S`, not about **how many**.

### 5.1 Why the growth mechanism you suspected cannot exist — two independent reasons

1. **There is nothing to grow against.** Every condition is `c·P | Expr` with the quotient handle
   `u = Expr/(c·P)` a *free, unbounded* variable, and by (H4) no variable in the data path is
   bounded. `Expr` does grow with the number of active leaves — the wire values are curve
   coordinates plus arbitrary `P`-multiples — but `c·P` is a fixed constant per condition
   (`≤ 16718259·P ≈ 2^280`) and `u` simply grows with `Expr`. **A growth/size argument needs a
   bounded variable somewhere in the system. There is not one.** This is a stronger and simpler
   statement than anything about wraparound.
2. **No-wraparound makes the ledger uniform.** `A ≢ 0 (mod P)` everywhere ⇒
   `det ∂(N1,N2)/∂(i5,i6) = A³` is a unit mod `P` at every live block ⇒ the two free chord wires
   solve the two mod-`P` congruences *exactly and identically at every block, for every `S`*.
   **So: no wraparound ⇒ the fold is uniform ⇒ the per-block ledger is `|S|`-blind.** That is your
   first branch, and it is the one that holds. It does not merely *bound* the growth mechanism; it
   removes the only place where block-to-block variation could have accumulated.

### 5.2 The counterfactual, and the bound it *would* have given
Block-locality means that if some set `Bad` of merge blocks were unsatisfiable in live mode, the
feasible selector sets would be exactly `{S : ∀v ∈ Bad, S∩I_v = ∅ or S∩J_v = ∅}` — a **downward-
closed set condition**, whose maximum weight is a tree DP. Computed on my recovered tree
(`af23_bound.py`):

| hypothetical `Bad` | resulting bound `w ≤ …` | status |
|---|---|---|
| all 255 merge blocks | **1** | refuted by the closures at `\|S\|` = 2 … 128 |
| the 189 blocks with a `c>1` congruence | **17** | refuted by `\|S\|` = 32, 64, 128 |
| the 121 blocks with a `c>1` off-pin | **44** | refuted by `\|S\|` = 64, 128 |
| the root block only (split 78\|178) | **178** | consistent with every closure, but far above §8's payoff band and above the null mean 128 |
| `Bad = ∅` (Theorem 5) | **256** — no bound | — |

> **THEOREM 4 (the ceiling on §8, independent of Theorem 5).** Any bound this mechanism can deliver
> has the form `w ≤ maxS(Bad)`. Because `close_T128s59.json` and `close_T128s7fix.json` are
> checker-verified closures at `|S| = 128` with **0 violated `c>1` conditions** (§7), `Bad` is
> avoided by a 128-set, hence **`maxS(Bad) ≥ 128`**.
> **§8 therefore cannot produce any bound below `w ≤ 128` — i.e. below the null mean — let alone
> reach its payoff band (`w ≲ 56` to beat rho, `w ≲ 24` to be actionable).**
>
> This holds *even if Theorem 5 is wrong*. It depends only on block-locality (§4.1, measured) and on
> two checker-verified artefacts.

---

## 6. Reconciliation with agent Z — no contradiction, and a strengthening

Z established: **zero linear constraints on the selectors, zero adder-shaped atoms, configuration
space exactly `2^256`.**

* My gate functions `L_v = OR(I_v) ∧ OR(J_v)` are **strictly non-linear** (monotone AND-of-ORs, up
  to depth 9). The identity `Σ_v L_v(s) = |s| − 1` is a *non-linear polynomial identity* on
  `{0,1}^256`, **not** a linear constraint. Z's finding predicts exactly this: any `|S|`-dependence
  must be non-linear, and mine is.
* My conclusion is that **no condition excludes any selector vector**, so the selector configuration
  space remains exactly `2^256` — Z's count, unchanged.
* **Strengthening.** Z ruled out constraints of degree 1. Theorem 5 rules out constraints of *any*
  degree that are functions of `|S|`, for the 927 + 766 families. The two results agree and mine is
  strictly the stronger statement on this sub-question.

**If anyone finds a contradiction here, it is mine to lose:** Z's `2^256` count is a census, mine is
a derivation with five hypotheses.

---

## 7. Empirical confirmation on artefacts that already exist (no probes run)

I evaluated **all 3707 conditions** `c·P | Expr` with my own parse on every `agentT_work/close_*.json`
and on the deliverable (`af24_check.py`). This is a measurement on existing files, not a probe.

| | |
|---|---|
| artefacts checked | **20** closures + the 39,026 deliverable |
| `#live merge blocks == \|S\| − 1` | **20/20**, `\|S\|` from 2 to 128 |
| artefacts with **0** violated `c>1` conditions | **19/20** |
| the one exception | **`close_T8.json`, at `\|S\| = 8`** — one `c = 15 333 171` congruence. `close_T8pair.json` and `close_T8w.json`, same `\|S\| = 8`, are clean |
| violated conditions, all artefacts | **exactly 1** in every case — the *same* ungated `c = 1` difference (the root-vs-target residual), present already at `\|S\| = 2` |

Two readings worth stating plainly:

1. **The lift closes at `|S| = 128` with zero `c>1` violations**, and `close_T128s59.json` /
   `close_T128s7fix.json` score **39,018/39,033** with the *byte-identical* 15-line footprint
   `[4573, 7123, 7469, 9648, 11854, 16622, 17726, 21382, 25539, 28653, 29437, 31061, 32894, 32916,
   34517]` as `close_T64.json`. Re-verified with `checker.py` here.
2. **The only `c>1` lift violation anywhere in the fleet's record is at `|S| = 8` — the low end —
   and is not reproduced by two sibling artefacts at the same `|S|`.** Whatever it is, it is a
   property of that closure run, not of `|S|`; and it is the wrong sign for a monotone story.

---

## 8. Scope, and what I did **not** establish

| claim | status |
|---|---|
| `927 / 2780` split; the five families; all counts in §1.1 | **exhaustive**, my own parse |
| gate law `L_v = OR(I_v)∧OR(J_v)`, 255/255 read-once, `I∩J=∅` | **exhaustive**, 0 exceptions |
| 510 supports laminar, 178\|78, 511-node binary tree | **exhaustive**, third independent reproduction |
| Theorem 1 (`Σ L_v = \|S\|−1`) | **proved** for any binary tree; verified on 20/20 artefacts |
| Theorem 2 (support table) | **exhaustive** over all 927 |
| Jacobian `[[A²,0],[B,A]]` | **re-derived** here by exact monomial expansion |
| rank 2 of the 3×2 matrix | **attributed to W**, not re-derived (my per-block difference pairing failed) |
| `A ≢ 0 (mod P)` everywhere | **attributed to U**, not re-derived |
| ≤1 `c>1` row per block per side; `gcd(c,α)=gcd(c,β)=1` | **exhaustive**, mine |
| no bounded data variable (H4) | **exhaustive** over all atoms, mine |
| **Theorem 5** | **a derivation, not a construction** — see below |
| **Theorem 4** (`maxS(Bad) ≥ 128`) | **solid**: block-locality (measured) + 2 checker-verified artefacts |

**The one real gap in Theorem 5.** I prove the ledger balances *locally* — every block has a spare
knob, and the knob a parent needs (its inputs' residues mod `c`) is exactly the spare its children
have. I do **not** exhibit a global bottom-up construction and run it through `checker.py`; that
would be a probe, and you asked for a derivation. The composition argument is CRT-clean
(`gcd(c,P)=1` separates the mod-`P` layer from the mod-`c` layer completely), the tree is laminar so
no wire is shared between two parents, and the recursion bottoms out at leaf `Y`-coordinates which
are free mod `P` and hence free mod `c`. I believe it, but it is a chain of five hypotheses, not a
certificate. **§7 is the strongest available substitute and it agrees at `|S| = 128`.**

**Direction matters.** Everything above classifies solutions of `atoms = 0`, which is a *subset* of
`equations = 0` (W's §5 boundary). My theorem is a **positive** ("this is satisfiable"), so the
inclusion runs the safe way: an atom-level solution is an equation-level solution. Had I concluded
"infeasible for large `|S|`", the same boundary would have made it unsound. Agent K's two wrong
verdicts were both negatives; this is not one.

**What would settle Theorem 5 completely, and what it costs.** Drive the bottom-up construction of
§4.4 explicitly at `|S| = 192` and `|S| = 250` on independent seeds, dump the full assignment, and
score it with `checker.py`; expect 39,018/39,033 with the same 15-line footprint. Agent T's
`t_close2wj.py` already does exactly this and took 1,173 s at `|S| = 128` (seed 59) and 2,289 s
(seed 7). **Cost: ~1–2 CPU-hours for the pair.** T's logs `t_c2wj_T192s47.out` /
`t_c2wj_T250s31.out` exist but the corresponding `close_*.json` do **not**, so those two points are
currently unverified — that is the single cheapest remaining experiment on this question.

---

## 9. Files

| file | what |
|---|---|
| `af1_parse.py` → `af_ast.pkl` | tokeniser + recursive-descent parser; the corrected additive peel |
| `af2_atoms.py` → `af_atoms.pkl` | 50,917 atoms in 3,028 shapes |
| `af3_defs.py`, `af4_P.py` → `af_defs.pkl`, `af_P.pkl` | definition DAG, alias closure, `P` and its 220 aliases |
| `af5_lift.py`, `af8_cond.py` → `af_lift.pkl`, `af_cond.pkl` | the 3707 conditions and the `927/2780` split |
| `af9_tree.py`, `af10_blocks.py` → `af_tree.pkl`, `af_blocks.pkl` | 256 selectors, 383 gate/complement pairs |
| `af11_gates.py` → `af_gates.pkl` | the 383 gates as boolean formulas |
| `af12_lam.py` → `af_lam.pkl` | read-once check, laminarity, the live-count law |
| `af13_map.py`, `af17_arith.py` → `af_map.pkl`, `af_rows.pkl` | condition→block map, the arithmetic census |
| `af14`–`af16`, `af18`, `af20` | residual anatomy, the per-block ledger, the Jacobian |
| `af19`, `af21`, `af22` | tree verification, exact supports, bounded-variable census |
| `af23_bound.py` | deliverable cross-check + the counterfactual tree DP (Theorem 4) |
| `af24_check.py` | all 3707 conditions evaluated on all 20 existing closures |
