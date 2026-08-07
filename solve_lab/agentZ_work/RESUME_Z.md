# RESUME_Z — Agent Z: can `w = |S|` be bounded directly from the instance?

**Verdict: NO. There is no constraint on `|S|` anywhere in the instance. The configuration
space really is the full `2^256` (minus the trivial `w ≥ 1`).** This is a clean confirmation,
established by a method the campaign had not run, not a restatement of prior belief.

Everything below is from my own parse of `EQUATIONS.txt`. No file from `agentP_work/`,
`agentS_work/`, `agentT_work/` or any other agent dir was imported or executed. I read
`RESUME_P.md`, `RESUME_S.md` and `selcouple*.py` for the *statements* under audit only.

## 0. Baseline re-verified
`python3 solve_lab/checker.py solve_lab/best/new_instance_partial_39026.json`
→ **39026/39033**, failing `[12231, 12270, 12350, 14584, 18673, 22044, 29125]`. Confirmed.

## 1. My pipeline (run in this order)
| file | what it does |
|---|---|
| `zparse.py` | recursive-descent parser; reduces each eq `scalar·L^k = 0` to its linear form `L`, splits `L` into additive atoms |
| `zsel.py` → `zsel.json` | selector identification by pure regex on the raw file |
| `zatoms.py` → `zatoms.pkl` | polynomial normalisation of every atom, **each expansion validated at 3 random points against direct tree evaluation** |
| `zclass.py` → `zclass.pkl` | classification of every atom that mentions a selector |
| `zaffine_a.py` → `zlin.pkl` | booleanity reduction (`x²→x` for all boolean vars), keeps the deg≤1 equations |
| `zelim.py` / `zelimQ.py` | sparse Markowitz elimination of every non-selector variable — mod `2^61−1` **and** exactly over ℚ |
| `zwire_full.py` | same system solved **symbolically**: every determined wire as an affine function of the 256 selectors |
| `zprop.py`, `zcount.py` | boolean/liveness layer; adder-vs-OR census |
| `zsep.py` | global monomial separability + **checker.py-anchored** weight sweep |

My parse: 39,033 equations, every one `scalar·L^k = 0`, **563,307 atom occurrences**,
**93,515 distinct atoms** (polynomial-normalised), degrees `{1: 58050, 2: 35464, 0: 1}`.

## 2. Selectors — independently re-derived
Regex `\(x_S\)\*\(\(x_C\)-\(K\)\)` over the raw file: **10,231 hits, 512 distinct triples,
256 distinct `S`, 512 distinct `C`, 512 distinct `K`** (bit lengths 287–296), plus exactly
**2** 256-bit literals (`P`, `Q`). So: **256 selectors, 2 coordinates each — confirmed.**

## 3. TASK 1 — multi-selector atoms: my count, and the reconciliation

| parse | atoms touching ≥2 distinct selectors |
|---|---|
| agent P | 0 |
| agent S | 48 (1 booleanity certificate + 47 "bundled") |
| **agent Z (this)** | **0** |
| **agent Z, coarsest possible unit (= a whole equation)** | **2,490** |

**The three numbers are not in conflict; they measure different granularities**, and the
ordering is `Z-atoms (finest) ⊂ P-atoms ⊂ S-atoms ⊂ equations`. A count at any of those
levels is a fact about that decomposition, so I computed the **granularity-free invariant**
instead — the expanded polynomial of an equation is unique, so its *monomials* are
parse-independent:

> **Over all 39,033 equations there are 819,975 monomials. 798,787 contain no selector,
> 21,188 contain exactly one, and 0 contain two or more.**
> (`zsep.py`; with multiplicity: 17,795 monomials of selector-degree 1 and 3,393 of
> selector-degree 2, the latter all `s²` from booleanity.)

**Consequence — the strongest parse-independent statement available:
every equation in the instance is AFFINE in the selector vector.** Each is
`Σ_i s_i·F_i(wires) + G(wires) = 0`. No two selectors ever multiply, anywhere. This is why
P's 0 and S's 48 are both correct: S's coarser atoms bundle several selectors, but each sits
in its own additive term, exactly as S said.

## 4. TASK 2 — classification of every selector-touching atom
2,387 distinct atoms mention a selector. Complete classification (`zclass.py`):

| class | count | shape | meaning |
|---|---|---|---|
| `LOAD_ON` | 512 | `s·(x − K)`, `K` a 287–296-bit constant | when `s=1`, leaf coordinate `x` takes `K` |
| `MIXED_OTHER` | 1,024 | `x_a − s`, `x_b − (1−s)`, and `LOAD_ON` bundled with one neighbouring wire | **alias definitions** — each introduces a *fresh* wire |
| `BOOLEANITY` | 256 | `c·(s − s²)` | `s ∈ {0,1}` |
| `PIN0` | 256 | `s` | the `−s` half of a booleanity/alias atom split at a window boundary |
| `PIN_CONST` | 256 | `1 − s` | the `(1−s)` half of an alias atom split at a window boundary |
| `QUAD_OTHER` | 83 | `s²` | the `s²` half of a split booleanity atom |

Every selector has exactly one `PIN0` and one `PIN_CONST` atom. **These are not pins.**
They are halves of `x_a − s` / `x_b − (1 − s)` that fall across an SLP-window boundary:
e.g. eq 4689 carries the single atom `[−1 + x4805 + s30207]` while eq 4022 carries the same
content split as `[x4805]` and `−[1 − s30207]`. Verified by direct comparison of the two
equations' expanded polynomials. **This is the trap in this measurement, and it is the reason
a naive "is there a `s − 1` atom?" grep says all 256 selectors are pinned when none is.**

**Equations whose variables are ALL selectors: exactly 13** — `7153, 11456, 12810, 13027,
13494, 13905, 16174, 16757, 18154, 19501, 20780, 30248, 33304`. Every one is
`Σ_i c_i·(s_i − s_i²) = 0` with **no constant term and no pure-linear term**: a booleanity
certificate. Verified through `checker.py`'s own compiled evaluator at **every weight
`w = 0 … 256`** with a random `σ` at each: all 13 satisfied at all 257 weights.

More generally **373 equations vanish identically on the boolean locus** (they are pure
booleanity certificates). Checker-verified satisfied at `w ∈ {0,1,2,3,8,17,32,64,128,192,
250,255,256}`: 373/373 every time.

## 5. TASK 3 — THE DECISIVE QUESTION: does anything constrain `|S|`?  **No.**

Method (name it if you cite it): **booleanity-reduced affine elimination.**

1. 3,484 variables carry a booleanity atom (256 selectors + 3,228 wires). Reduce `x²→x` for
   all of them — exact on the boolean locus.
2. After reduction the degree histogram is `{0: 373, 1: 10809, 2: 27851}`. The **10,809
   degree-≤1 equations are exact linear constraints valid on the boolean locus** — this is
   precisely where a cardinality constraint `Σ c_i s_i = const`, a parity constraint, or any
   one-hot/at-most-k condition would have to live.
3. Sparse Markowitz-eliminate all 11,707 non-selector columns from those 10,809 rows.

**Result: 6,829 pivots, 3,980 rows survive, and every surviving row is identically `0 = 0`.**

- **Genuine linear constraints on the selector vector: 0.**
- Inconsistent rows (`c = 0`, `c ≠ 0`): **0** — so the linear layer is satisfiable for
  *every* `σ`, not merely every boolean `σ`.
- Run **twice**: mod `2^61−1` and again in exact `Fraction` arithmetic over ℚ. Identical.
  Not a modular artefact.

Solved symbolically (`zwire_full.py`), the same system gives every determined wire as an
affine function of the selectors:

> **selector-support size of each solved wire: 2,550 wires depend on 0 selectors,
> 149 depend on exactly 1. Zero wires depend on 2 or more.**

**There is no wire anywhere in the linear layer that carries a sum of two selectors.** A
cardinality constraint needs a wire holding `Σ s_i`; no such wire is constructed.

Corroborating census (`zcount.py`): over the 9,527 distinct all-boolean atoms there are
**0 adder-shaped atoms `y = a + b`** and 0 three-variable atoms with exactly one cross
product. Nothing in the instance adds boolean values.

## 6. TASK 4 — the liveness layer
- The subsystem of equations all of whose variables are boolean is **4,763 equations, and
  after `x²→x` it is entirely LINEAR** (`{0: 373, 1: 4390}`) — no quadratic OR gate lives
  inside the pure-boolean layer.
- Solved symbolically: **2,201 of those wires are pinned to constants independent of every
  selector**; none has any selector support. These are P's "127 liveness slots pinned to 0"
  and their relatives.
- The `s → wire` couplings (`x_a = s`, `x_b = 1 − s`, `s·(x − K)`) live in the *mixed* layer,
  which the full elimination of §5 covers. It produces no constraint.
- **No combination of liveness constraints bounds how many leaves are simultaneously live.**
  Liveness composes by OR, which is monotone and saturating; there is no counting network,
  and §5 shows nothing downstream equates a selector sum to anything.

**The one real bound, and it is trivial:** `w ≥ 1`. The empty set leaves the root not live,
so the two target congruences cannot be met (P measured all-selectors-off at 36,388/39,033).
There is **no upper bound and no non-trivial lower bound.**

## 7. TASK 5 — verification by construction
No constraint was found, so there is no violating configuration to exhibit. What I did
instead, through `checker.py`'s own evaluator (not my parser):

| construction | result |
|---|---|
| 13 pure-selector equations, random `σ` at every `w = 0…256` | **satisfied 257/257 weights** |
| 373 identically-vanishing equations, deliverable assignment with selectors overwritten at `w ∈ {0,…,256}` | **373/373 at every weight** |

Both are exactly what "no constraint on `w`" predicts, and both would have failed loudly had
a cardinality or parity condition been hiding in the selector-only equations.

## 8. Clean statement of the negative result — with scope
> **The 256 leaf selectors are unconstrained as a set. The configuration space is the full
> `2^256` (`2^256 − 1` if you exclude the empty set). Established by booleanity-reduced
> affine elimination over the complete instance, run both mod `2^61−1` and exactly over ℚ,
> and corroborated by the parse-independent monomial invariant: not one of the 819,975
> monomials in the instance contains two distinct selector variables.**

**Scope, stated honestly.** The elimination is exact and complete for constraints that are
*linear on the boolean locus* — which is the entire class a cardinality, parity, one-hot or
at-most-`k` constraint belongs to. It does **not** by itself cover a constraint arising
nonlinearly through the 27,851 quadratic equations. Two facts close most of that gap: every
equation is affine in `s` (§3), and no wire in the linear layer depends on more than one
selector (§5), so no `Σ s_i` is ever formed. What remains uncovered by my method is a
constraint mediated by the *group-law* layer — and that is a constraint on *which* subset,
never on *how many*, because P's law is a commutative group operation whose output does not
see the cardinality of its input set.

**Knob set for every claim above: the 256 selector variables, everything else free.
Configuration: arbitrary — every statement in §5–§6 is a statement about all `σ`
simultaneously, not about one configuration.**

## 9. What this means for the campaign
- `MINIMUM_COST_SEARCH.md` §7's "Upper bound: none exists and none is obtainable from the
  instance" is **now measured rather than asserted**. The Hamming-weight angle is closed:
  no shortcut to `w` exists, so MITM bounds must be *assumed*, never derived.
- P's "no genuine cross-selector coupling" and S's "48 atoms, all bundled" are **both
  correct** and now reconciled at the monomial level.
- The `PIN0`/`PIN_CONST` window-boundary artefact (§4) is a live trap for any future agent
  grepping for pinned selectors. 256 false positives.

## 10. Highest-value next experiment
Not another selector-space measurement — that space is now provably featureless. The one
place a genuine constraint on *which* subsets are reachable still lives is the **infeasible
intermediate** noted by P: a merge whose two live inputs share an `x` but differ in `y`
gives `N1 = −B² ≠ 0`, which is unsatisfiable. That is a hard, checkable exclusion on `S`
that no search in this campaign has exploited. Measuring its density (what fraction of
subsets `S` hit such a pair at some merge) is cheap — the leaf coordinates are the 512
extracted constants, all in `zsel.json` — and if the density is non-negligible it prunes the
MITM tree rather than merely bounding it.

## 11. Best verified score
**39,026 / 39,033** — the existing deliverable, re-verified. I did not beat it, and my angle
was never going to: it was a bounding question, and the bound does not exist.

---

# 12. FOLLOW-UP (coordinator check-in 102): density of P's infeasible intermediate

**Verdict: the exclusion is real, exact, and prunes NOTHING. Density `2^-256` — one single
subset, and that subset was already a non-solution. The lead dies, cheaply and completely.**
Files: `zleaf.py`, `zdouble.py`, `zinfeas.py`, outputs `zleaves.json`, `zexpo.json`.

## 12.1 Leaves — measured directly, not assumed (`zleaf.py`)
Reduced the 512 extracted constants mod `p = 2^256 − 2^32 − 977`: **512 distinct residues**,
2 per selector. Orientation into `(x, y)` was **unambiguous, 256/256, 0 ambiguous cases**,
against the relation `y² − (x + Q/3)³ ≡ b`, and the recovered
`b = 64019533680030876408443198762210829058751700634554282185987325820393598524794`
**reproduces P's constant exactly** — an independent confirmation of P §3.

> **Leaf pairs sharing an x-coordinate: 0. Distinct leaf x-coordinates: 256/256.
> Leaves that are ± of another leaf: 0.**

## 12.2 Leaf → exponent map, recovered by doubling (`zdouble.py`)
All 256 leaves verified on `y² = X³ + b`. Doubling each: **255 of 256 doubles land on
another leaf**, exactly one leaf (selector `x_2779`) is nobody's double, and the successor
relation is a **single chain of length 256**. So the leaf set is exactly `{2^i·L0}`,
`i = 0..255`, **established by curve arithmetic rather than taken from the reduction**. Map
saved to `zexpo.json`.

## 12.3 Step 1 — the `±1` question, settled outright
`N − 1 = 2^6 · 3 · 149 · 631 · 107361793816595537 · 174723607534414371449 ·
341948486974166000522343609283189` (product re-verified `== N−1`).

> **`ord_N(2) = 1809251394333065553493296640760748560200586941860545380978205674086221273349`
> (250 bits), and it is ODD.**

- `+1` case needs `ord_N(2) | (i−j)`: impossible, `ord ≈ 2^250 > 255`.
- `−1` case needs `ord_N(2)` even: **it is not**, so `2^d ≡ −1 (mod N)` has **no solution
  for any `d` at all** — not merely none in `|i−j| ≤ 255`.
- Direct check `pow(2,d,N) ∈ {1, N−1}` for `d = 1..255`: **empty**.

Theory and the direct 12.1 measurement agree: **no two leaves can ever share an `x`.**

## 12.4 Step 2 — intermediates: the exact condition, not an estimate
A merge node `v` combines `A = Σ_{i ∈ S∩L_v} 2^i·G` and `B = Σ_{i ∈ S∩R_v} 2^i·G` with
`L_v, R_v` **disjoint** (each leaf lies in one subtree). Then

> `x(A) = x(B)` with `y(A) ≠ y(B)` ⟺ `A = −B` ⟺ **`k(S ∩ T_v) := Σ_{i ∈ S∩T_v} 2^i ≡ 0 (mod N)`**, `T_v = L_v ∪ R_v`.

Disjointness is what makes this exact: the two children's integer values simply **add**.
And because `0 ≤ k < 2^256` while **`2^256 < 2N`** (verified), the congruence has exactly one
nonzero solution:

> **`k(S ∩ T_v) = N` exactly. `popcount(N) = 192`, bit positions spanning 0 … 255.**

So a merge at `v` is infeasible **iff `S ∩ T_v = bits(N)`**, which requires
`bits(N) ⊆ T_v`, hence `|T_v| ≥ 192`.

**This is where the structure decides it.** The coordinator's ~`2/N`-per-pair heuristic is
the right null; the tree structure makes the event *rarer and rigid*, not commoner, because
the partial sums are integers below `2^256 < 2N` so the mod-`N` wrap can happen **at most
once**, collapsing a `~2^-255` random event into a single exact value.

With P's measured root split **178 | 78**, the root (`|T| = 256`) is the **only** node with
`|T_v| ≥ 192`, and there the condition reads `S = bits(N)`:

> **Exactly ONE configuration out of `2^256` is excluded, density `2^-256`. And it is
> already a non-solution: `k = N ⇒ k·G = O ≠ T`. The rule removes nothing the target
> congruence had not removed.**

Worst case over **all** binary trees on 256 leaves (not just this one): nodes with
`|T_v| ≥ 192` form a root path of at most 65 nodes, so the excluded fraction is at most
`65 · 2^-192 < 2^-186`. **Negligible under every tree shape.**

## 12.5 Verified by construction (`zinfeas.py`)
| construction | result |
|---|---|
| `Σ_{i ∈ bits(N)} 2^i·L0` | **= identity** (independently confirms `N·L0 = O`, i.e. `N` is the right modulus) |
| 12 random splits of `bits(N)` into two nonempty disjoint parts | **all 12**: `x(A) = x(B)`, `y(A) = −y(B)`, `y(A) ≠ y(B)` — infeasible, as predicted |
| control: 200 random weight-192 subsets, random split | **0 x-collisions** |

The predicted configuration fires, and only it fires.

## 12.6 Routing (superseded by §13 only in priority, not in content)
**Nothing to route to X, Y or AA.** The infeasible-intermediate exclusion does not prune the
MITM tree — it removes one leaf of it, which was already dead. Combined with §5–§8 and AB's
Theorem B, the instance now has **no measured structure in selector space at all**: no bound
on `w`, no cardinality constraint, and no subset exclusion of any usable density. The
`2^126.5` rho estimate stands unimproved, and every remaining lever is a *prior* on `k`, not
a fact about the equations.

---

# 13. AUDIT of agent AB's corrected cost model and Theorem D (check-in 104)

**Both survive. Neither survives unamended.** Scripts: `zaudit_cost.py` (first pass, which hit
AB's floor bug and is retained as the evidence for it), `zaudit_cost2.py`, `zaudit_thmD.py`.
AB's files were read for the *claims* only; nothing was imported. All combinatorics in exact
integers (`math.comb`), `log2` only at print time.

## 13.1 The corrected covering cost — **re-derived independently, and it is right**
`cost(W) = poly(W) · Vol₁₂₈(⌈W/2⌉)`. Splitting the 256 positions into halves of 128, a fixed
partition does **not** guarantee both parts are ≤ W/2, so a splitting system is needed
(Coppersmith / Stinson) at a `poly(W)` factor; per partition each half enumerates every subset
of size ≤ W/2, i.e. the **cumulative** `Vol₁₂₈`, not a single binomial.

Round-1's `C(256, W/2)` is wrong in two independent ways — wrong ground set (256 not 128) and
single binomial not cumulative — and I can add a third disqualifier AB did not state:
**`C(256,W/2)` is not monotone in `W`** (`C(256,128) = 2^251.7`, `C(256,127) = 2^251.7`, and it
falls thereafter). A ball cost must be monotone. AB's retraction is correct.

## 13.2 Boundary checks in both directions
| check | required | AB's model | verdict |
|---|---|---|---|
| `W = 0` (ball = one point) | cost/ball 1, #balls = whole region ⇒ full enumeration | 1, `Vol(0)=1` | **passes** — but AB never tested it; its scans start at `W = 2` |
| `W = 256` (ball = whole space) | **2^128** | **2^132.0** | **fails by 2^4** |
| monotone in `W` | yes | yes (with ⌈·⌉) | passes |
| never below the generic bound | yes | yes | passes, and I can *prove* it |

**The `W = 256` check is the one AB used to certify its own correction, and the corrected model
does not pass it.** `√W` is spurious exactly there: at `W = 256`, `W/2 = 128` equals the half
size, so every split is automatically balanced and no splitting system is needed. Without the
poly factor the model returns **2^128.0 exactly**. AB printed `2^132.0` and moved on without
saying the test was met only to within the poly factor.

**Never dips below the generic bound — now proved, not spot-checked.**
`(Vol₁₂₈(w/2))² = #{S : |S∩A| ≤ w/2 and |S∩B| ≤ w/2} ≤ Vol₂₅₆(w)` (restricted Vandermonde), so
`Vol₁₂₈(w/2) ≤ √(Vol₂₅₆(w))` **always**; verified exactly for every even `w ∈ [0,256]`. The
corrected MITM therefore sits above the `√(class size)` floor by at most the poly factor and
below it without one — **optimal to within `√W ≤ 2^4`.**

## 13.3 BUG FOUND: floor vs ceiling at odd radius
`ab_costfix.py` computes the half-list as `V128[W//2]`. A weight-`W` set with `W` **odd** cannot
split more evenly than `(⌈W/2⌉, ⌊W/2⌋)`, so the half-list must be `Vol₁₂₈(⌈W/2⌉)`. Using floor
**underprices every odd radius**:

| W | AB (floor) | correct (ceil) | understated by |
|---|---|---|---|
| 9 | 2^25.0 | 2^29.6 | **4.6 bits** |
| 19 | 2^46.4 | 2^49.9 | 3.6 bits |
| 55 | 2^94.9 | 2^96.8 | 1.9 bits |
| 107 | 2^126.4 | 2^126.9 | 0.6 bits |

**AB's published numbers are unaffected**, because every scan in `ab_costfix.py` steps
`range(2,257,2)` — even `W` only, where floor = ceil. **The bug fires the moment anyone reuses
the function on an odd radius.** My first pass did exactly that and produced spurious covering
optima at `W = 9, 15, 35, 55, 75, 103, 127` and a spurious `+1` on every budget row. Recorded
so the next agent does not rediscover it as a "finding".

## 13.4 The corrected budget table — **CONFIRMED exactly**
| budget | AB claims | I get | verdict |
|---|---|---|---|
| 2^47 | `w ≤ 18` | `w ≤ 18` | **CONFIRMED** |
| 2^58 | `w ≤ 24` | `w ≤ 24` | **CONFIRMED** |
| 2^80 | `w ≤ 40` | `w ≤ 40` | **CONFIRMED** |

Filled in: 2^30 → 10, 2^40 → 14, 2^70 → 32, 2^90 → 48, 2^126 → 104. Half-list memory at 2^47 is
2^44.2 entries — **AB's memory caveat is the binding one and it is understated**: on this box
(~2^30 entries) the reachable bound is `w ≤ 10`, not 18.

## 13.5 The crossover — AB quotes two different numbers for one quantity
The "largest affordable complement radius `W`" and the "rho crossover `w`" are **the same
function at the same budget**. AB's RESUME says `w ≈ 104`; AB's script prints `W = 107`.
With the floor bug fixed the reproducible value is:

> **crossover `w = 106`, break-even ceiling `B = 255 − 106 = 149`** (AB says 148).

Immaterial to every conclusion, but it is an internal inconsistency in a headline number.

## 13.6 Theorem B's qualitative claim — **CONFIRMED independently**
Minimising `(|{wt>B}|/Vol(W)) · cost(W)` over all `W ∈ [0,256]` (not AB's even grid), the
optimum is the degenerate **one-ball cover for every `B` below ≈247**: `B = 200 → W = 55`
(2^96.8), `B = 152 → W = 102` (2^125.7), `B = 128 → W = 126` (2^130.5), with `#balls = 2^0`
throughout. **The cheapest way to prove any nontrivial ceiling is to solve the instance.**

---

## 13.7 THEOREM D — the load-bearing barrier

### (1) Is `min(|D₀|,|D₁|)` the right normaliser? **Yes — the constant is not.**
Shoup simulation: `Bad = {k : the real view differs from the k-free simulation}` is fixed by the
simulation and has `|Bad| ≤ C(m,2)`. Then for `D_b` uniform on its class,
`|Pr_{D_b}[A=1] − Pr[A_sim=1]| ≤ |Bad|/|D_b|`, hence

> `Adv ≤ |Bad|·(1/|D₀| + 1/|D₁|) ≤ (m²/2)(1/|D₀| + 1/|D₁|) ≤ **m² / min(|D₀|,|D₁|)**`

**`min` is correct** — it is the dominant term of the harmonic sum. **AB's `m²/(2·min)` drops
the second side and is a factor 2 tighter than the argument supports.** The error direction
**overstates the barrier**, by exactly **0.5 bits** in `m`.

### (2) Does the single-root argument survive the automorphism group? **Yes — the encoding is the caveat.**
- **The affine form is untouched.** `λ` and negation are multiplication by fixed scalars of
  `Z_N`; applying either to `σ(α + βk)` gives `σ(λα + λβk)`, still `σ(α' + β'k)`. Every pairwise
  collision remains **one** affine equation over the field `Z_N`, hence **exactly one root**
  (`N` prime). **AB is right here, and this is the load-bearing step.**
- **The encoding is not.** Under `x`-coordinate encoding `σ(x)` and `σ(−x)` are one string; with
  GLV the whole order-6 orbit `{±1, ±λ, ±λ²}` collapses. A pair then collides if
  `α_i + β_i k = u(α_j + β_j k)` for any `u` in that group — **`AUT` equations per pair, each
  still with one root**, so `|Bad| ≤ AUT·C(m,2)` and the bound degrades by `√AUT`:
  `2^0.5` for negation, **`2^1.29` for the full order-6 group**. AB's knob line
  ("coordinate encoding excluded") is doing real work and should be promoted to the statement.

### (3) Does `B = 128 ⇒ m ≥ 2^127.5` follow? **Yes, from AB's inequality.**
Recomputed with an **exact digit-DP over `k ∈ [0,N)`** (self-checked: `#{k<N} = N`), not over
`[0,2^256)` as `ab_barrier.py` does: `|{w ≤ 128}| = 2^255.07`, `|{w > 128}| = 2^254.93`,
`min = 2^254.93`, `√min = 2^127.46`. AB's arithmetic reproduces. The universe choice is
immaterial (`N` vs `2^256` differ by `2^-128`). The DP also reproduces AB's free unconditional
result: **max popcount over `k < N` is 255.**

### (4) THE ONE SUBSTANTIVE PROBLEM — AB's headline compares across models
AB writes: *"B=128 needs `m ≥ 2^127.5` (solving costs 2^126.5) — deciding the weight predicate
is as hard as solving."* **Taken literally the numbers say deciding is *harder* than solving,
which cannot be true: any algorithm that solves also decides.** The contradiction is not in the
mathematics — it is a **model mismatch**. `2^127.5` is a generic-group bound with the
automorphisms **excluded**; `2^126.5` is a concrete rho cost with them **included**. In one
model:

| quantity | value |
|---|---|
| GGM LB, AB's constant, no automorphisms | 2^127.46 |
| GGM LB, corrected constant, `AUT = 6` | **2^125.67** |
| concrete rho + negation + GLV | 2^126.5 |

**`2^125.7 ≤ 2^126.5` — the inequality points the right way once the models agree.**
AB's qualitative conclusion is **unharmed**: the weight predicate admits **no generic shortcut**
and costs the same as solving to within the same `√6`. **What does not survive is the strict
claim that deciding is harder than solving**, and that phrasing should be struck before the
fleet plans against it.

### (5) Cross-check: does the corrected MITM meet the generic bound?
| B | GGM LB (corrected) | corrected MITM | gap |
|---|---|---|---|
| 10 | 2^28.5 | 2^29.7 | +1.19 |
| 20 | 2^48.5 | 2^50.0 | +1.50 |
| 24 | 2^55.3 | 2^56.8 | +1.59 |
| 40 | 2^77.8 | 2^79.6 | +1.87 |
| 56 | 2^94.7 | 2^96.8 | +2.10 |

**AB's "within 2^1 of optimal at B=20" is CONFIRMED** (1.50 bits). Across the whole actionable
range the corrected MITM is within ~2 bits of the generic floor, so **there is no room left in
the algorithm** — the class size is the whole story, exactly as AB argues.

## 13.8 Verdict, plainly
> **The corrected cost model survives. Its published budget table (2^47 → `w ≤ 18`,
> 2^58 → `w ≤ 24`, 2^80 → `w ≤ 40`) is confirmed exactly, and Theorem B's qualitative claim —
> that the cheapest proof of any nontrivial ceiling is to solve the instance — is confirmed
> independently. Theorem D's normaliser is right and its single-root argument survives the
> automorphism group.**
>
> **Three amendments, none fatal:** the `W = 256` sanity test is met only to within the `√W`
> poly factor (2^132 vs 2^128); `ball_cost` uses floor where it needs ceiling, a latent 4.6-bit
> underprice on odd radii that AB's even-only scans hide; and Theorem D's constant is a factor
> 2 too tight, which with the automorphism group moves `B = 128` from `m ≥ 2^127.5` to
> `m ≥ 2^125.7` and removes the claim that deciding is *harder* than solving.
>
> **Two numbers to fix in the record:** crossover `w = 106` and break-even `B = 149`
> (AB quotes 104 and 107 for the first, 148 for the second).
