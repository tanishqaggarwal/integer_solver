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
