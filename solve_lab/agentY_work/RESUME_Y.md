# Agent Y — the complement identity: the first UPPER bound on the Hamming weight of `k`

**Angle.** Every search in this campaign bounds `w = w(k)` from **below** (exhaust low weight,
find nothing, conclude `w > W`). The complement identity turns exactly the same machinery into an
**upper** bound.

    fold(S) + fold(S-bar) = sum_{i=0}^{255} 2^i G = (2^256 - 1) G  =:  A

so if `k·G = T` with ON-set `S`, then `k'·G = T'` with `T' = A − T`, ON-set `S-bar`, and

    k' = (2^256 - 1) - k        (bitwise complement of k as a 256-bit word)
    w(k') = 256 - w(k)

**Exhausting `w' ≤ W` on `T'` with no hit ⟹ `w' ≥ W+1` ⟹ `w ≤ 255 − W`.**

**Deliverable re-verified at start of this thread:**
`solve_lab/best/new_instance_partial_39026.json` → **39,026 / 39,033**, failing
`[12231, 12270, 12350, 14584, 18673, 22044, 29125]`. I produced no better assignment.

---

## 1. `T'` — computed three ways, all agreeing (`ysetup.py` → `ydata.json`)

Curve, ladder and `T` were re-derived from agent Q's raw instance-derived JSON (read-only) and
re-verified here, not inherited:

| check | result |
|---|---|
| `p == 2^256 − 2^32 − 977` (secp256k1 prime) | **True** |
| `N ==` secp256k1 group order | **True** |
| `a = 0`, 256/256 leaves on `y² = x³ + b` | **True** |
| `L_i == 2^i·G` by independent repeated doubling | **256/256, 0 bad** |
| `N·G == O`, `T` on curve, `N·T == O` | **True** |
| `G`, `T` identical to the values Q and X searched | **True** |

`A = (2^256 − 1)·G` computed **three independent ways**:

1. double-and-add with the scalar reduced, `(2^256−1) mod N`
2. **folding all 256 ladder leaves** `L_0 + L_1 + … + L_255`
3. double-and-add with the raw unreduced 256-bit scalar

```
A  = (56062832069265557652411710355921883695174853966112376871108859166282132096323,
      8778476077520406680454125479664240365770597536396101819536168999524180264869)
way1 == way2 == way3 : True     A on curve : True     N·A == O : True
```

**The complement target:**

```
T'.x = 34393883340176920870250405895813312293662019612542515124677151347541157631736
T'.y = 113211777249390963896039371927650543610689519859184504098294945786752299991906
```

| check | result |
|---|---|
| `T'` on curve | **True** |
| `N·T' == O` | **True** |
| `T + T' == A` | **True** |
| `T' ≠ T` | **True** |
| 12 random `S`: `fold(S)+fold(S̄)==A`, `k+k̄ == 2^256−1`, `fold==[k]G` | **12/12** |

## 2. The engine, and why X's table is reusable

The MITM table holds the low 64 bits of `x(Σ_{i∈A} 2^i G)` for every `|A| ∈ [1..4]`.
**It depends only on the ladder, not on the target** — so agent X's `tbl4s.bin` (1.42 GB) and
`bm4.bin` (512 MB) work unchanged for `T'`. I re-validated them independently
(`ytable_check.py`) before reuse:

| check | result |
|---|---|
| file size == `8·(C(256,1)+…+C(256,4))` = 177,589,056 keys | **True** |
| sorted ascending (full scan of all 177.6 M keys) | **True** |
| bitmap is exactly 2^32 bits | **True** |
| 60 random `|A|∈1..4` subsets, key computed here in Python bignum, present in table | **60/60** |
| …and passing the bitmap prefilter | **60/60** |
| 2,000 random 64-bit words present (negative control) | **0/2000** |

Engine: byte-identical copy of X's `xmitm.c` (`md5 1fec0554f59179c92e68eec7971193bd`) as
`ymitm.c`, recompiled. Its 256-bit field arithmetic was re-cross-checked against Python bignum
**for the complement base** (`ycheck.py`): `x(L5+L7)`, `x(T'−L3)`, `y(T'−L3)` all agree
limb-for-limb.

## 3. VALIDATION — planted answers, before any negative

`yplant.py` builds a fake `T'` from a chosen set; `ycheckplant.py` requires **every** exact split
of that set (within the scanned `i0` range) to appear as a HIT.

| planted `S'` | `|S'|` | scan size (α size) | `i0` range | exact splits expected | found |
|---|---|---|---|---|---|
| `{11,73,140,201,244}` | **5** | 2 (α=3) | full `[0,256)` | 10 | **10 — PASS** |
| `{17,64,111,158,199,222,240,253}` | 8 | 4 (α=4) | `[17,18)` | 35 | **35 — PASS** |
| `{6,41,88,129,163,190,214,233,250}` | 9 | 5 (α=4) | `[6,7)` | 70 | **70 — PASS** |
| `{150,170,185,200,215,228,239,247,252}` | 9 | 5 (α=4) | `[150,151)` | 70 | **70 — PASS** |
| `{22,60,101,144,177,203,246}` | 7 | 5 (α=2) | `[22,23)` | 15 | **15 — PASS** |

The required weight-5 plant passed at **full range**, and the expensive size-5 code path passed at
both a low and a high `i0` range and at two different α sizes. **The negatives below are negatives
of a search that has been shown to find what it is looking for.**

## 4. RESULTS on the real `T'`

Coverage: table `|α| ∈ [1..4]` × scan `|β| = b` covers `|S'| = b+1 … b+4`; `|S'| ≤ 4` is closed by
a direct table probe.

| `|S'|` | how covered | result |
|---|---|---|
| 0 | `T'` is an affine point, `≠ O` | no |
| 1 | full point comparison against all 256 `2^i·G` (exact, not truncated) | **0/256 match** |
| 2,3,4 | `key(T')` probed in the `|A|≤4` table (truncation gives false positives only, never false negatives — a miss is exact) | **absent** |
| 3–6 | scan `b=2` | see table below |
| 4–7 | scan `b=3` | " |
| 5–8 | scan `b=4` | " |
| 6–9 | scan `b=5` | " |

| scan size | candidates | `== C(256,b)` | hits | degenerate `dx=0` | time |
|---|---|---|---|---|---|
| 2 | 32,640 | ✔ | **0** | 0 | 0.0 s |
| 3 | 2,763,520 | ✔ | **0** | 0 | 0.3 s |
| 4 | 174,792,640 | ✔ | **0** | 0 | 23.5 s |
| 5 | 8,809,549,056 | ✔ | **0** | 0 | ~20 min |

> **COMPLEMENT WEIGHT `w' ≤ 9` IS EXHAUSTED. No hit.**

## 5. THE BOUND

    w' >= 10   =>   256 - w >= 10   =>   w <= 246

> ### `w < 247`, i.e. `w ≤ 246`.

**This is weak and should not be sold as anything else.** Under a uniform-`k` null
`w ~ Binomial(256, ½)`: mean 128, σ = 8, so `w ≤ 246` excludes essentially nothing —
`P(w > 246)` under the null is about `2^-190`. The bound is not evidence; it is a *certificate*.

**What it is worth is that it exists.** Before this thread the campaign had *no* upper bound of any
kind on `w` and had recorded that none was obtainable. One does exist, it comes from the same
machinery that produced the lower bounds, and it scales with budget at exactly the same rate. Every
weight level anyone ever adds to a forward sweep can be spent on the complement instead, and buys
one unit off the top.

### 5.1 The two-sided bracket

Agent X exhausted forward weight `w ≤ 9` (2^33.1 candidates, `tbl4s.bin` + scan sizes 2–5).
I exhausted complement weight `w' ≤ 9` (2^33.1 candidates, the same table, the same scan sizes,
against `T'`). Together:

> **`10 ≤ w ≤ 246`**, from `2^33.1` forward candidates (agent X, ~30 min on this box) plus
> `2^33.1` complement candidates (agent Y, ~21 min on this box reusing X's table).
> The interval has width 237 out of 257 possible values; the null distribution `Binomial(256,½)`
> puts `1 − 2^-190` of its mass inside it.

Sharpening either end costs the same: level `W → W+1` multiplies the scan by `(256−W+1)/(W+1)`,
about **42×** at `W = 9`. Both ends are 2^58 away from where the null lives.

## 5.2 The complement is the ONLY mask that works — and the family it sits in

Two structural facts, both worth carrying, because they say exactly how far this mechanism goes.

**(a) The complement is the unique nontrivial XOR mask.** For a mask `m`,
`k XOR m = k + m − 2·(k AND m)`. That is affine in `k` — i.e. `(k XOR m)·G` is a *fixed*
translate/reflection of `k·G`, computable without knowing `k` — only when `k AND m` is affine in
`k`, which happens exactly for `m = 0` (gives `0`) and `m = 2^256 − 1` (gives `k`). **So the
2^256 possible masks yield exactly two usable targets: `T` and `T'`.** There is no family of
cheap re-targetings hiding here; the upper bound has one and only one source.

**(b) It *is* the `D = {0..255}` member of a distance family.** For any center `D ⊆ {0..255}`,

    k = fold(D)-scalar + sum_{i in S\D} 2^i - sum_{i in D\S} 2^i

so a **signed-digit** MITM against the shifted target `T − fold(D)` at `m ≤ M` proves
`hamming_distance(S, D) > M`. `D = ∅` is the forward weight bound; `D = {0..255}` is this
thread's complement bound. The machinery costs the same for every center. What is missing is a
*prior* for choosing a center — a random `D` gives distance ≈ 128 just as the null gives
`w ≈ 128`, so only `∅` and the full set are motivated by anything.

**(c) Consequence for agent X's signed-digit sweep — a concrete fix.** Signed-digit search
subsumes the complement case, but **only if the digit alphabet includes exponent 256**:

    near-all-ones k  =  (2^256 - 1) - sum_{j=1..w'} 2^{e_j}  =  2^256 - 2^0 - sum_j 2^{e_j}

which is `w' + 2` signed terms with one exponent equal to **256**. With `e ∈ [0,255]` only —
the natural reading of "the 256 ladder points" — the whole near-all-ones family is **outside** a
signed-digit sweep entirely. Adding one extra digit point, `2^256·G = 2·L_255`, makes signed-digit
`m ≤ 11` cover complement weight `≤ 9`. Note the direct complement search is still cheaper for that
job: `2^33.1` here versus roughly `2^38`–`2^41` for signed-digit `m ≤ 11`.

## 6. Endomorphism orbit — derived, verified, edge-probed

`ylam.py` derives `beta` (cube root of 1 mod `p`), `lambda` (cube root of 1 mod `N`), and verifies
`phi(x,y) = (beta·x, y) == [lambda]` on **8/8 random points**:

```
beta   = 55594575648329892869085402983802832744385952214688224221778511981742606582254
lambda = 37718080363155996902926221483475020450927657555482586988616620542887997980018
```

Twelve targets, all distinct, all on the curve, all of order `N`: `±T, ±phi(T), ±phi²(T)`
(scalars `±k, ±λk, ±λ²k`) and the complement `A − X` of each. A hit on any recovers `k` by a known
automorphism. Note `c_T == T'`, so this thread's sweep is one of the twelve and X's is another.

**Already done at zero cost** (`yorbit_edge.py`): the `|S| ≤ 4` table probe on **all 12** — no hit
on any, and no target equals any `2^i·G` (full point comparison, 0/256 each).

`ylam.py` derives `beta` (cube root of 1 mod `p`), `lambda` (cube root of 1 mod `N`), verifies
`phi(x,y) = (beta·x, y) == [lambda]` on random points, and emits engine inputs for **12 targets**:
`±T, ±phi(T), ±phi²(T)` (scalars `±k, ±λk, ±λ²k`) and the complement `A − X` of each.
A hit on any one recovers `k` by a known automorphism.

**Coverage per unit cost, computed:**

| option | cost multiple | scalars covered multiple | scalars per unit cost |
|---|---|---|---|
| one more weight level, `w' ≤ 9 → ≤ 10` | `C(256,6)/C(256,5)` = **41.8×** | `C(256,10)/C(256,9)` = **24.7×** | **0.59** |
| the 12-target orbit at `w ≤ 9` | **12×** | **12×** | **1.00** |

Raw scalar count favours the orbit by **1.7×**. That comparison is honest only about arithmetic:
under the one hypothesis anyone actually holds ("the designer picked a low-weight `k`"), ten of the
twelve targets carry ~zero prior, so the orbit is a hedge against a *different* designer choice
(low-weight `λk`, or low-weight `N−k`), not more evidence about the one at hand. The deeper level
is worth less per candidate but is drawn from the hypothesis class that matters.

## 7. Files

| file | what |
|---|---|
| `ysetup.py` → `ydata.json` | re-derives + verifies curve, ladder, `T`, `A`, `T'` |
| `ytable_check.py` | independent validation of X's reused table + bitmap |
| `ygen.py` → `data_comp.txt`, `data_fwd.txt` | engine inputs (complement / forward control) |
| `ymitm.c` / `ymitm` | byte-identical copy of X's engine, recompiled |
| `ycheck.py` | C-vs-Python arithmetic cross-check on the complement base |
| `yedge.py` | `|S'| = 0, 1` and the `|S'| ≤ 4` table probe |
| `yplant.py`, `ycheckplant.py` | planted-answer validation |
| `yrun.sh`, `yrun.pid`, `rep_comp.txt`, `yrun_*.log` | the real sweep |
| `ylam.py` → `yorbit.json`, `data_*.txt` | the 12 endomorphism-orbit targets |

**Restart:** `python3 ysetup.py; python3 ygen.py; gcc -O3 -march=native -fopenmp -o ymitm ymitm.c`,
then `./yrun.sh`. If `solve_lab/agentX_work/tbl4s.bin` is gone, rebuild it with
`./ymitm table data_comp.txt 4 tbl4.bin`, `python3 ../agentX_work/xsort.py`,
`./ymitm bitmap data_comp.txt tbl4s.bin bm4.bin` (the table ignores the target, so any data file
works). Nothing here is a `*.pkl`, so a container restart does not wipe it; the 1.9 GB of `.bin`
lives in `agentX_work/` and is git-ignored.
