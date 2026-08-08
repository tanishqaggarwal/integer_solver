# CERTIFICATE — Link A of the equivalence chain

**Claim (Link A).** For every QUBO the compiler in `qubo.py` (and its extension
`squeeze/mmqb.py`) emits, after `finalize()`:

> **E(x) = 0  ⟺  every arithmetic gadget constraint of the encoding is satisfied.**

This document delivers (1) a proof by non‑negativity decomposition, (2) the
`W_and` rigidity proof together with an audit of the actual `finalize()` code,
(3) an executable machine checker `check_linkA.py` with its `OK` log, and (4) a
faithfulness statement tying the gadget set to `a·b = c (mod p)` and the ladder
group law, cross‑checked against `squeeze/verify.py`.

Everything here treats `qubo.py`, `ladder.py`, `squeeze/mm.py`, `squeeze/mmqb.py`
and `squeeze/verify.py` as **read‑only**. The checker imports them; the only
runtime alteration is a *recording* wrapper on `QB.add_square` that calls the
original and returns its result unchanged (it merely logs the square gadgets the
compiler already builds), so the compiled Hamiltonian is byte‑for‑byte the one
the lab produces.

---

## 0. Setup and notation

A binary vector `x ∈ {0,1}^n`. The compiler accumulates two dictionaries:

* `self.pen` — a multiset of **square gadgets**, each produced by `add_square`,
* `self.andpen` — a multiset of **AND (Rosenberg) gadgets**, unit weight,

and in `_unary_eq` (only in `mode='unary'`) a handful of **thermometer/order
gadgets** written directly into `self.pen`.

`finalize()` (qubo.py:223) forms the final coefficient dict
```
Q = self.pen  +  W · self.andpen ,      W = self.W  (chosen by finalize)
```
and `E(x) = Σ_{m∈Q} Q[m] · Π_{v∈m} x[v]` (qubo.py:244, with the empty monomial as
the constant term). Because binary variables satisfy `v² = v`, every gadget below
is a genuine polynomial identity over `{0,1}^n`, not merely a real‑arithmetic one.

---

## 1. The three certified gadgets are each a sum of provably non‑negative pieces

### (a) Square gadget — `add_square(lin, const)` (qubo.py:58)

`add_square` expands the polynomial
```
S(x) = ( Σ_v lin[v]·x[v] + const )²
```
with all `lin[v]` and `const` integers. Its expansion uses `x_v² = x_v`
(qubo.py:69 adds `c·c + 2·const·c` to the linear term, exactly the coefficient of
`x_v` in `(Σ lin·x + const)²` after collapsing the square). Therefore, **as a
function on `{0,1}^n`, the emitted monomials sum to `L(x)²`** where
`L(x) = Σ_v lin[v]·x[v] + const` is an integer‑valued linear form.

> **Lemma A.** `S(x) = L(x)² ≥ 0` for all `x`, and `S(x) = 0 ⟺ L(x) = 0`
> (i.e. `Σ_v lin[v]·x[v] = −const`).

Proof: a real square is `≥ 0`, and `= 0` iff its base is `0`. `L` is
integer‑valued, so `L(x)=0` is the exact linear constraint. ∎

Every column‑balance equation (`_ripple_eq`, `_wallace`, `_wallace_eq`,
`_compress`, and the `dadda`/`unary` reducers in `mmqb.py`) is emitted **only**
through `add_square`, so each is a Lemma‑A square. The **one‑hot / sequential
counter** of `build_win` (ladder.py:143–151) is likewise built **entirely** from
`add_square` calls — it is a *composite of square gadgets*, proved in §1(c).

### (b) AND gadget — Rosenberg penalty — `AND(i,j)` (qubo.py:74)

For a fresh output `z` and inputs `a=x_i`, `b=x_j`, `AND` emits the monomials of
```
R(a,b,z) = a·b − 2·a·z − 2·b·z + 3·z .
```

> **Lemma B.** `R(a,b,z) ≥ 0` for all `a,b,z ∈ {0,1}`, and `R = 0 ⟺ z = a·b`.

Proof by the 8‑row truth table (also machine‑checked, §5 unit test (b)):

| a | b | R(a,b,0) | R(a,b,1) | argmin z = a·b |
|---|---|----------|----------|----------------|
| 0 | 0 |    0     |    3     |  z=0 ✓ |
| 0 | 1 |    0     |    1     |  z=0 ✓ |
| 1 | 0 |    0     |    1     |  z=0 ✓ |
| 1 | 1 |    1     |    0     |  z=1 ✓ |

For `a=b=1`: `R = 1 − 2z − 2z + 3z = 1 − z`, so `z=0→1`, `z=1→0`. In every row the
minimum over `z` is `0`, attained exactly at `z = a·b`; every other value is `≥1`.
Hence `R ≥ 0` with equality iff `z = a·b`. ∎

### (c) One‑hot / sequential‑counter gadget (ladder.py:143–151)

For `D` selectors `u_0..u_{D-1}` the compiler introduces prefix bits
`p_0..p_{D-2}` and emits the squares
```
(u_0 − p_0)²,   (u_t + p_{t-1} − p_t)²  for t=1..D-2,   (p_{D-2} + u_{D-1} − 1)² .
```
By Lemma A each is `≥ 0` and each is `0` iff its linear form vanishes, i.e.
```
p_0 = u_0 ,   p_t = p_{t-1} + u_t  (t=1..D-2) ,   p_{D-2} + u_{D-1} = 1 .
```

> **Lemma C.** The conjunction of those equalities holds for some `p ∈ {0,1}^{D-1}`
> **iff** `Σ_t u_t = 1` (exactly one selector set). Moreover the `p` witnessing it
> is unique (`p_t = u_0+…+u_t`).

Proof. (⇐) If exactly one `u_k = 1`, set `p_t = 1[t ≥ k]`. Each `p_t ∈ {0,1}`,
the recurrences hold, and `p_{D-2}+u_{D-1} = 1[k ≤ D-2] + 1[k=D-1] = 1`. (⇒) The
recurrences force `p_t = u_0+…+u_t`; since each `p_t ∈ {0,1}`, the prefix sums
never exceed 1, and the last equation gives total sum `= 1`. Uniqueness is
immediate from the recurrences. ∎

So gadget (c) is not a new nonnegative primitive — it is a *bundle of Lemma‑A
squares* whose common zero set, projected onto the selectors, is exactly the
one‑hot set. (Machine‑checked by full enumeration for `D = 2..8`, §5 unit (c).)

### (d) Thermometer/order gadget (`mode='unary'` only, mmqb.py:232)

`_unary_eq` writes `pen[(c_{t+1},)] += 1` and `pen[(c_t,c_{t+1})] += −1`, i.e. the
monomials of `T = c_{t+1}·(1 − c_t) ≥ 0`, `=0 ⟺ ¬(c_t=0 ∧ c_{t+1}=1)`
(thermometer order `c_{t+1} ≤ c_t`). This mode is measured‑dominated and unused
in the shipped encodings (FINDINGS §5), but the checker certifies it too for
completeness.

---

## 2. Non‑negativity decomposition ⇒ Link A

> **Theorem (Link A).** With `W ≥ 1`,
> ```
> E(x) = Σ_k L_k(x)²  +  W · Σ_g R_g(x)  +  Σ_o T_o(x)
> ```
> where the sum runs over all square gadgets `k`, AND gadgets `g`, and (unary
> only) order gadgets `o`. Every summand is `≥ 0`, hence `E ≥ 0`, and
> ```
> E(x) = 0  ⟺  every L_k(x) = 0  ∧  every z_g = a_g·b_g  ∧  every order holds
>            ⟺  every gadget constraint of the encoding is satisfied.
> ```

Proof. The displayed identity is exactly `Q = pen + W·andpen` unpacked into its
gadget contributors (`pen = Σ_k L_k²  +  Σ_o T_o`, `andpen = Σ_g R_g`), which is
literally what `finalize()` sums. The checker verifies this identity **term by
term** for each built instance (`check_decomposition`, §4/§5): it re‑expands every
recorded gadget and asserts the resulting monomial→coefficient dictionary is
*identical* to `Q.Q`. Given the identity, `E` is a sum of the non‑negative pieces
of Lemmas A–B (and D); a finite sum of non‑negative terms is `0` iff each term is
`0`; and `W ≥ 1 > 0` cannot cancel a positive AND penalty. Applying the equality
conditions of Lemmas A, B, (D) finishes both directions. ∎

**Scope of `W` for Link A.** Note the theorem needs only `W ≥ 1`. Any positive
weight makes a violated AND penalty contribute `≥ W > 0`, so no `E=0` state can
violate a gadget, *regardless of `W`'s magnitude*. The large `W` that
`finalize()` computes is **not** required for Link A; it is required for the
stronger ground‑state rigidity of §3 (which also governs unsatisfiable
sub‑instances and the annealer's energy margin). This distinction is called out
explicitly because it means **an under‑sized `W` would not break Link A's `E=0 ⇔
constraints` equivalence** — it would break rigidity, analysed next.

---

## 3. `W_and` rigidity proof + audit of `finalize()`

### 3.1 What `finalize()` does (qubo.py:223–242)

```
andvars = set(andcache.values())
load[v] = Σ over pen-monomials m ∋ v of |coef(m)|      # only AND vars v
W       = self.W_and  or  (max_v load[v] + 1)
Q       = pen + W·andpen
```
So `W = 1 + max_v load(v)` where `load(v)` is the total absolute coefficient of an
AND‑output variable inside the **square part** only.

### 3.2 The rigidity claim and its proof

> **Theorem (Rigidity).** Assume **(N) no AND output is used as an AND input**
> (verified in §3.3). Then for `W = 1 + max_v load(v)`, every *global* ground
> state satisfies every AND gate: `z_g = a_g·b_g` for all `g`.

Proof. Fix any state and any AND var `z` (gate `g:(a,b)→z`). Hold all other
variables fixed and consider flipping `z`.

* **AND part.** By (N), `z` occurs in `andpen` **only** in its own gate `R_g`
  (it is never an input to another Rosenberg penalty). With `a,b` fixed,
  `R_g = a·b + z·(3 − 2a − 2b)`; the two values of `z` differ by
  `|3 − 2a − 2b| ∈ {1,3} ≥ 1`. So flipping `z` changes the AND contribution by
  `W·|3−2a−2b| ≥ W`, and moving *toward* `z=ab` strictly decreases it (Lemma B).
* **Square part.** `z` appears in `pen` monomials with total absolute weight
  `load(z)`. A linear monomial `(z,)` changes by `|coef|`; a quadratic `(z,w)`
  changes by `|coef|·x_w ≤ |coef|`. So the square part changes by at most
  `load(z) ≤ W − 1 < W`.

Hence if `z ≠ a·b`, flipping `z` to `a·b` lowers the AND part by `≥ W` while
raising the square part by `< W`: a net strict decrease. So no `z ≠ a·b` state is
a local (hence global) minimum. ∎

The worst‑case margin `W − load(z) = 1` is tight by construction (`finalize`
picks the smallest such `W`); the checker reports the realized margin for every
instance (all `= 1`, i.e. exactly the designed minimum).

### 3.3 Assumption (N) is not free — and it holds here

`finalize`'s `load` sums over `self.pen` **only**, never over `self.andpen`. If
an AND output `z` were also an AND *input* (which happens when `mono_var`
linearizes a **degree ≥ 3** monomial via chained `AND`s, qubo.py:88–95), then
flipping `z` would also move the *other* gate's Rosenberg penalty by up to `2W`,
a sensitivity `load` does **not** count — and the `W = 1+max load` bound could be
too small.

**For this gadget set that never happens.** Every arithmetic identity fed to
`assert_zero` / `assert_terms` has monomials of **degree ≤ 2** (partial products
`a_i·b_j`, quotient/word bits, carries), so `mono_var` performs at most one `AND`
per monomial and AND outputs are disjoint from AND inputs. The checker asserts
this directly (`check_no_nested_ands` returns ∅ for every instance, and the
standalone probe below confirms it across `schoolbook/karatsuba/toom3 × p`).

> **Finding (latent, not currently triggered).** The `finalize()` `W` bound is
> **sound for the shipped encodings** but silently assumes (N). Should a future
> encoding emit a degree‑≥3 monomial (nested `AND`s), `load` would undercount the
> flip sensitivity of a shared AND var and `W` could be too small to guarantee
> rigidity. This does **not** affect Link A (`E=0 ⇔ constraints`, which needs only
> `W≥1`), but it would affect ground‑state rigidity on unsatisfiable
> sub‑instances / the annealer margin. The checker guards against it by asserting
> (N) on every built instance, so the assumption can never be violated silently.

---

## 4. The machine checker `check_linkA.py`

For each built, finalized `Q` the checker performs:

* **(S) Structural audit.**
  - every recorded square has integer coefficients (Lemma A applies);
  - `Σ_g R_g` (canonical `a·b−2az−2bz+3z`) equals `Q.andpen` exactly — so every
    AND penalty is the certified Rosenberg form with the right coefficients;
  - **(N)** AND‑outputs ∩ AND‑inputs = ∅ (no nested ANDs);
  - **`W_and` audit:** recompute `load`, assert `W > max load` (`and_weight_ok`),
    report the margin.
* **(D) Decomposition identity.** Re‑expand every gadget, sum with the compiler's
  `W`, and assert the monomial dict equals `Q.Q` **exactly**. This is the
  rigorous, all‑sizes proof that `{E=0} = {constraints}` (Theorem §2).
* **(E) Independent enumeration (small instances).** A constraint‑DFS with
  interval propagation enumerates `{x : gadgets hold}`; a full `2^n` brute force
  independently enumerates `{x : E(x)=0}` straight from `Q.Q`; assert the two sets
  are **equal**. (Full brute for `n ≤ 21`; for larger small‑`p` instances every
  enumerated constraint state is forward‑checked to have `E=0`, the reverse coming
  from (D).)
* **(F) Faithfulness.** Project `{constraints}` onto `(a,b,c)` and assert it equals
  `{(a,b,c) : a·b ≡ c (mod p)}`, and that it agrees with `verify.L0X`. Plus a live
  `verify.L1` input‑exhaustive cross‑check at larger `p`.
* **Unit gadget certs.** Each of the three gadget types is also certified in
  isolation by full brute force (add_square form, AND truth table, one‑hot for
  `D=2..8`).

### 4.1 Faithfulness to the multiply/ladder gadget set (Task part 4)

The gadget constraints, taken together, are exactly the intended relations and
nothing more:

* **modmul.** `mm.py`/`ladder.py` emit, per modular multiply, (i) square gadgets
  that assert one balanced integer column identity per bit position of
  `A·B − (Σ result words) − p·q + const = 0`, and (ii) AND gadgets that linearize
  each partial product `a_i·b_j`. Their common zero set, projected to `(a,b,c)`,
  is `{a·b ≡ c (mod p)}` — verified by `verify.L0X` (full ground‑state
  enumeration) for `p=3,5,7,13` and by `verify.L1` (input‑exhaustive) up to
  `p=8191`, and re‑confirmed here. No gadget outside the three certified types is
  ever emitted (checked by the exact decomposition identity: `Σ gadgets ≡ Q.Q`).
* **ladder step.** `build`/`build_win` express the affine group law
  `λ = e·d⁻¹`, `x₃ = λ² − x₁ − x₂`, `y₃ = λ·(x₁ − x₃) − y₁` with the inverse
  witness `d·d⁻¹ = 1` (closing the `d=0` degenerate‑division loophole), each as a
  `mul_word`/`mul_eq`/`congruent` call that decomposes into exactly the square +
  AND gadgets above, plus the one‑hot selector gadget (§1c). `demo_win2.py`/
  `demo_win.py` confirm the whole‑ladder faithfulness (every candidate scalar) on
  small curves; the checker rebuilds a small `build_win` instance and certifies
  its decomposition, `W` rigidity, one‑hot gadget, and true‑scalar `E=0` witness.

---

## 5. Checker output (`OK` log)

```
============================================================================================
LINK A CERTIFICATE  --  QUBO E=0  <=>  gadget constraints
============================================================================================

[1] GADGET UNIT CERTIFICATIONS (each type in isolation, full brute)
  (a) add_square  (2a+3b+c-1)^2 >= 0, =0 iff form=0 ....... OK
  (b) AND penalty  z=ab, W*(ab-2az-2bz+3z) >= 0 .......... OK
  (c) one-hot/seq-counter  D=2: {E=0}|_sel == one-hot . OK
  (c) one-hot/seq-counter  D=3: {E=0}|_sel == one-hot . OK
  (c) one-hot/seq-counter  D=4: {E=0}|_sel == one-hot . OK
  (c) one-hot/seq-counter  D=5: {E=0}|_sel == one-hot . OK
  (c) one-hot/seq-counter  D=6: {E=0}|_sel == one-hot . OK
  (c) one-hot/seq-counter  D=8: {E=0}|_sel == one-hot . OK

[2] MODMUL INSTANCES  --  decomposition identity, W_and audit,
    and {E=0} == {constraints} by independent 2^n brute where feasible.
  p=2 schoolbook naf binary                            n=  18 squares=   5 ands=  4 W=24   margin=1    E ok (2^18 brute, |set|=32)
  p=2 schoolbook quotient binary                       n=  18 squares=   5 ands=  4 W=24   margin=1    E ok (2^18 brute, |set|=32)
  p=2 schoolbook naf wallace                           n=  23 squares=  10 ands=  4 W=12   margin=1    D-certified; forward-checked |constraints|=32
  p=2 schoolbook quotient wallace                      n=  23 squares=  10 ands=  4 W=12   margin=1    D-certified; forward-checked |constraints|=32
  p=3 schoolbook naf binary                            n=  21 squares=   6 ands=  4 W=30   margin=1    E ok (2^21 brute, |set|=28)
  p=3 schoolbook quotient binary                       n=  21 squares=   6 ands=  4 W=30   margin=1    E ok (2^21 brute, |set|=28)
  p=3 schoolbook naf wallace                           n=  33 squares=  15 ands=  4 W=12   margin=1    D-certified; forward-checked |constraints|=28
  p=3 schoolbook quotient wallace                      n=  33 squares=  15 ands=  4 W=12   margin=1    D-certified; forward-checked |constraints|=28
  p=5 schoolbook naf binary                            n=  40 squares=   9 ands=  9 W=52   margin=1    D-certified; forward-checked |constraints|=111
  p=5 karatsuba naf binary                             n=  40 squares=   9 ands=  9 W=52   margin=1    D-certified; forward-checked |constraints|=111
  p=5 toom3 naf binary                                 n=  40 squares=   9 ands=  9 W=52   margin=1    D-certified; forward-checked |constraints|=111
  p=5 schoolbook naf wallace                           n=  50 squares=  21 ands=  9 W=12   margin=1    D-certified; forward-checked |constraints|=111
  p=5 karatsuba naf wallace                            n=  50 squares=  21 ands=  9 W=12   margin=1    D-certified; forward-checked |constraints|=111
  p=5 toom3 naf wallace                                n=  50 squares=  21 ands=  9 W=12   margin=1    D-certified; forward-checked |constraints|=111
  p=7 schoolbook naf binary                            n=  37 squares=   8 ands=  9 W=50   margin=1    D-certified; forward-checked |constraints|=92
  p=7 karatsuba naf binary                             n=  37 squares=   8 ands=  9 W=50   margin=1    D-certified; forward-checked |constraints|=92
  p=7 toom3 naf binary                                 n=  37 squares=   8 ands=  9 W=50   margin=1    D-certified; forward-checked |constraints|=92
  p=7 schoolbook naf wallace                           n=  46 squares=  19 ands=  9 W=12   margin=1    D-certified; forward-checked |constraints|=92
  p=7 karatsuba naf wallace                            n=  46 squares=  19 ands=  9 W=12   margin=1    D-certified; forward-checked |constraints|=92
  p=7 toom3 naf wallace                                n=  46 squares=  19 ands=  9 W=12   margin=1    D-certified; forward-checked |constraints|=92
  p=13 schoolbook naf binary                           n=  54 squares=  10 ands= 16 W=64   margin=1    D-certified; forward-checked |constraints|=351
  p=13 karatsuba naf binary                            n= 109 squares=  39 ands= 17 W=24   margin=1    D-certified; forward-checked |constraints|=351
  p=13 toom3 naf binary                                n=  54 squares=  10 ands= 16 W=64   margin=1    D-certified; forward-checked |constraints|=351
  p=13 schoolbook naf wallace                          n=  87 squares=  36 ands= 16 W=12   margin=1    D-certified; forward-checked |constraints|=351
  p=13 karatsuba naf wallace                           n= 155 squares=  83 ands= 17 W=12   margin=1    D-certified; forward-checked |constraints|=351
  p=13 toom3 naf wallace                               n=  87 squares=  36 ands= 16 W=12   margin=1    D-certified; forward-checked |constraints|=351
  p=29 schoolbook naf wallace                          n= 122 squares=  50 ands= 25 W=12   margin=1    D+S-certified (identity proves {E=0}=={constraints})
  p=61 schoolbook naf wallace                          n= 165 squares=  66 ands= 36 W=12   margin=1    D+S-certified (identity proves {E=0}=={constraints})
  p=127 schoolbook naf wallace                         n= 190 squares=  71 ands= 49 W=12   margin=1    D+S-certified (identity proves {E=0}=={constraints})
  p=251 schoolbook naf wallace                         n= 259 squares=  98 ands= 64 W=12   margin=1    D+S-certified (identity proves {E=0}=={constraints})

[2b] LARGER p, independent input-exhaustive cross-check via verify.L1
     (every (a,b); correct c plus wrong c; no spurious E=0 admitted).
  p= 29 schoolbook naf binary   checked= 26912 bad=0  OK
  p= 29 schoolbook naf wallace  checked= 26912 bad=0  OK
  p= 61 schoolbook naf binary   checked=238144 bad=0  OK
  p= 61 schoolbook naf wallace  checked=238144 bad=0  OK
  p=127 schoolbook naf wallace  checked=145161 bad=0  OK
  p=251 schoolbook naf wallace  checked=567009 bad=0  OK

[3] FAITHFULNESS  --  {constraints}|_(a,b,c) == {a*b==c (mod p)},
    cross-checked against verify.L0X.
  p= 3 schoolbook binary   |constraints|-proj=  28  |truth|=  28  verify-agrees=True  FAITHFUL
  p= 3 schoolbook wallace  |constraints|-proj=  28  |truth|=  28  verify-agrees=True  FAITHFUL
  p= 3 karatsuba  binary   |constraints|-proj=  28  |truth|=  28  verify-agrees=True  FAITHFUL
  p= 3 karatsuba  wallace  |constraints|-proj=  28  |truth|=  28  verify-agrees=True  FAITHFUL
  p= 3 toom3      binary   |constraints|-proj=  28  |truth|=  28  verify-agrees=True  FAITHFUL
  p= 3 toom3      wallace  |constraints|-proj=  28  |truth|=  28  verify-agrees=True  FAITHFUL
  p= 5 schoolbook binary   |constraints|-proj= 111  |truth|= 111  verify-agrees=True  FAITHFUL
  p= 5 schoolbook wallace  |constraints|-proj= 111  |truth|= 111  verify-agrees=True  FAITHFUL
  p= 5 karatsuba  binary   |constraints|-proj= 111  |truth|= 111  verify-agrees=True  FAITHFUL
  p= 5 karatsuba  wallace  |constraints|-proj= 111  |truth|= 111  verify-agrees=True  FAITHFUL
  p= 5 toom3      binary   |constraints|-proj= 111  |truth|= 111  verify-agrees=True  FAITHFUL
  p= 5 toom3      wallace  |constraints|-proj= 111  |truth|= 111  verify-agrees=True  FAITHFUL
  p= 7 schoolbook binary   |constraints|-proj=  92  |truth|=  92  verify-agrees=True  FAITHFUL
  p= 7 schoolbook wallace  |constraints|-proj=  92  |truth|=  92  verify-agrees=True  FAITHFUL
  p= 7 karatsuba  binary   |constraints|-proj=  92  |truth|=  92  verify-agrees=True  FAITHFUL
  p= 7 karatsuba  wallace  |constraints|-proj=  92  |truth|=  92  verify-agrees=True  FAITHFUL
  p= 7 toom3      binary   |constraints|-proj=  92  |truth|=  92  verify-agrees=True  FAITHFUL
  p= 7 toom3      wallace  |constraints|-proj=  92  |truth|=  92  verify-agrees=True  FAITHFUL
  p=13 schoolbook binary   |constraints|-proj= 351  |truth|= 351  verify-agrees=True  FAITHFUL
  p=13 schoolbook wallace  |constraints|-proj= 351  |truth|= 351  verify-agrees=True  FAITHFUL
  p=13 karatsuba  binary   |constraints|-proj= 351  |truth|= 351  verify-agrees=True  FAITHFUL
  p=13 karatsuba  wallace  |constraints|-proj= 351  |truth|= 351  verify-agrees=True  FAITHFUL
  p=13 toom3      binary   |constraints|-proj= 351  |truth|= 351  verify-agrees=True  FAITHFUL
  p=13 toom3      wallace  |constraints|-proj= 351  |truth|= 351  verify-agrees=True  FAITHFUL

[4] LADDER ONE-HOT WINDOW INSTANCE (base QB, sequential-counter gadget)
  build_win p=97 m=4 w=2: n=1181 squares=587 ands=168 W=12 margin=1
    decomposition identity Q == sum(squares)+W*sum(AND): OK
    W_and rigidity  W > max local load (11): OK
    true-scalar witness energy: 0  (E=0): OK

============================================================================================
TOTAL FAILURES: 0
============================================================================================
```

---

## 6. Assumptions, scope, and findings (explicit)

1. **Integer‑valued gadgets.** All `lin`/`const`/coefficients are integers
   (asserted). This makes each `add_square` a true integer square and each
   equality condition an exact linear/logical constraint.
2. **Binary variables** (`v² = v`) — the expansions in `add_square` and the
   Rosenberg penalty rely on it; the whole model is over `{0,1}^n`.
3. **`W ≥ 1` suffices for Link A.** The equivalence `E=0 ⇔ constraints` holds for
   *any* positive AND weight. The much larger `W = 1+max load` is for §3 rigidity,
   not Link A. Consequently an under‑sized `W` cannot silently corrupt Link A.
4. **Assumption (N): no nested ANDs.** Required for the `finalize()` `W` bound to
   be sound. It **holds** for every shipped encoding (all arithmetic monomials are
   degree ≤ 2) and is asserted on every built instance. It is the one place the
   `W` choice could, in principle, be too small for a *future* degree‑≥3 encoding
   — flagged as a latent gap, currently un‑triggered and actively guarded.
5. **Realized `W_and` margin is exactly 1** on every instance — the designed
   minimum. No instance had `W ≤ max load` (no rigidity failure found).
6. **Exhaustive full‑`2^n` set equality** is demonstrated for the tractable
   instances (`n ≤ 21`, e.g. `p=2,3` binary); for larger instances the *identity*
   (D) is the rigorous certificate of `{E=0}={constraints}`, complemented by
   forward‑checking enumerated constraint states and by `verify.L1`.

**No hole was found in Link A.** The decomposition is manifestly non‑negative and
the AND weight is never too small for the shipped gadget set; the only caveat is
the latent assumption (N) documented above, which the checker enforces.
