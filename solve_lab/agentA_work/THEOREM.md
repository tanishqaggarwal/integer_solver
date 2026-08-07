# A lower bound on the number of unsatisfiable equations in `EQUATIONS.txt`

Agent A, self-contained. Everything below is exact integer or finite-field arithmetic;
no floating point is used anywhere, including in the verification of the deliverable.

---

## 0. The object

`EQUATIONS.txt` is a list of **39,033** polynomial equations over **38,748** integer
unknowns `x_0 … x_38747`. Parsing (`solve_lab/s9/atomize.py`, `poly.py`, `gates.py`,
`fwd.py`) factors every equation into a linear combination of shared subexpressions
called **atoms**: there are **42,267** atoms, each a polynomial of degree ≤ 4 in the
unknowns, and every equation has the form

    eq_e(x)  =  m_e * ( Σ_a c_{e,a} · β_a(x) )        or       m_e * ( Σ_a c_{e,a} · β_a(x) )²

with small integer coefficients `c_{e,a}` and an outer multiplier `m_e`. Since `m_e ≠ 0`,
equation `e` holds **iff** the inner linear form vanishes. Write `β(x) ∈ ℤ^42267` for the
vector of atom values.

The best assignment known to this campaign, `solve_lab/best/new_instance_partial_39026.json`
(copy: `agentA_work/A_best_39026.json`), satisfies **39,026 of 39,033** equations. Verify:

    python3 solve_lab/checker.py solve_lab/agentA_work/A_best_39026.json
    -> satisfied 39026/39033 (7 failing)
    -> failing line indices [12231, 12270, 12350, 14584, 18673, 22044, 29125]

Call this assignment `S`. This document proves that, within an explicitly delimited
neighbourhood of `S`, **7 is optimal** — and explains *why* it is 7.

---

## 1. Windows, knobs, and the code

### 1.1 Equation-closure windows

Fix a level `L ≥ 0`.

* `A_0` := the set of atoms occurring in the 7 equations that `S` fails.
* `A_{k+1}` := the set of atoms occurring in **any** equation that contains an atom of `A_k`.
* `R_L` := the set of equations containing at least one atom of `A_L`.
* `K_L` := the set of variables `u` such that **every** atom containing `u` lies in `A_L`.

The step `A_{k+1} = atoms(eqs(A_k))` is the essential one: it promotes every *foreign*
atom of every modelled equation into the model, where it is free to take any reachable
value and to **cancel** against the others. Nothing in this construction requires any atom
to be zero.

Measured (`eqwin.py`, `eqwin2.py`, `boundary.py`):

| L | atoms | equations | nontrivial rows | knobs | rank(N) | ℚ-inconsistent | violated by S |
|---|---|---|---|---|---|---|---|
| 0 | 24 | 27 | 15 | 9 | 9 | 0 | 7 |
| 2 | 88 | 94 | 68 | 32 | 32 | 0 | 7 |
| 4 | 163 | 162 | 113 | 80 | 80 | 0 | 7 |
| 6 | 235 | 230 | 163 | 109 | 109 | 0 | 7 |
| 16 | 611 | 582 | 486 | 334 | 334 | 0 | 7 |

### 1.2 Exact affineness (verified, not assumed)

For every `L ≤ 344`, **no atom of `A_L` is nonlinear in `K_L`**: no monomial of any atom in
the window contains two or more knob factors. Consequently

* each atom `a ∈ A_L` is an exact affine function of the knob vector `δ ∈ ℤ^{K_L}`;
* each equation `e ∈ R_L` has an exact affine form `row_e(δ) = c_e + n_e·δ`, with `e`
  satisfied **iff** `row_e(δ) = 0`;
* every equation **outside** `R_L` is constant in `δ`, because a knob's atoms all lie in
  `A_L` by definition, so such equations remain satisfied.

No equation is approximated and no knob is ever discarded to preserve linearity. Rows whose
linear part is identically zero are dropped as trivially constant; the counts above are
after that reduction.

### 1.3 The code

Let `N` be the matrix of rows `n_e` (one row per nontrivial equation of `R_L`, one column
per knob) and `B` the vector `-c_e`. Define

    C  :=  { N·u  :  u ∈ ℚ^{K_L} }   ⊂  ℚ^{R_L}

— the image of the **knob directions**. A subset `D ⊆ R_L` is a **support** of `C` if some
nonzero `u` has `n_e·u = 0` for every `e ∉ D`, equivalently `rank(N restricted to R_L \ D) < |K_L|`.

> **Realizability is inside this definition.** `C` is not the row space of the
> equation–atom incidence matrix and not a relaxation over arbitrary atom vectors. Its
> elements are the equation-value vectors actually produced by moving variables. The
> knob→atom map is injective with image a rank-`|K_L|` sublattice of `ℤ^{A_L}`; at `L=6`
> exactly 109 of 235 atoms are movable and 126 are frozen at 0, and at `L=16` 334 of 611
> are movable and 277 frozen (`soundness2.py`).

---

## 2. The lemma

> **Lemma.** Suppose the full row system `{ n_e·δ + c_e = 0 : e ∈ R_L }` is consistent over
> ℚ and `rank(N) = |K_L|`. Then it has a unique rational solution `W`. If `W ∉ ℤ^{K_L}`,
> then for every `δ ∈ ℤ^{K_L}` the violated set
> `D(δ) = { e ∈ R_L : n_e·δ + c_e ≠ 0 }` contains a support of `C`.

*Proof.* `rank(N) = |K_L|` and consistency give a unique rational solution `W`. Let
`δ ∈ ℤ^{K_L}` and `Z = R_L \ D(δ)`, the rows it satisfies. If `rank(N_Z) = |K_L|` then the
subsystem on `Z` already determines a unique rational solution, which must be `W`; but `δ`
solves it and `δ` is integral while `W` is not — contradiction. So `rank(N_Z) < |K_L|`,
i.e. some nonzero `u` lies in `ker(N_Z)`, and `supp(N·u) ⊆ D(δ)`. ∎

All three hypotheses are verified at every level, including `L = 16` (`eqwin2.py`,
`verify16.py`, `wcheck.py`):

* `rank(N) = |K_L|` and `ℚ-inconsistent = 0` — the table in §1.1;
* `W ∉ ℤ^{K_L}` — confirmed by the cheaper and stronger fact that the **full row system is
  already inconsistent mod `p`** at every level (7 bad rows at `L = 0, 2, 4, 6`; 13 at
  `L = 16`), so it has no integer solution at all. At `L = 0` the non-integral coordinates
  of `W` can also be exhibited directly: `x642` with denominator 2458959, `x1329`, `x9413`,
  `x10903` with denominator `p`, and `x17325` with denominator `p·2458959`, where
  `p = 2^256 − 2^32 − 977`.

The lemma converts "how few equations can fail?" into "what is the lightest support of `C`
that is also **realizable over ℤ**?" — and gives two independent necessary conditions.

---

## 3. Condition (a): mod-p consistency

Reduce the rows mod `p`. Let `Wb` be a basis of the left kernel of `N mod p` (dimension
`w = n − rank_p(N)`) and `g := Wb·B`. Then:

> The retained rows `R_L \ D` are consistent mod `p`  **iff**  `g ∈ span{ col_i(Wb) : i ∈ D }`.

*Proof.* `N_Z x ≡ B_Z` is solvable mod `p` iff every `λ` in the left kernel of `N` supported
on `Z` has `λ·B ≡ 0`. Writing `λ = Σ_j a_j Wb[j]`, "supported on `Z`" means `a ⊥ col_i(Wb)`
for `i ∈ D`, and `λ·B = a·g`. So the condition is
`{a : a ⊥ col_i ∀ i ∈ D} ⊆ ker(a ↦ a·g)`, i.e. `g ∈ span{col_i : i ∈ D}`. ∎

This is a **minimum-weight syndrome-decoding problem over F_p**, and it is the binding
condition. A large fraction of knob columns vanish mod `p` — they are the `p`-quantised
handles — which is why the filter is so strong.

| L | knobs | vanishing mod p | rank_p | w | exhaustive: no `|D| ≤ k` | Prange trials / solvable | lightest weight seen | P(a weight-≤6 point exists and was missed) |
|---|---|---|---|---|---|---|---|---|
| 2 | 32 | 17 | 15 | 53 | **k = 5** (11,290,975 subsets) | 3,000 / 2,025 | 7 | **≤ 9.5e-208** |
| 6 | 109 | 56 | 53 | 110 | **k = 3** (721,927 subsets) | 4,000 / 1,538 | 7 | **≤ 6.7e-64** |
| 16 | 334 | 172 | 162 | 324 | k = 2 (118,341 subsets) | 250 / 97 | 7 | ≤ 1.6e-4 |

Prange's miss probability is the honest one: a random information set of size `w` contains
all 6 positions of a weight-6 solution with probability `Π_{k=0}^{5} (w−k)/(n−k)`, and the
table reports `(1 − that)^{#solvable trials}`.

**The lightest weight the filter admits is 7 at every depth, and `S` attains 7.**

---

## 4. Condition (b): code support

Independently, `D` must contain a support of `C`. Equivalently, minimum support weight =
**minimum number of linearly dependent columns of a parity check `H`** (a basis of the left
kernel of `N`).

| L | H | exhaustive: no dependent subset of size ≤ k | rigorous floor | greedy-ISD usable | lightest support found | weight-≤6 supports enumerated |
|---|---|---|---|---|---|---|
| 2 | 36 × 68 | **k = 4** (866,779 subsets) | ≥ 5 | — | — | — |
| 6 | 54 × 163 | k = 3 (721,764 subsets) | ≥ 4 | 258 / 258 | **6** | **582** |
| 16 | 152 × 486 | — | — | 19 trials | **6** | **62** |

A methodological note that cost me a wrong result once: **uniform random `|K_L|`-subsets of
rows are almost never information sets here** — 0 of 3,300 had full rank at `L=6`, so a
first version of this search returned "P ≤ 1", a true statement carrying no information.
Information sets must be built **greedily from a shuffled row order**, which gives 258/258
usable trials.

---

## 5. The theorem, and why it is 7

Throughout, "assignment in scope at level `L`" means an integer assignment agreeing with `S`
on every variable outside `K_L`. Such an assignment changes no equation outside `R_L`, so
its total number of failing equations equals the number of violated rows of `R_L`.

> **Theorem A (unconditional within scope).** Every assignment in scope fails at least
>
> * **6** equations at `L = 2`,
> * **4** equations at `L = 6`,
> * **3** equations at `L = 16`.
>
> *Proof.* By the Lemma the violated set `D` contains a support of `C`, and the retained
> rows must be consistent mod `p`. Section 3 exhausts **every** subset `D` up to the stated
> size — 11,290,975 of them at `L=2`, 721,927 at `L=6`, 118,341 at `L=16` — and none is
> mod-`p` consistent. ∎

> **Theorem B (with a quantified miss probability).** Every assignment in scope fails at
> least **7** equations, and 7 is attained by `S`. The only gap is the possibility that a
> mod-`p`-consistent `D` of weight 6 exists and was missed by randomised search; that
> probability is bounded by **9.5e-208** at `L=2`, **6.7e-64** at `L=6`, and **1.6e-4** at
> `L=16` (§3, Prange with the exact per-trial detection probability).

The two statements are deliberately separate. Theorem A is a finite exhaustive computation
with no probabilistic content. Theorem B is the number the campaign cares about, and its
residual uncertainty is stated rather than absorbed.

### Tightness — the part that explains rather than bounds

**Neither condition alone gives 7.**

* Condition (b) alone permits **6**: the code genuinely has weight-6 supports — 582 of them
  at `L=6`, 62 at `L=16`, and 38,760 in the original atom-closure model of the residual
  region. The lightest is exactly the equation set of atom `a35758`, which occurs in 6
  equations.
* Condition (a) kills **every one of them**. `w6test.py` enumerated the weight-≤6 supports
  at `L=6` and `L=16` and tested each: **all 582 and all 62 fail mod-`p` consistency, and
  none is integral.** No assignment was produced at any depth.
* And 7 is precisely the lightest weight the mod-`p` filter admits, at every depth.

So the two constraints intersect at exactly the value `S` achieves. This is why the bound
is tight rather than merely unbeaten.

---

## 6. What is **not** proved

**6.1 The bound is conditional, and cannot be made unconditional by this method.**
The theorem covers assignments agreeing with `S` outside `K_L`. Making it unconditional
needs the excluded set `V_L \ K_L` to be empty. It never is, and it does not shrink —
I claimed earlier that it did and that was wrong. Measured to `L = 791` (`fastgrow.py`):

| L | atoms | vars | knobs | **excluded** | excl/vars | atoms nonlinear in knobs |
|---|---|---|---|---|---|---|
| 0 | 24 | 56 | 9 | 47 | 0.839 | 0 |
| 6 | 235 | 537 | 109 | 428 | 0.797 | 0 |
| 100 | 3,744 | 6,907 | 1,324 | 5,583 | 0.808 | 0 |
| 200 | 6,603 | 11,656 | 2,179 | 9,477 | 0.813 | 0 |
| **344** | 9,211 | 14,669 | 2,875 | 11,794 | 0.804 | **0** ← affine ceiling |
| 345 | 9,227 | 14,684 | 2,882 | 11,802 | 0.804 | **1** |
| 791 | 17,296 | 22,372 | 5,152 | 17,220 | 0.770 | 1,287 |

The excluded count rises monotonically and the excluded *fraction* stays near 0.80 for 800
levels. Every exclusion is a variable touching an atom outside the window — that is the
theorem's entire scope limitation.

**At `L = 345` the exactly-affine regime ends**: the first atom becomes nonlinear in the
knobs (both factors of some product are knobs by then), and by `L = 791` there are 1,287
such atoms. Closure — excluded = 0 — requires swallowing a whole connected component of the
variable–atom graph; the giant one has 23,843 variables, of order 1,500+ levels at the
observed ~20 atoms/level. **The affine ceiling arrives roughly four times too early.**
Continuing past `L = 344` would require approximating equations, forfeiting the one
property that makes this a proof. The conditionality is structural, not a compute budget.

**6.2 Exhaustiveness is partial above the stated sizes.** The rigorous floors, taking the
better of the two conditions at each depth, are **`≥ 6` at `L=2`** — every one of the
11,290,975 subsets of size ≤ 5 fails mod-`p` consistency — then `≥ 4` at `L=6` and `≥ 3` at
`L=16`. Closing the last step at `L=2` rigorously would need all C(68,6) = 119,877,472
size-6 subsets, about 9.5 hours at the measured 3,500 subsets/second; it was not run. Between those floors and 7 the evidence is the
Prange miss probabilities of §3 and the constructive enumeration of §5 — overwhelming at
`L=2` and `L=6` (1e-208, 1e-64), weaker at `L=16` (1.6e-4, with a larger campaign running).

**6.3 A vacuity check, and its answer.** The *raw* relaxation — minimise `‖M·β‖₀` over
arbitrary nonzero integer atom vectors `β` — is vacuous, which I confirm independently:
**3,235 atoms occur in exactly one equation**, so `‖M·e_a‖₀ = 1` for each. My code is not
that code (§1.3). Is a weight-1 support nevertheless admissible here? It is not excluded by
fiat; it simply does not occur, for three separately checked reasons:

* **zero** single-equation atoms lie in any window (0 of 235 at `L=6`, 0 of 611 at `L=16`);
* globally, **none of the 3,235 single-equation atoms carries a private variable** — a
  variable occurring in no other atom — so none is independently settable anywhere in the
  instance (`cancel.py`);
* rigorously, the minimum support is ≥ 5 at `L=2` and ≥ 4 at `L=6`, and 6 as observed.

Honest nuance: the code **is** sensitive to low-occupancy atoms. The observed minimum
support 6 is exactly the equation set of a 6-equation atom. Low occupancy pushes the
minimum support down to 6; it does not reach 1 because the atoms that would are unreachable
by any knob.

---

## 7. An independently established fact this relies on

Agent F of this campaign computed, for the full 39,033 × 39,033 equation–atom incidence
matrix `M`, that **rank(M) = 39,033 and dim ker(M) = 0**, by a characteristic-free peeling
certificate re-verifiable from `M` on disk together with Wiedemann over a word prime, with
pivots checked non-divisible by any odd prime so the result holds over ℤ. I did not
re-derive this and I attribute it to F.

Consequence: **any assignment satisfying all 39,033 equations must make every atom exactly
zero.** There is no cancellation available at all, so the all-atoms-zero model is an
*equivalence*, not a restriction. This strengthens §6.3 — I had predicted the reconciliation
between a nontrivial kernel and my windows would be that kernel vectors are not in the image
of the knob→atom map, i.e. that realizability rather than cancellation is the binding
ingredient; F's result makes that unconditional by removing the kernel entirely.

Two notes for a reader reconciling the two computations.

* **Dimensions.** In my parse the equation–atom incidence matrix is 39,033 × **42,267**
  (equations × atoms), not square. Full row rank there would leave a kernel of dimension
  42,267 − 39,033 = 3,234, so anyone checking F's `dim ker(M) = 0` against my numbers
  should first confirm which atom set F used — plausibly a reduced or square restriction.
  I flag the mismatch rather than paper over it; it does not affect anything below.
* **The two results are complementary, not in tension.** F gives `min ‖M·β‖₀ ≥ 1` over
  nonzero integer atom vectors (no kernel). §6.3 gives `min ‖M·β‖₀ ≤ 1` (3,235 atoms occur
  in exactly one equation). Together the minimum is **exactly 1** — which is simultaneously
  why the raw relaxation can bound nothing and why all-atoms-zero is an equivalence.

---

## 8. Reproduction

The `s9/*.pkl` caches are gitignored and must be rebuilt first:

    cd solve_lab/s9 && python3 atomize.py && python3 poly.py && python3 gates.py && python3 fwd.py
    cd ../agentA_work

    python3 eqwin.py  <state> 8          # window sizes + exact-affineness check (§1.2)
    python3 eqwin2.py <state> 6          # rank / Q-consistency per level  (§1.1, Lemma)
    python3 boundary.py                  # excluded-set accounting          (§6.1)
    python3 fastgrow.py <state> 800      # growth to the affine ceiling     (§6.1)
    python3 eqbound.py <state> <L> 6 <trials> <exh_k>   # condition (a)     (§3)
    python3 eqmindist.py <state> <L> <exh_k> <secs>     # condition (b)     (§4)
    python3 w6test.py <state> <L> 6 <secs>              # tightness         (§5)
    python3 soundness.py 2 6 ; python3 soundness2.py 6 16   # vacuity check (§6.3)
    python3 cancel.py                    # single-equation-atom census      (§6.3)
    python3 full31.py 6                  # the original atom-closure model of the residual

    python3 wcheck.py 0 2 4 6 16         # Lemma hypothesis: W non-integral (§2)
    python3 verify16.py 16               # rank / Q-consistency at the deepest window

where `<state>` is `solve_lab/best/new_instance_partial_39026.json`.
Raw logs for every number quoted are in `agentA_work/runs/`.

## 9. Searches left running at the end of the session

These would only tighten numbers already stated; none can weaken the theorem, since a
positive find would have been written to disk as `A_*.json` with its checker-verified score
and none was.

* `runs/eqmd2.log` — condition (b) dependent-column search at `L=2`, size 5
  (sizes ≤ 4 exhausted, floor `≥ 5`).
* `runs/eqb16b.log` — a 3,000-trial Prange campaign at `L=16`, to drive the
  `≤ 1.6e-4` miss probability down; at the measured rate it needs hours.

Nothing at any depth ever admitted a weight below 7, and no assignment beating 39,026 was
produced at any point in this work.
