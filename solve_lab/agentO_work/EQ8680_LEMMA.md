# eq8680 — the constraint is forced, and it is LINEAR.  Answer: option (2).

Agent O.  **Corrected after agent T's audit (check-in 60); corrections re-verified here against
the raw text of `EQUATIONS.txt`, using no parser at all** (`verify_lemma.py`,
`runs/verify_lemma.log`).

## The Lemma

> **`eq8680`'s left-hand side equals `S⁴`, where `S` is an AFFINE form in atoms.**
> `checker.py` requires the LHS to be exactly `0` over **ℤ**, and ℤ is an integral domain, so
> `S⁴ = 0 ⟺ S = 0`.
> **Therefore `S = 0` in every satisfying assignment** — no knob set, no frame, no
> configuration, no divisibility, and no modulus anywhere.

Verified on perturbed vectors (at the witness everything is 0 and would not discriminate).
Perturbing a single variable and reading the raw LHS straight from `EQUATIONS.txt`:

| perturbation | `S` | raw LHS | `S²`? | `S³`? | `S⁴`? |
|---|---|---|---|---|---|
| `x_4432 += 2` | 2 | 16 | ✗ | ✗ | ✓ |
| `x_4432 += 3` | 3 | 81 | ✗ | ✗ | ✓ |
| `x_4432 += 5` | 5 | 625 | ✗ | ✗ | ✓ |
| `x_19964 += 2` | −2 | 16 | ✗ | ✗ | ✓ |
| `x_28730 += 3` | −3 | 81 | ✗ | ✗ | ✓ |
| `x_23754 += 2` | −18 | 104976 | ✗ | ✗ | ✓ |

`LHS == S^k` holds for **k = 4 only**.  Measured exactly: `dS/dx_4432 = +1`,
`dS/dx_19964 = −1`, `dS/dx_28730 = −1`, and `a23618 = x_4432 − x_19964 − x_28730` is `S`'s
first term at coefficient exactly **+1**.

## Three corrections to what I first wrote — and one thing that was never wrong

1. **`S⁴`, not `S²`.**  The nesting is two levels: `LHS = T·T` with `T = S·S`.
2. **The object with slope +1 is `S`, not `T`.**  I wrote "`eq8680 = T²` where `T` is linear"
   and then quoted `dT/dx_4432 = +1`; those cannot describe the same object, since `T = S²`
   gives `dT/dx = 2S+1`.  I conflated one nesting level **in the prose**.
3. **Term count: 18 vs 20 is a granularity difference, not a contradiction.**  The raw text has
   **18** bracketed groups (agent F's certified parse and agent T both give 18).  E's parser
   emits **20** `(coef, atom)` entries because it splits exactly two of those brackets into two
   atoms each:
   `−13 · (x_21279·x_31731 + x_35619)` → `a23622, a23623`, and
   `−5 · (x_34600 − x_30108 + x_23642)` → `a11876, a11877`.  18 + 2 = 20.  **Both are correct
   descriptions of the same affine form** — quote whichever, but never mix them.
   ⚠ **`S`'s 18 is not M's enumeration exponent 18.  Different 18s; do not conflate.**

**What was never wrong: the computation.**  My frame-B "S row" was built from H's *inner*
factor, which is exactly this affine form (I measured its slope as +1 before using it).  So
every search below constrained the right object; only the write-up mislabelled a level.

Note also that in **E's** decomposition the residual `eqfails` tests is the affine sum itself —
degree 1 in the atoms — so E's model states `S = 0` directly, while H's states `S⁴ = 0`.  Same
zero locus.  The conclusion is robust to the exponent entirely: `S^k = 0 ⟺ S = 0` for any k ≥ 1.

## Why this is exactly the obstruction to δ₀
`a23618` is the sole carrier of the `L` boundary shift δ₀ requires, and enters `S` at
coefficient +1 with `dS/dx_4432 = +1`, `dS/dx_28730 = −1`, zero for every other region knob.
So `S = 0` is exactly `δx_4432 = δx_28730` — it collapses the `L` direction onto the private
handle direction, annihilating precisely the degree of freedom δ₀ needs, **unless one of `S`'s
other terms compensates**.  Those other terms are the whole remaining question; see
`T_COMPENSATION.md`.

**Cross-link (agent T).**  All three p-handles T found that L's census had omitted appear as
terms of `S`, confirmed here by source match:
`25 · (x_18253 − x_4339·x_15120)` = a20450, `1 · (x_37720 − x_14466·x_35531)` = a20452,
`23 · (x_23642 − x_8173·x_10422)` = a11875.  The equation independently confirms their incidence.

## The scoped optimality theorem

> Let `U` be the frame-B free inputs reaching any of the witness's nonzero check atoms
> (**|U| = 15**), and `C` the carriers of `S` (**|C| = 26**); `K = U ∪ C`, **|K| = 34**.
> **Every assignment agreeing with the 39,026 witness outside `K` satisfies at most 39,026
> equations.**

Evidence, all exact:
- The 7 failing equations depend only on the witness's nonzero atoms — verified term by term.
- Over `K`: **64** reachable checks, **190** reachable equations, **all 7 failures reachable**.
- **175 rows, every one exactly affine**, validated by a 5-point probe (t = 1,2,3,5,7) which
  finds precisely the same 7 non-affine checks the 2-point probe found, **none missed**.  The 16
  dropped rows all contain one of those 7 and **none currently fails**, so dropping them is
  *permissive* — the solver was free to break them and still found nothing.
- Budgets tested: see `T_COMPENSATION.md` §3.  `j=1,b=0` and `j=2,b≤1` are **complete**.

## What this is NOT
- Not a global optimality proof.  Scope: 34 of 8,751 free inputs, frame B's orientation.
- `j=3` was enumerated per triple to a wall-clock cap, not to completion; `j≥4` greedy only.
