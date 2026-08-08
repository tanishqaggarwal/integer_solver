# Convex / polyhedral structure in the modmul constraint system

Goal: find hidden convex polyhedral structure that lets the LP relaxation
*determine* variables exactly, so they can be dropped from the QUBO, with a
soundness gate (an LP-determined variable must take that value in every integer
solution). All numbers below are measured on real builds from
`../../squeeze/mm.py` + `../../squeeze/mmqb.py` (schoolbook, quotient/naf
reduction). scipy 1.17 (HiGHS / `milp`) does the LP and ILP work.

## TL;DR verdict

**The LP relaxation is weak exactly on the product/carry core, which is where
all the variables are.** The carry/summation network is **not** totally
unimodular, its LP relaxation is **not** integral, and **~0%** of carry
variables are LP-forced even when every partial product is fixed. The only
sound *convex* reduction is the linear-equality elimination (~34% of variables),
which is ordinary Gaussian presolve, not LP-face structure. Genuine LP-face
persistency removes **0-1 variables**. The arithmetic *is* fully determined by
the input bits (proven: unique integer witness), but that determinacy is
**nonconvex** and invisible to the LP.

---

## How the system is structured (why an LP is even well-defined)

Every relation in `mm.py` is asserted through `assert_terms`/`assert_zero`,
which (a) replaces each nonlinear monomial `a_i b_j` by a fresh AND variable
`z = a_i AND b_j`, and (b) balances the resulting **linear** identity
column-by-column with carries. So the E=0 feasible set is:

```
x in {0,1}^n
lin . x + const == 0      for every add_square(lin,const)   (linear equalities)
z == x_i AND x_j          for every AND cache entry           (nonconvex)
```

LP relaxation `P`: box `0<=x<=1`, keep the equalities exactly, and replace each
AND by its **McCormick hull** `z<=x_i, z<=x_j, z>=x_i+x_j-1, z>=0` (= conv hull
of `{z=x_i x_j}` on the cube). `P` contains the integer feasible set, so any
coordinate constant over `P` is constant over every integer solution -- the
soundness guarantee. Code: `lp_core.py`.

---

## 1. The carry/summation network is NOT totally unimodular (LP not integral)

`analysis.py` extracts the equality rows restricted to the carry/adder columns
(products/words moved to the RHS) and inspects them.

- **Matrix entries** are `{-2, -1, +1}`. The `-2` is the positional
  `... - s - 2*carry` term intrinsic to a full adder `a+b+c = s + 2d`.
  A matrix with an entry of magnitude 2 **cannot** be TU (its 1x1 minor is 2).
  A sampled 2x2 submatrix `[[-2,-1],[0,1]]` has determinant **-2** -- an explicit
  non-TU certificate.
- **Not a network matrix either** (network matrices are TU; entries are +-1).
- **LP is not integral.** With inputs fixed to a witness, the LP relaxation of
  the *whole* system has a positive-dimensional feasible polytope: for p=13,
  s=4, **20 coordinates are fractional** at an LP vertex, e.g. a carry bit at
  **0.5** and another at **0.75**. A single full adder already shows why:
  `a+b+c = s+2d` with `a,b,c` fixed leaves the 1-parameter family
  `s = k - 2d, d in [0,1/2]` in the box, and the network coupling does not close it.

**Carries LP-forced *given the products* (part-1 headline).** Fix the input
bits *and* all AND-product variables to a witness, then LP-min/max every carry
(`scaling.py`):

| instance | s | #carry vars | carries LP-pinned given products | fraction |
|---|---|---|---|---|
| p=13     | 4 | 54   | 2   | 3.7% |
| p=251    | 8 | 162  | 0   | 0.0% |
| p=65521  | 16| 600  | 2   | 0.33% |

The handful that pin are structurally-forced top carries (always 0). The
fraction goes to 0. **Projected to a 256-bit modmul: ~0 of its 134,138 carry
variables are LP-forced given the products.** The convex relaxation contributes
essentially nothing to determining the carry network.

## 2. LP-persistency of the full modmul

`scaling.py` (LP min/max per variable, HiGHS), `integer_gap.py` (ILP via
`scipy.milp` for the soundness gate).

- **Inputs free (removable-to-constant):** `0-1` variables across all sizes.
  The one occasional fixed variable is a top carry pinned to 0; verified to hold
  in **all 256 witnesses** of the p=13 instance (0 violations). The LP fixes
  nothing useful when inputs are free, because every internal variable is a
  *function* of the inputs, not a constant.
- **Inputs fixed:** the LP pins **exactly the AND products** (McCormick is exact
  at integer input vertices) and **almost no carries**:

| instance | internal vars | integer-determined | LP-determined | convex gap |
|---|---|---|---|---|
| p=13  | 79  | **79 (100%)** | 18 (all 16 products + 2 carries) | 61 |
| p=251 | 243 | **243 (100%)**| 64 (all 64 products, 0 carries)  | 179 |

  `integer_gap.py` proves via ILP that with inputs fixed the **witness is the
  unique integer solution** (int-min = int-max for every variable,
  `witness_unique=true`, 0 soundness mismatches). So integer-persistency is
  100% -- the whole circuit is determined by the 2s input bits -- but the LP
  captures only the McCormick-exact products. **The convex gap is precisely the
  carry network + output word.**

Soundness gate: **passed.** Every variable the LP claimed determined matched the
integer witness (`soundness_bad = 0` everywhere); the unique-witness ILP proves
no other integer solution exists to violate a claim.

## 3. Newton polytope / algebraic structure

`newton.py`. For the schoolbook product `A*B = sum 2^{i+j} a_i b_j`:

- **Newt(A*B) = Delta_{s-1} x Delta_{s-1}**, the product of two standard simplices:
  affine dimension `2(s-1)`, and all `s^2` product monomials are **vertices**
  (verified s=2..5). It equals the **Minkowski sum** Newt(A) (+) Newt(B)
  (verified: vertex set = Minkowski sum, exactly).
- A comb window's atoms (`sum_t u_t * table[t]`) form the **standard simplex
  Delta_{D-1}** (D=2^w): all D one-hot atoms are vertices; this polytope *is*
  integral (the assignment polytope), consistent with the ladder's prefix-counter
  one-hot being LP-exact.

**Reduction verdict:** because every one of the `s^2` product monomials is a
*vertex* of the Newton polytope, none is a convex combination of the others --
there is **no lower-dimensional face the solution lies on** and **no toric /
tropical collapse** that removes products. The Minkowski decomposition
`Newt(A*B)=Newt(A)(+)Newt(B)` is exactly the **bilinear / Karatsuba-Toom**
structure that `mm.py` *already* exploits to cut the number of *limb sub-products*
(from `n^2` toward `n^1.585`/`n^1.465`); it is an algebraic-identity saving,
not a convex-geometry one, and it does not shrink the `{0,1}` search set of a
fixed multiplier layout.

## 4. Effective dimension and the qubit-reduction bottom line

For the real secp256k1 256-bit modmul (`resources_256.py`, schoolbook/naf):

```
n (total vars)        200,699
 AND products          65,536   (= 256^2)
 carry/adder vars     134,138
 input bits (A,B)         512
 equalities (rank)     67,583   (rank == #equalities, full row rank, verified s<=16)
```

Layered picture of "must-search vs determined":

| layer | count (256-bit) | nature | sound convex reduction? |
|---|---|---|---|
| total QUBO variables | 200,699 | -- | -- |
| removable by **LP-face persistency** (const over P) | **~0** | convex | yes, but negligible |
| removable by **linear-equality elimination** | **67,583 (33.7%)** | convex (equalities) but = Gaussian presolve | yes; caveat below |
| affine hull after equality elimination | 133,116 | -- | -- |
| **carry core** (convexly irreducible) | 134,138 (67%) | nonconvex | **no** -- TU fails, LP not integral, 0% pinned |
| true intrinsic DOF (input bits) | 512 | nonconvex | unreachable by LP |

- **Genuine LP win:** ~0 variables. LP-face persistency does not remove anything
  the ordinary equality presolve doesn't already give.
- **Linear-elimination win:** 33.7% of variables are exact affine functions of
  the rest and can be substituted out. This is sound and real, **but it is
  Gaussian elimination on the equality rows, not LP-polytope structure**, and
  realizing it as a *smaller QUBO* is not free: substituting a variable that
  appears inside an AND (quadratic) penalty raises the degree, so the clean
  removable subset is the variables occurring only in linear-substitutable
  positions. Treat 33.7% as the affine-dimension upper bound on the reduction, to
  be realized by the presolver with possible densification.
- **The 67% carry core cannot be convexly removed at all.** This is the honest
  negative result: the LP relaxation is weak on the product/carry core, so
  convexity offers **no path** from 200k qubits toward the 512-bit intrinsic
  dimension. Everything below the equality-elimination layer stays nonconvex and
  must be handled by the annealer's integrality/AND machinery.

## Soundness summary

No unproven reductions are claimed. Every LP-determined value was checked against
the integer witness (`soundness_bad=0`); the free-input fixed variable was checked
against all 256 witnesses of the p=13 instance (0 violations); and `scipy.milp`
proves the witness is the unique integer solution given inputs, so
"LP-determined => same in every integer solution" holds by `P` containing the
integer hull.

## Files

- `lp_core.py` -- build a modmul instance, extract equalities + McCormick, LP
  persistency (random-objective filter + confirmed min/max).
- `analysis.py` -- equality rank (exact, over GF(2^61-1)), carry-subsystem
  extraction, TU necessary condition + sampled-determinant certificate,
  fractional-point demo.
- `scaling.py` -- persistency / carries-given-products across sizes; soundness
  across all witnesses.
- `integer_gap.py` -- ILP proof of unique witness + integer-vs-LP determinacy gap.
- `newton.py` -- Newton polytope = Delta x Delta, Minkowski check, comb-window simplex.
- `resources_256.py` -- real 256-bit secp256k1 counts.
- `scaling_results.json`, `integer_gap.json`, `newton.json` -- raw outputs.
