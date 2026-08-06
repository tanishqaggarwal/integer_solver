# Session 11, Part VI — the instance as an INTEGER PROGRAM

Parts I–V decoded the circuit and proved a deficit. This part drops the circuit reading
entirely and treats the problem as what it is:

> **38,748 integer variables, 39,033 polynomial equations; minimise the number of violated
> equations.**

Four IPs, all solved exactly (certificates, not heuristics). Scripts: `s11/ip1.py` … `s11/ip8.py`.

---

## IP #1 — minimum-cost defect placement  (`s11/ip1.py`)

    binary y_c = 1  iff constraint c is left violated
    feasible   :  the surviving constraints admit a perfect matching into their controls
    objective  :  | union of equations occupied by the violated constraints' atoms |

The live constraint set of a channel is small (14), so this is solved by **exact enumeration
over all 2^14 subsets**. For the (490,91) channel:

    full-set matching = 12 of 14   ->  deficit 2

    cheapest feasible placements
      15 equations  <- violate {mirror25118, mirror3719}      <== IP OPTIMUM
      20            <- {a29539, mirror25118}
      21            <- {a40065, a688}
      22            <- {a21050, a40065} ...

> **IP optimum = 15 failing equations, score 39,018** — exactly the score achieved and
> verified. The construction is certified optimal *for this channel*.

## IP #2 — a global lower bound  (`s11/ip2.py`)

A check owning a **private handle** (a free variable occurring in exactly one atom, linearly)
can always be closed by itself, so it can never be the thing you are forced to violate.
Absorbers must come from the rest. Over the whole instance:

    checks                         10,792
      self-closable (private handle)  121
      possible absorbers           10,671
    cheapest absorber                 1 equation   (3,234 checks live in a single equation)
    cheapest 2-absorber union         2 equations

> **Global lower bound: no defect placement anywhere can cost fewer than 2 failing equations,
> i.e. no score above 39,031.**

The same script exposed something the earlier cost model had missed: the checkpoint's seven
atoms `{22229, 22230, 35758..35762}` occupy **12** equations, yet only **7** fail. Five of the
twelve vanish because the nonzero atoms *cancel* there. Cancellation is real, and it is why
the checkpoint is cheap — so the right objective is over **equations**, not atoms.

## IP #3 / #4 — minimum-weight coset  (`s11/ip3.py`, `s11/ip4.py`)

At a point the reachable perturbations form an affine integer lattice. With `b` the vector of
equation inner-sums and `G` the matrix of handle effects, the objective is literally

        minimise   || b + G k ||_0      over integer k

a **minimum-weight coset problem**. Constrain the satisfied equations to stay satisfied (an
integer kernel via column HNF), then minimise nonzeros on the failing ones inside that kernel.
Because the handles interact, the exact map is polynomial, so IP #4 iterates the linearisation
(integer Newton on the lattice), accepting a step only if the true count does not increase:

    closehit2 (28 failing) -> model says 0, applied 27 -> re-linearise -> 15 -> 15   (fixed point)
    finish3   (15 failing) -> no improvement, kernel dim 10 over 158 kept equations

> Two structurally independent exact methods — matching (IP #1) and lattice coset (IP #3/#4) —
> both return **15** for this channel. That is a strong cross-check.

## IP #5 / #7 / #8 — the checkpoint, with no circuit orientation at all

The pipeline of Parts I–V forward-evaluates, which **destroys** the checkpoint: its score rests
on five deliberately broken gates, and `fw.forward` repairs them into a 37-failing state. The IP
formulation needs no orientation — the state is just the raw assignment — so it applies
directly.

At `best/new_instance_partial_39026.json` (7 failing):

    variables that can move the failing equations : 56   (19 with an EXACT LINEAR effect)
    collateral equations they disturb             : 123
    compensators that introduce no new collateral : 50
    reduced system                                : 130 x 69
    kernel dim preserving all satisfied equations : 2

Subset search with modular pre-screening (solvability mod two 61-bit primes is a cheap
necessary condition, run before the exact HNF solve):

    allow = 0 : infeasible
    allow = 1 : infeasible   (all 8 subsets)

Notably **every subset passes the modular screen but fails over ℤ** — the obstruction is pure
**divisibility**, not rank. The system is solvable mod p and unsolvable over the integers.

> So the checkpoint's 7 is rigid under every exact-linear move in its locality: one cannot even
> drop one failing equation and satisfy the rest.

---

## What the IP view settles

| question | answer | method |
|---|---|---|
| optimum for the (490,91) channel | **15** (score 39,018) | exact enumeration, 2^14 subsets |
| same, independently | **15** | minimum-weight coset + integer Newton |
| global lower bound | **2** (score ≤ 39,031) | absorber eligibility |
| is the checkpoint locally improvable? | **no** (allow 0 and 1 infeasible) | reduced IP + modular screening |
| nature of the obstruction | **divisibility over ℤ**, not rank mod p | modular screen passes, exact solve fails |

The last row is the sharpest new statement in this part: every previous session described the
wall in terms of rank, kernels and congruences mod p. In the integer program the mod-p system
is *solvable*; what fails is integrality. The gap between the global bound (2) and the best
known placement (7) is where any further progress must come from.

---

## IP #9 / #10 — the obstruction, reduced to ONE NUMBER

If the mod-p system is solvable and the integer one is not, the obstruction is an
invariant-factor condition. Solving the checkpoint's system over ℚ (`s11/ip9.py`):

    consistent over Q
    solution supported on exactly 7 of 69 variables:
        x_642, x_1329, x_9413, x_10903, x_17325, x_29854, x_31864
        (precisely the x_2099 ladder — the checkpoint's own defect)
    least common denominator D = 284727999958901794582548685725978209206968908223438855709900581314503500195397778817
                              = 2458959 * p          (84 digits;  2458959 = 3 * 819653, both prime factors)

The RREF denominator need not be minimal (the system has a kernel), so `s11/ip10.py` computes
the **true invariant**: the least `d` for which `M x = d·rhs` has an integer solution, testing
every divisor of D:

    d = 1        no        d = P          no
    d = 3        no        d = 3*P        no
    d = 819653   no        d = 819653*P   no
    d = 2458959  no        d = 2458959*P  SOLVABLE

> **The entire remaining obstruction at the 39,026 checkpoint is a single divisibility
> condition by `2458959·p`.** Not a rank deficiency, not an inconsistency — one number.
> Every proper divisor fails; that exact product succeeds.

And the failing values themselves are the reason (`s11/ip11.py`): at the checkpoint **0 of 7**
failing equation values are divisible by p (gcd of all seven = 1); at the session's 39,018
state, 0 of 15. The p-factor of the invariant is exactly the p-quantisation every previous
session measured indirectly — here it appears as one factor of one integer.

### The actionable consequence

The invariant factors as `3 · 819653 · p`. The p-part is the designed wall. **The other part is
seven digits.** So:

> If a state can be reached in which the failing equation values are ≡ 0 (mod p), the whole
> obstruction collapses to divisibility by **2458959** — a 7-digit modulus, attackable by
> exactly the CRT / quadratic-form method that closed the 8640431 condition in Part II
> (`s11/quad3.py`).

That is the single sharpest target this instance has produced, and it is stated entirely in
integer-programming terms: make the right-hand side p-divisible, then clear a 7-digit
invariant factor.

---

## IP #12 — the p-factor is UNIVERSAL

The invariant was `2458959·p` at the checkpoint and `8640431·p` at the session's 39,018 state.
Computing it at every saved state (`s11/ip12.py`), across constructions built by completely
different routes:

    state                        failing   D digits   p | D    D/p
    new_instance_partial_39026        7        84      yes     2458959
    finish3_named   (39,018)         15        85      yes     8640431
    closehit2       (39,005)         28        78      yes     1
    three                            39        85      yes     8640431
    tri7_best                        52        85      yes     8640431
    eqopt2_named                     15        85      yes     8640431
    cheapdefect_named                15        85      yes     8640431
    (quad3_hit, uv01_full: inconsistent over Q — genuinely rank-obstructed states)

    consistent over Q        : 7 / 7
    invariant divisible by p : 7 / 7
    cofactors D/p            : {1, 2458959, 8640431}

> **The residual integer program is always consistent over ℚ, and its invariant factor always
> contains p.** The small cofactors are state-dependent bookkeeping — 8640431 is the mirror
> handle's multiplier, 2458959 the ladder's — and both are clearable by CRT (Part II did exactly
> that for 8640431). The factor of **p** is the one that never leaves.

At `closehit2` the cofactor is **1**: the invariant there is *exactly* p, with no small part at
all. `s11/ip13.py` confirms it directly on the 357×190 system — `M x = d·rhs` is unsolvable for
d = 1, 2, 3 and solvable at d = p.

### The trapdoor, in one sentence

> Every state reachable in this instance leaves a residual integer program that is **solvable
> over the rationals** and whose **sole integrality obstruction is a single factor of
> p = 2²⁵⁶ − 2³² − 977**.

That is what ten sessions were circling under the names *p-quantisation*, *the conserved
obstruction*, *the deficit of 2*, and *"7 is an invariant"*. In integer-programming language it
is one invariant factor, and it is the same one everywhere.

It also says precisely why the score cannot be pushed by better search: the objective's floor is
set by how many equations must absorb that single factor of p, not by any failure of the search
to find a feasible point. A full solve requires *removing the p from the invariant* — i.e.
reaching a state whose failing right-hand side is already p-divisible — which is exactly the
problem the setter encoded.
