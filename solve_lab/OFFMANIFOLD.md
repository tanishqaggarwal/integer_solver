# The flawed assumption: forward-eval is NOT WLOG

The user confirmed a solution definitely exists. That means an assumption was wrong.
Found it:

## Every prior attempt searched only the forward-evaluation manifold
15+ agents parameterized the problem as "choose 8583 free inputs -> H.forward()
computes all gate outputs". That construction makes **every gate-definition atom
identically zero**. Measured at best_agentA: only **3 of 46,298 atoms are nonzero**.

But the real problem is "39033 polynomial equations in 38748 integer unknowns".
Gate outputs are just variables. Forward-eval is one *submanifold* of dimension 8583
inside a 38748-dimensional space — it discards ~30,165 degrees of freedom.

## Why that matters: equations are linear combinations of atoms
**32,443 of 39,033 equations contain >1 distinct atom**, so an equation can hold with
individual atoms NONZERO that cancel in the sum. Example, eq 22044 = `39649·F·F` with
```
F = (x_642 - x_28599*x_17325) - 24*G1 - 10*(x_28730 - x_17499*x_9413)
    - 29*G2 - 9*(x_23754 - x_26874*x_6947)
```
Forward-eval pins the three gate atoms to 0, collapsing F to `-24*G1 - 29*G2`. That
collapse is what produced the "G1 and G2 must both vanish" conclusion and the whole
dim-1 conserved-obstruction story. **Off-manifold, F=0 has vastly more solutions.**

## The gate-enforcing equations do NOT force atoms to zero
Critical check: eq 9123 (which enforces the x_642 gate) is a **single factor — a sum of
many atoms**, not `c·(x_642 - x_28599*x_17325)²`. So the gate atom need not vanish; it
only has to be balanced by the other atoms in that sum. Verified empirically: moving
x_642 off-manifold breaks exactly ONE equation (9123) — the coupling is local, not global.
And x_642 appears in 10 equations, **9 of which are among the 11 failures**.

## Consequence
All prior "infeasible / dim-1 conserved obstruction / defect=1 invariant" results are
statements about the 8583-dim submanifold, NOT about the actual problem. The correct
test is the Rouché-Capelli consistency of the repair system using ALL 38748 variables
(fullvar_jac.py).

## Verifier
`verify.py FILE.json` re-evaluates every equation in EQUATIONS.txt exactly over ℤ.
Accepts either a free-input file (completed via forward pass) or a FULL assignment —
the mission's success criteria explicitly allow a full assignment, so off-manifold
solutions are admissible. Confirmed independently: best_agentA_39022.json -> 39022/39033,
failing [2554,6816,8124,8680,9421,12231,12270,12350,14584,22044,29125].
(Root SOLUTION.json is a stale artifact of the previous instance: 12716/39033.)

---

# CORRECTION (this document's premise is REFUTED)

The claim above — that atoms may be nonzero and cancel — is **false for this instance**.

**The atom-coefficient matrix is square and nonsingular.** The "46,298 atoms, only 3
nonzero" count was wrong because `atoms/poly_atoms.jsonl` **mixes composite atoms with
primitives** (atom 42669 is the *entire* eq-22044 body; 20862/20864 are its constituents).
Under the primitive decomposition there are exactly **39,033 atoms and 39,033 equations**,
and the coefficient matrix C has rank **39033**, nullity **0** (verified over two independent
primes; full rank mod p ⟹ full rank over ℚ).

Therefore `C·t = 0 ⟹ t = 0`: **every atom must be exactly zero in any solution.** No
cancellation-in-the-sum is possible — the bundling is invertible by construction. So the
forward-evaluation manifold IS without loss of generality, and the earlier "prior analyses
were artifacts" conclusion is withdrawn.

**Empirical confirmation:** the ring-3 first-order direction (which is real — the linearized
ambient system is consistent) does NOT integrate to a solution. Applying the computed delta
took the system from 11 failures to **3,943** exact-ℤ failures (2,618 even mod p), because
the system is quadratic and the mod-p step size is ~2^256. Same divergence signature as
prior predictor-corrector attempts.

## What survives and is genuinely useful
- `checker.py` (at `solve_lab/checker.py`) accepts a FULL assignment and never re-derives
  gate outputs — verified by bumping x_642 by 1 and watching eq 9123 fail literally. So the
  answer format was never the constraint.
- **The clean problem statement**: the circuit DAG has 8,747 free inputs, 30,001 defined
  variables, and 9,032 check equations — of which **exactly 2 are violated**:
  ```
  α:  x_7068 ≡ x_2099  (mod 7376877·p)
  β:  x_4432 ≡ x_19964 (mod p)
  ```
  x_7068 and x_4432 are free inputs with tiny cones (31 vars, 7 checks). Local repair fails
  only because the ripple terminates at absorbers forcing δ ≡ 0 (mod p), colliding with α.
- Reaching 39033 therefore requires a **different global free-input assignment** satisfying
  those two congruences — not a perturbation of best_agentA.
