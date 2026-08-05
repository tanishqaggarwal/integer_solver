# METHOD SUMMARY

## Problem
`EQUATIONS.txt`: 39,031 polynomial equations `= 0` over integer unknowns
`x_0 … x_38747`. Find integers satisfying every equation exactly in ℤ.

## Result
**39,013 / 39,031 equations satisfied exactly** (verified, no floating point) by
`best/best_partial_39013.json`. 18 equations remain, all tied to **4 atoms** in a
single 256-bit densely-coupled core.

## What produced the solution

### 1. Structural reverse-engineering
Parsed every equation with Python `ast` (`x_N` are valid identifiers). Each LHS,
after stripping an outer scalar / square, is a **top-level `+`-chain** of terms
`coef * atom`. Expanding each atom to a canonical integer polynomial (gcd-reduced,
sign-normalized) shows the equations are **random linear combinations of 46,275
shared "atoms"**, each atom reused ~10× on average. Therefore an assignment that
makes **every atom vanish** satisfies every equation. (`poly_atoms.py`)

Atom taxonomy: 20,090 linear (add/sub/copy/scalar/NOT), 25,468 degree-2
(multiply/square/boolean), 717 degree-4. 3,484 boolean vars `x*(x−1)=0`; 1,103
unit pins `x=1`; 865 "huge" atoms `bit*(x_B − HUGE_290bit) − s*x_C` (bit-gated
constant/residue loads).

### 2. Exact checker  (`checker.py`)
Compiles each equation to a Python code object over an integer array and
re-evaluates with big integers. `python3 checker.py <assignment.json>` →
`satisfied k/39031`. Used to validate every candidate.

### 3. Integer propagation  (`propagate.py`, `solve_forward2.py`)
Direction-agnostic engine: whenever an atom reduces to a single unknown that is
linear (or a solvable quadratic), solve it; propagate to fixpoint.
- From the 1,103 unit pins → **5,897 variables forced** with zero choices.
- Setting the 1,156 free boolean inputs to 0 and re-propagating solves the value
  wires `x_B` through the huge atoms and computes all gates → **39,013/39,031**.

### 4. Core isolation  (`analyze_core.py`, `main_component.py`)
The residual splits into **one giant component (23,843 vars, 256 free bits)** plus
~297 tiny components that are homogeneous (satisfied by zeros). The 4 open atoms
are residue-consistency constraints (e.g. an add-gate whose two sides are pinned to
different 290-bit residues). Their backward cone touches ~255 of the 256 bits
(dense, cyclic) — no small local fix.

## What did not crack the core
- **z3 / SMT** on the component (23,843 vars, 28,386 nonlinear int constraints):
  `unknown` after >15 min. Conflict-cone z3 (radius 3) engulfs the whole component
  and also returns `unknown`. (`z3_main.py`, `cone_solve.py`)
- **Single-bit flip search**: individual bit flips reduce violated atoms 4→3 but
  introduce contradictions — bits are tightly coupled (need a consistent set).
- **Boolean-forcing propagation** (`bit=1` when `s*x_C≠0`): cannot bootstrap — the
  all-zeros state is itself a fixpoint; the witness is a non-trivial fixpoint.
- Modulus hunt: gcd of the 514 huge constants = 1; the power chains
  `x, x², x³` grow **unreduced** — it is not a single-modulus reduction.

## Assessment
The instance is an obfuscated arithmetic circuit whose 256-bit input is recoverable
only by inverting a densely-coupled selection/consistency kernel — the deliberately
hard core. 99.95% is solved deterministically and verifiably; the kernel resists
general-purpose SMT and local search. Remaining avenues (multi-bit combinatorial
search, exact GF(p) elimination + backtracking CP, or bit-blasted SAT) are listed
in `RESUME.md`.

## Session 5 addendum — improved to 39,019/39,031 (custom, no SAT/SMT)
- Reverse-engineered the 4 open atoms to 2 primitive sum-gates (27973, 27978) that
  propagation defined out-of-order from a 741-monomial combination atom. Setting
  x_9770/x_3183 to their gate-correct values → **39,019/39,031** (was 39,013).
- The remaining 4 atoms (1817, 30378, 40782, 44271) demand `x_18274=x_9770`,
  `x_17728=x_3183`. Because `x_8821=1`, these are identity-linked chains whose delta
  D=27766… must be absorbed at a **pinned** boundary (product 1816 forces x_26977=0;
  combos 44129/45064 pin x_15690/x_21092). Every local repair (matching, product-
  priority orientation, forward-cone recompute, value-driven augmentation, chain
  propagation, cluster solve — all in solve_lab/) plateaus at 8–12 violated atoms.
- Equivalent obstruction: two combo-computed subtrees must agree —
  `x_23268 = x_6616+x_21092` vs `x_18274 = x_15690−x_26870−x_34150` — differing by D.
  Both combos depend on ~all 256 bits; closing needs a global bit reconfiguration.
- `bit_flip_screen.py`: every single control-bit flip keeps violated atoms ≥6 (baseline
  4); propagation is non-confluent so single flips only relocate the twist. Confirms the
  kernel is quadratic-in-bits (needs a coordinated multi-bit set), matching prior sessions.

## Session 6 addendum — complete reduction of the obstruction
Re-verified the obstruction with an artifact-free model (`confluent_eval5`, mod-P) and
reduced it to a single clean statement.

- **The entire remaining obstruction = a twist match.** Mod-P, the ONLY atoms that ever
  float (for any bit setting) are the twist family; forcing `x_18274:=x_9770` and
  `x_17728:=x_3183` zeros all 4 open atoms including the 741-term cascade 40782
  (`test_40782.py`). So a solution ⇔ control bits with
  `x_9770(A)=x_18274(B)` and `x_3183(A)=x_17728(B)`.
- **Clean decoupling.** `x_9770,x_3183` depend ONLY on 22 control bits (`BITS22`);
  `x_18274,x_17728` depend ONLY on the other 233 bits (same 211-bit support). Disjoint.
- **Shared small denominator.** `x_18274=x_6773/x_8821`, `x_17728=x_17233/x_8821`.
  `x_8821 = 1 − Σ(18 specific bits)` is a *small integer*, exactly linear in the bits;
  the numerators `x_6773,x_17233` are high-degree (deep huge=huge−0 residue chains).
- **Disjoint residue pools.** The 22 bits gate 44 distinct ~290-bit residues; the 233
  bits gate 466; **zero overlap** — the two sides reach the common target through
  different residue arithmetic (the engineered coincidence).
- **Why it is hard.** 22-side is fully enumerable (2^22) and invertible; the 233-side is
  a high-degree map over 2^233 with no linear/low-degree inversion, and the witness must
  additionally lie in the 233-side's *integer-consistent variety* (all division wires
  exact — most B leave it, which is why every local move breaks 30+ atoms). Matching the
  two is a claw-find that needs the setter's trapdoor.

**Deliverable:** `best/best_partial_39019.json` — **39,019 / 39,031** exact in ℤ. The full
solve reduces to inverting the 233-side residue circuit onto a 22-side-reachable target;
that inversion is the open problem. Reusable model + probes in `solve_lab/` (confluent_eval5,
diag, test_40782, test_ratio, deg233, residue_pool, tab22, verify_frame).
