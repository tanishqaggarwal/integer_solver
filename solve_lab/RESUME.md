# RESUME — read this first

## Status: 39,019 / 39,031 equations satisfied exactly in Z (99.97%). 12 equations open.
## (was 39,013; Session 5 improved it by fixing gates 27973/27978 — see NOTEBOOK Session 5.)
## Best file: best/best_partial_39019.json. Remaining 4 atoms: 1817, 30378, 40782, 44271.
## Remaining obstruction is GLOBAL: two circuit subtrees must produce equal values
## (x_23268 = x_6616+x_21092  vs  x_18274 = x_15690−x_26870−x_34150, tied by check gate 1817);
## best differs by D=27766…; closing it needs a global bit reconfiguration, not a local repair.

## What this instance actually is
`EQUATIONS.txt` (39,031 equations, vars x_0..x_38747) is an **obfuscated arithmetic
circuit**. Every equation is `outer_scalar * (Σ c_k · atom_k)` or `(Σ c_k·atom_k)^2`,
set = 0, where the `atom_k` are **shared gate residuals** reused ~10× across equations
(46,275 distinct atoms). Any assignment that makes **every atom = 0** satisfies every
equation — that is the intended solution.

Atom vocabulary (all must vanish):
- add/sub/copy/scalar gates: `x_a-(x_b+x_c)`, `x_a-(x_b-x_c)`, `x_a-s*x_b`  (linear)
- multiply / square gates:   `x_a-x_b*x_c`, `x_a-x_b*x_b`
- NOT gates:                 `x_a-(1-x_b)`
- boolean gates:             `x_a*(x_a-1)` → 3,484 bit variables
- constant pins:             1,103 vars pinned to **1**
- huge (bit-gated) atoms:    `bit*(x_B - HUGE) - s*x_C`, HUGE ≈ 287-296 bits (514 distinct).
  bit=1 ⇒ x_B = HUGE + s*x_C (loads a residue); bit=0 ⇒ x_C = 0.

## How the current solution was built  (deterministic, ~6 s)
1. Integer propagation from the 1,103 pins → 5,897 vars forced (no choices).  `propagate.py`
2. Set the 1,156 free **boolean** inputs to 0; propagate (this SOLVES the value-wires
   x_B via the huge atoms).  Then zero-fill remaining. `solve_forward2.py`
3. Result `cand_forward2.json` == `best/best_partial_39013.json`: 39,013/39,031.

## The open core (why 18 remain)
Residual after propagation = **one giant component of 23,843 vars with 256 free bits**
(+ ~297 tiny components already satisfied by zeros). The 18 failing equations come from
**4 unsatisfied atoms** (27973, 27978, 41470, 45004) — residue-consistency constraints
(e.g. add-gate `x_9770 = x_35186 + x_3368` where the two sides are pinned to *different*
290-bit residues). The 4 atoms' backward cone touches ~255 of the 256 bits (dense, cyclic),
so there is no small local fix. z3 on the component (28k nonlinear int constraints) returns
`unknown`; conflict-cone z3 at radius 3 engulfs the component and also fails.

## Highest-EV next experiments (in order)
1. `flip_search.py` — single control-bit flips, re-propagate, count violations (running/last run
   → `flip_results.json`). If any single/￼pair flip drops violations to 0 → full solve. Extend to
   pair/triple flips over the ~255 control bits if singles don't close it.
2. Exact linear solve over GF(p) of all 20,090 linear atoms to see which bits become **forced**
   once residue-consistency is imposed (constraint propagation WITH the 4 consistency atoms as
   drivers + backtracking).
3. Long z3 with `qfnra-nlsat`/bounds as a lottery ticket (background).

## Verify anything
    python3 solve_lab/checker.py solve_lab/best/best_partial_39013.json
(loads assignment, exact integer re-eval of all 39,031 eqs; missing vars→0.)

## Do NOT redo
- Parsing / atom extraction (atoms/poly_atoms.jsonl is canonical, gcd/sign-normalized).
- The modulus hunt (gcd of the 514 huge constants = 1; the huge power-chains x^2,x^3 are
  NOT reduced — this is not a single-modulus reduction).
- Treating it as a named crypto instance — work the circuit.
