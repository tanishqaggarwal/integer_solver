# Lab Notebook — integer polynomial feasibility (EQUATIONS.txt)

## Session 1 (2026-08-02)

### Instance stats
- 39,031 equations, one per line; vars x_0..x_38747 (38,748 distinct, dense).
- Max line 6.5 KB, max paren-nesting depth 32 → parseable with Python `ast` (x_N are
  valid identifiers). Operators only + - *. 514 distinct huge literals (287-296 bit).

### Structure discovered
- Each equation = `outer * (Σ c_k atom_k)` or `(Σ c_k atom_k)^2`, = 0.
- Flattening the top-level **Add** chain (NOT Sub — Sub is inside an atom `x_t - rhs`)
  and stripping constant coefficients yields the atoms.
- Canonicalized atoms (gcd-reduced, sign-normalized) → **46,275 distinct**, each reused
  ~10× (max 38). ⇒ random-linear-combination design; zeroing every atom solves everything.
- Degrees: 20,090 linear, 25,468 deg-2, 717 deg-4. Booleans: 3,484. Pins: 1,103 (all =1).
- Huge atoms (865): `bit*(x_B-HUGE) - s*x_C`; bit=1 loads residue x_B=HUGE(+s x_C).

### Experiments
| run | method | result |
|-----|--------|--------|
| all-zeros | trivial | 27,352 / 39,031 |
| propagate (pins only) | `propagate.py` | 5,897 vars forced, 0 contradictions |
| propagate + naive zero-fill | forgot to re-propagate NOT gates | 32,316 |
| forward, zero all inputs | `solve_forward.py` | 38,990 (10 atoms viol) |
| **forward, zero only bits** | `solve_forward2.py` | **39,013 / 39,031 (4 atoms viol)** ✅ best |
| z3 on main component (23,843 var / 28,386 con) | `z3_main.py` | `unknown` (>15 min) |
| conflict-cone z3, r=3 | `cone_solve.py` | cone=27k vars → `unknown` |

### Core analysis
- Residual = 1 giant component (23,843 vars, **256 free bits**) + ~297 tiny comps (satisfiable
  by zeros because homogeneous). `analyze_core.py`, `main_component.py`.
- 4 violated atoms are residue-consistency (add-gates linking wires pinned to different
  290-bit residues). Values are computed through a dense cyclic web touching ~255 bits.
- Confirmed NOT a single-modulus reduction: power chains x_36614→x_36614^2 (592b)→^3 (888b)
  grow unreduced; HUGE constants act as **selectable loaded residues**, not a modulus.
- Near-solution has only **1 of 256** main bits set (x_24550, forced) — witness differs.

### Deliverables in place
- `checker.py` exact verifier (compile 10s, eval 0.1s).
- `best/best_partial_39013.json` (+ `failing_eqs.json`: 18 eq indices).
- Atom DB `atoms/poly_atoms.jsonl`; per-phase scripts; STATE.json/RESUME.md.

### Open / next
- `flip_search.py` single-bit flips (→ pair/triple if needed).
- GF(p) linear elimination to find forced bits under consistency; backtracking CP.
- background long z3.

### Session 1 — core probes (definitive)
- `flip_search.py`: 51 of 255 control bits improve violated-atoms 4→3; none reach <3.
  Flips *move* the violated set (e.g. +x_1263 → atoms [8523,25964,44093]) rather than
  removing them — rugged landscape.
- `greedy_combine.py`: greedy multi-flip plateaus at 3 (no 2nd bit improves).
- `enhanced_propagate.py` (Rule A: s·x_C≠0 ⇒ bit=1): cannot bootstrap; all-zeros is a
  fixpoint, witness is a non-trivial fixpoint.
- **Local-repair UNSAT proofs (z3):**
  - Only 4 vars (x_1642,x_4028,x_6236,x_10466) are unique to the 18 failing eqs →
    freeing just those: **UNSAT**.
  - Freeing all 55 vars in the failing eqs + constraining the 5,357 equations that touch
    them (rest fixed to near-solution): **UNSAT**.
  ⇒ The near-solution is a valid 39,013/39,031 assignment but is **not locally
  extendable** to a full solution; the 256-bit core must be solved globally.

### Conclusion (session 1)
99.95% solved deterministically and verified. The residual is a single 256-bit,
densely-coupled, cyclic circuit-inversion kernel with no exploitable modulus/structure;
it resists SMT (z3 unknown), conflict-cone decomposition, boolean-forcing propagation,
and single/greedy bit-flip local search, and local repair is provably UNSAT. A full
solve would require a purpose-built backtracking search / bit-blasted SAT (large,
uncertain) — documented in RESUME.md. Best verified artifact: best/best_partial_39013.json.
