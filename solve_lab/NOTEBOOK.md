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
