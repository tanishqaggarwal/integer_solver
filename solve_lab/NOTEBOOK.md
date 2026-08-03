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

## Session 2 (continued attack on the 256-bit core)
- **Multiplication density** (`main_component.py` follow-up): the main component has
  **91,080 distinct wire×wire products over 17,042 wires** (10,815 squares). ⇒ full
  bit-blast to SAT is infeasible (~10⁹ gates); the core is genuinely mult-heavy, not a
  small computation hidden behind decoys.
- **Linear consistency** (`linear_consistency.py`, exact GF(2⁶¹−1) Gaussian elim on all
  20,090 linear atoms): system is **consistent** (0 inconsistent rows), rank 19,381.
  199/256 core bits are pivots but **0 core bits are forced to a constant** by the linear
  part ⇒ a full zero-all-atoms solution very likely exists, but the bits are pinned only
  by the dense nonlinear coupling; linear elimination cannot shrink the core (product
  wires dominate the free set).
- **Residue structure**: only 4 residue-sized (~290-bit) distinct values in the near-
  solution; no a+b=c relations among them — no additive shortcut.
- **Pattern tests** (`pattern_test.py`): all-bits-0 (4 viol) remains best; all-bits-1 is
  computationally pathological (unbounded product blow-up, doesn't evaluate).
- Long complete-method z3 lottery (3600 s) left running in background.

### Verdict (session 2)
Full solution almost certainly exists (linear part consistent) but the 256-bit core is a
purpose-built, multiplication-dense obfuscated-circuit inversion with no exploitable
modulus/additive/linear structure. It is intractable by every general method tried
(SMT, propagation, local search, bit-blast, linear elimination). Best verified result
stands at 39,013/39,031. Solving the kernel would require the generator's trapdoor or a
problem-specific breakthrough.

## Session 3 (mod-p / GF(2) attack — "find heuristics")
- **All 256 core bits are genuinely free** (0 forced by pins; x_24550=1 is only a
  conditional artifact of the bits=0 completion, and is the locally-best single choice).
- **x_24550=1 loads TWO residues** (62388 into x_20659, 119182 into x_33718); the other
  active residues (91416, 125787) are *derived*. Violated atom 27973's value is exactly
  119182−91416 (two residues that must be made equal). No a+b=c or multiplicative relation
  among the 4 active residues → no modulus.
- **Boolean reduction** (x²→x for booleans, `bool_reduce.py`): forces 0 *main* bits
  (main-bit constraints stay multi-variable even after reduction).
- **GF(2) attack**: linearization rank 24515/34459 (`gf2_solve.py`); full GF(2) system as
  CNF is SAT in **4 s** (`gf2_sat.py`) but incremental SAT proves **all 256 core bits are
  FREE mod 2** (`gf2_forced.py`). GF(3) linearization even looser (88877 free monomials),
  0 forced. ⇒ every relaxation loses the multiplicative structure that pins the bits.
- **mod-P propagation** (`modp.py`, P=2⁶¹−1): faithful violated-atom proxy (5 vs Z's 4),
  avoids big-int blow-up so many-bit states stay cheap. Enables `modp_search.py`
  (complete pairs+triples over the 81 improving bits) — running.
- **Big-int blow-up** is the search bottleneck in Z (2-bit states can reach 2000+ bits);
  mod-P fixes that but per-eval is still ~1–4 s (Python loop over 46k atoms).

### Verdict (session 3)
Exhaustive: ~20 distinct methods (SMT, SAT, GF(2)/GF(3), linear algebra, boolean reduction,
propagation variants, local/greedy/SA/pairs search, bit-blast analysis, modulus hunts,
local-repair UNSAT). Every relaxation leaves the 256 core bits free; they are pinned only
by exact integer consistency across a 91k-multiplication circuit. No trapdoor found. The
core is, to the best of a broad standard toolkit, intractable without the generator's secret.
Best verified result: 39,013/39,031.

## Session 4 (2-hour deep push — fast eval, accurate search, cascade tests)
- **Fast DAG evaluator** (`fast_walk.py`/`fast_walk2.py`, order-replay, cached inv, 0.14–0.28s):
  ~20x faster but hits alternate all-zero fixpoints inside cyclic sub-components
  (baseline 20 vs propagation's 4) → not accurate enough to validate the witness.
- **Witness is a global set, not incremental**: dense random 76-bit assignment → 320 mod-P
  violations (vs 5 at all-0). Seeded Rule-A (`modp_ruleA.py`): a single seed bit forces
  ONLY itself — no cascade, because products need BOTH operands nonzero (chicken-and-egg).
- **mod-P propagation** (`modp.py`): accurate (5 vs 4), no blow-up; powers `modp_pairs2.py`
  (accurate complete pairs+triples over the 81 improving bits) — running.
- **Huge atoms are exact affine gates** `x_B = s·x_C + HUGE` (not modular reductions) →
  residues are random offsets, no modulus.
- **Bit-blast reconfirmed infeasible**: 122,221 ungated big×big (~290-bit) products in the
  main component.

### Verdict (session 4)
~25 distinct methods now tried. The 256-bit core is pinned only by exact integer consistency
across a mult-dense affine-gated circuit; it has no modulus/additive/linear/GF(p) structure,
no incremental/cascade path, and no tractable relaxation. Search is bounded by ~4s/eval and a
gradientless global landscape. This is an intentionally hard obfuscated-circuit inversion.
Deliverable stands at 39,013/39,031 verified. Accurate pair/triple searches continue as a
completeness exercise (rule out low-Hamming witnesses).

## Session 4 cont. — DPLL/backtracking ruled out
- **Partial assignments give NO early conflicts**: setting 5/20/50/128 (or all 256) control
  bits to decided values and propagating (without zero-filling value-wires) yields **0
  contradictions** every time (~13k–15k vars determined). Contradictions appear only after
  the ~25k free value-wires are all zero-filled. ⇒ conflicts are purely global; a
  CDCL/DPLL search over the bits gets no pruning and degenerates to 2^256. This closes the
  last complete-search avenue.
- Pairs-only accurate mod-P search over the 81 improving bits (3240 pairs) running as the
  final completeness sweep.

## Session 4 cont. — CRITICAL DOF CORRECTION
- The true degrees of freedom are NOT 256 bits. The main component has **256 bits PLUS
  ~4,945 genuinely-free value-inputs** (never a clean target, not huge-atom x_B, undetermined
  even with all bits set; 4,930 feed products). All prior bit-only searches (flip/greedy/SA/
  WalkSAT/pairs) were searching a **256-dim projection** of a ~5,200-dim problem → that is why
  they could not reach the witness.
- The 11 free value-inputs that ARE nonzero in the near-solution are **all residues** (HUGE
  constants): x_18274=91416.., x_3143=62388.., x_5528=119182.., x_12912=125787.., etc. So each
  free value-input takes a value in {0} ∪ {514 residues}.
- Component budget: 23,843 vars, 27,105 atoms → ~5,201 free vars (bits+value-inputs) pinned by
  ~8,463 consistency atoms ⇒ over-determined, unique witness. z3 returns "unknown" (not UNSAT)
  ⇒ likely satisfiable but undecided by SMT on the nonlinearity.
- Implication: the correct search space is a ~5,200-var integer/residue-selection problem, far
  larger than 256 bits. Still over-determined and mult-dense; no tractable method found, but the
  framing (residue-or-0 per value-input) is the accurate one for any future attack.
