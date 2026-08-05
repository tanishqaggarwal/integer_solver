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

## Session 4 cont. — UNSAT-CORE localization (valid progress) + z3 SAT-search blocker
- **z3 with all "free value-inputs" = 0 (bits free, everything else free) is UNSAT in 25-30s.**
  The **unsat core is only 4 value-inputs**: {x_13467, x_13780, x_29071, x_29561}. So the
  witness must set at least one of these nonzero. (Caveat: the "free value-inputs" set (4945)
  is overcounted — my clean-target detector missed sum-form definitions like
  `x_29071 = x_25916 + x_28986`, so several are determinable wires, not true inputs. But the
  UNSAT is real: z3 is free to adjust ALL other wires and still cannot satisfy the system with
  those vars pinned to 0.)
- **Iterative implicit-hitting-set** (`z3_iter.py`/`z3_iter2.py`/`z3_subst.py`): free the core
  vars, re-solve; grow via new cores until SAT. Sound in principle.
- **BLOCKER**: z3 finds UNSAT cores fast but **cannot construct a SAT model** — every SAT-check
  iteration (freeing even 4 value-inputs, residue domain or unbounded, substituted or not)
  runs >10 min without deciding. The 122k big×big products make model *construction*
  intractable for z3, even though refutation is fast. This is the crux limit: refutation-easy,
  satisfaction-hard.

### Standing conclusion
The instance is a large nonlinear integer system (obfuscated mult-dense arithmetic circuit).
z3 can refute constrained sub-instances quickly but cannot synthesize a satisfying assignment
for the residual; no relaxation, decomposition, incremental path, modulus, or search has closed
it across ~30 distinct methods. The unsat-core localization is the most promising lead for any
future work (it shrinks "which inputs matter"), but requires a SAT-model synthesizer that
handles the product structure (a purpose-built one, not general SMT). Verified deliverable
remains 39,013/39,031.

## Session 5 — custom orientation analysis → **39,019/39,031** (improved from 39,013)
Directive: no SAT/SMT; design own heuristics from the circuit structure.

### Root-cause of the 4 violated atoms (fully reverse-engineered)
- The 4 violated atoms reduce to **2 primitive gate residuals**: 27973 `x_9770 = x_35186 + x_3368`
  and 27978 `x_3183 = x_1642 + x_10466`. Atoms 41470 (= −1·Σ residuals, contains 27973) and
  45004 (= (Σ residuals)², contains 27978) are **redundant combinations** — they vanish once the
  primitives hold.
- Every *other* gate feeding these is already consistent (verified 26/28 local gates OK). Only
  `x_9770` and `x_3183` carry wrong values: `propagate` defined them from the 741-monomial
  **combination** atom 40782 (and 30378) in the wrong topological order, *before* their true
  sum-gates could fire — so 27973/27978 were left unmatched.
- **Fix that works**: set `x_9770 = x_35186+x_3368 = 119182…`, `x_3183 = x_1642+x_10466 = 62388…`.
  This satisfies 27973/27978/41470/45004. Net **39,019/39,031** (was 39,013). New best saved:
  `best/best_partial_39019.json`; failing eqs → `best/failing_eqs_39019.json` (12 lines).

### Why it can't (locally) go to 39,031 — the twist
- After the fix, exactly **4 atoms** remain: 1817 `x_18274 = x_9770`, 30378/44271 `x_17728 = x_3183`,
  and the compensation combo 40782. Because `x_8821 = 1`, the chains
  `x_9770 = x_18274 = x_6773` and `x_3183 = x_17728 = x_17233` are **identity-linked** and must
  all move together to the new values.
- But `x_18274`/`x_17728` are heavy fan-out inputs. Moving them forces `x_26517 = x_6773+x_34150`,
  `x_15690 = x_26517+x_26870`, … and the delta D = 27766… must be absorbed at a boundary that is
  **pinned**: `x_26977=0` (product 1816, since `x_20510=0`), `x_26870=x_6283·x_16160` (product),
  `x_15690`/`x_21092` (combos 44129/45064). D cannot be absorbed locally → it must propagate into
  those combos' inputs, i.e., a **global bit reconfiguration**.
- Equivalent framing: two subtrees must agree — `x_23268 = x_6616+x_21092` (=119182…) must equal
  `x_18274 = x_15690−x_26870−x_34150` (=91416…), mediated by check gate 1817. best makes them
  differ by D. Reaching 39,031 = finding bits that make the two subtrees equal.
- Confirmed local-optimum: every single identity-group move (`{18274,6773}`, `{17728,17233}`,
  both) *increases* failures (→38,990/38,996/39,000). 39,019 is a strict local max under local moves.

### Custom methods built this session (all no-SAT/SMT; scripts in solve_lab/)
- `prim_solve.py` primitives-only re-derivation; `match_solve.py` max bipartite matching
  (gate↔var); `full_forward.py`/`override_iterate.py` product-priority matching + iterate-to-
  fixpoint; `cone_fix.py`/`cone_fix2.py` surgical forward-cone recompute; `augment_repair.py`
  value-driven augmenting re-orientation; `chain_prop.py` delta chain-propagation;
  `cluster_analyze.py`/`cluster_solve.py` bounded-cluster extraction. All plateau at 8–12 violated
  atoms because the correction ripples past pinned product/combo boundaries — consistent with the
  "global bit reconfiguration required" conclusion. Verified deliverable now **39,019/39,031**.

## Session 5 cont. — CONFLUENT EVALUATOR built + validated (major tool)
- Built a deterministic forward evaluator (`confluent_eval4.py`) = best's prov orientation
  + override x_9770<-27973, x_3183<-27978, + **residue-load injection** for huge-atoms
  (x_B = bit*(HUGE + s*x_C); best's prov leaves x_B=0 so naive eval never loads). Uses best's
  acyclic eval-order with 9770/3183 moved to the end (they were computed BEFORE their input
  x_35186 — a real ordering bug). VALIDATED: `forward_Z([])` exactly reproduces 39,019
  (atoms 1817/30378/40782/44271, 0 non-divisible). This is confluent (no order-noise), unlike
  propagation — a bit's effect is now its TRUE semantic effect.
- With it, exactly **22 bits change the twist checks singly** — identical to the earlier BITS22.
- `confluent_enum.py`: vectorized (numpy, chunked, 31-bit prime) enumeration of all **2^22**
  patterns over those 22 bits. Result: **0 patterns zero both checks** (x_18274=x_9770 AND
  x_17728=x_3183). Trustworthy (evaluator validated). ⇒ the witness needs a bit OUTSIDE the 22.
- Backward cone of the check vars reaches ALL 255 control bits (via loads), but 233 of them have
  ONLY product (pair+) effects — quadratic-in-bits. So the twist = 2 multilinear equations over
  255 bits; 22 appear linearly, no linear-only solution. Next: search pairs/triples that activate
  the right product (best is all-0 bits & 39,019, so the true 1-bit set is likely small).

## Session 5 cont. — exhaustive bit-subset search + degree finding
Using the validated confluent evaluator (`pattern_search.py`, vectorized numpy over the
checks' 19,339-wire backward cone, two 31-bit primes), searched for a bit set zeroing both
twist checks:
- **2^22** over the 22 linear-effect bits → 0.  **all C(255,2) pairs** → 0.
- **triples with 2 of the 22** (55k) → 0.  **triples 1-of-22 + 2 outside** (595k) → 0.
- best has all 255 control bits = 0 (verified), so "set bits to 1" is the correct search.
- The free value-input DOF (e.g. x_31302 for x_26977 in check 1817) is blocked: x_20510 is
  pinned to 0 by ~65 atoms, so check 1817 genuinely requires x_18274 = x_9770.
- **Degree test**: the checks are NOT quadratic in the bits — some 22-bit triples deviate from
  the quadratic prediction ⇒ degree ≥3 multilinear. So neither brute-force (≤3 bits) nor a
  quadratic-form solve reaches the witness. This is high-degree circuit inversion = the intended
  hard core. Witness = a specific multi-bit (≥4) pattern; the search space over the 255 relevant
  bits with degree-≥3 coupling is not brute-forceable and not low-degree-solvable.

### Session 5 net
39,013 → **39,019 / 39,031** (verified in Z). Built + validated a confluent forward evaluator
(the correct circuit model incl. residue loading), pinned the obstruction to 2 high-degree
multilinear check equations over 255 bits, and exhausted brute-force up to triples. The witness
requires inverting the obfuscated selection kernel — unsolved, as across all prior sessions.
- **quad3in22** (3-of-22 + 1 outside, 359k) → 0 hits. All brute-force up to quads exhausted; witness is >=4 bits in a degree->=3 system over 255 bits (circuit inversion). Unsolved.

## Session 5 cont. — EVALUATOR BUG FOUND (v4→v5): x_18274 was frozen
- confluent_eval4 only solved LINEAR-output gates; it silently FROZE 1504 division-oriented
  wires whose defining gate puts the target inside a product (x_18274 <- 4954:
  x_6773 = x_8821*x_18274 ⇒ x_18274 = x_6773/x_8821). So x_18274/x_17728 appeared IMMOVABLE and
  the 2^22 / pair / triple 0-hit results were **artifacts** — the wrong function was searched.
- `confluent_eval5.py` adds division handling (v = -rest/(c*u)) + load + gate, validated
  (forward_Z([]) == 39,019). Now **x_18274 is moved by 232 single bits, x_17728 by 232**
  (was 0). The twist DECOUPLES: x_9770/x_3183 ← the 22 bits; x_18274/x_17728 ← ~232 other bits
  (largely disjoint). All prior bit-searches must be re-run with v5. This is the most promising
  lead of the session — the check x_9770 = x_18274 now has BOTH sides bit-movable.

## Session 5 cont. — v5 search: twist pairs found, but cascade obstruction
- With the fixed v5 evaluator, **all C(255,2) pairs** search: 3798 pairs zero the 3 twist checks
  (survive both primes) — of form (232-mover bit, 22-bit): one sets x_18274, the other x_9770 to
  match. So the twist IS bit-satisfiable (v4's 0-hits were the frozen-x_18274 artifact).
- BUT every such pair BREAKS 32–43 OTHER atoms (Z-verified): the 232-mover's residue load
  propagates and violates downstream gates. Via full propagation the best pair still gives 20
  violated atoms (vs baseline 4). So fixing the twist cascades into a new, larger twist.
- Key structural fact: **all-0 bits is the propagation MINIMUM (4 atoms = 39,019)**; every bit or
  pair set (via propagation or v5) strictly increases violations. The true 0-violation input is a
  specific bit pattern NOT reachable by local search/propagation from all-0 — this is the
  deliberately hard circuit-inversion core. Each activated bit must be balanced by its consistent
  partners across the whole 256-bit input; the cascade only terminates at the true input.
- Net: v5 corrected the model and proved the twist is bit-satisfiable in isolation, but a full
  witness needs the globally-consistent 256-bit input (circuit inversion), unsolved.

## Session 5 cont. — purpose-built searches (v5-based); witness confirmed isolated
- `greedy_v5.py`: parallel greedy constraint-repair over bits using the CORRECT v5 model
  (violations = real consistency-check failures). Seeded from a twist-pair (27 viol): **every**
  single-bit addition INCREASES violations (min 31) → strict local minimum. Confirms the witness
  is isolated: no descent path from the twist-fixing neighborhood.
- Decomposition-match tried: `x_18274` / `x_17728` are NOT linear in the 232 bits (specific pairs
  deviate) — quadratic+ coupling, so the twist-match is quadratic feasibility, not subset-sum.
- Orientation fix tried: `confluent_eval6.py` forces product-forward orientation to kill the
  cascade (products define outputs). It does NOT reproduce 39,019 (16 violations, worse than v5's
  4) and is cyclic — best's orientation is closer to the true circuit than naive product-forward,
  and no clean orientation is both all-0-consistent and cascade-free.
- Conclusion (with the corrected v5 model): the twist is bit-satisfiable but every fix cascades;
  the true 0-violation input is an ISOLATED point (all-0 is the local min, high-degree multilinear
  coupling, no linear/quadratic reduction, no descent path). Finding it = inverting the obfuscated
  selection kernel — the deliberately hard core. Deliverable stands at 39,019/39,031.

## Session 5 cont. — LINEAR-ALGEBRA BREAKTHROUGH: 232-part is slaved to the 22-part
- The twist residual is LINEAR in the 232 x_18274-side bits (verified: pair deltas add exactly).
- `linalg_attack.py`: computed each atom's GF(P) response to the 232 bits; 45,267/46,275 atoms are
  linear in them. Built the linear system from linear atoms only -> **rank 233/233, 0 inconsistent
  rows**. So over the 232 bits the consistency atoms UNIQUELY DETERMINE the 232-part given the
  22-part. At 22=0 the unique solution is 232=all-0 (=best, twist still violated).
- Consequence: the 232 control bits are NOT free — they are slaved (B(A) = M^{-1}·(-base(A))) to
  the 22 nonlinear bits. The real degrees of freedom = **just the 22 bits**. My earlier 2^22
  enumeration used 232=0 (wrong) -> found nothing; the correct enumeration must use 232=B(A).
- Next: enumerate 2^22 patterns A over the 22 bits; for each solve the linear system for B(A) and
  require B(A) in {0,1}^233 (strong filter — only the true A gives a 0/1 solution), then verify.

## Session 6 — CLEAN REDUCTION: entire obstruction = the twist match (2 equations)
Corrects several earlier lossy-forward-eval conclusions. All facts below re-verified
with confluent_eval5 (Z and mod-P). Scripts: diag.py, twist_struct.py, test_lin.py,
test_additive.py, extract_huge.py, test_40782.py, test_ratio.py, deg233.py, tab22.py.

### The 4 obstruction atoms (at best=all-0), exact structure
- atom 1817 : 6033033*x_9770 - 6033033*x_18274 + x_26977 = 0   (x_26977 identically 0)
- atom 44271: x_3183 - x_17728 = 0
- atom 30378: x_3183 - x_9982 - x_17728 = 0                      (x_9982 identically 0)
- atom 40782: big cascade check, 52 vars, deg 2-4, 741 terms, NO control bits.

### atom 40782 is IMPLIED by the twist (test_40782.py, decisive)
Forcing x_18274:=x_9770 and x_17728:=x_3183 in the all-0 state makes ALL FOUR
residuals (1817,30378,40782,44271) become exactly 0. So the ENTIRE remaining
obstruction is exactly:
        x_9770(A) = x_18274(B)   AND   x_3183(A) = x_17728(B)
Nothing else. (40782's 52 vars: 6 only-22-side, 2 only-233-side, 44 constants;
two of them, x_24252 & x_36641, are fixed ~250-digit constants, x_25471=1.)

### Clean decoupling (diag.py / test_lin.py)
- x_9770, x_3183 : moved ONLY by the 22 bits (19 and 21 of them). 0 of the 233.
- x_18274, x_17728: moved ONLY by the 233 bits (211 each; SAME 211-bit support). 0 of 22.
The two sides are functions of DISJOINT bit sets, coupled only by the 2 equations.
BITS22 = [1782,1858,2795,2800,3483,5443,10652,19520,21188,21588,23634,26947,27512,
          29682,30104,30596,30658,30792,33251,37748,37885,38116]; bits233 = control - BITS22.

### Denominator structure (test_ratio.py / deg233.py) — KEY
- x_18274 = x_6773 / x_8821 ;  x_17728 = x_17233 / x_8821   (SHARED denominator x_8821)
- x_8821 is EXACTLY LINEAR in the 233 bits (mod-P Mobius deg2=deg3=deg4 = 0). A subset-sum.
- x_6773, x_17233 (numerators) are HIGH-degree in the 233 bits.
- So the targets are: x_6773 = N1*x_8821 and x_17233 = N2*x_8821 with N1=x_9770(A),
  N2=x_3183(A). Numerators high-degree => not linearizable directly.

### Degrees
- x_9770, x_3183 : high-degree multilinear in the 22 bits (deg>=4). The 22-side is
  fully ENUMERABLE (2^22 = 4.2M) via its cone (7520 wires, 258 div).
- x_18274/x_17728: high-degree in the 233 bits; 2^233 NOT enumerable and no
  linear/low-degree inversion found.

### B=0 (all 233 bits off) essentially ruled out
best_partial_39019 already has ALL 255 control bits = 0. Solving would need
x_9770(A)=x_18274(0)=91416258160755509149180373473728639746431157665678710450404458852172057265575180278101002
x_3183(A)=x_17728(0)=125787314747601108116039725163361763116550465675981151838811516827327919228823597744635626
via the 22 bits alone. First 64k patterns: 0 single-coord matches mod two primes
(expected ~0 for a 290-bit target over 2^22). Full 2^22 scan launched (runs/tab22_full.log)
to settle this and to materialise the complete achievable set S = {(x_9770(A),x_3183(A))}.

### Correction to Session-5 "232-part slaved / rank 233" claim
That used a linearity filter built on the LOSSY integer forward-eval; the wires are in
fact genuinely nonlinear, so that reduction (and "B forced to 0") is NOT reliable. The
robust facts are those above, re-checked mod-P (artifact-free).

### Where it stands
A crisp 2-equation match between a fully-enumerable 22-bit side and a high-degree
233-bit side sharing linear denominator x_8821. The 233->target inversion is the
deliberately-hard trapdoor. Deliverable stands at 39,019/39,031 (verified exact in Z).
New scripts this session: diag.py twist_struct.py test_lin.py test_additive.py test_linP.py
trace_cascade.py extract_huge.py test_40782.py test_ratio.py deg233.py degree_probe.py tab22.py
Next tracks: (a) full S table + structure of S (moduli/factors); (b) residue-pool identity
between the two sides; (c) MITM/lattice via the x_8821 linear coordinate.

### Session 6 caveat (verify_frame.py) — reduction is mod-P / at consistent points
For RANDOM bit settings the integer forward-eval floats ~1132 different atoms (150-423 per
test), NOT just the twist. Reason: Z division wires x_v=num/den go non-integer for arbitrary
B, so forward-eval leaves a stale value and that atom breaks. Consequences:
- MOD-P the reduction is exact: mod-P forward-eval always fires (division invertible), so the
  ONLY floating atoms are the twist family — for ANY B. So "twist match mod P" <=> "all atoms 0
  mod P". The Session-6 structural facts (decoupling, x_8821 linear, denominators) are all mod-P
  and remain valid.
- In Z a full solution needs BOTH (i) the twist match AND (ii) every division wire on the 233-
  side exact (internal integer-consistency, which all-0 satisfies but most B break). test_40782's
  "twist ⇒ all 4 zero" was evaluated at the all-0-consistent state, so it shows the twist is the
  binding constraint GIVEN internal consistency — not that any twist-matching B solves in Z.
- Net: witness B* lies in the "integer-consistent variety" of the 233-side AND hits a 22-side
  target. That variety (+ high-degree map + disjoint 290-bit residue pools) is the trapdoor.
  This is why every local bit move breaks 30+ atoms (it leaves the consistent variety).

### Session 6 — residue/modulus probe (residue_probe.py): no algebraic shortcut
- gcd(all 510 residues)=1; gcd(D1, all residues)=1; gcd(D1,D2)=1 where D1=x_9770-x_18274
  (294 bit), D2=x_3183-x_17728 (295 bit) are the two twist gaps at all-0.
- No residue divides D1; largest gcd(D1, any residue)=39. Small-prime signatures of the
  four targets look unstructured. => NO hidden single modulus, no clean residue-lattice
  relation between the two sides. The 290-bit residues are effectively random with no shared
  modulus — a genuine (cryptographic-strength) trapdoor.
- DEFINITIVE (tab22 full, prime 2147483647, all 2^22): 0 patterns match x_18274(0) or
  x_17728(0) in even one coordinate => B=0 impossible; the witness needs the 233-side active.

### Session 6 — CORRECTED linear-algebra attack over all 255 bits (linalg255.py) — DECISIVE
Session-5's linalg fixed the 22 bits to 0 (=> rank 233, B=0). The corrected version treats
ALL 255 control bits as variables, uses the artifact-free mod-P forward-eval oracle (a valid
global acyclic orientation), keeps atoms linear in the bits (2-config filter), Gaussian-
eliminates over GF(P):
  atoms linear in 255 bits: 44980/46275; RANK **255/255**; free dim **0**; 0 inconsistent;
  unique solution = **all 255 bits = 0**.
Interpretation: the linearization of the circuit AROUND all-0 uniquely forces the trivial
point. But the true witness has bits != 0 (B=0 already ruled out by the full 2^22 scan). So the
witness must make some locally-linear atoms behave NONLINEARLY — i.e. it lies OUTSIDE the linear
neighborhood of all-0. **Linear algebra / lattice provably cannot reach the witness**; all-0 is
an isolated linear-consistent point and the solution is in the genuinely nonlinear region.
This is the clean, corrected version of (and supersedes) the Session-5 "slaved-B" claim.

## Session 6 — FINAL STATUS
Every attack class is now exhausted with a definitive negative:
 - enumeration: 22-side yes (2^22), 233-side no (2^233); claw-find needs ~2^558 evals.
 - B=0: ruled out (full 2^22 scan, 0 matches).
 - linear algebra over 255 bits: forces all-0 (rank 255) — witness is nonlinear-region.
 - lattice/subset-sum: numerators high-degree, no linear reduction; no residue-lattice relation.
 - modulus: gcd(all residues)=1, no hidden modulus.
 - slack: x_26977, x_9982 circuit-defined/pinned — twist gap cannot be absorbed.
 - local/greedy/SA/pairs/triples: all-0 is the isolated local min.
The instance is a genuine (cryptographic-strength) obfuscated-circuit trapdoor: a full witness
needs the setter's secret or a cryptanalytic break of the specific 233-side residue circuit.
DELIVERABLE: best/best_partial_39019.json = 39,019/39,031 exact in ℤ (checker-verified).

### Session 6 — IMPORTANT self-correction on forward-eval validity
Tested (mod-P): setting ANY control bit — 22-side OR 233-side — breaks 24-300 non-twist
atoms, not just the twist. all-0 is the ONLY point where the v5 orientation floats just the
4 twist checks. Consequence for honesty:
- The v5 confluent forward-eval is a valid oracle ONLY at (and infinitesimally near) all-0.
  For any nonzero bit pattern it satisfies its ORIENTED atoms but violates many CHECK atoms,
  i.e. it produces an assignment that is NOT a consistent circuit state.
- Therefore x_9770(A), x_18274(B) computed by forward-eval for nonzero A/B are HEURISTIC
  values of a local model, NOT the true circuit values. The Session-6 "enumerable 22-side",
  "decoupling", and "claw-find" framing describes this local model near all-0, and does NOT
  by itself give the global witness.
- What remains SOLID and verified: (i) best_partial_39019 satisfies 39,019/39,031 exactly;
  (ii) at all-0 exactly 4 atoms fail, with exact definitions
       1817: 6033033*(x_9770-x_18274)+x_26977;  44271: x_3183-x_17728;
       30378: x_3183-x_9982-x_17728;  40782: 741-term cascade (implied by the other 3 at all-0);
  (iii) x_26977=x_20510*x_31302 (atom 1816), x_9982 pinned 0 (atom 44272) — the obstruction is
       the rigid twist relation; (iv) all-0 is an isolated locally-consistent point (every bit
       move leaves it); (v) the 290-bit residues have no shared modulus.
- TRUE problem = full nonlinear CSP over the 46,275 atoms (obfuscated-circuit inversion). The
  forward-eval cannot represent the witness even at the correct bits (wrong orientation there),
  so the witness needs a real circuit solve, not this heuristic. A better global orientation
  (one that stays consistent under bit changes) would give a valid oracle but still leaves the
  2^233 claw-find. Net: genuine hard trapdoor; deliverable stands at 39,019/39,031.

### Session 6 — localizing the obstruction via full LINEAR elimination (linsolve_full/twist_core)
Corrected oracle understanding: `prop_oracle.py` shows CORRECT mod-P propagation leaves
x_9770/x_18274 UNDETERMINED (they are free vars of a ~23.8k-var coupled core); the v5
forward-eval was one heuristic filling. So the real object is that nonlinear core.
- Full GF(P) Gaussian on all 20090 linear atoms (`linsolve_full.py`): linear RANK 19381,
  0 inconsistent; ~19367 vars free w.r.t. the linear system (the residue VALUES, which the
  nonlinear load/product atoms determine). So the true DOF are the nonlinear atoms.
- `twist_core.py`: the twist vars reduce (mod P) to a FEW core wires:
    x_18274 ~ combo{x_31434, x_34236, x_35846};  x_9770 ~ same 3 + x_26977;
    x_17728 ~ combo{x_27912, x_28035};           x_3183 ~ same 2 + {x_6236, x_10466}.
  (mod-P coeffs carry 1/3 factors — GF(P) artifacts, not the Z relation; the SUPPORT is the
  takeaway.)
- Tracing those core wires to their defining atoms exposes the ESSENTIAL NONLINEARITY:
    x_35846^2 = x_3002   (square!)      x_28035^2 = x_36228  (square!)
    x_31434  = x_7101 / x_28035  (product/div)   x_26977 = x_20510*x_31302 (product)
    x_34236  = x_12293 + x_25804 (sum)   x_6236/x_10466/x_27912 = products.
  => the obstruction is a system of QUADRATIC / SQUARE-ROOT constraints over ~290-bit
  integers (x=+-sqrt(HUGE)), exactly the "huge power chains x,x^2 grow unreduced" note.
  This is the trapdoor kernel: matching two subtrees each built from squares/products of
  ~290-bit residues. Genuinely hard; no linear/lattice handle (confirmed) — needs solving
  the quadratic system over Z (setter's trapdoor or heavy cryptanalysis).

### Session 6 — perfect-square structure of the verifier checks (check_square.py)
- 530 of 5936 degree>=2 atoms are PERFECT SQUARES: 458 = (degree-2 form)^2, 72 = (linear form)^2.
  So those checks are Q^2=0 <=> Q=0, halving their degree. atom 40782 = Q_40782^2 (Q deg-2, 38
  terms) and 41285 = Q_41285^2 (Q deg-2, 42 terms, satisfied at best).
- Q_40782 combines BOTH twist cores: -6033033*(x_9770-x_18274) + 39*(x_3183-x_17728) - x_26977
  - 42*x_9982 + [18 product terms] + [linear terms]. It's ONE combined (weaker) equation, not a
  replacement for the two separate twist atoms 1817/44271.
- The 72 (linear)^2 hidden constraints add ZERO rank (redundant with the existing linear system);
  no new linkage between the 22-side and 233-side. (linsolve_plus.py)
- Net: nice structural insight (the verifier is ~9% squared forms) but no exploitable reduction of
  the core twist obstruction. Deliverable unchanged at 39,019.
- Background campaigns launched: sa_campaign.py (SA over 255 bits), mitm_lowB.py weight-3 (all-A x
  sparse-B meet-in-the-middle) — both running as good-faith long searches.

### Session 6 — MITM weight-3 result: DEFINITIVE negative
mitm_lowB.py weight-3: swept all 2,108,418 weight-<=3 233-side B patterns against the full
2^22 22-side hash (both x_9770,x_3183 coords, prime 2147483647), skipping degenerate zero
collisions. **0 nonzero hits.** So there is NO witness with a sparse (<=3-bit) 233-side residue
selection. Combined with the full-2^22 B=0 scan (0 hits) and prior <=3 total-weight searches,
the witness's 233-side selection is genuinely NON-SPARSE (>=4 residue bits) and unreachable by
feasible enumeration. SA campaign plateaued at 4 (all-0 is the isolated violation minimum).
FINAL: every search/structural avenue exhausted; genuine trapdoor; deliverable 39,019/39,031.

### Session 6 — TRAPDOOR MECHANISM fully reverse-engineered (the key result)
Why every forward-eval search (SA, mitm, greedy, pairs, enum22, local) was DOOMED, and what the
witness actually requires:

1. QUANTIZATION (codewords.py, quant_structure.py): under the confluent forward-eval, each twist
   wire is quantized to integer multiples of its all-0 value:
     x_9770 = m*g    (g=119182..., 296b, m in 27 vals {-9..27})
     x_3183 = m'*h   (h=62388...,  295b, m' in 45 vals)
     x_18274= m2*g2  (g2=91416..., 296b)
     x_17728= m2'*h2 (h2=125787..., 296b)
   with gcd(g,g2)=1 and gcd(h,h2)=2. The RIGID twist (x_9770=x_18274, x_3183=x_17728) then forces
   ALL FOUR to 0 (coprime units + small multipliers) -> only the degenerate (0,0,0,0), which fails
   ~300 other atoms. Since the instance is FEASIBLE, the forward-eval family CANNOT contain the
   witness. That is the root cause of every search plateau.

2. THE SLACKS. The wire DEFINITIONS carry product-slack terms the forward-eval zeros:
     a27973: x_9770  = x_35186 + x_3368 ,  x_3368  = x_12779 * x_24026   (a1660)
     a27978: x_3183  = x_1642  + x_10466,  x_10466 = x_12779 * x_27116   (a27976)
     a1817 : ... + x_26977 ,               x_26977 = x_20510 * x_31302   (a1816)
     a30378: ... - x_9982  ,               x_9982  = x_9897  * x_12518   (a1818)
   x_35186=m*g (the quantized part), x_3368 the slack. Both x_3368 and x_10466 are GATED by the
   single wire x_12779 = x_23380*x_36336 (a1652), which is 0 at best (x_36336=0).

3. THE WITNESS activates x_12779 (=> x_36336!=0 => a cascade down to specific control bits) so the
   slacks become nonzero and BRIDGE the coprime-quantization gap:
     x_3368  = x_18274 - x_35186 = m2*g2 - m*g
     x_10466 = x_17728 - x_1642  = m2'*h2 - m'*h
   i.e. the twist is satisfiable at a NONZERO value only by turning on the slack products. The
   forward-eval, which orients these products to 0, structurally cannot reach that state.

4. CONSEQUENCE / PATH: the true solve = activate the slack cascade (x_12779 -> x_36336 -> ... ->
   bits) AND make x_24026, x_27116 (etc.) hit the ~296-bit gap values (m2*g2-m*g)/const. This is a
   coupled product-chain inversion over the 255 bits -- the deliberately-hard trapdoor core. A
   solver must abandon the all-0 forward-eval orientation and drive these products nonzero.
   Deliverable remains 39,019/39,031; this is a complete mechanistic reverse-engineering of the kernel.

### Session 6 — slack activation is REACHABLE (resolves the feasibility paradox)
- x_12779 = x_23380*x_36336 is activated by 22-side bits: single flips of {1858,2795,5443,10652,
  19520,26947,27512,30104,...} move x_12779/x_36336/x_38073/x_14402; bit-PAIRS give x_12779=2.
- But the slack x_3368 = x_12779*x_24026 also needs x_24026 != 0 (moved by 0 single bits; activated
  deeper via x_38215 through a1813: x_14402*x_24026 = 321447*x_38215). x_27116 similarly.
- So the TRUE 22-side value is x_9770 = m*g + x_12779*x_24026 (quantized part + slack product); with
  the slack ON it ranges FAR beyond the 27 forward-eval values, and CAN equal x_18274=m2*g2. The
  coprime-quantization "impossibility" only holds with slacks OFF (the forward-eval regime).
- The full solve therefore = drive the 22-side bits to (a) activate x_12779 and x_24026/x_27116 and
  (b) make m*g + x_12779*x_24026 = x_18274(B) (and the h-version), i.e. a coupled bilinear/product
  match. This is the precise trapdoor kernel, now fully mapped. A purpose-built solver must search
  WITH the slacks active (the confluent forward-eval, which zeros them, cannot).

### Session 6 — the bridging slacks are BURIED DEEP (final piece)
x_24026 and x_38215 (the big slack wires that bridge the coprime-quantization gap) are activated by
ZERO 22-side pairs (0/231), zero 233-side pairs (sample 780), zero cross 22x233 pairs, and zero
22-side triples. They sit behind the full residue-load cascade and only turn nonzero for a deeply
coordinated (witness-level) bit pattern. So:
- slacks OFF  (any shallow/forward-eval config) => coprime quantization => twist degenerate only.
- slacks ON   requires the deep witness pattern (can't be reached by low-weight perturbation).
This is the crux of the trapdoor's hardness: the ONLY states that satisfy the twist non-trivially
require activating slack products buried behind the residue cascade, i.e. the setter's exact input.
No shallow search (any weight <=4, any single/pair/triple, any forward-eval point) can reach them.
COMPLETE mechanistic reverse-engineering; the residue-cascade inversion to activate x_24026 is the
irreducible trapdoor. Deliverable 39,019/39,031.

### Session 6 — solver attempts past the trapdoor (merge_solve.py)
- MERGE twist vars (x_18274->x_9770, x_17728->x_3183) bakes the twist into the atom set. The merged
  system is CONSISTENT under propagation (0 contradictions) - so the twist is not self-contradictory.
- But: unit propagation determines only the 5898 pin-forced vars, 0 control bits. Full GF(P) Gaussian
  on the merged linear atoms: rank 19380, 198 control bits are pivots but NONE determined to a constant
  (all depend on the 9447 free nonlinear-core vars). So linear methods + merge still cannot reach the bits.
- x_24026 / x_27116 (the bridging slacks): cone = 5523 wires spanning 243 of 255 control bits (dense) -
  no small activating set. x_38215 (numerator feeding a1813: x_14402*x_24026=321447*x_38215) has an
  11-wire, 0-control-bit cone: whether x_24026 is FREE (defined by a1660, bridging) or FORCED to 0
  (defined by a1813 with x_38215 const) is the orientation crux; the witness uses the a1660 orientation.
- CONCLUSION: the full solve = solve the densely-coupled nonlinear core (residue cascade + merged twist).
  It resists linear algebra, propagation, the merge, perfect-square reduction, and all search. This is
  the irreducible one-way trapdoor; inversion needs the setter's secret. Deliverable 39,019/39,031.

### Session 6 — slack activation traced (kernel_algebra.py, free_slack.py) — the coupling is exact
- x_38215=0 and x_29437=0 are CONSTANTS. So a1813 (x_14402*x_24026=321447*x_38215) => x_14402*x_24026=0,
  and a1815 => x_14402*x_27116=0. With a1657 (x_14402=1-x_12779): to get x_24026 != 0 (slack ON) need
  x_14402=0 i.e. **x_12779=1** (set by single 22-side bits: 19520,2795,1858,5443,26947,37748,27512,...).
- With x_12779=1, x_24026 is NOT a product but is LINEARLY pinned:
    a23394: x_24026 = x_12520 - x_29798 ;  a23402: x_24026 = x_29798 - x_1628
    a23391: x_29798 = x_23268 - x_36614 ;  a23395: x_12520 = x_24245 - x_23268
  so x_24026 = x_24245 - 2*x_23268 + x_36614  (a linear form in deeper residue-network wires).
- Setting bit 19520 alone: x_12779=1 but x_12520=x_29798=0 => x_24026=0 (slack still off), x_9770=0,
  viol=15. Activating the slack needs the deep wires x_12520/x_29798 != 0 (coordinated 233-side bits).
- COUPLING: x_12520 also appears in the degree-4 check a45004 (and x_29798 in a38195, both in a41470),
  so driving x_24026 to the bridge value m2*g2 - m*g necessarily perturbs those big checks. The witness
  must simultaneously (i) x_12779=1, (ii) x_24026 = x_18274 - x_35186 via the deep network, and (iii)
  keep a45004/a38195/a41470/a40782 satisfied. That simultaneous coupled system over the residue network
  IS the trapdoor. Complete to the finest level; the merged system is consistent (a solution exists) but
  the coupled nonlinear inversion resists linear algebra, propagation, merge, squares, lattice, search.

### Session 6 — the nonlinear core is 18,661 free vars (algebraically infeasible) — FINAL
Reducing ALL 26,185 nonlinear atoms through the linear network (build_pivots + reduce_var), the
distinct FREE variables they depend on = **18,661** (224 residue-loads, 57 control bits, 18,380
other free vars). So the nonlinear kernel is an 18,661-variable / 26,185-equation polynomial system
over GF(P) - far beyond Groebner (feasible only for ~20-50 vars) and every other algebraic method.
There is NO small separable core: the twist's free vars {24245,29798,36614,31434,34236,35846,26977,
35186} are each pinned by big degree-4 checks (a45004,a38195,a41470,a40782) that also mix in the
residue loads, so the whole thing couples. Combined with all prior results this is the definitive,
evidence-based conclusion: a correctly-built obfuscated-circuit trapdoor, not invertible by any
general or custom method available here. A witness exists (the merged system is consistent) but
recovering it is the designed one-way step. Deliverable: 39,019/39,031 (verified).

### Session 6 — extended exploration (3-hr campaign): more new methods, all negative
- GF(2)/2-adic: GF(2) linear rank on the atoms = 19381 (same as GF(P)); 198 control bits are pivots,
  0 determined to a constant. Mod 2 the residue products remain (boolean circuit = SAT); no help.
- Bridge linearity (bridge_linear.py): with the slack gate ON (x_12779=1 via a 22-bit), x_24026 and
  x_9770 are LINEAR in the 233 bits but their coefficients are all zero (x_24026 stays 0; needs
  weight>=4 coordinated activation). x_18274 stays nonlinear. So no subset-sum handle even slack-on.
- Factoring the quantization units: g=89*155682971*(243b), g2=2*79*1625329252399*(250b),
  h=2^2*3229*134807*(250b), h2=2*3^2*28843*14272028233*(240b), D1=3*13*3343*(250b), D2=2*19*71*(250b).
  Distinct small factors, distinct large cofactors, NO shared large prime -> no hidden modulus.
- Running long campaigns: sa_campaign2 (SA), ga_campaign (genetic algorithm over 255 bits, biased to
  slack gate). Both minimize mod-P violation count; checkpoint any state <4. Landscape min is 4
  (all-0); the witness is isolated, so these are low-odds but running per instruction.

## Session 7 — SLACK-ACTIVE BREAKTHROUGH: the twist IS satisfiable; obstruction reduced to a 9-atom frustrated core

Prior sessions concluded "forward-eval cannot represent the witness (slack identically 0)"
and treated the twist as a rigid coprime-quantization trapdoor. Both conclusions were
artifacts of ONE orientation choice. Re-examining the open atoms ORIENTATION-FREE changed
the picture (twist_neighborhood.py, slack_freedom.py):

The 4 open atoms are NOT symmetric:
  a44271: x_3183 - x_17728 = 0                          <- HARD (no slack): x_3183 = x_17728
  a30378: x_3183 - x_9982 - x_17728 = 0                 <- with a44272(x_9982=0) => same as a44271
  a1817 : 6033033*x_9770 - 6033033*x_18274 + x_26977=0  <- x_26977 = x_20510*x_31302 (a1816)!
  a40782: Q^2 = 0                                        <- perfect square (resid@best is exact square)
So the x_9770/x_18274 gap is bridged by the PRODUCT x_20510*x_31302 (a free slack), NOT forced
to 0.  x_9982 = x_9897*x_12518 (a1818) and a44272 forces it 0.  Only x_3183=x_17728 is rigid.

### The div-wire trap (why forward-eval forced the slack off)
x_24026 is a 'div' wire defined by a1813: x_14402*x_24026 = 321447*x_38215, i.e.
x_24026 = 321447*x_38215 / x_14402.  With x_12779=1 => x_14402=0 => division by zero =>
forward-eval's `elif den==0: val[v]=0` branch sets x_24026:=0.  THAT is the sole reason the
slack never activates.  x_38215 is const 0, so x_24026 is 0 either way in forward-eval.

### Slack-active evaluator (slack_active.py) — the twist now HOLDS
Two-pass: pass 1 computes side-values; then FREEZE x_24026 := x_18274-x_35186 and
x_27116 := x_17728-x_1642 as exogenous inputs and re-run forward-eval with x_12779=1.
Result: x_9770 = x_35186 + x_12779*x_24026 = x_18274  AND  x_3183 = x_1642 + x_12779*x_27116
= x_17728.  BOTH twist halves hold by construction.  Atoms 1817/30378/44271 are satisfied.
Confirmed for activating bits 1858/2795/5443/19520/26947/... (10 single 22-side bits set
x_12779=1).  This reaches the slack-active witness state plain forward-eval cannot represent.

### Cost: the ripple
Activating the slack perturbs ~46 downstream vars (via a23394 x_12520=x_24026+x_29798,
a23402, a23395 ...) and breaks 18-28 secondary CHECK atoms.  Setting the x_12779=1 bit itself
also breaks ~12-18 (it couples x_12779 to x_21941/x_36641/x_30323/... via a23149/a23150/...).

### Reduction to a 9-atom frustrated core (slack_repair.py, min-conflicts)
Seeded from bit 1858 (18 broken), integer min-conflicts (solve each broken atom for one of
its vars via divide/quadratic) drives 18 -> 9 and plateaus.  The persistent core:
  [2546, 7917, 23150, 23152, 23405, 39550, 40782, 42222, 44154]
  - 2546/7917/23150/23152: simple xA*xB = xC product-defs (trivially satisfiable alone)
  - 23405: linear x_15083 = 11855869*x_24852
  - 42222 (26 vars), 44154 (13 vars, involves x_12779): degree-2 checks
  - 40782, 39550: degree-4 PERFECT SQUARES (Q has 38 terms, deg 2 each; try_sqrt extracts).
    solve_for returns None for deg-4 => never repaired => the plateau. Replacing them by their
    roots Q=0 makes them solvable (slack_sa.py, SA over A' with roots).

### Status
The problem is no longer "unsatisfiable trapdoor" but "satisfy a 9-atom frustrated core in
slack-active space."  This is a genuine, large reduction. Local moves plateau because the
core's 129 vars are entangled with the 39k satisfied atoms (changing a core var breaks a
satisfied one). Active attack: SA with square-roots (slack_sa.py) + exact solve of the core
subsystem.  The witness genuinely differs from best_partial in x_12520/x_24245/x_23268/... and
requires a coordinated (not local) move — consistent with a hard-but-satisfiable circuit.

### Session 7 (cont.) — quantization is intrinsic; div-wire & slack escapes both bounded; SA grinds 18->6

Confirmed the coprime-quantization exactly (linearity_233.py, subset_sum_probe.py):
  x_9770 in g*Z, x_18274 in g2*Z, gcd(g,g2)=1     (single-bit deltas are all +/- base)
  x_3183 in h*Z, x_17728 in h2*Z, gcd(h,h2)=2
Single-bit deltas of every twist var equal +/-(its base), i.e. each var = base*(signed int).
So the twist x_9770=x_18274 needs g*m=g2*m2 with gcd(g,g2)=1 => m multiple of g2 (~10^88)
=> within the achievable small-|m| range, only m=0 (degenerate). Same for x_3183=x_17728
(gcd 2, needs m' multiple of h2/2 ~10^88 => 0). This is the trapdoor's core.

Div-wire ESCAPE (real but bounded): x_8821 = x_17810*x_27292 takes {-2,-1,0,1} (LINEAR in
18 bits, disjoint from the ~193 numerator bits). Since x_18274 = x_6773/x_8821 and x_17728 =
x_17233/x_8821 (div wires a4954/a13204), a nonunit x_8821 lets x_18274 escape g2*Z (24/300)
and x_17728 escape h2*Z (30/300).  BUT the escaped values are x_17233/(+/-2) = (h2/2)*Z, and
gcd(h, h2/2)=1 again => still degenerate for the twist.  x_17728 was a multiple of h (the
22-side quantum) in 0/500 samples.  The numerators x_17233/x_6773 are quantized to base*Z at
single-bit level and NONLINEAR for multi-bit (7/50 linear) => NOT a clean subset-sum, so the
density-0.8 LLL attack does not directly apply.

Re-orientation (reorient.py): forcing x_18274<-a11398 (linear gate x_6283+x_31434) and
x_17728<-a11388 instead of the div wires stays a VALID evaluator (forward_Z([]) still violates
only the twist-4) but makes x_18274 CONSTANT (1 value / 300).  The achievable set over VALID
(all-atom) assignments is intrinsic; re-orienting cannot break the quantization.  Strong
evidence the twist is a genuinely constructed one-way (knapsack/lattice) trapdoor.

The witness therefore uses the SLACK (x_3183 = x_1642 + x_12779*x_27116 with a huge x_27116),
and the difficulty is that a huge slack ripples into the verifier SQUARE checks (a40782=Q^2,
a39550=Q^2, ...).  Q40782 is degree-2; at the twist-satisfying point its twist terms
[39*(x_3183-x_17728), 6033033*(x_18274-x_9770), -42*x_9982, -x_26977] all vanish, leaving a
LINEAR equation in the ripple vars (28*x_10783 - x_30323 - 477654111*x_24252 + 4*x_32039 +
16*x_36641 + ...) = 0 that the ripple must satisfy exactly.

SA in slack-active space with squares replaced by roots Q=0 (slack_sa.py): reduces 18 -> 6
violated ([14440,23149,23151,23153,40782,44154] at last check; a39550 fixed).  Running a
4-way multi-start fleet (activators 1858/26947/27512/5443).  This is the first search
operating in the correct (slack-active, twist-holding) space; prior SA/GA were slack-OFF.

### Session 7 (cont.) — why the continuous slack knobs cannot fix the verifier square (crisp result)
x_24026 (9770-side slack) and x_31302 (df=None, a genuinely FREE variable) are the
only continuous knobs; x_27116 is pinned by the RIGID a44271 (x_3183=x_17728).
Write Q40782 = -6033033*x_9770 + 6033033*x_18274 - x_26977 + R, where R = 28*x_10783
+ (other ripple terms) is independent of x_9770 and x_26977. Then:
  a1817:      x_26977 = 6033033*(x_18274 - x_9770)
  Q40782=0:   x_26977 = 6033033*(x_18274 - x_9770) + R
Both hold  <=>  R = 0.  R is fixed entirely by the rigid 3183-side slack ripple
(x_10783 = x_16644*x_17301 ~ 10^272 at the slack-active state) and is UNAFFECTED by
x_24026 or x_31302 (verified: Q40782 is linear in x_24026 with slope 0 once a1817 is
maintained). So no continuous freedom can satisfy the verifier square -- only a
DISCRETE bit-setting whose ripple self-cancels (R=0). Finding such a setting is the
knapsack/lattice inversion the trapdoor protects; the 4-way slack-active SA fleet
searches exactly this (best 6/39031-atoms-in-A' so far). This is the tightest
statement of the obstruction: the witness = a 233/22-bit choice making the rigid
3183-slack ripple self-annihilate inside every verifier square.

### Session 7 (cont.) — CREATIVE structural attacks on the kernel internals
Non-standard probes (factor_squares/degree_test/constant_mining/x8821_zero):
1. Verifier square roots do NOT factor into linear forms: Q40782 quadratic-form rank
   30 (of 53), Q39550 rank 38 -> genuine irreducible quadrics, no linear disjunction.
2. Numerators x_6773(->x_18274)/x_17233(->x_17728) are >degree-2 in the 233 bits
   (deg-2 model fits only ~50-63%); low-order interactions PRESERVE quantization
   (all pair-interactions are multiples of the base g2/h2), so escape from g2*Z/h2*Z
   is only via HIGH-degree bit interactions + the div-wire (/x_8821).
3. R=0 reduction of the twist square a40782 (x_12779=1, twist held): using
   x_10783 = 14474367*x_12779*x_24252 (a23153) and x_30323 = x_12779*x_32039 (a23151),
   Q40782 residual R = -72372835*x_24252 + 3*x_32039 + 16*x_36641. With x_24252=a43825
   a perfect square (2025=45^2,729=27^2 => x_24252 = Q'^2 >= 0), R=0 demands
   72372835*Q'^2 = 3*x_32039 + 16*x_36641 -- a hard number-theoretic (square) condition
   on ripple vars. At slack-active state x_24252 ~ 10^250 but R=0 needs ~10^141.
4. x_16644, x_18950, x_31302, x_16644 are FREE vars (df=None) -> the broken product
   CHECK atoms a23152 (x_16644*x_17301=x_10783), a23150 (x_16320*x_18950=x_30323) are
   satisfiable by setting the free factor; they are NOT the real obstruction.
The real obstruction remains the verifier-square number-theoretic condition (item 3),
which is fixed by the bits (the knapsack). Testing the x_8821=0 escape regime next
(x8821_zero.py): when x_8821=0 the div wires collapse and x_18274/x_17728 are freed
onto their linear gates -- checking whether they can then hit the 22-side values.

### Session 7 (cont.) — DEFINITIVE: the collision gap is a CONSERVED INVARIANT
The 233-side twist activation (freeze x_18274:=x_9770, x_17728:=x_3183) gives a
SQUARE-FREE frustrated core; root-replacement (458 deg-4 squares -> deg-2 roots) lets
the repair reach a REPRODUCIBLE 4-atom core [8464, 19480, 41459, 44129] (multiple
independent seeds converge there). Structure: a44129 = irreducible rank-20 quadric
(contains all vars of the other 3); a8464 = product-def (x_10269=x_22895*x_24089,
x_24089 free); a41459/a19480 = linear (x_12390=x_14494, x_26526=0).

THE KEY OBSERVATION: at the 4-core, resid(a41459)=resid(a19480)=
63398753350954830538284979531311478224817569395477016427713014637060524103217265241016814
which is EXACTLY the twist gap G = x_17728-x_3183 (= resid of a44271/a30378 at best).
So the 233-side activation did NOT remove the collision gap -- it RELOCATED it from
(x_3183-x_17728) to (x_12390-x_14494). Across EVERY reformulation tried this session
-- twist atoms, slack-active verifier-square residual R, the x_8821 denominator web,
and now x_12390-x_14494 -- the SAME gap G=6.3e148 is conserved and merely moves to
different variables. The joint-solve via free x_24089 fails because a44129's
x_24089-dependence cancels (c1=c2=0): the free var cannot absorb the gap.

CONCLUSION: G is an INVARIANT of the system -- it can only be zeroed by a NATIVE
22/233 collision (x_3183=x_17728 with consistent bits), which the coprime
quantization (gcd(h,h2)=2, images tiny & coprime, 0 nonzero collisions in 3500
samples) forbids except at the degenerate 0. This is the tightest possible statement
of the trapdoor: the witness requires inverting the knapsack that produces the native
collision; no reformulation (slack, div-wire, re-orientation, root-reduction,
free-var joint-solve) can absorb the conserved gap. The tightest reduction achieved:
a single reproducible 4-atom core carrying G in one variable difference.

### Session 7 (cont.) — campaign round: x_24026-activation DEAD, LLL inapplicable
Exploration-queue round (activate_x12779_2.py, cond_linearity.py):
1. x_24026 ACTIVATION is structurally BLOCKED. The div-by-zero forcing x_24026=0
   occurs only at x_12779=1 (x_14402=0); at x_12779=2..5 forward-eval computes
   x_24026=321447*x_38215/x_14402 CORRECTLY -- BUT x_38215 = x_37917*x_30077 is
   IDENTICALLY 0 in forward-eval (image {0} over 400 high-weight samples; x_12779
   reaches {0,1,2,3,4,5}). So x_24026=0 regardless of x_12779. The 9770-side slack is
   un-activatable in forward-eval; the witness needs x_38215!=0 (an escape state).
2. LLL is INAPPLICABLE: x_17233 (num of x_17728) is nonlinear even CONDITIONAL on
   holding the 18 x_8821 bits fixed (9/40 linear; x_6773 5/40). Single-bit deltas are
   +/-base but multi-bit combinations are nonlinear (numerator loads interact via
   products). So the twist target is NOT a linear subset-sum -- it is base*f(bits)
   with f a nonlinear boolean fn of SMALL image (~12 values). No rich linear
   structure for lattice reduction; the density-0.65 subset-sum does not exist.
3. Homotopy/continuous relaxation: infeasible (38748 vars, ~10^250 values exceed
   float precision; exact "continuous" is meaningless for the discrete escape).
Every queue method reconfirms: the witness is a forward-eval-unreachable escape state
requiring the native collision / knapsack inversion the construction protects.

### Session 7 (cont.) — escape-source activation (x_12779=2) also blocked
NEW: at x_12779=2 (x_14402=-1) forward-eval computes x_24026=-321447*x_38215
CORRECTLY (no div-by-zero) -- so the slack CAN activate consistently IF x_38215!=0.
But x_38215=x_37917*x_30077 with x_37917==0 identically (image {0}); likewise x_29437
=x_7815*x_31807 with x_7815==0. x_37917/x_7815 are GATES (can't freeze without
breaking their atoms), and x_30077/x_31807 are free but multiply the 0-gates. So the
9770/3183-side slacks are un-activatable at EVERY x_12779 value in forward-eval. The
construction ensures each escape source is itself a 0-gate; the witness needs a
GLOBAL escape orientation (cascade of gates nonzero) that no local freeze reaches.
Note: x_18274!=0 IS compatible with x_12779=2 (39/61 samples); ~7053 free vars exist
but the load-active ones all multiply 0-gates.

### Session 7 (cont.) — CORE BREAKTHROUGH: escape cascade grounds at free var x_15
Tracing the "must be nonzero" requirement for x_24026 (slack) backward through the
confluent defining atoms (cascade_trace.py) reveals a LINEAR CHAIN of equalities:
  x_24026 <- x_38215 = x_30077*x_37917 ; x_37917=x_2524=x_9849=x_20564=x_3221=
  x_18850=x_4384=x_414 = x_15  (FREE var, df=None)
And ALL THREE slack cascades ground at the SAME free var x_15:
  9770 slack:  x_38215 = x_30077*x_15     (x_30077 free)
  3183 slack:  x_29437 = x_15*x_31807     (x_31807 free)
  9770 a1817:  x_26977 = x_20510*x_31302 = x_15*x_31302  (x_31302 free; x_20510=x_15)
So the slack IS forward-eval-activatable by SETTING x_15 (forward-eval kept x_37917
==0 only because x_15=0 at best). This overturns "forward-eval cannot represent the
witness" -- it CAN, via the free vars.

Consequence: at x_12779=2 (x_14402=-1) forward-eval computes x_24026=-321447*x_38215
=-321447*x_30077*x_15 correctly; x_3368=2*x_24026; x_9770=x_35186-642894*x_30077*x_15.
So the twist x_9770=x_18274 needs only  x_18274 == x_35186 (mod 642894), and
x_3183=x_17728 needs x_17728==x_1642 (mod 2). BOTH ARE REACHABLE (mod_reach.py:
61/387 x_12779=2 states satisfy both; 0 is in the image mod 642894 and mod each prime
factor 2,3,7,15307). The twist is a MODULAR match (~20-bit), NOT a 296-bit collision!
On a mod-satisfying state: freeze x_24026=(x_18274-x_35186)/2, x_27116=(x_17728-x_1642)
/2; set x_15=1, x_30077=(x_35186-x_18274)/642894, x_31807=(x_1642-x_17728)/2 -> twist
holds EXACTLY, a1813/a1815 SATISFIED, x_26977=0. Remaining breakage = pure ripple /
verifier-square residuals R_i (free_activate4.py searches mod-states for min |R|).
This is the real, tractable core: find a mod-satisfying state whose ripple self-
cancels, OR tune the ripple via additional free vars.

### Session 7 (cont.) — x_15 is a HUB; escape = ~15 chain vars all equal to one value V
Verified activation works (twist holds when run() is applied on top of solve(bestval,S)
so bits are set -- earlier free_activate3/4 had this bug). BUT the escape cascade vars
x_15,x_414,x_4384,x_18850,x_3221,x_20564,x_9849,x_2524,x_37917,x_22447,x_2322,x_7815,
x_21260,x_31902,x_20510 are ALL forced EQUAL (gate chain) and are HUBS: 49-88 atoms
each, 1033 atom-slots total. Setting x_15!=0 (required for slack, since x_38215=
x_30077*x_37917=x_30077*x_15) breaks ~350 atoms. So free-var activation is MORE
disruptive than the old freeze-x_24026 hack (18 broken -> 4-core). The escape is a
GLOBAL state: ~15 hub vars all = witness value V, consistent across 1033 atom-slots.
KEY next attack: MERGE the chain (x_414=x_4384=...=x_20510=x_15=V, one unknown) into
the system, substitute, and solve the reduced system for V + the twist -- this could
collapse the 350 ripple since the atoms then constrain V consistently rather than
appearing broken. Also: the mod-twist reduction stands (x_18274==x_35186 mod 642894
reachable per-prime via CRT; free_activate4/mod_reach). Old freeze approach (4-core,
invariant gap G) remains the tightest concrete state; reconcile G with the mod
constraint next.

### Session 7 (cont.) — the escape is a 220-var WIRE V, decoupled from best at 1st order
union_identities.py: union-find over all 3707 two-term identity/negation atoms
(c1*x_a+c2*x_b=0, |c1|=|c2|) finds ONE giant class of 220 variables all = x_15 (=V,
all +sign) -- a single obfuscated "wire" renamed 220 times. It is the ONLY large
class (all others size<=2). (x_3183 & x_17728 are merged too -- that's the twist atom
a44271.) The wire = escape master V; the slack needs V!=0 (x_38215=x_30077*x_15).
wire_analysis.py: the wire appears in 5233 atoms, up to degree 4 (V^4 in the deg-4
verifier squares). Substituting wire->V with all NON-wire vars at best_partial and
using the deg-2 square ROOTS:
  - 0 atoms FORCE V=0 (none is linear-only c1*V=0) -> V!=0 not trivially blocked
  - 0 atoms DEMAND a nonzero V either: every wire V-term is V*x_other or V^2*x_other
    where x_other=0 at best, so ALL V-coefficients vanish at best_partial.
=> The wire V is DECOUPLED from best_partial at first order. best_partial sits on the
V=0 branch (39019 atoms); the witness is on a SEPARATE V!=0 nonlinear branch where V
and its coupling-partner vars (currently 0) co-activate. This is exactly why every
local repair / perturbation of best_partial was UNSAT: the witness is not a
perturbation of best -- it is a different branch. The trapdoor = find the co-activated
(V, partners) branch; no first-order/linear handle from best reaches it.

### Session 7 (cont.) — the escape co-activates V + ~3707 FREE partner vars
partner analysis (squares->roots): the 220-var wire V couples (shares a monomial)
with 3785 distinct NON-wire vars; 3707 of them are FREE (df=None) and 3769/3785 are 0
at best_partial. So the witness is a SEPARATE branch that sets V!=0 AND ~3707 free
partner vars nonzero, co-activated so all 5233 wire atoms hold. best_partial (V=0, all
partners 0) is the degenerate branch. This is a large nonlinear CSP over ~3707 free
inputs -- the trapdoor. Total free vars = 7053; effective freedom ~256 bits (the
partners are cross-constrained by the ~11k check atoms). NEXT ATTACK: propagation-
solve seeded from V!=0 -- set the wire, then for each check atom with one unknown free
partner, solve+propagate; measure how many partners get pinned and whether it
converges (extends merge_solve which stalled at 5898 from pins alone). The mod-twist
reduction (x_18274==x_35186 mod 642894, reachable per-prime) is the other live handle.

### Session 7 (cont.) — propagation-solve stalls: escape is high-dimensional, not (bits,V)
propagate_wire.py: merged the 220-wire->x_15 and twist, seeded ALL 255 control bits +
x_15=V, propagated over GF(P) (linear + product/div, 1-unknown -> assign). Determined
only 14703/38748 with 0 contradictions and stalled. (Generic propagation is weaker
than the confluent orientation, which determines the ~31695 non-free vars; but the
point stands:) the ~7053 free vars are genuine inputs, NOT determined by (bits,V) via
propagation. The check atoms (~11k) constrain them to an effective ~256-dim solution
space, but NONLINEARLY -- so the witness is a high-dimensional co-activation of V +
~thousands of free partner vars, reachable only by solving the coupled nonlinear
system, not by propagation or perturbation from best_partial (V=0 branch).
FULL SESSION-7 characterization of the trapdoor (definitive):
 - obstruction = twist collision x_9770=x_18274 & x_3183=x_17728
 - the slack that bridges it grounds at a single free var x_15, aliased 220x (a "wire")
 - x_15=V is decoupled from best at 1st order (best is the V=0 branch); witness is a
   separate V!=0 branch co-activating ~3707 free partner vars
 - twist reduces to a modular match (x_18274==x_35186 mod 642894) reachable per-prime,
   but the verifier-square/ripple co-activation is the residual barrier
 - no first-order/linear/propagation handle reaches the branch => genuine one-way
   trapdoor; witness exists by construction but needs the coupled nonlinear solve.

================================================================================
SESSION 8 — ✅ SOLVED (all 39,031/39,031 equations, exact in ℤ)
================================================================================
The trapdoor fell by DROPPING the forward-eval / mod-P frame and working directly in
the RAW equation space. Full writeup in SOLVED.md; regenerator in build_solution.py.

1. At best_partial (39019/39031) exactly TWO shared atoms are nonzero:
     H = (x_17728 − x_3183) + x_9982     (x_9982=0 ⇒ H = G, the invariant gap)
     F = 6033033·(x_18274 − x_9770) − x_26977   (x_26977=0)
   Every one of the 12 failing equations is a linear combo of atoms that includes H or F.

2. The "rigid a44271: x_3183=x_17728" was a REFORMULATION ARTIFACT. In the raw text,
   (x_17728)-(x_3183) never appears alone — all 16 occurrences are paired with +x_9982.
   So the gap is a free PRODUCT SLACK, not a rigid equality:
     H=0 ⟺ x_9982  = −G   with atom 1818: x_9982  = x_12518·x_9897
     F=0 ⟺ x_26977 = F0   with atom 1816: x_26977 = x_20510·x_31302

3. The hub factors x_12518 (271 eqs) and x_20510 (237 eqs) both lie in the SAME 220-var
   identity class as x_15 (the "wire"). Setting the whole wire to sign·V (non-wire vars
   held at best) leaves V¹ and V²⁺ coefficients EXACTLY zero in all 5233 touched atoms —
   integer-exact, not just mod P. The wire is a genuinely free parameter.

4. Direct construction (V=1): wire member → its sign; the two rare partners (each in only
   2 atoms) → x_9897=−G, x_31302=F0; slack outputs x_9982=−G, x_26977=F0; rest = best.
   Satisfies 1818,1816,H,F and the verifier square a40782=Q² (Q→0). Checker: 39031/39031.
   x_12779 and x_24026 stay 0 — the div-wire / dirty-bits / x_12779≥2 saga was never needed.

WHY PRIOR SESSIONS MISSED IT: they searched the confluent forward-eval orientation (which
quantizes both twist sides to coprime units and zeros the slack products — literally cannot
represent the witness) and chased control-bit settings. The witness is not reachable by
propagation from the V=0 branch, but it IS a one-line algebraic construction once you read
the slacks off the raw equations and note the wire is quiet.

---

## Session 9 — independent re-derivation, a Jacobian bug, and a local-optimality proof

**Start:** ran the start-of-session ritual. `checker.py best/new_instance_partial_39022.json` →
`39022/39033 (11 failing)`, fails `[2554,6816,8124,8680,9421,12231,12270,12350,14584,22044,29125]`.
Deliverable confirmed before touching anything.

**Rebuilt the model from scratch** (`s9/`, independent of all prior tooling, whose `atoms/gates.jsonl`
intermediate no longer exists). `atomize.py` decomposes each equation into `(coeff, atom)` pairs and
**validates by exact integer re-evaluation on all 39,033 equations — 0 mismatches**. 42,267 atoms;
degree histogram {1: 19780, 2: 21788, 4: 699}. Exactly **3 atoms nonzero** at the partial:
22229, 22231, 37887 (the last a square whose root contains 22231) → the obstruction is 2 scalars.

**Canonical gate orientation from syntax** (`gates.py`): output = first bare-variable top-level
additive term. 31,475 gates / 7,273 free inputs / 10,792 checks. Validated: over all 7,273 unit
perturbations the forward ripple leaves **0** gates unsatisfied (`jcheck.py`).

**Identified the obstruction exactly.** `x_28599 = x_17499 = x_26064 = p = 2^256−2^32−977`
(secp256k1) exactly ⇒ `x_642 = p·x_17325`, `x_28730 = p·x_9413`, so the two atoms are the congruences
`x_7068 ≡ x_2099` and `x_4432 ≡ x_19964` mod p. Backward cones show `x_2099`/`x_19964` are 2-bit MUX
outputs (controls `x_4287=0`, `x_2081=1`) resolving to free inputs `x_6418`/`x_12553`, which load
gates 3576/3578 pin to constants K1/K2. So: **`x_7068 ≡ K1`, `x_4432 ≡ K2` (mod p)**.
`gcd(15804267, 7376877p) = 3` and `3 | D1`, so the pin side is Bézout-solvable — only p-divisibility bites.

**Mapped the conserved obstruction.** Canonical repair (`x_7068:=K1`, `x_4432:=K2`) ripples through
only **115 variables** and relocates the defect to mirror congruences 29539/7930
(`x_14853 ≡ x_1308` with `x_1308 = x_7068` identically; `x_24548 ≡ x_25442`). Repairing those with
the free leaves breaks the core M1/M2/M3 (19297/19299/30984). Derived the core exactly:
`L1 = 8646263S+1073965T`, `L2 = 10159099S+6926539T`, `L3 = 8272701S+5921311T`,
`S = A·u²−w²`, `T = B·u−w·c`. Core ⟺ `S ≡ T ≡ 0 mod p`.

**New: the core has a second branch.** Eliminating `w` from `T ≡ 0` and substituting into `S ≡ 0`
gives `u²·(A·c² − B²) ≡ 0`. So it is *not* only `u ≡ w ≡ 0` — `A·c² ≡ B² (mod p)` frees `u`
entirely. Blocked in practice by the mod-p pins on `x_22162` (atom 1618), `x_30213` (688),
`x_16742` (26731). Not previously recorded anywhere in this notebook; left open deliberately.

**Bug found in the Jacobian methodology.** First mod-p solve (all 7,273 free inputs, sparse
Gaussian elimination) reported INCONSISTENT at atom 42245 — which is a *perfect square*. At the
partial its root `E = 0`, so a finite-difference row for `E²` is `c²δ²`, quadratic, not linear.
783 check atoms are squares. `roots.py` extracts all their degree-2 roots (699 deg-2 + 84 deg-1,
none left over); `jac2.py` rebuilds on roots. The inconsistency then moves to atom **19297** — the
genuine core. Every earlier Jacobian/Newton/null-space conclusion in this notebook should be
re-checked against this: they may have been reading the artefact.

**Newton at a non-degenerate point.** At `u = 0` the core's gradient vanishes identically (S,T are
quadratic in u,w), so linearising at S0 is uninformative. Built S1 (both mirrors repaired, `u ≠ 0`)
and re-solved — still inconsistent at 19297, certificate now spanning 270 rows. Barrier is real and
correctly located.

**Bit scan.** 1,156 boolean free inputs (only 2 currently set). Flipped each (`bitscan.py`, ~1 s
total): exactly **two** — `x_2081`, `x_24601` — deactivate the core by forcing
`x_15298 = x_7715·x_34554 = 0`. Both are strictly worse locally (17 and 19 nonzero residuals vs 4).
These are the two known "quadrant" activators, now confirmed independently and exhaustively.

**Local-optimality proof (the session's main deliverable).** The handles absorb exactly the
p-multiples, so reachable defects are `A ∈ D1+pℤ`, `B ∈ D2+pℤ` — `x_7068` may also move by any
multiple of p (mirror handle `x_30163` absorbs it), non-multiples break a 12–16-equation atom.
Every failing equation is `m·(c₁A + c₂B)` with all other atoms zero, so it vanishes only if
`c₁·D1 + c₂·D2 ≡ 0 mod p`. Checked all nine linear ones — none. Eqs 8680 and 29125 need `B = 0`
exactly, impossible. Required ratios are `−c₁/c₂` with `|c| ≤ 40`; the actual
`D2/D1 mod p = 29086885819837044927165644576879200888791656444895144487473726012333934270214`
is a full 256-bit residue. And all six alternative placements cost 23–29 equations vs 11
(22229+22231+37887 = 11; 22229+3578 = 23; 29539+22231 = 24; 3576+22231 = 24; 22229+7930 = 25;
29539+7930 = 29). **Therefore 39,022/39,033 is a local optimum, and any improvement requires
cracking the core.** Strategies A/A1/A2/B1/B2 measured 39,002–39,013, consistent with this.

**Deliverable unchanged and re-verified: 39,022/39,033.** Writeup in `S9_STRUCTURE.md`.
