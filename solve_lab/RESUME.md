# RESUME — read me first

## STATUS (session 9): best verified **39,024 / 39,033**
Deliverable: `best/new_instance_partial_39024.json`
Verify: `python3 checker.py best/new_instance_partial_39024.json` → `satisfied 39024/39033 (9 failing)`
Independent check: `python3 s9/verify_ast.py best/new_instance_partial_39024.json` (AST walk, no eval/regex).
Failing lines: `[9123, 9421, 12231, 12270, 12350, 14584, 18673, 22044, 29125]`.

> **RETRACTED:** an earlier version of this file claimed 39,022 was a proved local optimum.
> That was wrong — see `S9_STRUCTURE.md` section 6. The proof assumed every atom outside the
> defect set must vanish. **It must not.** An equation is zero iff its *linear combination of
> atoms* is zero, so any atom whose whole equation footprint already lies inside the failing set
> is a FREE knob. Atom 22230 (`x_28730 − x_17499·x_9413`) is exactly that — `x_9413` and
> `x_28730` appear in no other atom — which frees `x_28730` from the lattice `p·ℤ`. Extending
> along the `35754…35762` ladder gives 5 knobs over 13 equations, 4 of them simultaneously
> zeroable: 11 + 2 − 4 = 9 failing. **Work in EQUATION space, not atom space.**
**Read `S9_STRUCTURE.md` first** — it supersedes the older analyses below on every point of conflict.

### The 60-second version
Exactly **3 atoms** are nonzero at the partial (22229, 22231, and the square 37887 whose root
contains 22231). Since `x_28599 = x_17499 = p = 2^256−2^32−977` (secp256k1 prime) exactly, the
entire residual system is two congruences:

    x_7068 ≡ K1 (mod p)      and      x_4432 ≡ K2 (mod p)

with `x_17325`, `x_9413` as free quotient handles and `K1, K2` the constants that atoms 3576/3578
pin `x_6418`, `x_12553` to. Everything else in the 39,033 equations is already exact in ℤ.

### Three results from session 9 that change the picture
1. **METHODOLOGICAL BUG FOUND — re-check anything built on a Jacobian.** 783 check atoms are
   perfect squares `E²`; at the partial `E = 0`, so a finite-difference row for `E²` is *quadratic*
   (`c²δ²`), not linear. Feeding squares into a linear/Newton/null-space solve solves the wrong
   system. `s9/roots.py` extracts all 783 roots. After the fix the reported mod-p inconsistency
   moved off a *spurious* square (42245) and onto the *genuine* core atom (19297). Earlier
   sessions' "inconsistent Jacobian" verdicts were reading the artefact.
2. **The core has a SECOND branch, never previously recorded.** Core ⟺ `S ≡ T ≡ 0 mod p`, and
   eliminating `w` gives `u²·(A·c² − B²) ≡ 0`. So it is *not* only `u ≡ w ≡ 0`: the alternative
   `A·c² ≡ B² (mod p)` frees `u` completely. Blocked only by the mod-p pins on `x_22162`,
   `x_30213`, `x_16742`. **This is the one door opened and not closed — start here.**
3. **39,022 is a local optimum, proved (not just observed).** Handles absorb exactly the
   p-multiples, so the reachable defects are `A ∈ D1+pℤ`, `B ∈ D2+pℤ`. Each failing equation is
   `m·(c₁A + c₂B)`, so it can vanish only if `c₁D1 + c₂D2 ≡ 0 mod p` — checked for all nine
   (none), and eqs 8680/29125 need `B = 0` exactly (impossible). All six alternative defect
   placements cost 23–29 failing equations vs 11. **Improving the score requires cracking the core.**

### BIGGEST NEW LEAD (session 9): the core is NOT a wall — both cores can be zeroed
`x_8599 = 1` (88 of the 1,156 boolean free inputs do it, keeping `x_21839 = 1`) reroutes
`x_12186` from the computed `x_30454` to the **free input `x_5096`**. Set `x_5096 := K1`,
`x_14853 := x_12186`, let the gates close C1/C2 → `u = w = S = T = L1 = 0` exactly: **the core
(19297/19299/30984) is satisfied simultaneously with C1 and C2 for the first time.**
It lights a second core of identical shape (26733/28438/32342, gated by `x_38170`), and *that one
is zeroable too*: ~90 boolean free inputs shift its controls by exact complements mod p
(`u' + δ = p` and `w' + δ' = p` exactly). Two-bit constructions (e.g. `x_2527` + `x_1502`,
`s9/construct3.py`) leave **both cores clean**, with the entire residual being activated load pins.
Closing those pins moves variables in the first core's cone and re-lights it (11→16→13, stalls).
**So the real invariant is the pin/mirror cascade, not the core.** Attack that next:
enumerate each bit's pin set and search for a bit subset whose pins are jointly satisfiable.

The core's other branch (`A·c² ≡ B² mod p`) is **CLOSED**: it needs `A` to be a quadratic residue;
`A` is fixed by the pins and is a non-residue for both reachable `x_12186` residues (a ~50% sample
of shifts *would* be residues — the setter's pin of `x_14853` to `K1` is what closes it).

### THE LOCK, IN CLOSED FORM (read S9_STRUCTURE.md section 12)
`x_12186 ≡ x_22649 (mod p)` exactly 1:1 and `x_22649` is FREE — it is the circuit's output wire
(found with `s9/modtrace.py`). Moving `{x_22649, x_22152, x_14853, x_7068}` together by
`δ = K1 − x_12186` **kills chain 1 outright** (atom 22229 = 0) while holding the core's `u = 0`.
What stops it is load pin **31670**, gated by bit `x_24601`, which pins the output wire to exactly
the original `x_12186` — a constant differing from `K1` by precisely `D1`.

Flipping the MUX control **`x_4287 = 1`** makes `x_2099`/`x_19964` read the FREE inputs
`x_9118`/`x_8731` instead of the pinned `x_6418`/`x_12553`, so **atoms 22229 AND 22231 are both
zero simultaneously** — the first time. Its price is three new loads that collapse (3-from-2) to
`x_4306 ≡ x_27177 ≡ 0 (mod p)`, which IS solvable — but `∂x_27177/∂x_8731 = 0`, so `x_9118` mod p
is uniquely forced to `33371159155735472537534252650716501592825364489306217536352743247010353604716`,
while the mirror+core chain needs it to equal the pinned `x_12186`
(`82007976…230357` if `x_24601=1`, or `0` if `x_24601=0`). **Neither matches — that mismatch IS the
trapdoor.** Any further attack should target that single scalar identity, not search assignments.

### Do NOT redo
- Greedy ripple-repair from the canonical orientation (converges to the core in 3 rounds, 39,013).
- Single-bit flips: all 1,156 boolean free inputs scanned; only `x_2081`, `x_24601` deactivate the
  core (`x_15298 = 0`), both strictly worse locally (17/19 residuals vs 4).
- Repairing via `x_7068`/`x_4432`/`x_6418`/`x_12553` (strategies A/A1/A2/B1/B2 in `s9/drive.py`):
  39,002–39,013, all worse than 39,022.
- Mod-p linear solve over all 7,273 free inputs at fixed bits — inconsistent at S0 *and* S1 even
  with the corrected root model (270-row certificate).

### Next experiments, in priority order
1. **The pin/mirror cascade** (see the lead above) — with both cores now zeroable this is the only
   remaining invariant. Per activated bit, enumerate its load pins `bit·(x_B − HUGE) = s·x_C`, then
   search for a bit subset whose pins are jointly satisfiable (exact-cover over the pinned free
   inputs). `s9/pinclose.py` is the naive greedy version and it stalls; the set-cover view is new.
2. `s9/construct3.py` reached both-cores-clean; extend it rather than restarting.
3. Full two-bit flip scan (~667k pairs). Session 9 only ran 88x14 = 1,232 targeted pairs.
4. Quadrant re-solve at `x_2081 = 0` / `x_24601 = 0` with the corrected root-based residual model.

### Toolchain (`s9/`, self-contained, caches regenerable)
`cd solve_lab/s9 && python3 atomize.py && python3 poly.py && python3 gates.py && python3 fwd.py`
then e.g. `python3 state0.py`, `python3 drive.py A`, `python3 bitscan.py`, `python3 newton.py`.
`atomize.py` validates the whole decomposition against the raw file (0 mismatches over 39,033 eqs).

### Git
Branch `claude/read-prompt-324ju4`.

---
# RESUME — read me first

## ⚠️ CURRENT (re-randomized) INSTANCE — 39,013/39,033; core REDUCED but not cracked
The EQUATIONS.txt in the repo is a NEW re-randomized instance (39,033 eqs). Full analysis in
**`NEW_INSTANCE_STATUS.md`** and **`CORE_REDUCTION.md`** (read both). Best verified partial:
**39,013 / 39,033** (`best/new_instance_partial_39013.json`, quadrant (1,1)).

### BREAKTHROUGH (latest session): the 20-equation core is fully reduced.
All 20 remaining verifier squares = integer combos of three monsters M1,M2,M3, wiring-defined by
two base gates S=x_35389, T=x_6671. Core ⟺ M1=M2=M3=0 ⟺ **S≡0 and T≡0 mod p** (quadrant 1,1),
then set private quotient handles x_30317,x_2936,x_5146. See CORE_REDUCTION.md for the full chain
down to control differences x_29322=x_14853-x_12186, x_3558=x_24908-x_16742.

### WIRE ESCAPE (the actionable path — read METHODOLOGY.md + CORE_REDUCTION.md's wire section)
Agent B PROVED the wire=p mod-p system is a rigid isolated point: rank(J_sat)=3035/3036 active
cols, null space dim 1, 19/20 core conditions directly contradict the wiring. The ~5547 "dead"
free inputs feed ONLY products against the p-wire (wire=p≡0 mod p → wire·handle≡0). This is WHY
the witness is unreachable on the wire=p branch.
THE ESCAPE: the 220-var identity wire (root 38100, forced to p by x_26064's single-var atom) is,
per METHODOLOGY.md, only "meant to vanish" — the witness lives on the wire≠p branch. Set the whole
wire = sign·1: then wire·handle = handle ≢ 0 mod p, ACTIVATING all ~5547 quotient handles → a huge
new null space. Core collapses: M1=L1+x_30317→x_30317=−L1 (trivial), M3→x_2936=537773·L3 (trivial),
M2→x_5146=L2/6672769 (needs 6672769|L2 — a 2^23 modulus; L2 mod 6672769 is message-controllable:
4239005 at 39013, 2032135 at 39018). Only ~13 "active unpackings" break (wire members as standalone
terms + (x_26064−p) checks + wire·x_31342 products): [8429,11166,11915,12594,23869,25313,26785,
31400,32300,36106,36767,37257]. NEXT: build the wire=1 global solve over the activated handles
(Dixon lift), heal the 13, set L2≡0 mod 6672769, verify. Agents B (wire=1 consistency via
tangent-linear) and E (wire=1 construction) are on it. Best partial: 39,018 (best_agentD_39018.json).

### OBSTRUCTION on the wire=p branch (superseded by the escape above): residues are pinned.
The sparse wiring solution is unique (only 30 slack inputs nonzero, rank 30). Sparse-witness
null-space solves up to the FULL closure (6,114 inputs, 7,119 constraints) are INCONSISTENT —
S,T residues are linearly pinned by the wiring; no local move reaches the core. Both quadrants
(1,1) and (0,0) have the SAME 20 hard squares via different branches; multi-role control
variables couple the core to the whole system. Remaining paths: a global nonlinear solver
(basin-hopping / large mod-p Gaussian) or the setter's witness. NEXT: try a global linear solve
from the all-zero (0,0) point (cleaner linearization: all products vanish → residue conditions
become linear), or re-attempt with a genuinely different activator/quadrant that unpins S,T.

Definitive findings (exhaustive):
- Gate DAG fully ACYCLIC; forward-eval from free inputs satisfies all wiring automatically.
- The ONE large identity wire (220 vars) is PINNED to **p = 2^256-2^32-977 = the field
  prime** (x_26064=p, appears 13x); it is the twist multiplier. No free wire exists (unlike the
  solved instance, whose wire was FREE — that was the whole solve).
- The remaining obstruction is a **256-bit boolean codeword message** (disjoint 178+78-bit cones
  of control bits x_7715,x_34554). x_9274=OR(controls)=1 is FORCED, so activation is mandatory.
- Each set bit triggers a huge additive load; the GF(p) load matrix over the 256 bits is FULL
  RANK (no bit self-cancellation). Data can absorb mod p but the ℤ-lift imposes p-divisibility
  carries -> every wrong-message data solve is ℤ-inconsistent (SNF pivot=p); iterative repair
  diverges (27->300+). Inhomogeneous GF(p) message solve inconsistent (rank 3). No small vinegar
  linearization (21,922 vars bilinear). => structured GF(p) codeword/MQ; needs the setter secret.
- Key tools: build_twist.py (activate+route MUX), newton2.py (simultaneous absorber solve),
  p_message.py (GF(p) message solve), scc.py, localize.py, scan_bits.py, p in huge_consts.json.
- Untried heavy attack: lattice/LLL on the mod-p codeword (standard, likely resisted by design).

---
# RESUME — read me first

## ✅ SOLVED — all 39,031 / 39,031 equations satisfied exactly in ℤ
Solution: `best/SOLUTION.json` (also repo root `SOLUTION.json`), 38,748 vars, 2,954 nonzero.
Verify: `python3 checker.py best/SOLUTION.json` → `RESULT: OK`.
Reproduce: `python3 build_solution.py`. Full writeup: **`SOLVED.md`**.

The twist reduced to activating two product-slacks (H: x_9982=x_12518·x_9897=−G;
F: x_26977=x_20510·x_31302=F0). The hubs x_12518/x_20510 sit in a 220-var identity wire
that is *exactly quiet*, so it is a free parameter; set wire=sign·1 and the two rare
partners x_9897=−G, x_31302=F0. Direct algebraic construction, no search. The earlier
forward-eval framing (control bits, div-wire, x_12779≥2, dirty bits, high-dim co-activation)
was an artifact of the orientation — none of it was needed.

---
## (historical) prior best partial
Best verified partial was **39,019 / 39,031** (`best/best_partial_39019.json`).

## SESSION 7 — SLACK-ACTIVE SOLVER BUILT; obstruction reduced to R=0 (read NOTEBOOK Session 7 tail)
The slack-active evaluator EXISTS now: `slack_active.py` (freeze x_24026:=x_18274-x_35186,
x_27116:=x_17728-x_1642 with x_12779=1 via a single 22-side bit e.g. 1858). It makes BOTH twist
halves hold by construction — the state plain forward-eval cannot represent. Activating the
slack ripples into ~18 verifier CHECK atoms; SA-with-square-roots (`slack_sa.py`, replaces the
deg-4 squares a40782/a39550 by their deg-2 roots Q=0 via `check_square.try_sqrt`) drives the
frustrated core 18 -> 6. Run the 4-way fleet: `python3 slack_sa.py <activator> <seed> <out.json>`
with activators in {1858,26947,27512,30104,5443,...}.
CRISP OBSTRUCTION: for verifier square a40782, satisfying it AND a1817 reduces to R=0 where
R = 28*x_10783 + (ripple terms), x_10783=x_16644*x_17301, all fixed by the RIGID 3183-slack
(a44271: x_3183=x_17728, so x_27116=x_17728-x_1642 is pinned). The continuous knobs x_24026 and
the FREE var x_31302 (df=None) CANNOT change R (Q40782 slope in x_24026 is 0 once a1817 held).
So the witness = a DISCRETE 233/22-bit choice whose rigid-slack ripple self-annihilates in every
verifier square (a knapsack). The div-wire escape (x_8821=x_17810*x_27292 in {-2,-1,0,1}) lets
x_18274/x_17728 leave their g2/h2 lattice but only onto (base/2)*Z, still coprime => degenerate.
NEXT: keep the slack-active SA fleet running; or attack R_i=0 across the 530 squares as a system
(linearize ripple monomials); or find the setter's 233-bit knapsack solution (LLL blocked by
numerator nonlinearity — 7/50 linear). NOTE: single 22-side bits give x_12779=1 (not 2).

## TRAPDOOR MECHANISM — fully reverse-engineered (Session 6, read NOTEBOOK Session 6 tail)
The obstruction (atoms 1817,30378,40782,44271) is the twist x_9770=x_18274 & x_3183=x_17728.
KEY: the confluent forward-eval QUANTIZES both sides to COPRIME units and ZEROS the slack
products, so it can NEVER represent the (feasible) witness — this is why every forward-eval
search (SA/mitm/greedy/pairs/enum/local) plateaued. Specifically:
- Under forward-eval: x_9770=m*g, x_3183=m'*h, x_18274=m2*g2, x_17728=m2'*h2 (g=119182..,
  g2=91416..; gcd(g,g2)=1, gcd(h,h2)=2). Rigid twist => degenerate 0 only. (codewords.py, quant_structure.py)
- BUT the wire DEFS carry product slacks: x_9770 = x_35186(=m*g) + x_3368, x_3368=x_12779*x_24026;
  x_3183 = x_1642(=m'*h) + x_10466, x_10466=x_12779*x_27116. Both gated by x_12779=x_23380*x_36336.
- forward-eval sets x_12779=0 (slacks off) -> quantization. The WITNESS activates x_12779 (22-side
  bit pairs give x_12779=2) AND x_24026/x_27116 (deeper, via x_38215) so
      x_9770 = m*g + x_12779*x_24026 = x_18274 = m2*g2   (bridges the coprime gap).
- So the TRUE solve = search WITH the slacks active. With slacks on, x_9770 is NOT limited to 27
  values and CAN equal x_18274; the decoupling (x_9770<-22 only) is a slacks-OFF artifact.

NEXT-STEP for a solver: build an evaluator/search that DRIVES x_12779, x_24026, x_27116 nonzero
(find their activating bit cascades: x_12779<-{1858,2795,5443,10652,19520,26947,27512,30104,...},
x_24026<-x_38215<-...), then solve m*g + x_12779*x_24026 = x_18274(B) (coupled product match).
Do NOT rely on the all-0 forward-eval regime — it structurally excludes the witness.

## Earlier (still true) reduction
- `A` = the 22 control bits `BITS22`; `B` = the other 233 bits.
- `x_18274 = x_6773/x_8821`, `x_17728 = x_17233/x_8821` (SHARED denominator x_8821).
- `x_8821` is **exactly linear** in the 233 bits; numerators are high-degree.
- best_partial_39019 sets ALL 255 control bits = 0.
- twist eqs: 1817 = 6033033*(x_9770-x_18274)+x_26977; 44271 = x_3183-x_17728;
  30378 = x_3183-x_9982-x_17728. (x_26977, x_9982 identically 0.)

## How to evaluate (the correct model)
`confluent_eval5.build5()` -> (A_atoms, kind, info, seq, bestval, ncyc). Build `seq`:
```python
order = json.load(open('eval_order.json'))['order']
defset = set(v for v in kind if kind[v] != 'const')
seq = [v for v in order if v in defset and v not in (9770,3183)]
seq += [v for v in (9770,3183) if v in defset]
seq += [v for v in defset if v not in set(order) and v not in (9770,3183)]
```
`make_forward(kind,info,seq,bestval)` -> Z solver `solve(list(bestval), setbits)`;
`make_forward(...,mod=P)` -> mod-P solver. forward_Z([]) violates exactly {1817,30378,40782,44271}.
The forward-eval satisfies every ORIENTED gate/load/div atom by construction for ANY bit set;
only the twist "check" atoms float — so it is a valid oracle for x_9770/x_3183/x_18274/x_17728.
NOTE: integer forward-eval is *lossy* (leaves a stale value when a division isn't exact) — use
the mod-P solver for any linearity/degree probing.

## Highest-EV next experiments
1. `runs/tab22_full.log` — full 2^22 (x_9770,x_3183) mod two 31-bit primes; saves
   tab22_9770_{p}.npy / tab22_3183_{p}.npy. When done: confirm B=0 fails; hash S and inspect
   structure (common factors, moduli). S then lets you INVERT the 22-side in O(1) (lookup).
2. Residue-pool identity: `extract_huge.py` -> huge_network.json (865 huge atoms; 512 simple
   loads bit*(x_B-HUGE)=s*x_C). Check whether x_9770(A) and x_18274(B) are combinations of the
   SAME HUGE residues => matching becomes combinatorial, not brute 2^233.
3. MITM/lattice via x_8821 (the linear coordinate on the 233 side) — see NOTEBOOK Session 6.

## Exhausted this session (do NOT redo)
- SAT/SMT (user directive: custom heuristics only; z3/cvc5 return unknown anyway).
- v4 evaluator / anything freezing x_18274 (fixed in v5).
- Linear algebra / lattice: `linalg255.py` (CORRECT, over all 255 bits, mod-P) has RANK 255/255
  and forces ALL bits = 0. The witness (!= all-0) is OUTSIDE the linear neighborhood of all-0, so
  linear/lattice/subset-sum provably cannot reach it. Supersedes the Session-5 "slaved-B" claim.
- B=0: ruled out (full 2^22 scan, 0 matches).
- Modulus (gcd residues=1), residue-lattice relation (none), slack vars (x_26977/x_9982 rigid).
- Local search / greedy / SA / pairs / triples from all-0 — all plateau (all-0 is the local min).

## The ONE remaining avenue (unimplemented)
A custom NONLINEAR solver / backward circuit-inversion: pin x_18274=N1, x_17728=N2 for a chosen
(N1,N2) from the 2^22 table S, and propagate/search backward through the 233-side acyclic circuit
(residue-load selects + product/sum gates) to determine the bits. Big build, uncertain (z3 failed
the analogous forward CSP). This is the only path not proven dead — everything else is exhausted.
The instance is a genuine obfuscated-circuit trapdoor; a full witness likely needs the setter's
secret or a cryptanalytic break of the specific 233-side residue circuit.

## Git
Branch `claude/read-prompt-5t2raw`. Commit+push after meaningful experiments.
