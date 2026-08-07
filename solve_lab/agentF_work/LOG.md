# Agent F log — multi-modular / p-adic angle on EQUATIONS.txt

All results below are derived from my own parse of the raw file; nothing imported from other agents.
Every score claim is verified by `solve_lab/checker.py`.

## 1. Baseline re-verification
`python3 solve_lab/checker.py solve_lab/best/new_instance_partial_39026.json`
-> `satisfied 39026/39033 (7 failing)`, failing lines `[12231,12270,12350,14584,18673,22044,29125]`. CONFIRMED.

## 2. Independent parse and the atom model  (`parse.py`,`parse2.py`,`parse3.py`,`core2.py`,`circ4.py`)
- Every equation LHS is `(integer scalar) * S^k` or `(c1+c2)*S`; so each equation holds iff its **core S** is 0.
- S is a left-nested spine `A0 + c1*A1 + c2*A2 + ...` (<=26 terms, coefficients in [-40,40]).
- There are exactly **39,033 distinct atoms** in only 31 syntactic shapes.
- 30,001 atoms are *definitions* `x_out - f(inputs)`; the definition graph is a **DAG** (`sched.py`: greedy
  topological schedule places all 30,001, zero cycles), leaving **8,747 free inputs**.
- Residual system after forward elimination: **9,032 atoms** (4,621 redundant definitions + 4,411 others).
  **All residual atoms zero  =>  all 39,033 equations hold.**  (One-directional; kept as stated.)
- Residual atom census: 3,958 booleans `x(x-1)`, 927 divisibility checks `A - M*h`, 509 conditional pins
  `b*(x-C) - M*h`, 2 unconditional constant pins, the rest wire copies/adders.
- `fwd.py` Engine: forward pass 16 ms; its equation score matches checker.py exactly (validated at 39,001,
  39,013, 39,022, 39,023, 39,024).

## 3. Where the difficulty sits, in purely integer terms
- Setting **all 8,747 free inputs to 0** gives **39,005** with only **THREE** nonzero residual atoms:
      (x24468 - K1) - x32989
      (x2300 - x9274)
      8863713*(x18956 - K2) - x14257
  with the two literal constants (296 bits each)
      K1 = 91416258160755509149180373473728639746431157665678710450404458852172057265575180278101002
      K2 = 125787314747601108116039725163361763116550465675981151838811516827327919228823597744635626
- The literal **p = 115792089237316195423570985008687907853269984665640564039457584007908834671663** occurs in
  the file as the constant defining x26064 (prime; 256 bits).  Every "handle" variable is p times a free
  input, so every check has the exact shape
        A - B = M * p * h,     h a free integer, M a ~24-bit literal.
  Hence each check is the congruence  A ≡ B (mod M*p).
- Decoding the three broken atoms: x9274 = a + b - a*b where a,b are outputs of two OR-trees over 178 and
  78 boolean free inputs; x2300 is the literal 1, so a=1 or b=1.  x24468 and x18956 are selected by
  a,b from three wire pairs, and the two unconditional pins then force
        selected_x ≡ K1 (mod p),   selected_y ≡ K2 (mod p).
- Turning on one boolean in each tree and setting the selected pair to (K1,K2) gives **39,013** with 4
  nonzero atoms, all of them the conditional pins `b*(w - C) - M*p*h` of the two ON booleans.
- Each ON boolean j forces, through a rigid chain of congruences mod M*p, two specific wires to two
  literal constants C of the file; those chains terminate at the four "selected" wires.  So a full
  all-atoms-zero solution requires the four selected wire values to be four of the 509 pin constants, tied
  together by one further degree-3 polynomial congruence mod p that the circuit checks.

## 4. Exact integer Jacobian and lattice tooling
- `jac.py` builds the exact integer Jacobian of all 9,032 residual atoms w.r.t. all 8,747 free inputs by
  probing (about 100 s).  At the 39,013 state: 5,747 affine knobs, 260 non-affine, 5,779 nonzeros -- the
  matrix is almost diagonal (5,715 columns of degree 1).
- `intsolve.py` solves small integer systems exactly by column HNF (unit-tested on 200 random systems).
- `lin.py`/`breaks.py`/`comp.py` decompose the residual system into bipartite components and solve them.

## 5. Where the lift obstruction actually concentrates  (the core multi-modular result)
`modq.py` builds a forward engine over Z/qZ and repairs residuals there.
**For every prime q != p tested, the FULL system of 39,033 equations is solved exactly mod q**, with
0 nonzero residual atoms and 0 failing equations:
    q = 4093, 7919, 65537, 104729, 1000003, 15485863, 999999937, 10^9+7, 2^31-1, 2^61-1, 10^12+39   (~4 s)
    q = 2^127-1  (4.2 s)   and   q = 2^255-19  (4.4 s)      <- size control, same bit-length as p
So the obstruction is not about modulus size: at any q coprime to M*p the term `M*p*h` is surjective, so
every congruence check is vacuous mod q and only the boolean constraints survive (satisfied by any
single-boolean-on configuration).  At q = p the handle term vanishes identically and the pins bite.
**The lift obstruction is concentrated at the single modulus p.**  The mod-p solve does not converge
(timed out at 900 s), consistent with all the content being there.

## 6. The lift above p^1 is unobstructed
`relaxed_pin.json`: if the two unconditional constant pins are relaxed (the selected pair set to the values
the degree-3 check demands instead of (K1,K2)) then everything else closes **exactly over Z** -- the free
handles absorb every quotient with no residue left over.  So there is no p-adic obstruction at p^2 or
higher: a solution of the single mod-p condition would lift to Z immediately.  The entire instance is
one equality mod p.

## 7. Exhaustive search of the all-atoms-zero model
- Two ON booleans in the same OR-tree are impossible: with two of them on, the exact integer Jacobian gives
  six broken atoms in six components, four of size 8 that are **integrally unsolvable** (verified with
  `intsolve.solve_int`), and only two selected-wire knobs are available to close four competing chains.
- With one ON boolean per tree, the four selected wire values are determined mod p by the pin constants
  (verified on 9 different boolean choices: the selected values equal the pin constants mod p exactly).
- Enumerating **all 178 x 78 = 13,884** combinations against the degree-3 congruence: **zero solutions**.
  Also zero for the degenerate cases.
=> No assignment makes all 39,033 residual atoms vanish.
**Caveat, explicit**: this does NOT prove the instance infeasible.  Equations are integer linear
combinations of atoms, so a full solution could carry nonzero atoms that cancel.  That would need a
nonzero realisable atom vector in ker(M), M = the 39,033 x 39,033 equation-atom incidence matrix
(525,982 nonzeros).  I did not compute rank(M).  **No infeasibility claim is made.**

## 8. Cancellation / frame analysis (how the 39,026 deliverable actually works)
- The 39,026 state is NOT forward-consistent: it carries **7 nonzero atoms**, three residual and four
  broken *definitions*, and only 12 equations touch them, of which 5 are satisfied by cancellation.
- Its reachable residual lattice has 7 private generators.  Enumerating all 2^12 row subsets with exact
  integer solving: **exactly one 5-subset is integrally reachable** (all 924 six-subsets are solvable over
  Q -- integrality is the obstruction).  So 39,026 is optimal *in that frame*.  Prior claim CONFIRMED.
- Every single-atom extension of that frame (7 candidates) still gives 7 failures.
- Forward-closing the 39,026 free inputs gives 39,020 with 4 broken atoms; fixing two of them individually
  gives **39,022 with only 2 broken atoms** (`best_F_39022.json`); relocating the break into a definition
  and re-solving the frame gives **39,023** (`F_frame.json`) and **39,024** (`F_best.json`).
  The coset of the residual value inside the lattice is what decides how many rows cancel.

## 9. Multi-modular sweep with checkpointing  (`modm.py`, results in `modm_results/`)
Every result below = "solve the FULL 39,033-equation system modulo m; count nonzero residual atoms and
failing equations".  Each prime/prime-power is checkpointed to its own JSON as it completes.
- **Primes**: all 25 primes < 110, plus 1009, 10007, 100003, 1000003, 10000019, 10^9+7, 2^31-1, 2^61-1,
  2^89-1, 2^127-1, 2^255-19  ->  **all SOLVED, 0 nonzero atoms, 0 failing equations**, ~3-4 s each.
- **Prime powers**: 2^64, 2^100, 3^40, 3^80, 5^10, 5^20, 5^30, 5^40, 7^25, 11^20, 13^20, 1009^8, 65537^4,
  1000003^3, (2^31-1)^2  ->  **all SOLVED, 0 failing equations**.
  (5^k initially failed with one choice of ON booleans -- a pivot in my greedy repair happened to be
  divisible by 5.  With ON booleans (24601, 2081) every 5^k solves cleanly.  Not an obstruction.)
- **Moduli equal to the handle multipliers themselves** (M = 9027329, 13921279, 8422691, 2818759,
  13818563, 1001745, 7376877, 3*7376877): **all SOLVED, 0 failing equations**.  So even the moduli that
  are NOT coprime to the handle coefficient are unobstructed -- the extra factor p in `M*p*h` still makes
  the term surjective there.
=> **The integer lift is obstructed at exactly ONE modulus, p, and only at level p^1.**
   Everywhere else the system is not merely solvable but solvable in ~3 seconds by naive propagation.

## 10. Handle structure validated at scale
Restricting to the genuine divisibility shapes among the 9,032 residual atoms -- 927 of the form
`A - (M * x_h)` and 3,201 of the form `A - x_h` -- and evaluating x_h under 4 independent random draws of
all 8,747 free inputs: **3,173 of them have x_h ≡ 0 (mod p) in every draw**.  The 955 exceptions are all
mis-classified boolean atoms `(X*X)-X` (788) and wire copies `(X-X)` (167), not handles.  In particular
both unconditional constant pins and every chain row have p-multiple handles.
=> The "every check is a congruence mod M*p with a free quotient" reading is confirmed structurally,
   not just on the handful of atoms visible in the Jacobian.

## 11. CRT reconstruction  (`crt.py`, per-prime vectors checkpointed in `crt_sols/`)
Solved the full system independently modulo 20 distinct 7-digit primes (each: 0 nonzero atoms, 0 failing
equations, ~3.5 s), then CRT-combined the 20 solution vectors coordinatewise into a single assignment
modulo Q = product of the 20 primes, **Q of 399 bits**.
Verification with the mod-Q engine: **0 nonzero residual atoms and 0 of 39,033 equations failing mod Q.**
(That is exactly what the ring isomorphism Z/Q = prod Z/q_i predicts, and it confirms the per-prime solves
are genuine.)  So the system has solutions modulo arbitrarily large integers coprime to p.
Taking balanced representatives of the mod-Q solution gives an integer assignment scoring **38,991/39,033**
under checker.py (`crt_balanced.json`) -- i.e. CRT reconstruction alone does NOT lift, because the mod-q
solution spaces are enormous and the 20 chosen solutions are mutually unrelated.  A lift would need the
per-prime solves steered toward a common integer candidate, and the only thing that pins such a candidate
is the mod-p condition, which is exactly the obstructed one.

## 12. Summary of the multi-modular picture
- Solvable, in seconds, modulo: every prime tested from 3 to 2^255-19; every prime power tested
  (2^100, 3^80, 5^40, 7^25, 11^20, 13^20, 1009^8, 65537^4, 1000003^3, (2^31-1)^2); every handle
  multiplier M tested; and modulo a 399-bit composite built by CRT.
- Obstructed at exactly one modulus, p, at level 1 only (a mod-p solution lifts to Z with no further
  work -- demonstrated by `relaxed_pin.json`).
- The mod-p condition itself is a single degree-3 congruence tying four of the file's 509 literal pin
  constants to the two 296-bit literals K1, K2; all 13,884 admissible constant combinations were checked
  and none satisfies it.  No infeasibility of the instance is claimed -- see the ker(M) caveat in §7.

## 13. Deliverables (all checker-verified)
| file | score | note |
|---|---|---|
| `solve_lab/best/new_instance_partial_39026.json` | 39,026 | shared baseline, re-verified by me |
| `agentF_work/best_F_39024.json` | 39,024 | best produced by my own pipeline end-to-end |
| `agentF_work/F_frame.json` | 39,023 | 2 nonzero atoms + break relocated into a definition |
| `agentF_work/best_F_39022.json` | 39,022 | minimal residual: only 2 nonzero atoms in the whole system |
| `agentF_work/crt_balanced.json` | 38,991 | balanced representative of the 399-bit mod-Q solution |
| `agentF_work/relaxed_pin.json` | 39,014 | everything closes over Z once the two constant pins are relaxed |

## 14. mod p and mod p^2 (the obstructed modulus)
Same greedy repair, same code path, ON booleans (24601, 2081):
    m = p    -> NOT solved: 21 nonzero atoms, 124 of 39,033 equations failing (33 s, 6 sweeps)
    m = p^2  -> NOT solved: 15 nonzero atoms, 102 failing (29 s, 4 sweeps)
versus 0 / 0 in ~3.5 s for EVERY other modulus tried (60 of them, up to 255 bits).  Checkpointed as
`modm_results/P_p.json`, `modm_results/P_pp2.json`.  This is the sharpest single measurement in the lab:
the difficulty of EQUATIONS.txt is localised at one prime literal.

## 15. **ker(M) = 0 — the campaign's gate, answered exactly**
`buildM.py` assembles M (rows = equations, cols = atoms, coefficients from the spine decomposition, with
repeated atoms inside one equation merged): **39,033 x 39,033, 525,982 nonzeros**, row degree 1..24
(mean 13.48), column degree 1..22 (mean 13.48), max |coefficient| 80, no zero row or column.
Exactly ONE row and ONE column have degree 1.

`peel.py` / `peel_cert.py`: a **characteristic-free peeling argument**.  A row whose surviving support is a
single atom j with nonzero coefficient forces r_j = 0 over Z.  Starting from the single degree-1 row, the
cascade **forces all 39,033 atoms to zero** in 2 s.

    forced 39033 atoms of 39033      ->      rank(M) = 39033 ,   dim ker(M) = 0

The elimination order is saved as a **checkable certificate** (`peel_order.npy`, 39,033 pairs
(atom j, row i)) and re-verified by an independent pass that reloads M from disk and checks, for every
step, that row i contains atom j with a nonzero coefficient and that every OTHER atom of row i was
already forced zero at an earlier step.  **Certificate verified: True.**
This is exact and holds over Z and over every field of characteristic > 80 (the largest |coefficient|).

### Consequence
**Any integer assignment satisfying all 39,033 equations must make all 39,033 atoms exactly zero.**
So the "all-atoms-zero" model is not a restriction at all -- it is equivalent to a full solve.  Every
optimality statement in this lab that was conditional on that model is therefore unconditional in that
respect, and the cancelling-residual route to a full solution does not exist.
Combined with the exhaustion in section 7 (at most one ON boolean per OR-tree; the selected wire values
are then pin constants mod p; all 178 x 78 = 13,884 combinations checked against the degree-3 congruence,
zero solutions) this says the instance is infeasible -- **conditional only on the two measured (not
symbolically proved) links**: (i) that two ON booleans in the same OR-tree are always contradictory,
verified by exact integer Jacobian on one such pair, and (ii) that an ON boolean forces the selected wire
to its pin constant mod p, verified on 9 different boolean choices.  I state those as measurements.

## 16. Attempt to make the two measured links unconditional (partial)
`modp_uf2.py`: computed the closure Z of wires that are **provably == 0 (mod p) for every assignment**
(seeded by the literal p and closed through the definition DAG): **7,202 of 38,748 wires**.  Using Z, built
a weighted union-find over wires from every definition and residual atom that reduces mod p to
`alpha*u = beta*v` with alpha,beta units: 4,223 links from definitions + 352 from residual atoms,
**0 conflicts**.  Result: the 509 conditional-pin wires do NOT join the coordinate classes statically.
Reason (measured, not a bug): the path from a pin wire to the selected wire runs through the selector
products `x34606*x1 + x5647*x2 + x15298*x3`, which are bilinear and only become linear once the booleans
are fixed.  So links (i) and (ii) of section 15 are intrinsically configuration-dependent and cannot be
discharged by a static rigidity argument; they stay measurements (9 boolean choices for (ii), one exact
integer-Jacobian pair for (i)).

## 17. Wiedemann cross-check (independent algorithm)
`wiedemann.py`: Krylov sequence a_i = u^T M^i v over GF(q) for i < 2n, then Berlekamp-Massey; M is
singular over GF(q) iff x divides the sequence minimal polynomial, i.e. iff its trailing coefficient is 0.
**Implementation validated first on controls**: 4/4 randomly generated nonsingular 70x70 matrices reported
NONSINGULAR (BM degree 70, nonzero trailing coefficient) and 3/3 rank-69 matrices reported SINGULAR
(trailing coefficient exactly 0), with the true rank computed independently by dense elimination mod q.
Then run on the real M for q = 2^31-1 and q = 2147483629 (see `wiedemann.log`).

## 18. Sharpening the ker(M)=0 certificate
- The 39,033 pivot coefficients used by the cascade take only **2 distinct values, all of absolute value
  1 or 2** (1,144 of them are the even one).  So the whole argument is a unit/2-pivot cascade.
- Divisibility of the pivots: 0 of 39,033 for every odd prime tested (3, 5, 7, 11, 13, 2^31-1); 1,144 for 2.
  Therefore **rank(M) = 39,033 over Z and over every field of characteristic != 2** (the argument is
  simply silent in characteristic 2, which is irrelevant here).
- The cascade is order-independent: three randomized peeling disciplines each force 39,033/39,033 atoms.

## 19. Exhaustive check of measured link (ii)  (`link2_sweep.py`, checkpointed to `link2_results.json`)
For every one of the 256 conditional-pin booleans: turn it on (with a fixed partner in the other tree),
repair the chains, and compare the resulting selected wire pair with that boolean's two pin constants
mod p.  Runs ~65 s per boolean and checkpoints after each; every boolean completed so far matches.
Resume by re-running the script -- it skips booleans already recorded.

### Wiedemann result (independent confirmation)
    q = 2147483647 (2^31-1):  sequence 853 s, Berlekamp-Massey 690 s
    minimal-polynomial degree L = 39033 = n,  trailing coefficient = 268435456 != 0
    =>  M NONSINGULAR over GF(2^31-1),  rank = 39033,  dim ker = 0
Two independent methods (characteristic-free peeling certificate, and Wiedemann over a word prime) agree.
A second prime q = 2147483629 is running for redundancy (`wiedemann.log`).

## 20. **CORRECTION — link (i) is NOT established; treat the infeasibility argument as OPEN**
Wiedemann finished on the second prime as well:
    q = 2147483629 : minpoly degree 39033 = n, trailing coefficient 11716781 != 0 -> M NONSINGULAR
So **rank(M) = 39,033 and dim ker(M) = 0** is confirmed by three independent computations
(peeling certificate over Z; Wiedemann over 2^31-1; Wiedemann over 2147483629).  That part stands.

What does NOT stand is link (i).  I built a *rigorous* configuration-conditional mod-p rigidity engine
(`cfg_rigid2.py`): fix a boolean configuration, note that wires whose free-input support lies inside the
boolean free inputs then have configuration-determined values, so selector products linearise; then run an
AFFINE weighted union-find (value[x] = A*value[root] + B) over every definition and residual atom, with a
ZERO-closure fixed point.  On a single-boolean configuration it derives, with **0 conflicts**, that all six
selected wires are pinned to explicit constants mod p -- reproducing my earlier empirical measurements
digit for digit.  That makes link (ii) a derivation rather than a measurement.

But on SAME-TREE PAIRS the engine derives **no contradiction and no forcing**:
  - both booleans' pin wires are still forced to their own constants (4/4 checked), conflicts = 0;
  - x1 and y1 are simply NOT derived.
Diagnosis: exactly 108 wires lose their forcing when the second same-tree boolean turns on, and among the
23 whose definition inputs are still forced sits `x11317 := (x11532 + x14681)` -- an ADDER that carried
x1's whole value in the single-boolean configuration, with x14681 holding boolean 47's contribution and
x11532 provably zero.  With two booleans on, x11532 is no longer provably zero.
**That is the signature of an accumulator, not of a conflict.**  If the selected wire is a SUM of the
per-boolean contributions, then many booleans may be on simultaneously, the reachable set of the selected
pair is far larger than the 13,884 one-per-tree combinations I enumerated, and **my exhaustion in section 7
is invalid**.  The earlier exact-integer-Jacobian evidence for link (i) rested on one pair and excluded the
quadratic coordinate knobs -- the same weakness I criticised elsewhere.
**Therefore: no infeasibility is claimed.  The instance is OPEN.**  What is proved is only ker(M)=0.

## 21. **THE ACCUMULATOR IS REAL — the instance is a 96-STAGE COMBINATION TREE**
Census of residual atoms of the "gated check" shapes `(g*w)+h`, `(C*(g*w))-h`, `(g*w)-(M*h)`, grouped by
the gate g:  **96 gates carry exactly 3 checks each, and each gate's three checked wires depend on its OWN
distinct six-tuple of free inputs.**  96 gates, 96 distinct six-tuples, no overlap.
Gate-support sizes (how many of the 256 conditional-pin booleans each gate sees):
    0:7  1:20  2:32  3:9  4:11  6:2  7:4  8:1  9:1  10:1  11:1  14:1  21:2  22:1  50:1  88:1  256:1
That is a **binary combination tree**: the 256 pin constants are the leaves, 96 internal stages each
combine two values into one by the instance's own degree-3 law, and the root stage (gate x15298, support
256) produces the value that the two unconditional constant pins compare against (K1, K2).
Concretely, the second stage found earlier, gate x24533 (support 50, all tree-A booleans), owns the
six-tuple {736, 5186, 11532, 14681, 25591, 38551}; and the wires that lost their forcing when a second
same-tree boolean turned on -- `x11317 := x11532 + x14681`, `x751 := x736 - x38551`,
`x1622 := x25591 - x14681` -- are exactly that stage's own coordinates feeding the root stage's input x1.
So a second ON boolean does not contradict the first: it activates a subordinate stage whose OUTPUT
becomes the root stage's input.

### Consequence
**Section 7 of this log is invalid.**  My enumeration of "178 x 78 = 13,884 combinations" covered only
configurations in which a single stage is active with two LEAF inputs -- depth 0 of a 96-stage tree.
The reachable value of the selected pair is everything obtainable by composing the degree-3 law over
subtrees of the 96 stages with leaves drawn from the 256 pin constants.  That is exponentially large.
**No infeasibility is claimed; the instance is wide open, and the search should be run over stage
configurations, not over leaf pairs.**

## 22. **The composition law is CONFIRMED EXACTLY — stages compose with the same degree-3 law**
Decoded stage x24533 completely (mod-p rigidity, `cfg_rigid2.build`):
    input A = (x14681, x38551)   input B = (x25591, x736)   output = (x11532, x5186)
    gate x24533 = 1 exactly when BOTH inputs are live (two ON booleans below it); with one ON boolean the
    gate is 0, the stage is inert, and the single leaf value passes straight through to the parent.
Its three checks are LINEAR in the output pair; the 3x2 system has rank 2 and is CONSISTENT (third check
residual exactly 0), and the output it demands is
    x = 50819865790474622290929616283831020419846599045035819092742156686783125696627
    y = 1516899452477486833305684599119574148292933010271216764472923894680745180820
which equals **chordK(A,B) digit for digit**, with the SAME universal constant
K = 97553848499418123410591666447050222001188385549510401465815187079080512838891
that governs the root stage.  MATCH: True.

### Statement of the instance, corrected
EQUATIONS.txt encodes a **96-stage binary combination tree** over one degree-3 law
    out_x = l^2 - a_x - b_x - K ,  out_y = l*(a_x - out_x) - a_y ,  l = (b_y-a_y)/(b_x-a_x)   (mod p)
with the 256 conditional-pin constants as leaves, each stage gated by an AND of two OR-groups over its own
leaf subset, and the root (gate x15298, support 256) required to produce (K1 mod p, K2 mod p).
Everything else -- 30,001 definitions, 9,032 residual atoms, all the handles -- is bookkeeping that closes
exactly over Z once the mod-p tree evaluates correctly (sections 5, 6, 11).
**The accumulator composes: the reachable root value is exponential in the number of active leaves.**
My 13,884-combination enumeration was over depth-0 configurations only and says nothing about the instance.

## 23. Both sweeps FINISHED — final numbers
**Link (ii), EXHAUSTIVE (256/256 booleans, `sweep_ii.json`)**
  - contradictions derived: **0** (no boolean disagrees with anything);
  - **248 / 256** fully confirm: both selected wires derived and both equal that boolean's pin constants;
  - **8** only partially derived (7 give 1 of 2 coordinates, 1 gives 0) -- and every value that WAS derived
    is a pin constant of that boolean, so none of the 8 disagrees; 3 of them (18184, 22579, 33434) are the
    booleans that carry only ONE pin, so a second coordinate cannot be forced by design.
  Link (ii) therefore holds with **zero counterexamples**; the 8 are engine incompleteness, not conflict.

**Link (i): REFUTED as an argument (`sweep_i.json`)**  Same-tree pairs produce no contradiction; they
activate a subordinate stage.  The sweep's coverage number is "fraction tested", never "fraction excluded".

## 24. The stage law is UNIFORM across the whole tree
`stage_law2.py` searches, for every stage, all role partitions AND all coordinate orderings (which member
of a pair is x, which is y, which input comes first), solving the 3x2 linear system once per partition and
testing the 16 orderings arithmetically; 2 independent random input draws must agree.
    stages with a full six-tuple of free inputs analysed : 72
      admit a chord-with-offset law                      : 72   (100%)
      and the offset is the SAME UNIVERSAL K             : 72   (100%, zero exceptions)
      offset different from K                            : 0
      no consistent law                                  : 0
    stages skipped (six-tuple has only 4 free inputs -- leaf-adjacent stages, one input is a literal): 24
So all 96 stages run ONE law, with ONE constant
K = 97553848499418123410591666447050222001188385549510401465815187079080512838891.

## 25. FINAL STATEMENT OF THE INSTANCE
EQUATIONS.txt is a **96-stage binary combination tree of depth 6** over the single degree-3 law
    out_x = l^2 - a_x - b_x - K,  out_y = l*(a_x - out_x) - a_y,  l = (b_y - a_y)/(b_x - a_x)   (mod p)
with the **256 conditional-pin constants as leaves**, each stage gated by an AND of two OR-groups over its
own leaf subset (gate = 1 iff both of its inputs are live, otherwise the live input passes through), and
the root stage (gate x15298, support 256) required to produce (K1 mod p, K2 mod p).
Everything else closes exactly over Z: ker(M) = 0 forces all atoms to vanish, the definition graph is a
DAG, and the handles absorb every quotient (sections 5, 6, 11, 15).
**The task is therefore: choose a subset S of the 256 leaves so that folding the fixed tree gives the
target.  The reachable set is exponential in |S|.  Every enumeration in this campaign, mine included,
covered only depth-0 configurations.**

## 26. **THE LAW IS INVERTIBLE — YES.**  Meet-in-the-middle is available.
200/200 random triples, both directions, exact:
    given (input A, output O) -> input B :  l = (o_y + a_y)/(a_x - o_x);  b_x = l^2 - a_x - o_x - K;
                                            b_y = a_y + l*(b_x - a_x)
    given (input B, output O) -> input A :  l = (o_y + b_y)/(b_x - o_x);  a_x = l^2 - b_x - o_x - K;
                                            a_y = b_y + l*(a_x - b_x)
Each inversion is O(1) (one modular inverse).  So the target can be pushed DOWN the tree from the root as
easily as leaf values are pushed up, and a meet-in-the-middle over the root's two input slots replaces a
flat search: enumerate ON-subsets of one slot's leaf support forward, invert the target through the other
slot, and match in a hash table.  Because the inversion is exact and O(1) at EVERY stage, the same trick
applies at internal nodes, so the meeting point can be pushed deeper than the root.
**This single fact is what makes the tree attackable rather than merely understood.**

## 27. The 24 leaf-adjacent stages
17 have a six-tuple of 4 free inputs and 7 have only 2, with boolean support 0.  Their missing inputs are
hard-wired literals inside the circuit (the checked wires unfold to `(C1*u + C2*v)` forms over constants),
so those stages take one or two constant inputs rather than a live subtree.  They must be resolved before
the fold evaluator is complete; the unfolding is mechanical from the definition DAG.

## 28. What is NOT done (state honestly)
The **fold evaluator is not built and therefore not validated**, and no subset search was run.  What exists
is: the 96-stage tree (`tree96.json`), the uniform law and its universal constant, exact invertibility,
and the per-stage role assignments (`stage_roles.json`, produced by `stage_law2.py`).  The missing piece is
the WIRING: each stage input is a gated multiplexer (a sum of selector*value terms, as in
`x11317 := x11532 + x14681`), so the map "which stage output feeds which stage input slot" still has to be
read off the definition DAG.  Until that is done and the evaluator reproduces a known state, **no search
result from it should be believed.**

## 29. The wiring is a 3-WAY MULTIPLEXER per input slot
Each stage input slot is a free variable w constrained by a residual atom `((w - z) - handle)`, and z
unfolds to a sum of exactly **three** `selector * value` terms.  Decoded 39 of the 144 slots with the
current pattern (the other 105 use atom shapes the regex does not yet match).  Examples:
    stage x15298 slot inA (wires 12186,16742) -> (x24673*x10261) + (x2754*x30454) + (x38170*x5096)
    stage x30973 slot inA (wires  5460,38101) -> (x24454*x4622)  + (x31584*x7256) + (x30699*x27616)
    stage x30973 slot inB (wires 23971,28486) -> (x57*x29578)    + (x22550*x26655)+ (x11778*x20720)
This is the same shape as the root selection `x13682 = x34606*x1 + x5647*x2 + x15298*x3` found at the very
start of this lab: 3 gated candidates per slot, at most one live.
So the complete picture is: 96 stages, each with 2 input slots, each slot a 3-way gated mux over child
outputs or leaf constants, one uniform invertible degree-3 law with one universal constant K, root gate
x15298 required to produce (K1 mod p, K2 mod p).
Artifacts for the next session: `tree96.json` (stages, six-tuples, gate supports, containment depth),
`stage_roles.json` (per-stage input/output role assignment and coordinate ordering for all 72 full stages),
`sweep_ii.json`, `sweep_i.json`, `M.npz`, `peel_order.npy`.
