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
