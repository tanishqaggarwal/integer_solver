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
