# Agent F log — p-adic / multi-modular lifting angle

## 2026-08-07 start
- Verified baseline: `python3 solve_lab/checker.py solve_lab/best/new_instance_partial_39026.json`
  -> `satisfied 39026/39033 (7 failing)`, failing `[12231,12270,12350,14584,18673,22044,29125]`. CONFIRMED.
- Built independent parser (`parse.py`,`parse2.py`,`parse3.py`) from raw text; no reliance on prior lab code.
  - 39,033 equations. Each LHS is a scalar multiple / square / (c1+c2) multiple of a single core S.
  - S is a left-nested spine `A0 + c1*A1 + c2*A2 + ...`, up to 26 atoms, small coefs.
  - 96,883 distinct atoms; 31 atom SHAPES only. Dominant shapes:
    (X-(X*X)) 160719, (X-(X+X)) 56409, (X-X) 45987, X 43773, (X-(C*X)) 28396,
    (X-(X-X)) 28082, (X-0) 20417, ((X*X)-X) 20242, (X-(1-X)) 17207, (X*X) 16413,
    (X*(1-X)) 14964, (X*(X-1)) 14950, (X-1) 13369, ...
  - So the instance is a straight-line circuit: most atoms are *definitions* x_out - f(inputs);
    plus boolean atoms X*(X-1) and constant pins (X-0),(X-1).

## Structural decomposition (independent, verified)
- Each equation LHS = (scalar) * S^k  or  (c1+c2)*S; every equation <=> core S = 0.
- S = left-nested spine of atoms with small integer coefficients.
- 39,033 distinct atoms. 30,001 are *definitions* x_out = f(inputs); the definition graph is a **DAG**
  (greedy topological schedule succeeds, 0 cycles). 8,747 free inputs.
- Residual system after forward elimination: **9,032 atoms** (4,621 redundant defs + 4,411 non-def atoms).
  All atoms zero => all 39,033 equations hold.
- `fwd.py` Engine: forward pass = 16 ms, equation score matches checker.py exactly (validated at 39,001).
- **All free inputs = 0 gives score 39,005 with only THREE nonzero residual atoms.**
    (x24468-K1)-x32989 ; (x2300-x9274) ; 8863713*(x18956-K2)-x14257
- p = 115792089237316195423570985008687907853269984665640564039457584007908834671663 (secp256k1 field prime)
  appears as literal x26064; handle vars enter as p*h. CONFIRMS prior lab's p.
- The three atoms decode to: OR(a,b)=1 with a=x7715,b=x34554 selectors, and the SELECTED coordinate
  pair must be (K1,K2) mod p. Selecting (a,b)=(1,1) picks the free pair (x22162,x30213).
- Doing that (turn one bit in each OR tree + set x22162=K1,x30213=K2,x24468=K1,x18956=K2) -> **39,013**,
  4 nonzero atoms, all "conditional constant pins" b*(x-C)-M*h of the two on-bits.
- Jacobian (exact, integer, by probing): 5,747 affine knobs, 260 NON-affine (the 254 OR-tree bits + 6 EC
  coordinates, which are quadratic). Residual system is near-diagonal: 5,715 columns of degree 1.
  Bipartite components: 3,669 singletons + four components of size 8-9 containing the broken atoms.
- Those components are INCONSISTENT in the affine linearization *because the coordinate knob is excluded*
  (quadratic). This is the real coupling.
- Prior 39,026 deliverable: its free inputs forward-close to score 39,020 (4 vars differ => they use a
  non-circuit-consistent state where broken atoms cancel inside equations). Its 4 broken atoms are
  3130,3131,3132,7251 -- a DIFFERENT set from mine.

## The instance decoded (agent F, independent)
Branch structure: a=x7715=OR(178 boolean free bits), b=x34554=OR(78 boolean free bits); all 256 bits are
boolean-constrained.  x9274=OR(a,b) must be 1.  Selected coordinate pair:
  (a,b)=(1,0) -> (x1,y1)=(x12186,x16742); (0,1) -> (x2,y2)=(x14853,x24908); (1,1) -> (x3,y3)=(x22162,x30213).
Unconditional pins force  selected_x = K1 (mod p),  selected_y = K2 (mod p).
Each ON bit j forces, through a rigid mod-p chain, TWO coordinate values: tree-a bit -> (x1,y1),
tree-b bit -> (x2,y2), and the values are LITERALLY the bit's two pin constants mod p (verified on 9 bits).
Two bits on in the same tree => contradiction (verified empirically).
EC layer (only active when x15298=ab=1): three residual atoms vanish iff
    A = (x2-x1)^2*(x3+x1+x2+K) - (y2-y1)^2 = 0   and   B = (y3+y1)(x2-x1) - (y2-y1)(x1-x3) = 0
with K = 97553848499418123410591666447050222001188385549510401465815187079080512838891 (measured, universal;
verified constant across 8 bit pairs and random coordinate perturbations; y-law offset exactly 0).

### Exhaustive enumeration (all-atoms-zero model)
- (1,0)/(0,1) branch: need a bit whose pin constants are (K1 mod p, K2 mod p). NONE of the 509 pin constants
  equals K1 or K2 mod p.
- (1,1) branch: need chord_K(P1,P2) = (K1,K2) over 178x78 bit pairs, both pin orderings (350x156 ordered
  candidates). ZERO hits -- not even the x-equation alone.
- degenerate P1=P2: requires a tree-a bit and a tree-b bit with identical pin constants. The 353 tree-a
  constants and 156 tree-b constants are DISJOINT.
=> Under "all residual atoms zero" the instance is INFEASIBLE.  (Not a claim about the instance itself:
   equations are linear combinations of atoms, so nonzero atoms may still cancel.)

### Reproduced/beat prior states
- Prior 39,026 uses bits x24601 (tree a) + x2081 (tree b), branch (1,1), and sets (x2,y2)=(x1,y1) by force,
  breaking bit 2081's two chain atoms (3132, 7251).
- Forward-closing their free inputs: 39,020 / 4 broken atoms. Fixing atoms 3130 and 3131 (each individually
  solvable: knob x8731; and 5113045*d(x9118) - p*d(x1329) = rhs) gives **39,022 with only 2 broken atoms**
  (`T2.json`, checker-verified 39022).  Their extra 2 broken atoms buy 4 more equations by cancellation.
- With only atoms {3132,7251} nonzero, 11 equations touch them; 2 contain 3132 alone, and the 9 mixed rows
  have 9 distinct coefficient ratios => at most 1 can be cancelled => max 39,023 for that residual pattern.

## Cancellation layer
- Prior 39,026 is NOT circuit-consistent: it has **7 nonzero atoms**, three residual and four DEFINITION
  atoms (x28730-(x17499*x9413), x29854-(x22665*x1329), x31864-(x28961*x10903), x642-(x28599*x17325)).
  So it deliberately breaks definitions to buy cancellation. My all-atoms-forward model is a restriction.
- Break relocation experiment: 64 placements of the two chain breaks (`breaks.py` + exact integer HNF solver
  `intsolve.py`); best is the original position (39,022). Others: 39,021, 39,021, 39,018, ...
- No atom's equation footprint is contained in the 11 equations touched by the two forced atoms, so no
  free dial helps without a joint solve.
- Therefore beating 39,026 requires a *syndrome-decoding* search over realizable atom vectors, not more
  circuit reasoning.

## Multi-modular / p-adic results (agent F core deliverable)
1. **Curve identification (independent).** Shift X = x + K/3 mod p with the measured K.  Then 253 of the 256
   conditional-pin bits give a point (X_j, y_j) with y^2 - X^3 EQUAL to the same constant
   b = 64019533680030876408443198762210829058751700634554282185987325820393598524794,
   and the TARGET point (K1 + K/3, K2) lies on the SAME curve.  So the object is
      E_b : y^2 = x^3 + b  over  F_p,  p = secp256k1 field prime  (a sextic twist of secp256k1, j=0).
   The correct pin ordering (which constant is x, which is y) is thereby determined for every bit.
2. **Group order, computed and verified independently.**  From 4p = t^2 + 27 m^2 with t = p+1-n_secp256k1
   the six j=0 twist orders are enumerated; exactly one annihilates points of E_b:
      N = 115792089237316195423570985008687907852837564279074904382605163141518161494337   (256 bits)
   Verified: annihilates 5 random points of E_b; no other candidate order does.
   **N is PRIME**: no factor below 10^6, sympy.isprime = True, Miller-Rabin bases 2..100 all pass.
   => Pohlig-Hellman is worthless on E_b; there is no small-order or smooth structure to exploit.
   Quadratic twist order = 3^2 * 13^2 * 3319 * 22639 * (211-bit prime) -- also useless (and irrelevant:
   the target is on E_b, and a static instance admits no invalid-curve interaction).
   N-1 = 2^6 * 3 * 149 * 631 * (231-bit prime).
3. **The obstruction is EXACTLY at p — proved computationally.**  `modq.py` builds a forward engine over
   Z/qZ and repairs residuals mod q.  For q in {4093, 7919, 65537, 104729, 1000003, 15485863, 999999937,
   1000000007, 2147483647, 2^61-1, 10^12+39} the FULL system is solved: **0 nonzero residual atoms and 0 of
   39,033 equations failing mod q**, in ~4 s each.
   Structural reason: every pin/chain check has the shape  A - B = M*p*h  with h a free variable.  Mod any
   q coprime to M*p the handle term M*p*h is surjective, so every such check is VACUOUS mod q; only the
   boolean constraints survive and they are satisfied by any single-bit-on configuration.
   At q = p the handle term dies identically and the pins bite -- that is where all the content is.
   => I's prediction "no mod-p freedom except the selector bits; the lift obstruction concentrates at p"
      is CONFIRMED by direct construction, not by inference.
4. **My 39,022 residual is the same point-addition object.**  Atoms 3132 `(x4432-x19964)-x28730` and
   7251 `(x7068-x2099)-(7376877*x642)` are the two chain rows that force (x2,y2) to bit x2081's pinned
   point; the deliverable state instead sets (x2,y2)=(x1,y1) (the degenerate doubling that makes
   A = B = 0 identically) and pays with those same two chain atoms.
