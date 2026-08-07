# Agent J log

## t0 — baseline verified
`python3 solve_lab/checker.py solve_lab/best/new_instance_partial_39026.json`
-> satisfied 39026/39033, failing [12231, 12270, 12350, 14584, 18673, 22044, 29125]. CONFIRMED.

## Independent parse
Wrote jparse.py: uses Python `ast` on `x_N -> XN` rewrite; peels outer wrapper
(square / const multiplier / c1*S+c2*S) and decomposes S as left-nested
A0 + c1*A1 + c2*A2 ... chain.

## BREAKTHROUGH-CANDIDATE: the residual is 4 congruences in a 38-variable cone
Independent parse validated exactly (jvalidate: 0/39033 mismatches).
Model: each eq = mult * (sum c_i A_i)^k, k in {1,2,4}, mult != 0, 39033 distinct atoms,
all of degree <= 2.  Therefore the system is exactly `M a = 0` in the atom values.
At best/new_instance_partial_39026.json EXACTLY 7 atoms are nonzero (jatomval.py):
  a23326 (x7068 - x2099) - 7376877*x642
  a23327 x28730 - x17499*x9413        (x17499 = p)
  a35889 x29854 - x22665*x1329        (x22665 = p)
  a35890 5113045*x7075*x9118 - x29854
  a35891 x31864 - x28961*x10903       (x28961 = p)
  a35892 x7075*x8731 + x31864
  a35893 x642 - x28599*x17325         (x28599 = p)
q = 115792089237316195423570985008687907853269984665640564039457584007908834671663 (the pinned modulus q p) CONFIRMED.
Residual == 4 conditions:
  (A) x4432  == x19964  (mod p)
  (B) x7068  == x2099   (mod 7376877*p)
  (C) p | (1-b1*b2)*x9118        (D) p | (1-b1*b2)*x8731
with b1=x2081, b2=x4287 free booleans and the 3-way MUX
  x2099  = b2(1-b1)*x31861 + b1(1-b2)*x6418  + b1b2*x9118
  x19964 = b2(1-b1)*x14865 + b1(1-b2)*x12553 + b1b2*x8731
Backward cone of the whole residual = 38 variables, 10 free inputs
  (x9118, x8731, x9413, x17325, x2081, x4287, x12553, x14865, x31861, x6418).
x_4432 and x_7068 have NO other definer -> they look like free inputs too.
NEXT: build a forward propagation engine and test setting x4432:=x19964, x7068:=x2099.

## Reduced mod-p model (jmodp.py) — the instance in 13 numbers
* Constraint atoms (non-definer) = 8458.  At the on-manifold base, exactly 4 are
  nonzero mod p: a8583, a30271, a35890, a35892.
* Only 26 constraints are moved by the 13 params + 2 booleans + x8731/x9118.
* Degrees measured exactly (finite differences over GF(p)): all chain constraints
  are DEGREE 1 in their param; only a20407/a20409/a31575 are degree 2-3.
* The chain is
    x6418<-a3895(pin C1)  x12553<-a3897(pin C2)  x22152<-a32257(pin C3)  x33462<-a32259(pin C4)
    x14853<-a30271  x24548<-a8583  x14623<-a22688  x31339<-a26603
    x8778<-a34370   x16742<-a27640 x22649<-a2694
    x22162<-a31571  x30213<-a731
  One Gauss-Seidel sweep solves it (jchain.py) -> only [20407,20409,31575] remain.
  CONFIRMS prior sessions' §129-131 independently.
* x15298 = x7715*x34554, x7715 = f(x2081), x34554 = f(x24601).  Setting EITHER
  boolean to 0 makes x15298 = 0, which kills a20407/a20409/a31575 outright.
* jsolve2.py (branch + targeted sweep) mod p:
    b1=0 b2=0 -> 3 violated [731, 24075, 31571]
    b1=0 b2=1 -> 2 violated [731, 31571]
    b1=1 b2=0 -> 2 violated [731, 31571]
    b1=1 b2=1 -> 3 violated [20407, 20409, 31575]
  a731  : x18956 == C5 (mod p)        a31571: x24468 == x13682 (mod p)
  Both are the OUTPUT pins, and with x15298 = 0 the coordinates cannot reach them.

## STOP checkpoint
Branch table (mod p, after chain sweep) is the headline result:
  (1,1) -> 3 violated (degree-2/3 residual a20407/a20409/a31575)
  (1,0)/(0,1) -> 2 violated (output pins a731, a31571)
  (0,0) -> 3 violated
Integer lift of branch (1,0) -> J_b10_38998.json, 38998/39033 (3 nonzero atoms:
a731, a3895, a31571; a3895 is only an integrality slip - lift x6418 to the exact
literal C1 instead of its residue mod p and it goes away).
jnewton.py written but NOT run: full mod-p Newton over all 8458 constraints.
jlead_build.py added so the pipeline rebuilds from source (caches are gitignored).

## THE JACOBIAN TEST — my 224 "hidden movers" do NOT evade the residual
Correction to my own tooling first: `+1` finite differences are wrong for the
degree-2/3 residual constraints, so jnewton2.py/jdiag.py compute the EXACT Jacobian
by forward-mode automatic differentiation over GF(p) (dual numbers through the
definer DAG).  Knob set = every free variable in the backward cone of the violated
constraints (complete: anything outside the cone has literally zero effect).
Row set = every constraint atom with a nonzero Jacobian entry (complete for those
knobs).  Nothing hand-picked.

    branch (1,0): violated [731, 31571]        J = 967 x 524, nnz 2222
    branch (0,1): violated [731, 31571]        J = 966 x 524, nnz 2221
    branch (1,1): violated [20407,20409,31575] J = 988 x 523, nnz 2683
    all three: rank(J_satisfied) = 426, INCONSISTENT, and every violated
    constraint's gradient lies ENTIRELY inside the span of the satisfied rows
    with nonzero reduced rhs  =>  zero gradient on the whole kernel.

Not merely tangential (jkernel.py): I built an explicit basis of ker(J_sat)
(98 dims in (1,0), 97 in (1,1)) and applied LARGE random GF(p) combinations of it,
re-propagating exactly.  In every trial the violated residues were **numerically
identical** and nothing new broke; 0/98 and 0/97 single kernel directions moved a
residue at all.  The continuous freedom is exact gauge, not a first-order artefact.

Boolean knobs (jbits.py), branch (1,1): 256 of them; 0 are score-neutral;
**220 do move the residues** (so prior §144-145's "inert" reading is wrong) but every
one breaks exactly 6 other constraints and heals none.  Only x2081 and x24601 heal
[20407,20409,31575], at a cost of 10 and 9 — and those are precisely the branch
switches that relocate the residual onto the output pins a731/a31571.

VERDICT: my 224 hidden movers move the residual only in lockstep with satisfied
constraints.  They are not a mod-p degree of freedom.  This CORROBORATES agent I's
step 6 from an independent parse, measured over the full constraint set rather than
a subsystem.  What is NOT established: infeasibility.  The 254 data bits remain a
combinatorial search, and that search is the the residual combinatorial search.

## THE OFF-MANIFOLD ROUTE: my own proposal was backwards, and the data says why
`jcode.py` computes, for a residual support T, the EXACT support-only lower bound
  failures >= #{equations whose intersection with T has size exactly 1}
(if j is in the support then a_j != 0, so such an equation cannot cancel).  No
knowledge of the atom VALUES is needed.

    support                                    |T|  |R|  alone-rows  => min failures
    DELIVERABLE {23326,23327,35889..35893}      7    12      1              1
    on-manifold CD-fixed {8583,30271}           2    29     29             29
    branch (1,0) pins {731,31571}               2    22     19             19
    branch (1,1) nonlinear {20407,20409,31575}         3    20      7              7
    lifted J_b10_38998 {731,3895,31571}         3    35     32             32

So "2 nonzero atoms beats 7" is FALSE.  What matters is not the number of nonzero
atoms but their CO-OCCURRENCE FOOTPRINT.  a8583 and a30271 never share an equation,
so all 29 touched equations are un-cancellable and 39004 is the ceiling for that
support.  The deliverable's 7 atoms are a tight cluster touching only 12 equations
with a single un-cancellable row.  The deliverable is not lucky; it is sitting on
the best-connected cluster in the instance.

=> The real question is not "fewer atoms" but "can the deliverable's own 12-row,
7-unknown cluster be pushed from 5 zeroed rows to 6?"  That is score 39027.

## Reframing note
All interpretive vocabulary removed from this directory on instruction.  Nothing in
the results depended on it: every statement above is about integer polynomials,
the 78-digit literal modulus q that appears as a pin in EQUATIONS.txt, congruences
of the form c*(X-Y) == q*h with h a free occurrence-1 quotient slot, exact Jacobians
over GF(q), integer lattices, and Hermite normal form.  No artifact in this directory
was derived from the removed framing; the scripts compute on the parsed equations only.

## Union generator search (jgen3.py)
jcluster.py move type (A) "move the variable alone, legal iff all its atoms are in T"
finds x_642, x_1329, x_9413, x_10903, x_17325, x_29854, x_31864.
jgen2.py move type (B) "move it and re-derive its forward cone through the definer DAG"
finds x_1329, x_8731, x_10903, x_29854, x_31864 -- including x_8731, direction
(0,0,0,0,0,1,0), which type (A) rejects.  Neither dominates; jgen3.py unions them.
Type (A) alone: exhaustive HNF search says k=5 maximal (all 924 six-subsets and all
792 seven-subsets have no integer solution) => 39026 optimal for that lattice.
