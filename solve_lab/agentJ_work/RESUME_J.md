# RESUME_J — agent J

Everything below is a statement about integer polynomials in `EQUATIONS.txt`:
constants, congruences, linear relations. No interpretive framing is used or needed.

## Best verified
Baseline `solve_lab/best/new_instance_partial_39026.json` = **39026/39033**, re-verified
by me with `solve_lab/checker.py` (fails [12231,12270,12350,14584,18673,22044,29125]).
My own artifact `agentJ_work/J_b10_38998.json` checker-verifies at **38998**.
**No improvement on the deliverable.**

## Independence
Whole pipeline rebuilt from `EQUATIONS.txt`; no prior-session or other-agent artifact
was read into any computation.  Gate on everything: `jvalidate.py` reconstructs all
39,033 raw equations from my model and compares at a random point — **0 mismatches**.

## 1. The model (exact)
 * eq_i = mult_i * (sum_j c_ij A_j)^k, k in {1,2,4}, mult_i != 0 always; 39033 distinct
   atoms, every one of degree <= 2.  So the instance IS exactly `M a = 0` in the atom
   values, and **score = 39033 - ||M a||_0**.
 * A single 78-digit constant
   q = 115792089237316195423570985008687907853269984665640564039457584007908834671663
   appears as a literal pin; dozens of variables are copies of it.  Nearly every
   constraint has the shape `c*(X - Y) - q*h` with `h` a free variable of occurrence 1,
   i.e. an assertion `X == Y (mod q)` with a free quotient absorbing the integer part.
 * Definer DAG acyclic: 30575 defined vars, 8173 free, 8458 constraint atoms.
   Propagating the deliverable's free inputs leaves only **4** constraints violated.
 * Free inputs: 3747 have occurrence 1 (pure quotient slots), only 21 are nonzero, and
   only 17 are non-slot: 2 zero/one-valued (x2081, x24601), x8731/x9118 (freely
   zeroable — verified, costs nothing), and **13 integers of 295-296 bits**.
   Exactly **4 of the 13 are pinned to literal constants** (a3895, a3897, a32257,
   a32259), gated by x2081 / x24601.  The other 9 form a chain that is **degree 1** in
   each knob (verified by exact finite differences over GF(q)); one Gauss-Seidel sweep
   solves it, leaving only a20407/a20409/a31575 (degrees 3, 3, 2).

## 2. Branch structure
x15298 = x7715(x24601) * x34554(x2081).  Setting EITHER zero/one variable to 0 makes
x15298 = 0, which annihilates a20407/a20409/a31575 outright.  Mod q, after the sweep:
      (x2081,x24601) = (1,1) -> 3 violated [20407, 20409, 31575]
      (1,0) or (0,1)         -> **2 violated [731, 31571]**
      (0,0)                  -> 3 violated [731, 24075, 31571]
Mechanism: with x15298=1 the two pin constraints a731/a31571 read the FREE x30213 /
x22162 (so they are satisfiable) but the degree-2/3 residual is live; with x15298=0
they instead read x14853 / x12186, which the chain has already pinned.

## 3. The continuous freedom is EXACT GAUGE (not a tangent-space claim)
Tooling correction I had to make first: `+1` finite differences are wrong for the
degree-2/3 constraints, so `jnewton2.py`/`jdiag.py` compute the EXACT Jacobian by
forward-mode automatic differentiation over GF(q).  Knobs = every free variable in the
backward cone of the violated constraints (complete: anything outside the cone has
literally zero effect on them).  Rows = every constraint with a nonzero entry
(complete).  Nothing hand-picked.

    branch (1,0): violated [731,31571]         J = 967 x 524, nnz 2222
    branch (0,1): violated [731,31571]         J = 966 x 524, nnz 2221
    branch (1,1): violated [20407,20409,31575] J = 988 x 523, nnz 2683
    all three: rank(J_satisfied) = 426; INCONSISTENT; every violated constraint's
    gradient lies ENTIRELY in the span of the satisfied rows with nonzero reduced
    rhs => exactly zero gradient on the whole kernel.

Then the step this lab had never taken (`jkernel.py`): explicit basis of ker(J_sat) —
98 dims in (1,0), 97 in (1,1) — with LARGE random GF(q) combinations applied and
re-propagated exactly.  Every trial left the violated residues **numerically
identical** with nothing new broken; 0/98 and 0/97 single kernel directions moved a
residue at all.  So the continuous freedom is exact gauge, not a first-order artefact.

## 4. The discrete knobs have a FLAT cost profile
256 free variables carry a constraint `a*z^2 + b*z = 0`, i.e. z in {0, -b/a}.
Measured all 256 by hand (`jbits.py`): **0 are score-neutral**, and **220 do move the
residues** — which refutes the lab's long-standing "these are inert" reading.  But in
branch (1,1) every one of them breaks exactly 6 constraints and heals none, except
x2081/x24601 which heal all three at cost 10/9 and are simply the branch switches.
In branch (1,0) **no single one heals anything at all**.
Stated as landscape: the discrete knobs have a flat cost profile at 5-6 with no
gradient for a search to follow.  That is a fact about this system, independent of any
interpretation of what the system computes.

## 5. Off-manifold coding: my own "fewer atoms is better" proposal was BACKWARDS
`jcode.py` gives an exact support-only lower bound needing no atom values:
  failures >= #{equations meeting the support in exactly ONE atom}
(if j is in the support then a_j != 0, so such an equation cannot cancel).

    support                                    |T|  |R|  alone-rows  min failures
    DELIVERABLE {23326,23327,35889..35893}      7    12      1             1
    on-manifold 2-atom {8583,30271}             2    29     29            29
    branch (1,0) pins {731,31571}               2    22     19            19
    branch (1,1) residual {20407,20409,31575}   3    20      7             7
    lifted J_b10_38998 {731,3895,31571}         3    35     32            32

"2 nonzero atoms beats 7" is FALSE.  The binding quantity is not the number of nonzero
atoms but their CO-OCCURRENCE FOOTPRINT.  a8583 and a30271 never share an equation, so
all 29 touched equations are un-cancellable and 39004 is that support's ceiling.  The
deliverable's 7 atoms are a tight cluster touching only 12 equations with a single
un-cancellable row.  The deliverable is not lucky — it sits on the best-connected
cluster in the instance.

## 6. The real question, and my exhaustive answer so far
Within the deliverable's cluster: 12 equations, 7 atom values, score = 39033-(12-k).
Deliverable achieves k=5.  k=6 is 39027.  Reachable atom vectors = a* + L where L is
generated by the moves that wake no atom outside T.  Two independent move types
(move alone / move and re-derive the forward cone) find DIFFERENT generators, so they
must be unioned — `jcluster.py` (7 gens) and `jgen2.py` (5 gens, incl. x_8731 which
the strict test rejects), unioned in `jgen3.py`.
Exhaustive integral subset search (Hermite normal form) over the strict 7 generators:
**all 924 six-subsets and all 792 seven-subsets have NO integer solution; k=5 is
maximal.**  Confirms 39026 optimal *for that generator set* — the union search is the
live computation.

## Re-enter (caches gitignored; rebuild from source, ~3 min)
    cd /home/user/integer_solver/solve_lab/agentJ_work
    python3 jparse2.py && python3 jpoly.py && python3 jlead_build.py && python3 jfit.py
    python3 jvalidate.py    # MUST print 0 mismatches before trusting anything
    python3 jcode.py        # support-only lower bounds
    python3 jgen3.py        # union generators + exhaustive subset search
    python3 jsolve2.py ; python3 jdiag.py 1 1 ; python3 jkernel.py 1 1 ; python3 jbits.py 1 1

## Single highest-value next experiment
Enlarge L honestly rather than enlarging the search over a fixed L.  Both move types
are single-variable; the untested class is COMPENSATED PAIRS — move u and v together so
their effects on every atom outside T cancel while their effect on the 7 cluster values
does not.  That is a kernel computation on the Jacobian of the outside-T atoms
restricted to the cluster's variable cone, and it is small.  If the enlarged lattice
still fails every six-subset, 39026 is optimal for this cluster in a strong and
checkable sense, and the remaining hope is a DIFFERENT tight cluster: scan all atoms
for support sets with |R| small and few alone-rows (jcode.py already computes that
statistic for any candidate T).
