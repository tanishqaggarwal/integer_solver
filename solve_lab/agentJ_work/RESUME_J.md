# RESUME_J — agent J (reduced-parameterization attack)

## Best verified
Baseline `solve_lab/best/new_instance_partial_39026.json` = **39026/39033**, re-verified
myself with `solve_lab/checker.py` (fails [12231,12270,12350,14584,18673,22044,29125]).
My own best states are worse: 39004 (on-manifold, mod-p clean) and
`agentJ_work/J_b10_38998.json` = 38998. No improvement on the deliverable.

## VERDICT: the reduced parameterization is REAL, and I re-derived it from scratch
Independent pipeline, no prior-session code reused:
  jparse2.py  -> jmodel2.pkl   (paren-exact text parse of EQUATIONS.txt)
  jvalidate.py                 **0 / 39033 mismatches** vs the raw equations (mod 2^127-1)
  jpoly.py    -> jpoly.pkl     (exact monomial expansion of every atom)
  jlead_build.py -> jlead.pkl  (syntactic definer per atom)
  jengine.py / jman.py         definer DAG (acyclic; 30575 defined vars, 8173 free)
  jmodp.py                     GF(p) forward model + all 8458 constraint residues
  jchain.py / jsolve2.py       chain Gauss-Seidel; jlift.py integer lift; jmovers.py

Established facts (each measured, not assumed):
 * eq_i = mult_i*(sum_j c_ij A_j)^k, k in {1,2,4}, mult_i != 0 always; 39033 distinct
   atoms, every one of degree <= 2.  So the instance IS exactly "M a = 0" over atoms.
 * p = 115792089237316195423570985008687907853269984665640564039457584007908834671663
   (secp256k1 prime) is a literal pin; dozens of variables are copies of it.
 * 8458 constraint atoms (non-definer).  On-manifold from the deliverable's free
   inputs only **4** are violated: a8583, a30271, a35890, a35892.
 * Free inputs 8173; 3747 are occ==1 handles; only 21 nonzero; only 17 non-handle:
   2 booleans (x2081, x24601), x8731/x9118 (freely zeroable - verified, costs nothing),
   and **13 numbers of 295-296 bits**: x6418 x8778 x12553 x14623 x14853 x16742 x22152
   x22162 x22649 x24548 x30213 x31339 x33462.
 * Exactly **4 of the 13 are pinned to literal constants** (a3895, a3897, a32257,
   a32259), gated by x2081 / x24601.  The other 9 form a chain that is **degree 1**
   in each knob (verified by exact finite differences over GF(p)); one Gauss-Seidel
   sweep solves it and leaves only a20407/a20409/a31575 (degrees 3,3,2 - the
   point-addition residual).  Reproduces prior §129-131 by a different route.
 * NEW, sharpest statement of the obstruction:
   x15298 = x7715(x24601) * x34554(x2081).  Turning EITHER boolean off sets x15298=0,
   which annihilates a20407/a20409/a31575 outright.  Branch table (mod p, after sweep):
      (x2081,x24601) = (1,1) -> 3 violated [20407, 20409, 31575]  (EC residual)
      (1,0) or (0,1)         -> **2 violated [731, 31571]**       (the two OUTPUT pins)
      (0,0)                  -> 3 violated [731, 24075, 31571]
   Mechanism: with x15298=1 the output pins a731/a31571 consume the FREE x30213/x22162
   (so they are satisfiable) but the EC residual is live; with x15298=0 they instead
   consume x14853 / x12186, which the chain already pins.  That is the trade.

## Re-enter (caches are gitignored; rebuild from source, ~3 min)
    cd /home/user/integer_solver/solve_lab/agentJ_work
    python3 jparse2.py && python3 jpoly.py && python3 jlead_build.py && python3 jfit.py
    python3 jmodp.py        # base residues + per-parameter sensitivity table
    python3 jsolve2.py      # branch table
    python3 jlift.py 1 0    # integer lift + score, writes J_b10_<score>.json

## Single highest-value next experiment
`jnewton.py <b1> <b2>` (written, NOT yet run): full mod-p Newton — knobs = every free
variable in the cone of the violated constraints (~500), rows = every constraint any
knob moves, exact Jacobian by +1 probes, Gaussian elimination over GF(p).  Run it in
branch (1,0).  If the tangent system is CONSISTENT the two output pins are reachable
and the instance closes; if inconsistent, that inconsistency — measured on the FULL
constraint set rather than a hand-picked subsystem — is the first honest infeasibility
signal this lab has produced.
