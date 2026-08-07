# RESUME_J — agent J (reduced-parameterization attack)

## Best verified
Baseline `solve_lab/best/new_instance_partial_39026.json` = **39026/39033**, re-verified
by me with `solve_lab/checker.py` (fails [12231,12270,12350,14584,18673,22044,29125]).
My own best states are worse: 39004 (on-manifold, mod-p clean) and
`agentJ_work/J_b10_38998.json` = 38998.  **No improvement on the deliverable.**

## Independence
Whole pipeline rebuilt from `EQUATIONS.txt`; no prior-session or other-agent artifact
was read into any computation.  Gate on everything: `jvalidate.py` reconstructs all
39,033 raw equations from my model and compares at a random point — **0 mismatches**.

## VERDICT 1 — the reduced parameterization is REAL (re-derived independently)
 * eq_i = mult_i*(sum_j c_ij A_j)^k, k in {1,2,4}, mult_i != 0 always; 39033 distinct
   atoms, all of degree <= 2.  The instance IS exactly "M a = 0" over the atoms.
 * p = 115792089237316195423570985008687907853269984665640564039457584007908834671663
   (secp256k1 prime) is a literal pin; dozens of variables are copies of it.
 * Definer DAG acyclic: 30575 defined vars, 8173 free, 8458 constraint atoms.
   On-manifold from the deliverable's free inputs only **4** constraints are violated.
 * Free inputs: 3747 are occ==1 handles, only 21 nonzero, only 17 non-handle:
   2 booleans (x2081, x24601), x8731/x9118 (freely zeroable — verified, costs nothing),
   and **13 numbers of 295-296 bits**.  Exactly **4 of the 13 are pinned to literal
   constants** (a3895, a3897, a32257, a32259), gated by x2081 / x24601; the other 9 form
   a chain that is **degree 1** in each knob (exact finite differences over GF(p)).
   One Gauss-Seidel sweep solves it, leaving only a20407/a20409/a31575 (deg 3,3,2).

## VERDICT 2 — the branch structure (new; nobody had worked here)
x15298 = x7715(x24601) * x34554(x2081).  Turning EITHER boolean off sets x15298 = 0 and
annihilates a20407/a20409/a31575 outright.  Mod p, after the chain sweep:
      (x2081,x24601) = (1,1) -> 3 violated [20407, 20409, 31575]   (EC residual)
      (1,0) or (0,1)         -> **2 violated [731, 31571]**        (the output pins)
      (0,0)                  -> 3 violated [731, 24075, 31571]
Mechanism: with x15298=1 the output pins consume the FREE x30213/x22162 (satisfiable)
but the EC residual is live; with x15298=0 they instead consume x14853 / x12186, which
the chain already pins.  That is the trade the instance is built around.

## VERDICT 3 — my "224 hidden movers" do NOT give a mod-p degree of freedom
I flagged 224 currently-zero free variables that move the residual, as a challenge to
the claim that only the selector bits are free.  The full Jacobian test settles it
against me.  (Tooling correction first: `+1` finite differences are wrong for the
degree-2/3 residual constraints, so `jnewton2.py` / `jdiag.py` compute the EXACT
Jacobian by forward-mode AD over GF(p).)  Knobs = every free variable in the backward
cone of the violated constraints (complete — anything outside has literally zero
effect).  Rows = every constraint atom with a nonzero entry (complete).  Nothing
hand-picked, unlike prior §135.

    branch (1,0): violated [731,31571]         J = 967 x 524, nnz 2222
    branch (0,1): violated [731,31571]         J = 966 x 524, nnz 2221
    branch (1,1): violated [20407,20409,31575] J = 988 x 523, nnz 2683
    all three: rank(J_satisfied) = 426; INCONSISTENT; and every violated
    constraint's gradient lies ENTIRELY in the span of the satisfied rows with
    nonzero reduced rhs => exactly zero gradient on the whole kernel.

Not merely tangential (`jkernel.py`): explicit basis of ker(J_sat) — 98 dims in (1,0),
97 in (1,1) — with LARGE random GF(p) combinations applied and re-propagated exactly.
Every trial left the violated residues **numerically identical** with nothing new
broken; 0/98 and 0/97 single kernel directions moved a residue at all.  The continuous
freedom is **exact gauge**, not a first-order artefact.

Boolean knobs (`jbits.py`), the one thing a tangent test is structurally blind to:
256 of them, **0 score-neutral**, **220 do move the residues** (so prior §144-145's
"inert" reading is wrong), but in branch (1,1) every one breaks exactly 6 constraints
and heals none except x2081/x24601, which heal all three EC constraints at a cost of
10/9 — i.e. they are just the branch switches, relocating the residual to the output
pins.  In branch (1,0), **no bit heals anything at all**.

**Conclusion: corroborates agent I's step 6 from an independent parse, measured over
the full constraint set.**  What is NOT established: infeasibility.  The 254 data bits
remain a combinatorial search, and that search is the ECDLP.

## Re-enter (caches are gitignored; rebuild from source, ~3 min)
    cd /home/user/integer_solver/solve_lab/agentJ_work
    python3 jparse2.py && python3 jpoly.py && python3 jlead_build.py && python3 jfit.py
    python3 jvalidate.py    # MUST print 0 mismatches before trusting anything
    python3 jsolve2.py      # branch table
    python3 jdiag.py 1 1    # exact Jacobian + frozen-residual diagnostic (caches jjac_11.pkl)
    python3 jkernel.py 1 1  # large exact kernel moves
    python3 jbits.py 1 1    # boolean knob census

## Single highest-value next experiment
Stop looking for continuous freedom — it is provably gauge in all three branches.
The only untested object is the 254-bit combinatorial space, and single-bit costs are
uniform (5-6), which is the signature of a ladder rather than a searchable landscape.
The one experiment that could still change the picture cheaply: take the deliverable's
39,026 state and ask whether a SECOND off-manifold coding exists — the deliverable
achieves 7 failing equations from 7 nonzero atoms via cancellation, while my cleanest
on-manifold state has only 2 nonzero atoms but 29 failing equations.  Minimising
failing equations over the 2-atom residual (choose the handle lifts so the two atom
values land in the kernel of as many equation rows as possible) is a small exact
integer-linear problem and is the only route to >39,026 that does not require the DL.
