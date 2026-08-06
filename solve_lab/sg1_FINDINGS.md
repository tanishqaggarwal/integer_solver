# sg1 — Damped global solve on wire=p branch: findings

Best reached: **39022/39033 (11 fails)** = sg1_best.json (= agentA). Did NOT beat it.
core+G1/G2 are **NOT jointly satisfiable by local/damped methods** on the wire=p branch.

## Methods tried (all confirm the wall)
- Full-slave forward (sg1_fullslave.py, sg1_slavemap.json): 2074 checked free inputs identified;
  iterative slaving keeps every equality-check satisfied. KEY: under any unchecked-free perturbation
  the slave DIVERGES (integer-divisibility/carry breaks, baddiv grows 5->31) — the checks are
  p-adically rigid, not merely mod-p rigid.
- Damped coordinate descent / Gauss-Seidel (sg1_solve.py, sg1_heal.py): from 39013 -> plateau 20;
  from agentA -> plateau 11 (immediate). Healing the residue-fixed state (41 fails) only reverts
  x_14853/x_16742, snapping back to 20. Demonstrates the conservation law directly.
- Damped/LM Newton on core residual (sg1_newton.py): on the plain-forward manifold the mod-p
  Jacobian of (x_29322%p, x_3558%p) wrt all non-control unchecked frees is ZERO — the residue is
  a function of the free inputs x_14853,x_12186,x_16742 ONLY. Moving them zeroes the residue but
  breaks 23 equality-checks (13 from x_14853, ~10 from x_16742).

## The obstruction (precise, my characterization)
CONSERVATION LAW: the two sub-p residues r1,r2 must "live" somewhere:
- 39013 basin: r1,r2 sit in x_29322=x_14853-x_12186 and x_3558=x_24908-x_16742 -> core fails (20).
- agentA basin: core solved (r's absorbed) but r1,r2 reappear as G1=x_2099-x_7068 and
  G2=x_4432-x_19964 residues -> 11 fails. VERIFIED x_29322%p(39013) = p - r1(agentA).
- Joint satisfaction requires r1=r2=0. Both are pinned nonzero:
  r1 = 61705020361863629770768910187978745858728889529652486596432934143473517757811
  r2 = 33310166114805471624282140578459083391052142224394967852279417483154815501175

WHY pinned: the only unchecked handles inside G1,G2 are x_17325 (via x_642=x_28599*x_17325) and
x_9413 (via x_28730=x_17499*x_9413). Confirmed x_28599=x_17499=p (wire), so both products are
**p-quantized** (d(x_642)/d(x_17325)=p exactly). The gaps r1,r2 are sub-p -> no integer handle
absorbs them. Zeroing r1,r2 is a CVP/codeword problem on the 256-bit message = the setter's secret.

## Correction to prior docs
CORE_FINAL/RESUME say the residues are "linearly pinned". More precisely: x_29322%p IS reachable
by moving the free inputs x_14853/x_12186 directly, but that breaks 23 equality-checks whose
re-satisfaction (slaving) is integer-divisibility-INCONSISTENT (the carries don't close). The
rigidity is p-adic (carry-level), which is stronger/cleaner than "linear pinning".
