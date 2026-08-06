# Channel (1,0) is rigorously ℤ-infeasible (conserved obstruction Q)

Established by a full global mod-p feasibility analysis over all 8583 free-input
columns (mod-p forward-tangent Jacobian J_sat of the 39022 satisfied equations).

## Result
The 11 fails ⟺ {G1=0, G2=0} (10 are integer combos a·G1+b·G2; **eq 8680 = G2³
exactly**). Mod-p Rouché–Capelli test of closing both gaps while keeping J_sat=0:

```
rank(J_sat)                    = 3035   (nullity 5548)
rank([J_sat; gradG1])          = 3036   gradG1 INDEPENDENT  -> closing G1 alone: CONSISTENT
rank([J_sat; gradG1; gradG2])  = 3036   gradG2 IN rowspace  -> G2 locked to G1
rank([coeff | rhs])            = 3037   -> INCONSISTENT (both gaps together)
```
Robust across reversed gap order and random row shuffles.

## The certificate
On null(J_sat), any move changing G1 by δ changes G2 by exactly r1·δ (locked 1-dim
coupling, `r1·r2 ≡ 1 mod p`). Closing both needs −G2 = r1·(−G1), i.e. the conserved
quantity
```
Q = (G2 − r1·G1) mod p = 32805349568647205012024350336202625691777084639260590353243920465632505312816  ≠ 0
```
**Beyond first order:** an actual finite mod-p move that closes G1 (verified by full
circuit re-eval) leaves G2 landing on exactly Q, still unclosable. A direct encoding
(fix x_7068→M1, x_4432→M2 as mod-p targets) gives the identical verdict. Consistent
with the prior exact-ℤ z3 UNSAT over all local selector combos.

## Consequence
No perturbation/CRT-lift of best_agentA (channel (1,0)) can reach 39033. Since every
absorber shift is a multiple of p, mod p is the only binding modulus and it is locked.
The other *local* MUX channels are worse and independently UNSAT. Reaching 39033 there-
fore requires a **non-local message/selector reconfiguration** (a globally different
assignment), not a repair of agentA. Channel (1,1) has a *different* obstruction
(branch-B gadget on x_21279 = x_9062·x_20434) under separate analysis.

Scripts: modp_ranks.py (ranks), modp_extract_Q.py (r1, Q), modp_finite_move.py
(beyond-first-order confirmation).

## Independent confirmation (second method) + the mechanism
A separate deflation / second-order-Newton agent reached the identical verdict by a
different route, and pinned the mechanism:
- Rank-defect **exactly 1 at every scale**: full system (3036 vs 3037), 11-fail
  window (9 vs 10), 27-eq window (19 vs 20), 785-eq closure (129 vs 130). Robust; it
  *relocates* under moves but never vanishes.
- Obstruction covector c supported on 9 eqs [2043,2554,6494,8124,8680,8687,9421,12231,
  12270]; c·J=0 over all 39033 rows, c·R=1. The obstruction IS movable by off-manifold
  moves (30/30 distinct values) but c·H[v,v]=0 for v∈null(J) (why first-order calls it
  "conserved").
- **Why every local knob is impotent:** null(J) is spanned by the **5547 "dead" free
  inputs**, which are mod-p inert — they feed products against the hardcoded constant
  p-wire x_26064 (220-member identity class all = p), so wire·handle ≡ 0 mod p. Moving
  any of them leaves all 39022 satisfied, the 11 still failing, c·R still = 1. The 3036
  active inputs are mutually independent (no constraint-preserving move among them). So
  the obstruction can only be moved through active/range directions, which are pinned by
  the other equations — killing it needs a mod-p "teleport" to a random ~p-sized value,
  which lands far from any solution (predictor-corrector/window solvers diverge
  11→124→228). The wire-escape used to solve the *previous* instance is structurally
  closed here (constant-p literal, not a free wire).

CONVERGENT CONCLUSION (mod-p feasibility agent + deflation agent + direct ranks.py):
best_agentA is a rigorously-characterized local maximum at 39022/39033. Reaching 39033
requires a non-local reconfiguration (the setter's global 256-bit witness), not any local
repair. Channel (1,1) branch-B and a direct reduced-core exact solve are under test as the
remaining non-local routes.
