# Massively-underdetermined MQ analysis of EQUATIONS.txt — regime determination

Best verified unchanged: **39022/39033** (`best/new_instance_partial_39022.json`, 11 fail, all fail mod p).
No full solution found: the massively-underdetermined regime does NOT hold for this instance.

## Formulation
All atoms are degree 1, 2, or 4; the 694 degree-4 atoms are perfect squares Q^2 (=0 <=> Q=0),
so the flat system (all 38748 vars, mul-gates as quadratic constraints) is a genuine MQ system.

## d and m (three formulations, all consistent)
FLAT (cleanest, everything degree<=2):
  d = 16641  free params after full linear solve (union-find + Gaussian, rank as computed)
  m = 16771  independent quadratic forms (exact GF(p) rank of substituted deg-2 atoms) [+694 sq]
  => m > d : the quadratic core is OVER-DETERMINED.

FREE-INPUT (mission's "#free-inputs - rank"):
  free inputs = 8583 ; only 386 atoms are affine (linear), rank 301
  d = 8583 - 301 = 8282  ;  m = 6213 nonlinear atoms (value*value core ~4351)
  BUT residuals in the free/knob inputs are degree 3 and 6 (cascaded products), NOT quadratic.

LOCAL at 39022:
  Jacobian rank 3036 = active free-input cols 3036  => continuous null space 0
  exactly 1 independent obstruction (the G1/G2 conservation law; "cert size 71" = atoms in it)
  the 11 failing eqs are INDEPENDENT of all 5547 inactive knobs.

## Regime verdict
KPG/Thomae-Wolf/MHT need d >= m(m+1) (or the reduction gain floor(d/m)-1 to be large).
Here in EVERY formulation:
  - flat: m(m+1) ~ 2.8e8  vs d=16641   (fails ~2e4x) and m>d anyway
  - free: m(m+1) ~ 3.9e7  vs d<=8282    (fails ~5e3x); d/m ~ 1.33 => Thomae-Wolf gain 0
The system is DETERMINED-to-OVER-DETERMINED (omega = d/m in [0.98, 1.33]), never the
omega ~ m (thousands) required for polynomial-time MQ. NOT massively underdetermined.

## Why the approach cannot solve it (mechanistic)
1. Obstruction unreachable: the 11 failing eqs depend only on the fully-determined active
   inputs (null space 0); they are independent of every underdetermined ("knob") direction.
2. Not quadratic: the knob residuals are degree 3 and 6, so no MQ change-of-variables applies.
3. Wire lock (mod p): atom 37110 (x_26064 - p) forces the 220-var wire = 0 mod p, so all
   wire-gated handles (incl. the G1/G2 partners x_17325, x_9413) are inert mod p.
Consistent with prior 4-agent conclusion: a hardened 256-bit-message trapdoor.

## Scripts (in scratchpad): mp.py, param.py, compute_m.py, rank_m.py, dfree.py, crux.py, obs22.py
