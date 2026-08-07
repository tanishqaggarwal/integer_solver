# Agent G — RESUME (relaxation & rounding angle)

## Status
- Best verified: **39,026 / 39,033** = the inherited deliverable
  `solve_lab/best/new_instance_partial_39026.json` (re-verified by me with
  `python3 solve_lab/checker.py solve_lab/best/new_instance_partial_39026.json`
  -> `satisfied 39026/39033 (7 failing)`, failing `[12231,12270,12350,14584,18673,22044,29125]`).
- No improvement of my own yet.

## Confirmed so far
- 39,026 claim CONFIRMED (checker, 16 s).
- Coordinate probe (`g01_probe.py`): at the 39,026 deliverable x1==x2 and y1==y2
  (P1 == P2, a doubling, not a generic addition). At AG_39013/EC_39014/PF_39015 they differ.
  None of P1,P2,P3 lies on y^2=x^3+7 mod p at any state examined -> the
  "secp256k1 point addition" reading of Part XXVI is at best a formal identity,
  NOT points on the secp256k1 curve. (p IS the secp256k1 field prime.)

## Re-enter
```
cd /home/user/integer_solver/solve_lab/agentG_work
python3 g01_probe.py [state.json]     # coordinates / curve test on any state
```
Shared machinery: `solve_lab/s9/eff/lib.py` (L.load, L.all_atom_values, L.failing_eqs,
L.polys, L.atom_out, L.topo), `solve_lab/s10/ad.py` (ad.fwd forward eval),
`solve_lab/s10/suppfree.py` (free-input support bitsets).

## Next experiment (highest value)
Exact **symbolic** forward evaluation over F_p with the unpinned advice values as
indeterminates (gate output coefficients are all +-1, so forward eval is an honest
polynomial map). Goal: write the residual as a SMALL explicit polynomial system over
F_p in <= 9 unknowns, then solve it globally (resultants / Groebner / root finding
with flint) instead of locally. That is the "global relaxation" my angle calls for.
