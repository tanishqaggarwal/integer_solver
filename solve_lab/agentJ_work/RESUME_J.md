# RESUME_J (agent J) — verdict on the reduced parameterization

## Best verified
`solve_lab/best/new_instance_partial_39026.json` = 39026/39033 (re-verified myself).
No improvement of my own yet.  My on-manifold states: 39009 (see below), still < 39026.

## VERDICT: the claimed reduction IS REAL (independently re-derived, not copied)
Pipeline (all mine, all exact, `jvalidate.py` = 0/39033 mismatches):
 jparse2.py -> jmodel2.pkl   each eq = mult*(sum c_i A_i)^k, k in {1,2,4}, mult != 0,
                             39033 distinct atoms, all degree <= 2.
 jpoly.py   -> jpoly.pkl     monomial expansion of every atom
 jengine.py + jman.py        definer DAG (acyclic, 30575 defined vars, 8173 free),
                             forward propagation, exact scoring.
Facts established:
 * p = 115792089237316195423570985008687907853269984665640564039457584007908834671663
   (secp256k1 prime) appears as a pinned literal; many vars are copies of it.
 * At the 39026 deliverable exactly 7 atoms are nonzero (off-manifold trick).
 * ON-MANIFOLD (propagate the deliverable's free inputs through the definer DAG):
   score 38995, and only **4** constraint atoms are violated:
     a8583  9367949*(x24548-x25442) - x7927   ->  x24548 == x25442 (mod p)
     a30271 12846437*(x14853-x1308) - x29967  ->  x14853 == x1308  (mod p)
     a35890 5113045*x7075*x9118 - x29854      ->  p | x7075*x9118
     a35892 x7075*x8731 + x31864              ->  p | x7075*x8731
 * FREE INPUTS: 8173 total, 3747 are occ==1 handles, only **21 are nonzero**, and of
   those only **17 are non-handle**: 2 booleans (x2081, x24601), x8731/x9118 (huge),
   and **13 numbers of 295-296 bits**:
     x6418 x8778 x12553 x14623 x14853 x16742 x22152 x22162 x22649 x24548 x30213 x31339 x33462
   -> the "thirteen 296-bit numbers" claim is CONFIRMED.
 * x9118 := 0, x8731 := 0 (with handles x1329, x10903 := 0) kills a35890/a35892 with NO
   new violation: on-manifold 39004.  Then x24548 := x25442 kills a8583: **39009**,
   but breaks a22688 (11436039*(x14623-x27522) - x36864).

## Re-enter
    cd /home/user/integer_solver/solve_lab/agentJ_work
    python3 jparse2.py && python3 jpoly.py && python3 jfit.py     # rebuild caches
    python3 -c "import jman as J; J.run({9118:0,1329:0,8731:0,10903:0}, tag='CD')"

## Next experiment (highest value)
Build the mod-p reduced system: forward-propagate the whole circuit in GF(p) as a
function of the 13 parameters, extract the ~192 modular assertions as polynomials
mod p, and solve that system directly (it is small).
