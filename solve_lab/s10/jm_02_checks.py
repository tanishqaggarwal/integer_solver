"""jm step 2: anatomy of every CHECK atom that a pin move breaks, and the
free inputs able to repair each of them."""
import os, sys, json
sys.path.insert(0, '/home/user/integer_solver/solve_lab/s10')
import jm_lib as J
import lib as L, tools as T, ad
P = J.P
w = J.base_state()

BROKEN = [3576, 29539, 40826, 7930, 41512, 21617, 37662, 3578, 19297, 19299,
          30984, 36185, 40812, 3568, 3570, 19088, 22233, 22235, 29465, 36602,
          37887, 688, 1618, 3580, 3582, 25676, 29466, 30976, 30978, 33796,
          40608, 3584, 3586, 7932, 7934, 7936, 29467, 22231]

for a in BROKEN:
    if a not in range(L.NA):
        continue
    neq = len(L.atom2eq.get(a, {}))
    isg = a in J.atom_out
    print(f'\na{a}  neq={neq} gate={isg}  vars={sorted(L.avars[a])}')
    print(f'   src: {L.atom_src[a][:230]}')
    for u in sorted(L.avars[a]):
        d = J.definer.get(u)
        c = J.cone(u)
        fr = sorted(x for x in c if x in J.FREESET)
        print(f'     x_{u:<6} free={d is None} val={str(w[u])[:26]:<28} '
              f'cone={len(c):<5} freeinputs={fr[:12]}{"..." if len(fr) > 12 else ""}')
