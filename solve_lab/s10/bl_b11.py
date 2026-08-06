"""bl_b11: branch (1,1) done properly.

With x_7075 = 0 the two surviving p-quantisation CHECKS in the canonical frame,
   a35759: 5113045*x_7075*x_9118 - x_29854 = -x_29854 = -P*x_1329
   a35760: x_31864 - P*x_10903,  x_31864 = -x_7075*x_8731 = 0
collapse to  x_1329 = 0  and  x_10903 = 0, and x_1329 / x_10903 occur in NO other
atom.  So the whole seven-atom cluster becomes free.  Also satisfy the two loads
that x_4287=1 switches on: a3568 (pins x_31861) and a3570 (pins x_14865).
"""
import os, sys, json, time
HERE='/home/user/integer_solver/solve_lab/s10'; LAB='/home/user/integer_solver/solve_lab'
sys.path.insert(0, os.path.join(LAB,'s9','eff')); sys.path.insert(0, HERE)
import lib as L, tools as T, ad
from bl_core import BOOL, CANON, F2, pot, engine, FORBID
P = 2**256-2**32-977

def rep(v, tag):
    p, av, nz = pot(v)
    print(f'{tag}: score {p[0]} nz {len(nz)} {nz[:24]}', flush=True)
    return p, av, nz

for src, F in (('mod9118_0.json', CANON), (os.path.join(LAB,'best','new_instance_partial_39026.json'), F2)):
    name = os.path.basename(src)
    v0 = L.load(src if os.path.isabs(src) else os.path.join(HERE, src)); F.fwd(v0, rounds=8)
    print(f'\n===== {name} / frame {"F2" if F is F2 else "CANON"} =====', flush=True)
    rep(v0, ' base')
    v = list(v0); v[2081] = 1; v[4287] = 1
    F.fwd(v, rounds=8); rep(v, ' (1,1)')
    for u in (1329, 10903, 29854, 31864):
        v[u] = 0
    F.fwd(v, rounds=8); p, av, nz = rep(v, ' (1,1)+zero p-wires')
    # satisfy the two newly-switched-on loads a3568 (x_31861) and a3570 (x_14865)
    for a, x in ((3568, 31861), (3570, 14865)):
        t = T.solve_lin(a, x, v)
        if t is not None: v[x] = t
    F.fwd(v, rounds=8); p, av, nz = rep(v, ' (1,1)+pwires+loads')
    # branch (1,1) switches ON three x_21279-gated quantisation checks:
    #   a19088 : x_21279*x_9106  = 13523997*x_9629 = 13523997*P*x_950
    #   a22233 : 6122989*x_21279*x_2239 = x_23754  = P*x_6947
    #   a22235 : x_21279*x_31731 = -x_35619        = -P*x_33168
    # zero the free multipliers so each collapses to "the sum wire must be 0"
    for u in (950, 6947, 33168):
        if u in F.FREE: v[u] = 0
    F.fwd(v, rounds=8); p, av, nz = rep(v, ' (1,1)+zero gated multipliers')
    for a in (19088, 22233, 22235):
        for u in sorted(L.avars[a]):
            t = T.solve_lin(a, u, v)
            if t is None or t == v[u]: continue
            tr = list(v); tr[u] = t; F.fwd(tr, rounds=8)
            p2, av2, nz2 = pot(tr)
            if p2 > p: v, p, av, nz = tr, p2, av2, nz2
            d = F.definer.get(u)
            if d is None: continue
            vv2 = list(v); vv2[u] = t
            for z in sorted(L.avars[d]):
                if z == u or z not in F.FREE: continue
                t2 = T.solve_lin(d, z, vv2)
                if t2 is None: continue
                tr = list(v); tr[z] = t2; F.fwd(tr, rounds=8)
                p2, av2, nz2 = pot(tr)
                if p2 > p: v, p, av, nz = tr, p2, av2, nz2
    rep(v, ' (1,1)+targeted a19088/a22233/a22235')
    T.save(v, os.path.join(HERE, f'bl_b11_{"f2" if F is F2 else "canon"}_seed.json'))
    t0 = time.time()
    cur, vv, nz2 = engine(v, F, iters=80, budget=600, verbose=True, tag='b11')
    print(f' ENGINE -> {cur[0]} nz {len(nz2)} {nz2}  ({time.time()-t0:.0f}s)', flush=True)
    T.save(vv, os.path.join(HERE, f'bl_b11_{"f2" if F is F2 else "canon"}_{cur[0]}.json'))
    if cur[0] > 39026:
        T.save(vv, os.path.join(HERE, 'bl_best.json')); print(' *** SAVED bl_best.json', flush=True)
