"""bl_branch: the four MUX branches, evaluated in the WITNESS frame (frame2) and
canonically, with the free p-wires that the seven checks need.

Branch (1,1) makes x_2099 = x_9118 (a FREE 2405-bit p-wire) instead of the pinned
constant K1 -- that is the only branch in which the congruence C0 = x_7068-x_2099
is not pinned by a constant load.
"""
import os, sys, json, time
HERE='/home/user/integer_solver/solve_lab/s10'; LAB='/home/user/integer_solver/solve_lab'
sys.path.insert(0, os.path.join(LAB,'s9','eff')); sys.path.insert(0, HERE)
import lib as L, tools as T, ad
from bl_core import BOOL, CANON, F2, pot, engine, FORBID
P = 2**256-2**32-977
SEVEN=[22229,22230,35758,35759,35760,35761,35762]

w26 = L.load(os.path.join(LAB,'best','new_instance_partial_39026.json')); F2.fwd(w26)
v0  = L.load(os.path.join(HERE,'mod9118_0.json')); CANON.fwd(v0)

for name, F, base in (('F2/w39026', F2, w26), ('CANON/mod9118', CANON, v0)):
    print(f'\n===== {name} =====', flush=True)
    for b1 in (0, 1):
        for b2 in (0, 1):
            v = list(base); v[2081]=b1; v[4287]=b2
            F.fwd(v, rounds=8)
            p, av, nz = pot(v)
            print(f'  ({b1},{b2}) plain : score {p[0]:>6}  nz {len(nz):>3} {nz[:10]}  '
                  f'x_7075={v[7075]} x_2099bits={v[2099].bit_length()}', flush=True)
            # zero the two p-wires the (1,1) branch wants to be zero
            if (b1, b2) == (1, 1):
                v2 = list(v); v2[1329]=0; v2[10903]=0
                F.fwd(v2, rounds=8)
                p2, av2, nz2 = pot(v2)
                print(f'         +x_1329=x_10903=0: score {p2[0]}  nz {len(nz2)} {nz2[:10]}', flush=True)
            t0=time.time()
            cur, vv, nz3 = engine(v, F, iters=40, budget=200, verbose=True, tag=f'{name}({b1},{b2})')
            print(f'  ({b1},{b2}) engine: score {cur[0]}  nz {len(nz3)} {nz3[:10]}  ({time.time()-t0:.0f}s)', flush=True)
            if cur[0] > 39026:
                T.save(vv, os.path.join(HERE, f'bl_best.json'))
                print('  *** SAVED bl_best.json', flush=True)
