"""S12 step 8: TARGETED closure rebuild.

bulk.py activated at random and rows/cols grew in lockstep.  Here we rebuild the
closure after the activations that the systematic sweep ranked best -- cheapest
per knob, and (the key criterion) free inputs that reach the cluster but few
other checks, so columns are added without adding rows.
Reports  rows x cols, rank, kernel, consistency defect, and the minimum-equation
-cost violated set.
"""
import os, sys, json, time, collections, random
HERE='/home/user/integer_solver/solve_lab/s10'; LAB='/home/user/integer_solver/solve_lab'
sys.path.insert(0, os.path.join(LAB,'s9','eff')); sys.path.insert(0, HERE)
import lib as L, tools as T, ad, ac_lib as A
from ac_sacr import closure_of, decode
P = ad.P
BAD = [21617, 29539]
B = A.Base(os.path.join(HERE,'mod9118_0.json'))
supp0 = A.grad_supp(B.v0, BAD)

def report(tag, act):
    v = list(B.v0)
    for u, val in act.items(): v[u] = val
    if act: A.fwd_local(v, list(act))
    changed = {w for w in range(L.NVARS) if w < L.NVARS and v[w] != B.v0[w]}
    sc, newnz, newchk, lost, gained, av, nz = B.cost(v, changed)
    supp = A.grad_supp(v, BAD)
    t0 = time.time()
    rows, Us, cols = closure_of(v, BAD)
    w = {c: len(L.atom2eq.get(c, {})) for c in rows}
    best = None
    orders = {'nat': rows, 'cost-desc': sorted(rows, key=lambda c: (-w[c], c))}
    random.seed(3)
    for k in range(3):
        o = list(rows); random.shuffle(o); orders[f'r{k}'] = o
    rank = nbad = None
    for name, o in orders.items():
        x, D, rank, nbad = decode(rows, Us, cols, av, o)
        eqs = L.eqs_of_atoms(D)
        if best is None or len(eqs) < best[0]: best = (len(eqs), name, D)
    print(f'{tag}: act={ {k:(v_ if abs(v_)<1000 else "big") for k,v_ in act.items()} }\n'
          f'   integer score {sc}  knobs +{len(supp-supp0)}  checks broken {len(newchk)}  '
          f'eqs lost {len(lost)}\n'
          f'   closure {len(rows)} x {len(Us)}  rank {rank}  KERNEL {len(Us)-rank}  '
          f'inconsistent {nbad}\n'
          f'   min-cost decode: {best[0]} equations (order {best[1]}), '
          f'|D|={len(best[2])} -> mod-p ceiling {L.NEQ-best[0]}  ({time.time()-t0:.0f}s)',
          flush=True)
    return dict(tag=tag, act={str(k): str(v_) for k, v_ in act.items()}, score=sc,
                knobs=len(supp-supp0), rows=len(rows), cols=len(Us), rank=rank,
                kernel=len(Us)-rank, nbad=nbad, cost=best[0], D=best[2])

if __name__ == '__main__':
    cands = json.load(open(os.path.join(HERE, sys.argv[1])))
    out = []
    for c in cands:
        out.append(report(c['tag'], {int(k): int(v_) for k, v_ in c['act'].items()}))
        json.dump(out, open(os.path.join(HERE, sys.argv[2]), 'w'))
