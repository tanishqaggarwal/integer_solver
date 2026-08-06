"""S12 step 20: STACKED frozen activation + repair.

The frozen repair buys +2 knobs for 6 equations (39009 -> 39003) instead of the
raw 19.  Stack several activations, freeze them all, repair, and measure the real
exchange rate; then rebuild the closure on the repaired state and re-run the
minimum-equation-cost decode.
"""
import os, sys, json, time, random
HERE='/home/user/integer_solver/solve_lab/s10'; LAB='/home/user/integer_solver/solve_lab'
sys.path.insert(0, os.path.join(LAB,'s9','eff')); sys.path.insert(0, HERE)
import lib as L, tools as T, ad, ac_lib as A
import engine as E
from ac_sacr import closure_of, decode
P = ad.P
B = A.Base(os.path.join(HERE,'mod9118_0.json'))
BAD = [21617, 29539]
supp0 = A.grad_supp(B.v0, BAD)
TOP = [24365, 12054, 16586, 17406, 11368, 23751, 22562]
random.seed(77)
res = []
for k in (1, 2, 3, 4, 7):
    sel = TOP[:k]
    v = list(B.v0)
    for z in sel: v[z] = 1
    A.fwd_local(v, sel)
    ch = {w for w in range(L.NVARS) if v[w] != B.v0[w]}
    sc0, newnz, newchk, l0, g0, av, nz = B.cost(v, ch)
    E.FORBID = {2081, 4287} | set(sel)
    v2, cur = E.run(v, f'stack{k}', iters=40, budget=300)
    alive = [z for z in sel if v2[z] != 0]
    k2 = len(A.grad_supp(v2, BAD) - supp0)
    av2 = L.all_atom_values(v2)
    nz2 = sorted(a for a in range(L.NA) if av2[a])
    print(f'STACK {k}: raw {sc0} -> repaired {cur[0]}  knobs +{k2}  alive {len(alive)}/{k}  '
          f'nonzero atoms {nz2}', flush=True)
    if cur[0] >= 39000 and k2 > 0:
        t0 = time.time()
        rows, Us, cols = closure_of(v2, [a for a in nz2 if a in A.CHECKSET] or BAD)
        w = {c: len(L.atom2eq.get(c, {})) for c in rows}
        best = None
        for name, o in (('cost-desc', sorted(rows, key=lambda c: (-w[c], c))), ('nat', rows)):
            x, Dv, rank, nbad = decode(rows, Us, cols, av2, o)
            e = len(L.eqs_of_atoms(Dv))
            if best is None or e < best[0]: best = (e, Dv, rank, len(Us))
        print(f'   closure {len(rows)} x {len(Us)}  rank {best[2]}  KERNEL {best[3]-best[2]}  '
              f'min-cost decode {best[0]} equations -> mod-p ceiling {L.NEQ-best[0]}  '
              f'({time.time()-t0:.0f}s)', flush=True)
    res.append((k, sc0, cur[0], k2, len(alive)))
    if cur[0] > 39026:
        T.save(v2, os.path.join(HERE,'ac_best.json')); print('*** BEAT 39026 ***')
print('\nk  raw  repaired  knobs  alive')
for r in res: print(' ', r)
