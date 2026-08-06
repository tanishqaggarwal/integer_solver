"""S12 step 11: the same systematic activation sweep, but in the frame that
actually holds the record (39,026).  fwd_local makes a FULL sweep of all 7,273
free inputs affordable, so this is exhaustive at first order.
"""
import os, sys, json, time, collections
HERE='/home/user/integer_solver/solve_lab/s10'; LAB='/home/user/integer_solver/solve_lab'
sys.path.insert(0, os.path.join(LAB,'s9','eff')); sys.path.insert(0, HERE)
import lib as L, tools as T, ad, ac_lib as A
P = ad.P; FREE = sorted(ad.FREE); FORBID = {2081, 4287}
B = A.Base(os.path.join(LAB,'best','new_instance_partial_39026.json'))
print(f'base score {B.score0}  nonzero atoms {sorted(B.nz0)}  failing {sorted(B.fail0)}', flush=True)
BADC = sorted(a for a in B.nz0 if a in A.CHECKSET)
supp0 = A.grad_supp(B.v0, BADC)
print(f'nonzero CHECK atoms {BADC}; their gradient support {len(supp0)} free inputs', flush=True)
t0 = time.time(); res = []
VALS = (1, -1)
for i, z in enumerate(FREE):
    if z in FORBID: continue
    for val in VALS:
        if B.v0[z] == val: continue
        v = list(B.v0); v[z] = val
        A.fwd_local(v, [z])
        changed = {w for w in range(L.NVARS) if v[w] != B.v0[w]}
        if not changed: continue
        sc, newnz, newchk, lost, gained, av, nz = B.cost(v, changed)
        res.append((sc, z, val, len(newnz), len(gained), len(lost)))
    if i % 1000 == 0: print(f'  {i}/{len(FREE)} ({time.time()-t0:.0f}s) best {max(res)[0] if res else 0}', flush=True)
res.sort(reverse=True)
print(f'\nswept {len(res)} single moves in {time.time()-t0:.0f}s')
print('top (score, x, val, atoms_broken, eqs_gained, eqs_lost):')
for r in res[:30]: print(f'  {r}')
json.dump(res[:4000], open(os.path.join(HERE,'ac_sweep26.json'),'w'))
neutral = [r for r in res if r[0] >= B.score0]
print(f'\nmoves that do not lose ground (score >= {B.score0}): {len(neutral)}')
