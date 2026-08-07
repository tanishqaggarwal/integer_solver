"""S12 step 15: the systematic activation study in the RECORD frame (39,026).

Cluster = the 7 nonzero atoms {22229,22230,35758..35762} touching 12 equations.
Enumerate its structural cone, every zero free input in it, every pair, and
measure knobs gained by the cluster's gradient support vs equations lost.
"""
import os, sys, json, time, itertools, collections, random
HERE='/home/user/integer_solver/solve_lab/s10'; LAB='/home/user/integer_solver/solve_lab'
sys.path.insert(0, os.path.join(LAB,'s9','eff')); sys.path.insert(0, HERE)
import lib as L, tools as T, ad, ac_lib as A
P = ad.P; definer = L.definer
FREE = set(ad.FREE); FORBID = {2081, 4287}
B = A.Base(os.path.join(LAB,'best','new_instance_partial_39026.json'))
BADA = sorted(B.nz0)
supp0 = A.grad_supp(B.v0, BADA)
suppC = A.grad_supp(B.v0, [a for a in BADA if a in A.CHECKSET])
print(f'base {B.score0}; cluster atoms {BADA}')
print(f'gradient support of ALL 7 nonzero atoms: {len(supp0)} free inputs; of the 2 checks: {len(suppC)}', flush=True)

def cone(seeds):
    c, st = set(), list(seeds)
    while st:
        t = st.pop()
        if t in c: continue
        c.add(t)
        a = definer.get(t)
        if a is None: continue
        for w in L.avars[a]:
            if w != t: st.append(w)
    return c
seeds = set()
for a in BADA: seeds |= set(L.avars[a])
CC = cone(seeds)
pool = sorted(u for u in CC if u in FREE and B.v0[u] == 0 and u not in FORBID)
print(f'cluster cone {len(CC)} vars; zero free inputs in it: {len(pool)} '
      f'(of which already in supp0: {len(set(pool)&supp0)})', flush=True)

t0 = time.time(); sing = {}
for z in pool:
    rec = []
    for val in (1, -1, 0x9e3779b97f4a7c15):
        v = list(B.v0); v[z] = val
        A.fwd_local(v, [z])
        changed = {w for w in range(L.NVARS) if v[w] != B.v0[w]}
        sc, newnz, newchk, lost, gained, av, nz = B.cost(v, changed)
        rec.append((len(A.grad_supp(v, BADA) - supp0), len(lost), sc, val))
    sing[z] = rec
print(f'singles done ({time.time()-t0:.0f}s)')
gk = sorted(((r[0], -r[1], z, r[3], r[2]) for z, rs in sing.items() for r in rs), reverse=True)
print('best singles (knobs, -eqs_lost, x, val, score):')
for g in gk[:15]: print(f'  {g}')
free_sing = [z for z in pool if sing[z][0][1] == 0]
print(f'zero-cost singles in the cone: {len(free_sing)}')

t0 = time.time(); res = []; mx = 0
pairs = list(itertools.combinations(pool, 2))
for i, (u, z) in enumerate(pairs):
    v = list(B.v0); v[u] = 1; v[z] = 1
    A.fwd_local(v, [u, z])
    k = len(A.grad_supp(v, BADA) - supp0)
    mx = max(mx, k)
    m = max(sing[u][0][0], sing[z][0][0])
    if k == 0: continue
    changed = {w for w in range(L.NVARS) if v[w] != B.v0[w]}
    sc, newnz, newchk, lost, gained, av, nz = B.cost(v, changed)
    res.append((k, k - m, len(lost), u, z, sc))
    if i % 5000 == 0: print(f'  {i}/{len(pairs)} hits {len(res)} max {mx} ({time.time()-t0:.0f}s)', flush=True)
print(f'\nswept {len(pairs)} pairs in {time.time()-t0:.0f}s; {len(res)} grew the support; max knobs {mx}')
syn = [r for r in res if r[1] > 0]
print(f'genuinely second-order pairs: {len(syn)}')
res.sort(key=lambda r: (-r[0], r[2]))
print('top pairs (knobs, extra_over_singles, eqs_lost, x_u, x_z, score):')
for r in res[:20]: print(f'  {r}')
freeg = [r for r in res if r[2] == 0]
print(f'\nPAIRS GROWING THE SUPPORT AT ZERO EQUATION COST: {len(freeg)}')
for r in freeg[:20]: print(f'  {r}')
json.dump({'pool': pool, 'sing': {str(k): v for k, v in sing.items()}, 'res': res},
          open(os.path.join(HERE,'ac_act26.json'),'w'))
