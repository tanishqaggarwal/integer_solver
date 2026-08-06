"""S12 step 6: EXHAUSTIVE pair sweep over the complete candidate set.

All 189 zero free inputs in the cluster cone, all 17,766 unordered pairs, at
value 1 (and value -1 for the pairs that grow the support).  Metrics per pair:
knobs gained by the cluster's gradient support, atoms broken, check atoms broken,
equations lost.  Ranked by knobs per equation lost.
"""
import os, sys, json, time, itertools, collections
HERE='/home/user/integer_solver/solve_lab/s10'; LAB='/home/user/integer_solver/solve_lab'
sys.path.insert(0, os.path.join(LAB,'s9','eff')); sys.path.insert(0, HERE)
import lib as L, tools as T, ad, ac_lib as A
P = ad.P
definer = L.definer; FREE = set(ad.FREE); FORBID = {2081, 4287}
B = A.Base(os.path.join(HERE,'mod9118_0.json'))
BAD = [21617, 29539]
supp0 = A.grad_supp(B.v0, BAD)
D = json.load(open(os.path.join(HERE,'ac_single.json')))
pool = D['pool']
single = {int(z): r for z, r in D['res'].items()}
zerocost = [z for z in pool if single[z][0]['lost'] == 0]
print(f'pool {len(pool)}  zero-cost singles {len(zerocost)}  '
      f'singles with knobs>0 {sum(1 for z in pool if single[z][0]["knobs"]>0)}', flush=True)

t0 = time.time(); res = []
pairs = list(itertools.combinations(pool, 2))
for i, (u, z) in enumerate(pairs):
    v = list(B.v0); v[u] = 1; v[z] = 1
    A.fwd_local(v, [u, z])
    changed = {w for w in range(L.NVARS) if v[w] != B.v0[w]}
    supp = A.grad_supp(v, BAD)
    k = len(supp - supp0)
    if k == 0: continue
    sc, newnz, newchk, lost, gained, av, nz = B.cost(v, changed)
    res.append((k, len(lost), len(newchk), len(newnz), u, z, sc))
    if i % 2000 == 0: print(f'  {i}/{len(pairs)}  hits {len(res)} ({time.time()-t0:.0f}s)', flush=True)
print(f'swept {len(pairs)} pairs in {time.time()-t0:.0f}s; {len(res)} grew the support')
mx = max(r[0] for r in res) if res else 0
print(f'max knobs gained by any pair: {mx}')
hist = collections.Counter(r[0] for r in res)
print(f'knob histogram: {dict(sorted(hist.items()))}')
res.sort(key=lambda r: (-r[0], r[1]))
print('\nTOP by knobs then cheapest (knobs, eqs_lost, checks_broken, atoms_broken, x_u, x_z, score):')
for r in res[:25]: print(f'  {r}')
byeff = sorted(res, key=lambda r: (-(r[0]/(r[1]+1)), r[1]))
print('\nTOP by knobs-per-equation-lost:')
for r in byeff[:25]: print(f'  knobs={r[0]} lost={r[1]} chk={r[2]} atoms={r[3]} x_{r[4]},x_{r[5]} score={r[6]}')
free_growth = [r for r in res if r[1] == 0]
print(f'\nPAIRS THAT GROW THE SUPPORT AT ZERO EQUATION COST: {len(free_growth)}')
for r in free_growth[:30]: print(f'  {r}')
json.dump({'res': res}, open(os.path.join(HERE,'ac_pairs.json'),'w'))
