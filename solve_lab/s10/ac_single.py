"""S12 step 5: exhaustive SINGLE-activation profile of every zero free input in
the cluster cone (the complete candidate set: 189)."""
import os, sys, json, time, random
HERE='/home/user/integer_solver/solve_lab/s10'; LAB='/home/user/integer_solver/solve_lab'
sys.path.insert(0, os.path.join(LAB,'s9','eff')); sys.path.insert(0, HERE)
import lib as L, tools as T, ad, ac_lib as A
P = ad.P
definer = L.definer
FREE = set(ad.FREE); FORBID = {2081, 4287}
random.seed(7)
B = A.Base(os.path.join(HERE,'mod9118_0.json'))
BAD = [21617, 29539]
supp0 = A.grad_supp(B.v0, BAD)
print(f'base score {B.score0}  cluster support {len(supp0)}', flush=True)

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
for a in BAD: seeds |= set(L.avars[a])
CC = cone(seeds)
pool = sorted(u for u in CC if u in FREE and B.v0[u] == 0 and u not in FORBID)
print(f'pool {len(pool)}   (of which in supp0: {len(set(pool)&supp0)})', flush=True)

t0 = time.time(); out = {}
VALS = [1, -1, 0x9e3779b97f4a7c15]
for i, z in enumerate(pool):
    rec = []
    for val in VALS:
        v = list(B.v0); v[z] = val
        A.fwd_local(v, [z])
        changed = {u for u in range(L.NVARS) if v[u] != B.v0[u]}
        sc, newnz, newchk, lost, gained, av, nz = B.cost(v, changed)
        supp = A.grad_supp(v, BAD)
        newvars = sorted(u for u in changed if B.v0[u] == 0 and v[u] != 0)
        rec.append({'val': val, 'score': sc, 'nchanged': len(changed),
                    'newnz': sorted(newnz), 'newchk': sorted(newchk),
                    'lost': len(lost), 'knobs': len(supp - supp0),
                    'supp': len(supp), 'newvars': newvars})
    out[str(z)] = rec
    if i % 25 == 0: print(f'  {i}/{len(pool)} ({time.time()-t0:.0f}s)', flush=True)
json.dump({'pool': pool, 'supp0': sorted(supp0), 'res': out},
          open(os.path.join(HERE,'ac_single.json'),'w'))
print(f'done ({time.time()-t0:.0f}s)')
best = []
for z, rec in out.items():
    for r in rec:
        if r['knobs'] > 0: best.append((r['knobs'], -r['lost'], int(z), r['val'], r['score'], len(r['newchk'])))
best.sort(reverse=True)
print(f'\nSINGLE activations that grew the cluster support: {len(best)}')
for b in best[:20]: print('  knobs=%d lost_eqs=%d  x_%d=%s  score %d  checks broken %d' % (b[0], -b[1], b[2], b[3], b[4], b[5]))
free_ones = sorted(((r['lost'], len(r['newchk']), len(r['newvars']), int(z), r['val'])
                    for z, rec in out.items() for r in rec), key=lambda t: (t[0], t[1]))
print('\ncheapest activations (lost_eqs, checks_broken, new_nonzero_vars, x, val):')
for f in free_ones[:20]: print(f'  {f}')
