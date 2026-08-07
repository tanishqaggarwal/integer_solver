"""S10 step 114: boolean scan INSIDE the cluster cones, with re-solving.

Part X scanned all 1,156 booleans in the WITNESS FRAME -- flip applied, nothing
re-solved -- and reported >= 7 for every one.  Every boolean in the cluster's
ancestor cones is a discrete freedom no linearisation can see, so rescan them
with forward re-evaluation plus a short enriched repair.
"""
import os, sys, time, collections
HERE = os.path.dirname(os.path.abspath(__file__))
LAB = os.path.dirname(HERE)
sys.path.insert(0, os.path.join(LAB, 's9', 'eff'))
sys.path.insert(0, HERE)
import lib as L, tools as T, ad
from engine import pot, moves
P = ad.P
definer, atom_out = L.definer, L.atom_out
FREE = set(ad.FREE)
BOOL = set()
for _a, _p in enumerate(L.polys):
    ks = list(_p.items())
    if len(ks) == 2:
        sq = [m for m, c in ks if len(m) == 2 and m[0] == m[1]]
        li = [m for m, c in ks if len(m) == 1]
        if sq and li and sq[0][0] == li[0][0]: BOOL.add(li[0][0])
v0 = L.load(os.path.join(HERE, 'mod9118_0.json'))
cone = set()
for tgt in (27522, 1308, 6858, 25442):
    st = [tgt]
    while st:
        t = st.pop()
        if t in cone: continue
        cone.add(t)
        a = definer.get(t)
        if a is None: continue
        for w in L.avars[a]:
            if w != t: st.append(w)
cands = sorted((cone & BOOL & FREE))
print(f'boolean free inputs inside the cluster cones: {len(cands)}', flush=True)
base, _, nz0 = pot(v0)
print(f'base {base[0]} nonzero {nz0}', flush=True)

def quick(v, iters=6):
    cur, av, nz = pot(v)
    for _ in range(iters):
        got = None
        for a in nz:
            for u, nv in moves(a, v, av, nnewton=8):
                tr = list(v); tr[u] = nv
                ad.fwd(tr, rounds=5)
                p2, av2, nz2 = pot(tr)
                if p2 > cur: got = (p2, tr, av2, nz2); break
            if got: break
        if not got: break
        cur, v, av, nz = got[0], got[1], got[2], got[3]
    return cur, v, nz

t0 = time.time()
best = (base, None, None)
res = []
for i, u in enumerate(cands):
    v = list(v0); v[u] = 1 - v[u]
    ad.fwd(v, rounds=6)
    cur, v, nz = quick(v)
    res.append((cur[0], u, len(nz)))
    if cur > best[0]:
        best = (cur, u, v)
        T.save(v, os.path.join(HERE, f'bool_{u}_{cur[0]}.json'))
        print(f'  *** x_{u}: score {cur[0]} nonzero {nz}', flush=True)
    if i % 10 == 0:
        print(f'  {i}/{len(cands)} ({time.time()-t0:.0f}s) best {best[0][0]}', flush=True)
res.sort(reverse=True)
print(f'\ntop 12 boolean flips: {res[:12]}')
print(f'BEST {best[0][0]} via x_{best[1]}')
