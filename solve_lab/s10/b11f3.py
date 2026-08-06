"""S11 step 34: frame 3 in branch (1,1), where x_21279 = 1 makes x_8731 a LIVE
driver of x_19964 -- exactly the knob the pair move needs, and one that costs
nothing in frames 2/3.
"""
import os, sys, time
HERE = os.path.dirname(os.path.abspath(__file__))
LAB = os.path.dirname(HERE)
sys.path.insert(0, os.path.join(LAB, 's9', 'eff'))
sys.path.insert(0, HERE)
import lib as L, tools as T, ad
from frame3 import DETACH, definer, ORDER, FREE, CHECKS, fwd, score, SSET
P = ad.P
FORBID = {2081, 4287}
base = L.load(os.path.join(LAB, 'best', 'new_instance_partial_39026.json'))
fwd(base)
v = list(base); v[2081] = 1; v[4287] = 1
fwd(v, rounds=10)
av = L.all_atom_values(v)
print(f'frame3 + branch(1,1): score {score(v)}  x_21279={v[21279]} x_7075={v[7075]}')
print(f'  nonzero {[a for a in range(L.NA) if av[a]]}', flush=True)
d0 = v[19964]
w = list(v); w[8731] += 1
fwd(w, rounds=10)
print(f'  d(x_19964)/d(x_8731) = {w[19964]-d0}   (was 0 in branch (1,0))')

def pot(vv):
    a2 = L.all_atom_values(vv)
    nz = [a for a in range(L.NA) if a2[a]]
    return (L.NEQ - len(L.failing_eqs(a2)), -len(nz)), a2, nz

def moves(a, vv, a2, nn=20):
    out = []
    for x in sorted(set(L.avars[a])):
        if x in FORBID: continue
        tgt = T.solve_lin(a, x, vv)
        if tgt is None or tgt == vv[x]: continue
        if x in FREE: out.append((x, tgt))
        else:
            dd = definer.get(x)
            if dd is None: continue
            vx = list(vv); vx[x] = tgt
            for u in sorted(set(L.avars[dd])):
                if u == x or u not in FREE or u in FORBID: continue
                nv = T.solve_lin(dd, u, vx)
                if nv is not None: out.append((u, nv))
    r = a2[a] % P
    if r:
        vm = [x % P for x in vv]
        lam = {}
        # reverse AD in this frame
        import collections
        lamd = collections.defaultdict(int)
        for x in L.avars[a]: lamd[x] = (lamd[x] + ad.dpart(a, x, vm)) % P
        for t in reversed(ORDER):
            lt = lamd.get(t, 0)
            if not lt: continue
            dd = definer[t]
            dz = ad.dpart(dd, t, vm)
            if dz % P == 0: continue
            f = -lt * pow(dz, -1, P) % P
            for x in L.avars[dd]:
                if x == t: continue
                dx = ad.dpart(dd, x, vm)
                if dx: lamd[x] = (lamd[x] + f * dx) % P
            lamd[t] = 0
        g = {u: lamd[u] % P for u in FREE if lamd.get(u, 0) % P}
        for k, u, dz in sorted((len(L.var_atoms[u]), u, dz) for u, dz in g.items()
                               if u not in FORBID)[:nn]:
            out.append((u, vv[u] + (-r * pow(dz, -1, P)) % P))
    return out

cur, av, nz = pot(v)
print(f'\nengine start {cur[0]} nonzero {len(nz)}', flush=True)
t0 = time.time()
for it in range(50):
    if time.time() - t0 > 1800: break
    got = None
    for a in nz:
        for u, nv in moves(a, v, av):
            tr = list(v); tr[u] = nv
            fwd(tr, rounds=8)
            p2, a2, n2 = pot(tr)
            if p2 > cur: got = (a, u, p2, tr, a2, n2); break
        if got: break
    if not got:
        print(f'  it{it}: stuck {cur[0]}  nonzero {nz}'); break
    a, u, p2, tr, a2, n2 = got
    print(f'  it{it}: a{a} via x_{u}  {cur[0]} -> {p2[0]}  nz {len(nz)}->{len(n2)}',
          flush=True)
    v, cur, av, nz = tr, p2, a2, n2
    if p2[0] > 39026:
        T.save(v, os.path.join(HERE, f'B11F3_{p2[0]}.json'))
        print(f'   *** BEATS THE DELIVERABLE {p2[0]}', flush=True)
print(f'FINAL {cur[0]}  nonzero {nz}')
