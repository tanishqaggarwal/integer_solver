"""CL step 3c: price a mod-p Newton move on EVERY free input in the gradient support
of a21617 / a29539, with handle absorption."""
import os, sys, json, collections, time
HERE='/home/user/integer_solver/solve_lab/s10'; LAB='/home/user/integer_solver/solve_lab'
sys.path.insert(0, os.path.join(LAB,'s9','eff')); sys.path.insert(0, HERE)
import lib as L, tools as T, ad
P = 2**256 - 2**32 - 977
definer, atom_out = L.definer, L.atom_out
FREE = set(ad.FREE)

v0 = L.load(os.path.join(HERE,'mod9118_0.json'))
av0 = L.all_atom_values(v0)
NZ0 = set(a for a in range(L.NA) if av0[a])
BAD0 = set(L.failing_eqs(av0))
S0 = L.NEQ - len(BAD0)
print(f'base score {S0}  nonzero {sorted(NZ0)}')

TARGETS = [(21617, 21616, 5040), (29539, 29538, 30163)]

def absorb(v, target, gate, handle):
    w = atom_out[gate][1]
    tgt = T.solve_lin(target, w, v)
    if tgt is None: return False
    vv = list(v); vv[w] = tgt
    nv = T.solve_lin(gate, handle, vv)
    if nv is None: return False
    v[handle] = nv
    ad.fwd(v, rounds=8)
    return True

def price(target, gate, handle, us=None, tag=''):
    vm = [x % P for x in v0]
    g = ad.grad(target, vm)
    r = av0[target] % P
    rows = []
    US = us if us is not None else sorted(g)
    for u in US:
        d = g.get(u, 0) % P
        if d == 0: continue
        delta = (-r * pow(d,-1,P)) % P
        v = list(v0); v[u] = v[u] + delta
        ad.fwd(v, rounds=8)
        av = L.all_atom_values(v)
        ok = av[target] % P == 0
        if ok: absorb(v, target, gate, handle)
        av = L.all_atom_values(v)
        nz = set(a for a in range(L.NA) if av[a])
        bad = set(L.failing_eqs(av))
        s = L.NEQ - len(bad)
        rows.append((s, u, len(L.var_atoms[u]), ok, sorted(nz), len(bad-BAD0), len(BAD0-bad)))
    rows.sort(key=lambda t: -t[0])
    print(f'\n===== a{target} Newton pricing {tag} ({len(rows)} moves) =====')
    for s,u,nc,ok,nz,br,fx in rows[:25]:
        print(f'  x_{u:<6} cons={nc:<3} modp0={ok}  score {s:>6} ({s-S0:+d})  broke {br:>3} eqs, fixed {fx:>3}  nz={nz[:9]}')
    return rows

R1 = price(21617, 21616, 5040)
R2 = price(29539, 29538, 30163)
json.dump({'a21617': [[r[0],r[1],r[2],r[3],r[4],r[5],r[6]] for r in R1],
           'a29539': [[r[0],r[1],r[2],r[3],r[4],r[5],r[6]] for r in R2]},
          open(os.path.join(HERE,'cl_newton.json'),'w'))
