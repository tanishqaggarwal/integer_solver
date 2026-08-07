"""S11 step 19: price every repair path for the A1 = 0 chain.

A1 = 0 buys eq 29125 (6 satisfied instead of 5) but moves x_28730, breaking
a7930.  Repairing a7930 through x_24548 then breaks a21617.  Enumerate ALL repair
paths at each link and find one that terminates.
"""
import os, sys, itertools
HERE = os.path.dirname(os.path.abspath(__file__))
LAB = os.path.dirname(HERE)
sys.path.insert(0, os.path.join(LAB, 's9', 'eff'))
sys.path.insert(0, HERE)
import lib as L, tools as T, ad
from frame2 import definer, ORDER, FREE, CHECKS, fwd, score, grad
P = ad.P
SSET = {22229, 22230, 35758, 35759, 35760, 35761, 35762}
FORBID = {2081, 4287}
base = L.load(os.path.join(LAB, 'best', 'new_instance_partial_39026.json'))
fwd(base)

def outside(v):
    av = L.all_atom_values(v)
    return [a for a in range(L.NA) if av[a] and a not in SSET], av

def paths(a, v, av):
    """all repair moves for atom a: direct free, 2-level, and mod-p Newton."""
    out = []
    for w in sorted(set(L.avars[a])):
        if w in FORBID: continue
        tgt = T.solve_lin(a, w, v)
        if tgt is not None and tgt != v[w]:
            if w in FREE: out.append(('direct', w, tgt))
            else:
                d = definer.get(w)
                if d is not None:
                    vv = list(v); vv[w] = tgt
                    for u in sorted(set(L.avars[d])):
                        if u == w or u not in FREE or u in FORBID: continue
                        nv = T.solve_lin(d, u, vv)
                        if nv is not None: out.append((f'via x_{w}', u, nv))
    r = av[a] % P
    if r:
        vm = [x % P for x in v]
        try: g = grad(a, vm)
        except Exception: g = {}
        for k, u, d in sorted((len(L.var_atoms[u]), u, d) for u, d in g.items()
                              if u not in FORBID and d % P)[:12]:
            out.append(('newton', u, v[u] + (-r * pow(d, -1, P)) % P))
    return out

v0 = list(base); v0[28730] = v0[9413] * P
fwd(v0, rounds=8)
nz0, av0 = outside(v0)
print(f'A1 = 0 -> broken outside the seven: {nz0}', flush=True)

results = []
for kind, u, nv in paths(7930, v0, av0):
    v1 = list(v0); v1[u] = nv
    fwd(v1, rounds=8)
    nz1, av1 = outside(v1)
    print(f'  a7930 {kind:<12} x_{u:<7} -> outside {nz1}', flush=True)
    if not nz1:
        print('    *** CHAIN TERMINATES HERE'); results.append((0, [(u, nv)])); continue
    for kind2, u2, nv2 in paths(nz1[0], v1, av1):
        if u2 == u: continue
        v2 = list(v1); v2[u2] = nv2
        fwd(v2, rounds=8)
        nz2, av2 = outside(v2)
        eqs = set()
        for a in nz2: eqs |= set(L.atom2eq[a])
        print(f'      a{nz1[0]} {kind2:<12} x_{u2:<7} -> outside {nz2} '
              f'({len(eqs)} eqs)', flush=True)
        results.append((len(eqs), [(u, nv), (u2, nv2)]))
        if not nz2:
            print('        *** CHAIN TERMINATES -- everything outside the seven is clean')
results.sort(key=lambda t: t[0])
print(f'\nbest chain leaves {results[0][0] if results else "n/a"} equations broken '
      f'outside the seven')
