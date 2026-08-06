"""CL engine: generic mod-p Newton cascade with automatic p-handle detection.

For a residual CHECK atom a:
  * find its p-handle:  a variable w of a defined by a gate  w - wire*handle
    with v[wire] == p and handle a free input with <=2 consumers.
  * Newton move on a free input u in grad(a) so that a == 0 (mod p),
  * then set the handle so a == 0 exactly over Z.
"""
import os, sys, json, collections, time, heapq
HERE='/home/user/integer_solver/solve_lab/s10'; LAB='/home/user/integer_solver/solve_lab'
sys.path.insert(0, os.path.join(LAB,'s9','eff')); sys.path.insert(0, HERE)
import lib as L, tools as T, ad
P = 2**256 - 2**32 - 977
definer, atom_out = L.definer, L.atom_out
FREE = set(ad.FREE)
FORBID = {2081, 4287}

_hcache = {}
def handles_of(a, v):
    """list of (wire, handle) pairs usable to absorb atom a over Z."""
    if a in _hcache: return _hcache[a]
    out = []
    for w in sorted(L.avars[a]):
        d = definer.get(w)
        if d is None: continue
        Pp = L.polys[d]
        # w - A*B  form
        for m, c in Pp.items():
            if len(m) == 2 and w not in m:
                x, y = m
                for wire, h in ((x, y), (y, x)):
                    if v[wire] % P == 0 and v[wire] != 0 and h in FREE and h not in FORBID:
                        out.append((w, d, wire, h))
    _hcache[a] = out
    return out

def absorb(v, a):
    """set a p-handle of check a so that a == 0 exactly.  True on success."""
    for w, d, wire, h in handles_of(a, v):
        tgt = T.solve_lin(a, w, v)
        if tgt is None: continue
        vv = list(v); vv[w] = tgt
        nv = T.solve_lin(d, h, vv)
        if nv is None: continue
        v[h] = nv
        ad.fwd(v, rounds=8)
        if L.evalpoly(L.polys[a], v) == 0: return True
        v[h] = 0
    return False

def newton(v, a, u, vm=None, g=None, av=None):
    """return a NEW state with free input u shifted so that atom a == 0 mod p, then absorbed."""
    if vm is None: vm = [x % P for x in v]
    if g is None: g = ad.grad(a, vm)
    if av is None: av = L.all_atom_values(v)
    d = g.get(u, 0) % P
    if d == 0: return None
    r = av[a] % P
    delta = (-r * pow(d, -1, P)) % P
    w = list(v); w[u] = w[u] + delta
    ad.fwd(w, rounds=8)
    if L.evalpoly(L.polys[a], w) % P != 0: return None
    absorb(w, a)
    return w

def stats(v):
    av = L.all_atom_values(v)
    nz = [a for a in range(L.NA) if av[a]]
    bad = L.failing_eqs(av)
    return av, nz, L.NEQ - len(bad), set(bad)

def is_check(a):
    return a not in atom_out
