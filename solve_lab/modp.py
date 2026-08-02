#!/usr/bin/env python3
"""Mod-P propagation engine: identical logic to the integer engine but all
arithmetic is done modulo a large prime P. Keeps values bounded (no blow-up),
so evaluation stays fast even with many bits set. An atom that is nonzero in Z
is nonzero mod P with overwhelming probability, so the mod-P violated-atom count
faithfully tracks the true count."""
import json, time
from collections import defaultdict, deque
from propagate import load_atoms, atom_vars

NVARS = 38748
P = (1 << 61) - 1

_invcache = {}
def inv(a):
    a %= P
    r = _invcache.get(a)
    if r is None:
        r = pow(a, P - 2, P); _invcache[a] = r
    return r

def substitute_modp(poly, val):
    out = defaultdict(int)
    for m, c in poly.items():
        coef = c % P; newm = []
        for v in m:
            if val[v] is not None:
                coef = (coef * val[v]) % P
            else:
                newm.append(v)
        out[tuple(sorted(newm))] = (out[tuple(sorted(newm))] + coef) % P
    return {m: c % P for m, c in out.items() if c % P != 0}

def solve_single_modp(poly):
    uv = set()
    for m in poly: uv.update(m)
    u = next(iter(uv))
    c0 = c1 = c2 = 0
    for m, c in poly.items():
        d = len(m)
        if d == 0: c0 = (c0 + c) % P
        elif d == 1: c1 = (c1 + c) % P
        elif d == 2: c2 = (c2 + c) % P
        else: return ('skip', None)
    if c2 % P == 0:
        if c1 % P == 0:
            return ('contra', None) if c0 % P != 0 else ('skip', None)
        return ('val', (u, (-c0 * inv(c1)) % P))
    # quadratic mod P: x*(x-1) style boolean -> domain {0,1}; else skip
    if c0 % P == 0 and (c1 + c2) % P == 0:
        return ('dom', (u, {0, 1}))
    return ('skip', None)

class ModPEngine:
    def __init__(self, atoms):
        self.atoms = atoms
        self.var_atoms = defaultdict(list)
        for ai, poly in enumerate(atoms):
            for v in atom_vars(poly):
                self.var_atoms[v].append(ai)
        self.val = [None] * NVARS
        self.contra = 0
        self.wl = deque(range(len(atoms)))
        self.inwl = [True] * len(atoms)

    def assign(self, v, x):
        x %= P
        if self.val[v] is not None:
            if self.val[v] != x: self.contra += 1
            return
        self.val[v] = x
        for ai in self.var_atoms[v]:
            if not self.inwl[ai]:
                self.inwl[ai] = True; self.wl.append(ai)

    def propagate(self):
        while self.wl:
            ai = self.wl.popleft(); self.inwl[ai] = False
            poly = substitute_modp(self.atoms[ai], self.val)
            uv = atom_vars(poly)
            if len(uv) == 0:
                if poly.get((), 0) % P != 0: self.contra += 1
                continue
            if len(uv) == 1:
                k, data = solve_single_modp(poly)
                if k == 'val': self.assign(*data)
                elif k == 'dom':
                    u, roots = data
                    if len(roots) == 1: self.assign(u, next(iter(roots)))
                elif k == 'contra': self.contra += 1

    def viol_count(self):
        val = [x if x is not None else 0 for x in self.val]
        vi = 0
        for poly in self.atoms:
            s = 0
            for m, c in poly.items():
                t = c % P
                for x in m: t = (t * val[x]) % P
                s = (s + t) % P
            if s % P != 0: vi += 1
        return vi

def make_base(atoms):
    eng = ModPEngine(atoms); eng.propagate()
    return list(eng.val)

def eval_bits(atoms, bset, base_val, ones):
    eng = ModPEngine(atoms)
    eng.val = list(base_val); eng.wl = deque(); eng.inwl = [False] * len(atoms)
    for b in ones:
        if eng.val[b] is None: eng.assign(b, 1)
    eng.propagate()
    for v in [b for b in bset if eng.val[b] is None]:
        if eng.val[v] is None: eng.assign(v, 0); eng.propagate()
    for v in range(NVARS):
        if eng.val[v] is None: eng.assign(v, 0); eng.propagate()
    return eng.viol_count(), eng

if __name__ == '__main__':
    from repair import boolean_vars
    t0 = time.time()
    atoms = load_atoms(); bset = boolean_vars(atoms)
    base = make_base(atoms)
    print(f"base propagation modp in {time.time()-t0:.1f}s")
    t1 = time.time()
    v0, _ = eval_bits(atoms, bset, base, [])
    print(f"bits=0: {v0} violated atoms mod P, eval={time.time()-t1:.2f}s  (Z gives 4)")
