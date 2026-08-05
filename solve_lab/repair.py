#!/usr/bin/env python3
"""Provenance-tracking forward solver + backward trace to find the free bits
controlling the residual violations, then z3-search over just those bits."""
import json, time, sys
from collections import defaultdict, deque
from propagate import load_atoms, atom_vars, substitute, solve_single

NVARS = 38748

class ProvEngine:
    def __init__(self, atoms):
        self.atoms = atoms
        self.var_atoms = defaultdict(list)
        for ai, poly in enumerate(atoms):
            for v in atom_vars(poly):
                self.var_atoms[v].append(ai)
        self.val = [None] * NVARS
        self.prov = [None] * NVARS      # (atom_index or 'free', deps tuple)
        self.domain = {}
        self.contra = []
        self.wl = deque(range(len(atoms)))
        self.inwl = [True] * len(atoms)

    def assign(self, v, x, prov):
        if self.val[v] is not None:
            if self.val[v] != x:
                self.contra.append((v, self.val[v], x, prov))
            return
        self.val[v] = x
        self.prov[v] = prov
        self.domain.pop(v, None)
        for ai in self.var_atoms[v]:
            if not self.inwl[ai]:
                self.inwl[ai] = True
                self.wl.append(ai)

    def propagate(self):
        val = self.val
        while self.wl:
            ai = self.wl.popleft()
            self.inwl[ai] = False
            orig = self.atoms[ai]
            poly = substitute(orig, val)
            uv = atom_vars(poly)
            if len(uv) == 0:
                if poly.get((), 0) != 0:
                    self.contra.append(('const', ai))
                continue
            if len(uv) == 1:
                kind, data = solve_single(poly)
                deps = tuple(sorted(atom_vars(orig) - uv))
                if kind == 'val':
                    self.assign(data[0], data[1], (ai, deps))
                elif kind == 'dom':
                    u, roots = data
                    if u in self.domain: roots = self.domain[u] & roots
                    self.domain[u] = roots
                    if len(roots) == 1:
                        self.assign(u, next(iter(roots)), (ai, deps))
                elif kind == 'contradiction':
                    self.contra.append(('atom', ai))

    def n_assigned(self):
        return sum(1 for x in self.val if x is not None)


def boolean_vars(atoms):
    bset = set()
    for poly in atoms:
        vs = atom_vars(poly)
        if len(vs) == 1 and len(poly) == 2:
            v = next(iter(vs))
            if (poly.get((v,)) == 1 and poly.get((v, v)) == -1) or (poly.get((v,)) == -1 and poly.get((v, v)) == 1):
                bset.add(v)
    return bset


def main():
    t0 = time.time()
    atoms = load_atoms()
    bset = boolean_vars(atoms)
    eng = ProvEngine(atoms)
    eng.propagate()
    n_pin = eng.n_assigned()

    free_bits = [v for v in bset if eng.val[v] is None]
    for v in free_bits:
        if eng.val[v] is None:
            eng.assign(v, 0, ('free', ()))
            eng.propagate()
    remaining = [v for v in range(NVARS) if eng.val[v] is None]
    for v in remaining:
        if eng.val[v] is None:
            eng.assign(v, 0, ('free', ()))
            eng.propagate()
    print(f"assigned {eng.n_assigned()}, contra {len(eng.contra)}, time {time.time()-t0:.1f}s")

    v = [x if x is not None else 0 for x in eng.val]
    def ev(poly):
        s = 0
        for m, c in poly.items():
            t = c
            for var in m: t *= v[var]
            s += t
        return s
    violated = [ai for ai, poly in enumerate(atoms) if ev(poly) != 0]
    print(f"violated atoms: {len(violated)}")

    # backward trace from violated atoms to free bits
    control = set()
    freevars_seen = set()
    seen = set()
    stack = []
    for ai in violated:
        stack.extend(atom_vars(atoms[ai]))
    while stack:
        w = stack.pop()
        if w in seen: continue
        seen.add(w)
        p = eng.prov[w]
        if p is None: continue
        tag, deps = p
        if tag == 'free':
            freevars_seen.add(w)
            if w in bset:
                control.add(w)
            continue
        for d in deps:
            if d not in seen:
                stack.append(d)
    print(f"backward cone: {len(seen)} vars, free vars {len(freevars_seen)}, control BITS {len(control)}")
    print(f"control bits sample: {sorted(control)[:40]}")

    json.dump({f"x_{i}": v[i] for i in range(NVARS)}, open('solve_lab/cand_repair0.json', 'w'))
    json.dump(sorted(control), open('solve_lab/control_bits.json', 'w'))
    json.dump(violated, open('solve_lab/violated_repair0.json', 'w'))
    # also save the non-bit free vars in the cone (value inputs that might need setting)
    freeval = sorted(freevars_seen - control)
    json.dump(freeval, open('solve_lab/control_freevals.json', 'w'))
    print(f"non-bit free vars in cone: {len(freeval)} sample {freeval[:20]}")
    print("wrote control_bits.json, control_freevals.json, cand_repair0.json")

if __name__ == '__main__':
    main()
