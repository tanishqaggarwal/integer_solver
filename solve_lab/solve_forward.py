#!/usr/bin/env python3
"""Forward-evaluation solver: seed pins, propagate, then set free inputs to a
default (0) in a leaf-preferring order and propagate, computing all derived
wires correctly. Detects contradictions (constraints that force a nonzero bit)."""
import json, time, sys
from collections import defaultdict, deque
from propagate import load_atoms, atom_vars, substitute, solve_single

NVARS = 38748

class Engine:
    def __init__(self, atoms):
        self.atoms = atoms
        self.var_atoms = defaultdict(list)
        for ai, poly in enumerate(atoms):
            for v in atom_vars(poly):
                self.var_atoms[v].append(ai)
        self.val = [None] * NVARS
        self.domain = {}
        self.contra = []
        self.wl = deque(range(len(atoms)))
        self.inwl = [True] * len(atoms)

    def assign(self, v, x):
        if self.val[v] is not None:
            if self.val[v] != x:
                self.contra.append((v, self.val[v], x))
            return False
        self.val[v] = x
        self.domain.pop(v, None)
        for ai in self.var_atoms[v]:
            if not self.inwl[ai]:
                self.inwl[ai] = True
                self.wl.append(ai)
        return True

    def propagate(self):
        val = self.val
        while self.wl:
            ai = self.wl.popleft()
            self.inwl[ai] = False
            poly = substitute(self.atoms[ai], val)
            uv = atom_vars(poly)
            if len(uv) == 0:
                if poly.get((), 0) != 0:
                    self.contra.append(('const', ai))
                continue
            if len(uv) == 1:
                kind, data = solve_single(poly)
                if kind == 'val':
                    self.assign(*data)
                elif kind == 'dom':
                    u, roots = data
                    if u in self.domain:
                        roots = self.domain[u] & roots
                    self.domain[u] = roots
                    if len(roots) == 1:
                        self.assign(u, next(iter(roots)))
                elif kind == 'contradiction':
                    self.contra.append(('atom', ai))

    def n_assigned(self):
        return sum(1 for x in self.val if x is not None)


def main():
    t0 = time.time()
    atoms = load_atoms()
    eng = Engine(atoms)
    # seed with pins already known: run initial propagation
    eng.propagate()
    print(f"after initial propagation: {eng.n_assigned()} assigned, contra={len(eng.contra)}")

    # inputs list (never-target vars) from summary
    summ = json.load(open('solve_lab/atoms/summary.json'))
    inputs = set(summ['inputs'])

    # Phase 1: set unassigned INPUT vars to 0 (leaf-first), propagate
    order = [v for v in range(NVARS) if v in inputs and eng.val[v] is None]
    for v in order:
        if eng.val[v] is None:
            eng.assign(v, 0)
            eng.propagate()
    print(f"after zeroing inputs: {eng.n_assigned()} assigned, contra={len(eng.contra)}")

    # Phase 2: anything still unassigned -> 0, propagate
    for v in range(NVARS):
        if eng.val[v] is None:
            eng.assign(v, 0)
    eng.propagate()
    print(f"after zero-fill all: {eng.n_assigned()} assigned, contra={len(eng.contra)}")
    print(f"time {time.time()-t0:.1f}s")

    # build assignment and check violated atoms
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
    from collections import Counter
    degc = Counter(max((len(m) for m in atoms[ai]), default=0) for ai in violated)
    bigc = sum(1 for ai in violated if any(abs(c) >= 10**20 for c in atoms[ai].values()))
    print(f"violated by degree: {dict(sorted(degc.items()))}, with-big-const: {bigc}")

    out = {f"x_{i}": v[i] for i in range(NVARS)}
    json.dump(out, open('solve_lab/cand_forward0.json', 'w'))
    json.dump(violated, open('solve_lab/violated_forward0.json', 'w'))
    print("wrote cand_forward0.json, violated_forward0.json")
    print(f"contradictions sample: {eng.contra[:10]}")

if __name__ == '__main__':
    main()
