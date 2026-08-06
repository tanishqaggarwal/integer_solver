#!/usr/bin/env python3
"""Enhanced propagation with a boolean-forcing rule (case split):
for an atom with a boolean unknown b, if substituting b=0 yields a nonzero
constant (impossible) then force b=1, and vice-versa. This lets a bit-gated
huge atom  bit*(x_B-HUGE)=s*x_C  force bit=1 whenever s*x_C is known nonzero,
cascading reduction/selection decisions instead of defaulting bits to 0."""
import json, time, sys
from collections import defaultdict, deque
from propagate import load_atoms, atom_vars, substitute, solve_single

NVARS = 38748

def boolean_vars(atoms):
    bset = set()
    for poly in atoms:
        vs = atom_vars(poly)
        if len(vs) == 1 and len(poly) == 2:
            v = next(iter(vs))
            if (poly.get((v,)) == 1 and poly.get((v, v)) == -1) or (poly.get((v,)) == -1 and poly.get((v, v)) == 1):
                bset.add(v)
    return bset

def subst_one(poly, var, value):
    out = defaultdict(int)
    for m, c in poly.items():
        coef = c; newm = []
        for v in m:
            if v == var: coef *= value
            else: newm.append(v)
        out[tuple(newm)] += coef
    return {m: c for m, c in out.items() if c != 0}

class Engine:
    def __init__(self, atoms, bset):
        self.atoms = atoms
        self.bset = bset
        self.var_atoms = defaultdict(list)
        for ai, poly in enumerate(atoms):
            for v in atom_vars(poly):
                self.var_atoms[v].append(ai)
        self.val = [None] * NVARS
        self.contra = []
        self.wl = deque(range(len(atoms)))
        self.inwl = [True] * len(atoms)

    def assign(self, v, x):
        if self.val[v] is not None:
            if self.val[v] != x: self.contra.append((v, self.val[v], x))
            return
        self.val[v] = x
        for ai in self.var_atoms[v]:
            if not self.inwl[ai]:
                self.inwl[ai] = True; self.wl.append(ai)

    def propagate(self, use_rule_a=True, max_unknowns_ruleA=4):
        val = self.val
        while self.wl:
            ai = self.wl.popleft(); self.inwl[ai] = False
            poly = substitute(self.atoms[ai], val)
            uv = atom_vars(poly)
            if len(uv) == 0:
                if poly.get((), 0) != 0: self.contra.append(('const', ai))
                continue
            if len(uv) == 1:
                kind, data = solve_single(poly)
                if kind == 'val': self.assign(*data)
                elif kind == 'dom':
                    u, roots = data
                    if len(roots) == 1: self.assign(u, next(iter(roots)))
                elif kind == 'contradiction': self.contra.append(('atom', ai))
                continue
            # Rule A: boolean forcing via case split
            if use_rule_a and len(uv) <= max_unknowns_ruleA:
                bins = [b for b in uv if b in self.bset]
                for b in bins:
                    p0 = subst_one(poly, b, 0)
                    p1 = subst_one(poly, b, 1)
                    bad0 = (len(atom_vars(p0)) == 0 and p0.get((), 0) != 0)
                    bad1 = (len(atom_vars(p1)) == 0 and p1.get((), 0) != 0)
                    if bad0 and not bad1:
                        self.assign(b, 1); break
                    if bad1 and not bad0:
                        self.assign(b, 0); break

    def n_assigned(self):
        return sum(1 for x in self.val if x is not None)

def main():
    t0 = time.time()
    atoms = load_atoms()
    bset = boolean_vars(atoms)
    eng = Engine(atoms, bset)
    eng.propagate()
    print(f"enhanced propagation fixpoint: {eng.n_assigned()} assigned, contra {len(eng.contra)}, {time.time()-t0:.1f}s")
    mainbits = set(json.load(open('main_comp.json'))['main_bits'])
    ones = [b for b in mainbits if eng.val[b] == 1]
    zeros = [b for b in mainbits if eng.val[b] == 0]
    unk = [b for b in mainbits if eng.val[b] is None]
    print(f"main bits: =1 {len(ones)}, =0 {len(zeros)}, unknown {len(unk)}")

    # complete: zero remaining free bits then value wires, propagate
    for v in [b for b in bset if eng.val[b] is None]:
        if eng.val[v] is None: eng.assign(v, 0); eng.propagate()
    for v in range(NVARS):
        if eng.val[v] is None: eng.assign(v, 0); eng.propagate()
    print(f"after completion: {eng.n_assigned()} assigned, contra {len(eng.contra)}")

    v = [x if x is not None else 0 for x in eng.val]
    def ev(poly):
        s = 0
        for m, c in poly.items():
            t = c
            for var in m: t *= v[var]
            s += t
        return s
    violated = [ai for ai, poly in enumerate(atoms) if ev(poly) != 0]
    print(f"violated atoms: {len(violated)}  -> {violated[:12]}")
    json.dump({f"x_{i}": v[i] for i in range(NVARS)}, open('cand_enhanced.json', 'w'))
    print("wrote cand_enhanced.json")

if __name__ == '__main__':
    main()
