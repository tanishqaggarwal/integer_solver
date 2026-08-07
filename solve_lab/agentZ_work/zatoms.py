#!/usr/bin/env python3
"""Agent Z: polynomial normalisation of atoms + the multi-selector census.

Every expansion is validated against a direct numeric evaluation of the parse
tree at random points (standing rule: never trust a symbolic expansion).
"""
import os, sys, json, pickle, random, collections
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from zparse import parse, key, varset, flat_add, flat_mul, split_coeff, reduce_L, atoms_of, ev

HERE = os.path.dirname(os.path.abspath(__file__))
EQ = os.path.join(HERE, '..', '..', 'EQUATIONS.txt')

# ---------------------------------------------------------------- polynomial
def pmul(a, b):
    o = {}
    for ka, va in a.items():
        for kb, vb in b.items():
            k = tuple(sorted(ka + kb))
            o[k] = o.get(k, 0) + va * vb
    return {k: v for k, v in o.items() if v}

def padd(a, b):
    o = dict(a)
    for k, v in b.items():
        o[k] = o.get(k, 0) + v
        if o[k] == 0:
            del o[k]
    return o

def poly(t):
    if t[0] == 'n':
        return {(): t[1]} if t[1] else {}
    if t[0] == 'v':
        return {(t[1],): 1}
    if t[0] == '+':
        o = {}
        for u in t[1]:
            o = padd(o, poly(u))
        return o
    o = {(): 1}
    for u in t[1]:
        o = pmul(o, poly(u))
    return o

def pev(p, env):
    s = 0
    for mon, c in p.items():
        t = c
        for v in mon:
            t *= env.get(v, 0)
        s += t
    return s

def pkey(p):
    return tuple(sorted(p.items()))

# ---------------------------------------------------------------- main
def main():
    sys.setrecursionlimit(100000)
    sel = set(json.load(open(os.path.join(HERE, 'zsel.json')))['selectors'])
    print("selectors:", len(sel))

    lines = []
    with open(EQ) as f:
        for line in f:
            line = line.strip()
            if line:
                lines.append(line.rsplit('=', 1)[0])

    rng = random.Random(20260807)
    atom_poly = {}        # pkey -> poly
    atom_eqs = collections.defaultdict(list)
    eq_L_vars = []        # per eq: (selvars, othervars)
    multi_atoms = {}      # pkey -> poly, atoms with >=2 distinct selectors
    eq_atomkeys = []
    nval = 0
    ndeg = collections.Counter()

    for i, lhs in enumerate(lines):
        E = parse(lhs)
        L, note = reduce_L(E)
        vs = varset(L)
        eq_L_vars.append((frozenset(vs & sel), frozenset(vs - sel)))
        aks = []
        for c, a in atoms_of(L):
            p = poly(a)
            k = pkey(p)
            if k not in atom_poly:
                atom_poly[k] = p
                # validate expansion numerically, 3 random points
                av = sorted(varset(a))
                for _ in range(3):
                    env = {v: rng.randrange(-50, 50) for v in av}
                    assert ev(a, env) == pev(p, env), ("expansion mismatch", i)
                nval += 1
                ndeg[max((len(m) for m in p), default=0)] += 1
                sv = set()
                for m in p:
                    sv |= (set(m) & sel)
                if len(sv) >= 2:
                    multi_atoms[k] = p
            aks.append((c, k))
            atom_eqs[k].append(i)
        eq_atomkeys.append(aks)
        if i % 5000 == 0:
            print("  ...", i, flush=True)

    print("distinct atoms (polynomial-normalised):", len(atom_poly))
    print("expansions validated:", nval)
    print("atom degree histogram:", sorted(ndeg.items()))
    print("ATOMS TOUCHING >=2 DISTINCT SELECTORS:", len(multi_atoms))

    # equation-level census
    ge2 = [i for i, (s, o) in enumerate(eq_L_vars) if len(s) >= 2]
    onlysel = [i for i, (s, o) in enumerate(eq_L_vars) if s and not o]
    print("EQUATIONS whose linear form mentions >=2 distinct selectors:", len(ge2))
    print("EQUATIONS whose variables are ALL selectors:", len(onlysel), onlysel[:50])
    hist = collections.Counter(len(s) for s, o in eq_L_vars)
    print("per-equation selector-count histogram:", sorted(hist.items()))

    pickle.dump({'atom_poly': atom_poly, 'multi_atoms': multi_atoms,
                 'eq_L_vars': eq_L_vars, 'eq_atomkeys': eq_atomkeys,
                 'atom_eqs': dict(atom_eqs)},
                open(os.path.join(HERE, 'zatoms.pkl'), 'wb'))
    print("saved zatoms.pkl")

if __name__ == '__main__':
    main()
