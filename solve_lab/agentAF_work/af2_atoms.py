#!/usr/bin/env python3
"""agent AF, step 2: atom inventory + definition DAG from my own parse."""
import sys, os, pickle, time
from collections import Counter, defaultdict
sys.setrecursionlimit(100000)
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from af1_parse import peel_sum, is_const, strip_outer

def shape_of(n, d=0):
    """structural shape string, variables/constants abstracted"""
    if n[0] == 'v':
        return 'V'
    if n[0] == 'c':
        return 'K'
    if n[0] == 'neg':
        return '-(%s)' % shape_of(n[1], d+1)
    return '(%s %s %s)' % (shape_of(n[1], d+1), n[0], shape_of(n[2], d+1))

def varsof(n, out):
    if n[0] == 'v':
        out.add(n[1])
    elif n[0] == 'c':
        pass
    elif n[0] == 'neg':
        varsof(n[1], out)
    else:
        varsof(n[1], out); varsof(n[2], out)
    return out

def main():
    d = pickle.load(open(os.path.join(HERE, 'af_ast.pkl'), 'rb'))
    bodies = d['bodies']
    atom_id = {}
    atoms = []
    eq_atoms = []           # eq -> list of (coef, atom_id)
    for i, bl in enumerate(bodies):
        body = bl[0]
        terms = peel_sum(body)
        la = []
        for c, a in terms:
            if a not in atom_id:
                atom_id[a] = len(atoms); atoms.append(a)
            la.append((c, atom_id[a]))
        eq_atoms.append(la)
    print('distinct atoms: %d' % len(atoms))
    sh = Counter(shape_of(a) for a in atoms)
    print('distinct atom shapes: %d' % len(sh))
    for s, k in sh.most_common(40):
        print('  %8d  %s' % (k, s))
    # atoms per equation histogram
    h = Counter(len(x) for x in eq_atoms)
    print('atoms/eq histogram:', sorted(h.items())[:30])
    pickle.dump({'atoms': atoms, 'eq_atoms': eq_atoms}, open(os.path.join(HERE, 'af_atoms.pkl'), 'wb'), 2)

if __name__ == '__main__':
    main()
