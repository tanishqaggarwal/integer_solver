#!/usr/bin/env python3
"""Normalize each atom to a polynomial dict {monomial(tuple sorted vars): coef}."""
import pickle, os, re, ast, sys, collections

sys.setrecursionlimit(100000)
HERE = os.path.dirname(os.path.abspath(__file__))


def polymul(a, b):
    r = collections.defaultdict(int)
    for ma, ca in a.items():
        for mb, cb in b.items():
            r[tuple(sorted(ma + mb))] += ca * cb
    return {m: c for m, c in r.items() if c}


def polyadd(a, b, s=1):
    r = dict(a)
    for m, c in b.items():
        r[m] = r.get(m, 0) + s * c
        if r[m] == 0:
            del r[m]
    return r


def to_poly(node):
    if isinstance(node, ast.Constant):
        return {(): node.value} if node.value else {}
    if isinstance(node, ast.Name):
        return {(int(node.id[1:]),): 1}
    if isinstance(node, ast.UnaryOp):
        p = to_poly(node.operand)
        if isinstance(node.op, ast.USub):
            return {m: -c for m, c in p.items()}
        return p
    if isinstance(node, ast.BinOp):
        l = to_poly(node.left); r = to_poly(node.right)
        if isinstance(node.op, ast.Mult):
            return polymul(l, r)
        if isinstance(node.op, ast.Add):
            return polyadd(l, r, 1)
        if isinstance(node.op, ast.Sub):
            return polyadd(l, r, -1)
    raise ValueError(ast.dump(node))


def build():
    D = pickle.load(open(os.path.join(HERE, 'atoms.pkl'), 'rb'))
    polys = []
    for s in D['atom_src']:
        polys.append(to_poly(ast.parse(s, mode='eval').body))
    return D, polys


if __name__ == '__main__':
    D, polys = build()
    with open(os.path.join(HERE, 'polys.pkl'), 'wb') as f:
        pickle.dump(polys, f)
    # classification
    kinds = collections.Counter()
    for p in polys:
        deg = max((len(m) for m in p), default=0)
        nt = len(p)
        kinds[(deg, nt)] += 1
    print("(degree, nterms) histogram:")
    for k, v in sorted(kinds.items()):
        print("   ", k, v)
    print("total", len(polys))
