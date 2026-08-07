#!/usr/bin/env python3
"""Expand every atom into a monomial dict, classify, and build the definition DAG."""
import ast, re, os, pickle, sys, time
from collections import Counter, defaultdict

HERE = os.path.dirname(os.path.abspath(__file__))
VAR = re.compile(r'x_(\d+)')


def expand(node):
    """ast node -> dict {sorted tuple of var ids : int coef}"""
    if isinstance(node, ast.Constant):
        return {(): int(node.value)}
    if isinstance(node, ast.Name):
        return {(int(node.id[1:]),): 1}
    if isinstance(node, ast.UnaryOp):
        if isinstance(node.op, ast.USub):
            return {k: -v for k, v in expand(node.operand).items()}
        if isinstance(node.op, ast.UAdd):
            return expand(node.operand)
        raise ValueError('unop')
    if isinstance(node, ast.BinOp):
        a = expand(node.left); b = expand(node.right)
        if isinstance(node.op, ast.Add):
            r = dict(a)
            for k, v in b.items():
                r[k] = r.get(k, 0) + v
        elif isinstance(node.op, ast.Sub):
            r = dict(a)
            for k, v in b.items():
                r[k] = r.get(k, 0) - v
        elif isinstance(node.op, ast.Mult):
            r = {}
            for k1, v1 in a.items():
                for k2, v2 in b.items():
                    k = tuple(sorted(k1 + k2))
                    r[k] = r.get(k, 0) + v1 * v2
        else:
            raise ValueError('binop')
        return {k: v for k, v in r.items() if v}
    raise ValueError('node ' + ast.dump(node))


def main():
    M = pickle.load(open(os.path.join(HERE, 'jmodel2.pkl'), 'rb'))
    atoms = M['atoms']
    t0 = time.time()
    polys = []
    for a in atoms:
        src = VAR.sub(r'X\1', a)
        polys.append(expand(ast.parse(src, mode='eval').body))
    print(f"expanded {len(polys)} atoms in {time.time()-t0:.1f}s")

    degc = Counter()
    for p in polys:
        degc[max((len(k) for k in p), default=0)] += 1
    print("degree histogram:", dict(degc))

    # classification: which variables occur linearly-alone with coef +-1
    shapes = Counter()
    defs = []          # (atom_id, defined_var) candidates
    for i, p in enumerate(polys):
        lin = {k[0]: v for k, v in p.items() if len(k) == 1}
        higher = set()
        for k in p:
            if len(k) >= 2:
                higher.update(k)
        cands = [v for v, c in lin.items() if abs(c) == 1 and v not in higher]
        defs.append(cands)
        shapes[(len(p), max((len(k) for k in p), default=0), len(cands))] += 1
    print("(#monomials, deg, #defcands) top:", shapes.most_common(20))
    print("atoms with 0 def candidates:", sum(1 for d in defs if not d))
    pickle.dump({'polys': polys, 'defcands': defs},
                open(os.path.join(HERE, 'jpoly.pkl'), 'wb'))


if __name__ == '__main__':
    main()
