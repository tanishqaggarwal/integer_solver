#!/usr/bin/env python3
"""Expand every leaf 'atom' term into a canonical integer polynomial, dedup,
and measure reuse across equations. Determines whether the system is a set of
shared residual atoms (each meant to be 0) combined randomly."""
import ast, re, json, time, sys
from math import gcd
from functools import reduce
from collections import defaultdict, Counter

EQ_PATH = __file__.rsplit('/', 1)[0] + '/../EQUATIONS.txt'
OUT = __file__.rsplit('/', 1)[0] + '/atoms'

def const_val(node):
    if isinstance(node, ast.Constant):
        return node.value
    if isinstance(node, ast.UnaryOp) and isinstance(node.op, ast.USub) and isinstance(node.operand, ast.Constant):
        return -node.operand.value
    return None

# ---- polynomial as dict: monomial(tuple sorted var ids) -> int coef ----
def p_const(c):
    return {(): c} if c != 0 else {}
def p_var(i):
    return {(i,): 1}
def p_add(a, b):
    r = dict(a)
    for m, c in b.items():
        r[m] = r.get(m, 0) + c
        if r[m] == 0: del r[m]
    return r
def p_neg(a):
    return {m: -c for m, c in a.items()}
def p_mul(a, b):
    r = {}
    for m1, c1 in a.items():
        for m2, c2 in b.items():
            m = tuple(sorted(m1 + m2))
            r[m] = r.get(m, 0) + c1 * c2
            if r[m] == 0: del r[m]
    return r

def expand(node):
    if isinstance(node, ast.Constant):
        return p_const(node.value)
    if isinstance(node, ast.Name):
        return p_var(int(node.id[2:]))
    if isinstance(node, ast.UnaryOp) and isinstance(node.op, ast.USub):
        return p_neg(expand(node.operand))
    if isinstance(node, ast.BinOp):
        a = expand(node.left); b = expand(node.right)
        if isinstance(node.op, ast.Add): return p_add(a, b)
        if isinstance(node.op, ast.Sub): return p_add(a, p_neg(b))
        if isinstance(node.op, ast.Mult): return p_mul(a, b)
    raise ValueError(ast.dump(node))

def canon(poly):
    """Return canonical hashable key: gcd-reduced, sign-normalized sorted items."""
    if not poly:
        return ()
    g = reduce(gcd, (abs(c) for c in poly.values()))
    items = sorted(poly.items())
    # sign: make coef of smallest monomial positive
    lead = items[0][1] // g
    sign = -1 if lead < 0 else 1
    return tuple((m, sign * c // g) for m, c in items)

def strip_outer(node):
    while True:
        if isinstance(node, ast.BinOp) and isinstance(node.op, ast.Mult):
            a, b = node.left, node.right
            ca, cb = const_val(a), const_val(b)
            if ca is not None and cb is not None: return node
            if ca is not None: node = b; continue
            if cb is not None: node = a; continue
            if ast.unparse(a) == ast.unparse(b): node = a; continue
            return node
        if isinstance(node, ast.UnaryOp) and isinstance(node.op, ast.USub):
            node = node.operand; continue
        return node

def flatten_add(node):
    terms = []
    def rec(n):
        if isinstance(n, ast.BinOp) and isinstance(n.op, ast.Add):
            rec(n.left); rec(n.right)
        else:
            terms.append(n)
    rec(node)
    return terms

def strip_coef(node):
    if isinstance(node, ast.BinOp) and isinstance(node.op, ast.Mult):
        ca, cb = const_val(node.left), const_val(node.right)
        if ca is not None and cb is None: return node.right
        if cb is not None and ca is None: return node.left
    return node

def main():
    t0 = time.time()
    atom_eqs = defaultdict(list)     # canon_key -> list of eq indices
    atom_repr = {}                   # canon_key -> human string (first seen)
    n = 0
    with open(EQ_PATH) as f:
        for i, line in enumerate(f):
            line = line.strip()
            if not line: continue
            n += 1
            lhs = line.rsplit('=', 1)[0]
            core = strip_outer(ast.parse(lhs, mode='eval').body)
            for term in flatten_add(core):
                atom = strip_coef(term)
                poly = expand(atom)
                key = canon(poly)
                if not key:
                    continue
                atom_eqs[key].append(i)
                if key not in atom_repr:
                    atom_repr[key] = ast.unparse(atom)
    print(f"parsed {n} eqs in {time.time()-t0:.1f}s")
    print(f"distinct atoms: {len(atom_eqs)}")
    reuse = Counter(len(v) for v in atom_eqs.values())
    print("reuse distribution (times_atom_appears -> #atoms):")
    for k in sorted(reuse)[:20]:
        print(f"   appears {k}x : {reuse[k]} atoms")
    print(f"   max reuse: {max(reuse)}")
    # classify by degree / shape
    deg_hist = Counter()
    nvars_hist = Counter()
    has_bigconst = 0
    bool_atoms = []
    for key in atom_eqs:
        mons = [m for m, c in key]
        deg = max(len(m) for m in mons)
        deg_hist[deg] += 1
        allv = set()
        for m in mons: allv.update(m)
        nvars_hist[len(allv)] += 1
        # big constant?
        for m, c in key:
            if abs(c) >= 10**20:
                has_bigconst += 1
                break
        # boolean atom: {(a,a):1,(a,):-1}
        if len(key) == 2 and len(allv) == 1:
            d = dict(key)
            a = allv.pop()
            if d.get((a, a)) == 1 and d.get((a,)) == -1:
                bool_atoms.append(a)
    print(f"atom degree histogram: {dict(sorted(deg_hist.items()))}")
    print(f"atom #distinct-vars histogram: {dict(sorted(nvars_hist.items()))}")
    print(f"atoms containing a >=10^20 constant: {has_bigconst}")
    print(f"boolean atoms x*(x-1): {len(bool_atoms)}  sample: {sorted(bool_atoms)[:20]}")

    import os
    os.makedirs(OUT, exist_ok=True)
    # serialize atoms
    with open(OUT + '/poly_atoms.jsonl', 'w') as g:
        for key, eqs in atom_eqs.items():
            g.write(json.dumps({
                "poly": [[list(m), c] for m, c in key],
                "repr": atom_repr[key],
                "n_eq": len(eqs),
                "eqs": eqs[:50],
            }) + "\n")
    with open(OUT + '/bool_vars.json', 'w') as g:
        json.dump(sorted(bool_atoms), g)
    print("wrote atoms/poly_atoms.jsonl, bool_vars.json")

if __name__ == '__main__':
    main()
