#!/usr/bin/env python3
"""Pass 1: parse EQUATIONS.txt exactly with `ast`, expand to canonical integer
polynomials, and cache to disk.

For every equation line i we record:
  * outer:  ('sq', ) if the LHS was  CORE*CORE  (identical unparse),
            ('const', c) if the LHS was  c*CORE  with c a nonzero int constant,
            ('id',)  otherwise  (nothing stripped)
    -- in all these cases  LHS == 0  <=>  CORE == 0  (c != 0).
  * terms:  the top-level '+'-chain of CORE, each as (coef, atom_id)
            where atom_id indexes a table of canonical atom polynomials and
            coef is the literal integer scalar that multiplied it.
  * poly:   the FULL expanded polynomial of CORE (exact, sum of coef*atom).

Cache format: pickle at bitconstraints/cache.pkl
  {'atoms': [poly,...], 'eq_terms': [[(coef,aid),...],...],
   'eq_poly': [poly,...], 'eq_outer': [...]}
A poly is a tuple of (monomial, coef) with monomial a sorted tuple of var ids
(repeats = powers), sorted by monomial.
"""
import ast, sys, time, pickle, os
from math import gcd
from functools import reduce

HERE = os.path.dirname(os.path.abspath(__file__))
EQ_PATH = os.path.join(HERE, '..', '..', 'EQUATIONS.txt')
CACHE = os.path.join(HERE, 'cache.pkl')

sys.setrecursionlimit(100000)


def const_val(node):
    if isinstance(node, ast.Constant):
        return node.value
    if isinstance(node, ast.UnaryOp) and isinstance(node.op, ast.USub):
        v = const_val(node.operand)
        return None if v is None else -v
    return None


def p_const(c):
    return {(): c} if c != 0 else {}


def p_add(a, b):
    r = dict(a)
    for m, c in b.items():
        v = r.get(m, 0) + c
        if v == 0:
            r.pop(m, None)
        else:
            r[m] = v
    return r


def p_neg(a):
    return {m: -c for m, c in a.items()}


def p_mul(a, b):
    r = {}
    for m1, c1 in a.items():
        for m2, c2 in b.items():
            m = tuple(sorted(m1 + m2))
            v = r.get(m, 0) + c1 * c2
            if v == 0:
                r.pop(m, None)
            else:
                r[m] = v
    return r


def expand(node):
    if isinstance(node, ast.Constant):
        return p_const(node.value)
    if isinstance(node, ast.Name):
        return {(int(node.id[2:]),): 1}
    if isinstance(node, ast.UnaryOp) and isinstance(node.op, ast.USub):
        return p_neg(expand(node.operand))
    if isinstance(node, ast.UnaryOp) and isinstance(node.op, ast.UAdd):
        return expand(node.operand)
    if isinstance(node, ast.BinOp):
        a = expand(node.left)
        b = expand(node.right)
        if isinstance(node.op, ast.Add):
            return p_add(a, b)
        if isinstance(node.op, ast.Sub):
            return p_add(a, p_neg(b))
        if isinstance(node.op, ast.Mult):
            return p_mul(a, b)
    raise ValueError(ast.dump(node)[:200])


def freeze(poly):
    return tuple(sorted(poly.items()))


def canon(poly):
    """gcd-reduced, sign-normalised canonical key of a polynomial."""
    if not poly:
        return ()
    g = reduce(gcd, (abs(c) for c in poly.values()))
    items = sorted(poly.items())
    lead = items[0][1] // g
    sign = -1 if lead < 0 else 1
    return tuple((m, sign * c // g) for m, c in items)


def strip_outer(node):
    """Peel nonzero-constant factors, unary minus, and X*X. Returns (node, tag)."""
    tags = []
    while True:
        if isinstance(node, ast.BinOp) and isinstance(node.op, ast.Mult):
            a, b = node.left, node.right
            ca, cb = const_val(a), const_val(b)
            if ca is not None and cb is not None:
                return node, tags
            if ca is not None:
                if ca == 0:
                    return node, tags
                tags.append(('const', ca))
                node = b
                continue
            if cb is not None:
                if cb == 0:
                    return node, tags
                tags.append(('const', cb))
                node = a
                continue
            if ast.dump(a) == ast.dump(b):
                tags.append(('sq',))
                node = a
                continue
            return node, tags
        if isinstance(node, ast.UnaryOp) and isinstance(node.op, ast.USub):
            tags.append(('const', -1))
            node = node.operand
            continue
        return node, tags


def flatten_add(node):
    terms = []
    stack = [(node, 1)]
    while stack:
        n, s = stack.pop()
        if isinstance(n, ast.BinOp) and isinstance(n.op, ast.Add):
            stack.append((n.left, s))
            stack.append((n.right, s))
        else:
            terms.append((n, s))
    return terms


def strip_coef(node):
    """Return (coef, atomnode) peeling literal scalar factors / unary minus."""
    c = 1
    while True:
        if isinstance(node, ast.UnaryOp) and isinstance(node.op, ast.USub):
            c = -c
            node = node.operand
            continue
        if isinstance(node, ast.BinOp) and isinstance(node.op, ast.Mult):
            ca, cb = const_val(node.left), const_val(node.right)
            if ca is not None and cb is None:
                c *= ca
                node = node.right
                continue
            if cb is not None and ca is None:
                c *= cb
                node = node.left
                continue
        return c, node


def main():
    t0 = time.time()
    atom_key_to_id = {}
    atoms = []
    eq_terms = []
    eq_poly = []
    eq_outer = []
    with open(EQ_PATH) as f:
        for i, line in enumerate(f):
            line = line.strip()
            if not line:
                eq_terms.append([])
                eq_poly.append(())
                eq_outer.append([('blank',)])
                continue
            lhs, rhs = line.rsplit('=', 1)
            assert rhs.strip() == '0', (i, rhs[:40])
            core, tags = strip_outer(ast.parse(lhs, mode='eval').body)
            terms = []
            total = {}
            for tnode, sign in flatten_add(core):
                c, atom = strip_coef(tnode)
                c *= sign
                ap = expand(atom)
                key = canon(ap)
                aid = atom_key_to_id.get(key)
                if aid is None:
                    aid = len(atoms)
                    atom_key_to_id[key] = aid
                    atoms.append(key)
                terms.append((c, aid))
                # exact contribution uses the *raw* atom poly, not the canon key
                total = p_add(total, {m: c * v for m, v in ap.items()})
            eq_terms.append(terms)
            eq_poly.append(freeze(total))
            eq_outer.append(tags)
            if (i + 1) % 5000 == 0:
                print(f"  {i+1} eqs  {time.time()-t0:.1f}s  atoms={len(atoms)}",
                      flush=True)
    print(f"parsed {len(eq_terms)} eqs in {time.time()-t0:.1f}s, "
          f"{len(atoms)} distinct atoms")
    with open(CACHE, 'wb') as g:
        pickle.dump({'atoms': atoms, 'eq_terms': eq_terms,
                     'eq_poly': eq_poly, 'eq_outer': eq_outer}, g, 4)
    print("wrote", CACHE, os.path.getsize(CACHE))


if __name__ == '__main__':
    main()
