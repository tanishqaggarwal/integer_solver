#!/usr/bin/env python3
"""Independent parse of EQUATIONS.txt into atoms (agent J).

Rewrite x_123 -> X123, use Python's `ast`.  Peel outer wrapper
   (C)*(S) | (C1)*(S)+(C2)*(S) | (S)*(S) | (S) | (C)*((-1)*(S))
then decompose S as the LEFT-NESTED chain  A0 + c1*A1 + c2*A2 + ...
(each appended term is literally `(int)*(atom)`).
"""
import ast, re, sys, pickle, time, os
from collections import Counter

HERE = os.path.dirname(os.path.abspath(__file__))
EQ = os.path.join(HERE, '..', '..', 'EQUATIONS.txt')
VAR = re.compile(r'x_(\d+)')


def norm(node):
    return ast.dump(node)


def is_const_int(node):
    if isinstance(node, ast.Constant) and isinstance(node.value, int):
        return node.value
    if isinstance(node, ast.UnaryOp) and isinstance(node.op, ast.USub):
        v = is_const_int(node.operand)
        return None if v is None else -v
    return None


def split_coef_mult(node):
    """If node is (int)*(expr) or (expr)*(int) return (c, expr) else None."""
    if isinstance(node, ast.BinOp) and isinstance(node.op, ast.Mult):
        cl = is_const_int(node.left)
        cr = is_const_int(node.right)
        if cl is not None and cr is None:
            return cl, node.right
        if cr is not None and cl is None:
            return cr, node.left
    return None


def core_of(tree):
    if isinstance(tree, ast.BinOp) and isinstance(tree.op, ast.Mult):
        l, r = tree.left, tree.right
        cl, cr = is_const_int(l), is_const_int(r)
        if cl is not None and cr is None:
            k, m, c = core_of(r); return k, m * cl, c
        if cr is not None and cl is None:
            k, m, c = core_of(l); return k, m * cr, c
        if norm(l) == norm(r):
            k, m, c = core_of(l); return 'sq', m, c
    if isinstance(tree, ast.BinOp) and isinstance(tree.op, ast.Add):
        a = split_coef_mult(tree.left); b = split_coef_mult(tree.right)
        if a and b and norm(a[1]) == norm(b[1]):
            k, m, c = core_of(a[1]); return k, m * (a[0] + b[0]), c
    if isinstance(tree, ast.UnaryOp) and isinstance(tree.op, ast.USub):
        k, m, c = core_of(tree.operand); return k, -m, c
    return 'lin', 1, tree


def decompose(node):
    """S -> [(coef, atomnode)] using the left-nested chain rule."""
    terms = []
    cur = node
    while isinstance(cur, ast.BinOp) and isinstance(cur.op, (ast.Add, ast.Sub)):
        cm = split_coef_mult(cur.right)
        if cm is None:
            break
        sgn = 1 if isinstance(cur.op, ast.Add) else -1
        terms.append((sgn * cm[0], cm[1]))
        cur = cur.left
    # cur is A0 (possibly with its own leading integer factor)
    cm = split_coef_mult(cur)
    if cm is not None:
        terms.append((cm[0], cm[1]))
    else:
        terms.append((1, cur))
    terms.reverse()
    return terms


def main():
    t0 = time.time()
    eqs = []
    atom_ids = {}
    atoms = []
    shapes = Counter()
    with open(EQ) as f:
        for idx, line in enumerate(f):
            line = line.strip()
            if not line:
                continue
            lhs = line.rsplit('=', 1)[0].strip()
            src = VAR.sub(r'X\1', lhs)
            tree = ast.parse(src, mode='eval').body
            kind, mult, core = core_of(tree)
            lst = []
            for c, a in decompose(core):
                key = ast.unparse(a)
                aid = atom_ids.get(key)
                if aid is None:
                    aid = len(atoms); atom_ids[key] = aid; atoms.append(key)
                lst.append((c, aid))
            eqs.append({'i': idx, 'kind': kind, 'mult': mult, 'terms': lst})
            shapes[kind] += 1
            if idx % 10000 == 0:
                print(f"  {idx} ... {time.time()-t0:.1f}s", file=sys.stderr)
    print(f"parsed {len(eqs)} eqs, {len(atoms)} distinct atoms in {time.time()-t0:.1f}s")
    print("kinds:", dict(shapes))
    print("terms/eq histogram:", Counter(len(e['terms']) for e in eqs).most_common(15))
    with open(os.path.join(HERE, 'jmodel.pkl'), 'wb') as f:
        pickle.dump({'eqs': eqs, 'atoms': atoms}, f)
    # sample atoms
    for a in atoms[:20]:
        print("ATOM:", a)


if __name__ == '__main__':
    main()
