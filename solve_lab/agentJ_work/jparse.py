#!/usr/bin/env python3
"""Independent parse of EQUATIONS.txt into atoms (agent J).

Strategy: rewrite x_123 -> X123 and use Python's own `ast` module to get a
faithful expression tree.  Then peel outer structure:
   LHS forms observed: (C)*(S), (C1)*(S)+(C2)*(S), (S)*(S), (S), (C)*((-1)*(S))
and decompose the core S as a left-nested sum  A0 + c1*A1 + c2*A2 + ...
Each Ai is an "atom".

Outputs a pickle:  jmodel.pkl  with
   eqs   : list of dicts {kind, mult, atoms:[(coef, atom_id)]}
   atoms : list of ast-normalised atom source strings
"""
import ast, re, sys, pickle, time, os

HERE = os.path.dirname(os.path.abspath(__file__))
EQ = os.path.join(HERE, '..', '..', 'EQUATIONS.txt')
VAR = re.compile(r'x_(\d+)')


def norm(node):
    """Canonical string for an ast node."""
    return ast.dump(node)


def unparse(node):
    return ast.unparse(node)


def is_const_int(node):
    """Return int value if node is an integer literal (possibly negated)."""
    if isinstance(node, ast.Constant) and isinstance(node.value, int):
        return node.value
    if isinstance(node, ast.UnaryOp) and isinstance(node.op, ast.USub):
        v = is_const_int(node.operand)
        return None if v is None else -v
    return None


def split_sum(node):
    """Flatten a left-nested Add/Sub chain into list of (sign, term)."""
    out = []
    def rec(n, s):
        if isinstance(n, ast.BinOp) and isinstance(n.op, ast.Add):
            rec(n.left, s); rec(n.right, s)
        elif isinstance(n, ast.BinOp) and isinstance(n.op, ast.Sub):
            rec(n.left, s); rec(n.right, -s)
        else:
            out.append((s, n))
    rec(node, 1)
    return out


def peel_coef(node):
    """term -> (coef, atomnode).  Handles (c)*(A) and (A)."""
    if isinstance(node, ast.BinOp) and isinstance(node.op, ast.Mult):
        cl = is_const_int(node.left)
        cr = is_const_int(node.right)
        if cl is not None and cr is None:
            return cl, node.right
        if cr is not None and cl is None:
            return cr, node.left
    return 1, node


def parse_line(line):
    lhs = line.rsplit('=', 1)[0].strip()
    src = VAR.sub(r'X\1', lhs)
    tree = ast.parse(src, mode='eval').body
    return tree


def core_of(tree):
    """Peel the outer wrapper. Returns (kind, mult, corenode)."""
    # (S)*(S) square, or (C)*(S)
    if isinstance(tree, ast.BinOp) and isinstance(tree.op, ast.Mult):
        l, r = tree.left, tree.right
        cl, cr = is_const_int(l), is_const_int(r)
        if cl is not None and cr is None:
            k, m, c = core_of(r)
            return k, m * cl, c
        if cr is not None and cl is None:
            k, m, c = core_of(l)
            return k, m * cr, c
        if norm(l) == norm(r):
            k, m, c = core_of(l)
            return 'sq', m, c
    if isinstance(tree, ast.BinOp) and isinstance(tree.op, ast.Add):
        # (C1)*(S)+(C2)*(S)
        c1, a1 = peel_coef(tree.left)
        c2, a2 = peel_coef(tree.right)
        if norm(a1) == norm(a2):
            k, m, c = core_of(a1)
            return k, m * (c1 + c2), c
    if isinstance(tree, ast.UnaryOp) and isinstance(tree.op, ast.USub):
        k, m, c = core_of(tree.operand)
        return k, -m, c
    return 'lin', 1, tree


def main():
    t0 = time.time()
    eqs = []
    atom_ids = {}
    atoms = []
    shapes = {}
    with open(EQ) as f:
        for idx, line in enumerate(f):
            line = line.strip()
            if not line:
                continue
            tree = parse_line(line)
            kind, mult, core = core_of(tree)
            terms = split_sum(core)
            lst = []
            for s, t in terms:
                c, a = peel_coef(t)
                key = unparse(a)
                aid = atom_ids.get(key)
                if aid is None:
                    aid = len(atoms); atom_ids[key] = aid; atoms.append(key)
                lst.append((s * c, aid))
            eqs.append({'i': idx, 'kind': kind, 'mult': mult, 'terms': lst})
            shapes[kind] = shapes.get(kind, 0) + 1
            if idx % 5000 == 0:
                print(f"  {idx} ... {time.time()-t0:.1f}s", file=sys.stderr)
    print(f"parsed {len(eqs)} eqs, {len(atoms)} distinct atoms in {time.time()-t0:.1f}s")
    print("kinds:", shapes)
    with open(os.path.join(HERE, 'jmodel.pkl'), 'wb') as f:
        pickle.dump({'eqs': eqs, 'atoms': atoms}, f)


if __name__ == '__main__':
    main()
