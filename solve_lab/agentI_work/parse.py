#!/usr/bin/env python3
"""Independent parser for EQUATIONS.txt (agent I).

Produces:
  atoms.pkl : {
    'atom_src'   : list[str]           canonical source of each distinct atom
    'atom_ast'   : list[ast node]      (not pickled; rebuilt on demand)
    'eq_terms'   : list[list[(coef,int atom_id)]]   core of each equation
    'eq_outer'   : list[str]           description of outer wrapper
    'atom_vars'  : list[frozenset[int]]
  }
"""
import sys, ast, re, pickle, time, os

sys.setrecursionlimit(100000)

HERE = os.path.dirname(os.path.abspath(__file__))
EQ = os.path.join(HERE, '..', '..', 'EQUATIONS.txt')
VAR_RE = re.compile(r'x_(\d+)')


def dump(node):
    """Canonical source string of an AST node (structural key)."""
    return ast.dump(node)


def as_int(node):
    if isinstance(node, ast.Constant) and isinstance(node.value, int):
        return node.value
    if isinstance(node, ast.UnaryOp) and isinstance(node.op, ast.USub):
        v = as_int(node.operand)
        return None if v is None else -v
    return None


def strip_outer(node):
    """Reduce  c*P, P*P, c1*P+c2*P, (-1)*P, P*P*P ... down to core P."""
    outer = []
    while True:
        changed = False
        if isinstance(node, ast.BinOp) and isinstance(node.op, ast.Mult):
            l, r = node.left, node.right
            li, ri = as_int(l), as_int(r)
            if li is not None and ri is None:
                outer.append(('mul', li)); node = r; changed = True
            elif ri is not None and li is None:
                outer.append(('mul', ri)); node = l; changed = True
            elif li is None and ri is None and dump(l) == dump(r):
                outer.append(('sq', 0)); node = l; changed = True
            elif li is None and ri is None:
                # A*A*A parses as (A*A)*A ; check
                if isinstance(l, ast.BinOp) and isinstance(l.op, ast.Mult) \
                   and dump(l.left) == dump(l.right) == dump(r):
                    outer.append(('cube', 0)); node = r; changed = True
        elif isinstance(node, ast.BinOp) and isinstance(node.op, ast.Add):
            # c1*P + c2*P
            l, r = node.left, node.right
            def split(n):
                if isinstance(n, ast.BinOp) and isinstance(n.op, ast.Mult):
                    li, ri = as_int(n.left), as_int(n.right)
                    if li is not None: return li, n.right
                    if ri is not None: return ri, n.left
                return None, None
            c1, p1 = split(l); c2, p2 = split(r)
            if p1 is not None and p2 is not None and dump(p1) == dump(p2):
                outer.append(('lin', c1 + c2)); node = p1; changed = True
        if not changed:
            break
    return node, outer


def flatten_core(node, out):
    """core is left-nested: (((a0) + c1*a1) + c2*a2) ... -> list of (coef, atomnode)."""
    if isinstance(node, ast.BinOp) and isinstance(node.op, ast.Add):
        l, r = node.left, node.right
        # right must be coef*atom or bare atom
        if isinstance(r, ast.BinOp) and isinstance(r.op, ast.Mult):
            ri = as_int(r.left)
            if ri is not None:
                flatten_core(l, out)
                out.append((ri, r.right))
                return
            ri = as_int(r.right)
            if ri is not None:
                flatten_core(l, out)
                out.append((ri, r.left))
                return
        flatten_core(l, out)
        out.append((1, r))
        return
    if isinstance(node, ast.BinOp) and isinstance(node.op, ast.Sub):
        l, r = node.left, node.right
        # only split if it's clearly the sum-chain; heuristic: don't split Sub
        out.append((1, node)); return
    out.append((1, node))


def main():
    t0 = time.time()
    atom_key = {}
    atom_src = []
    atom_vars = []
    eq_terms = []
    eq_outer = []
    n = 0
    with open(EQ) as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            lhs = line.rsplit('=', 1)[0]
            src = VAR_RE.sub(r'X\1', lhs)
            tree = ast.parse(src, mode='eval').body
            core, outer = strip_outer(tree)
            terms = []
            flatten_core(core, terms)
            tl = []
            for c, a in terms:
                k = dump(a)
                aid = atom_key.get(k)
                if aid is None:
                    aid = len(atom_src)
                    atom_key[k] = aid
                    s = ast.unparse(a)
                    atom_src.append(s)
                    atom_vars.append(frozenset(int(m) for m in re.findall(r'X(\d+)', s)))
                tl.append((c, aid))
            eq_terms.append(tl)
            eq_outer.append(outer)
            n += 1
            if n % 5000 == 0:
                print(f"  {n} eqs, {len(atom_src)} atoms, {time.time()-t0:.0f}s", flush=True)
    print(f"parsed {n} equations, {len(atom_src)} distinct atoms in {time.time()-t0:.0f}s")
    with open(os.path.join(HERE, 'atoms.pkl'), 'wb') as f:
        pickle.dump({'atom_src': atom_src, 'atom_vars': atom_vars,
                     'eq_terms': eq_terms, 'eq_outer': eq_outer}, f)
    print("wrote atoms.pkl")


if __name__ == '__main__':
    main()
