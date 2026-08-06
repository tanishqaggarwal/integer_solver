"""Decompose every equation into (coeff, atom) pairs. Validated by exact re-evaluation."""
import ast, re, json, pickle, sys, time, os
import harness as H

def unfold(node):
    """Return list of (coeff, atom_node) for the top-level additive atom chain."""
    out = []
    while isinstance(node, ast.BinOp) and isinstance(node.op, ast.Add):
        r = node.right
        c, a = split_scalar(r)
        if c is None:
            break
        out.append((c, a))
        node = node.left
    out.append((1, node))
    out.reverse()
    return out

def const_of(node):
    if isinstance(node, ast.Constant) and isinstance(node.value, int):
        return node.value
    if isinstance(node, ast.UnaryOp) and isinstance(node.op, ast.USub):
        c = const_of(node.operand)
        return None if c is None else -c
    return None

def split_scalar(node):
    """If node is (int)*(X) return (int, X) else (None,None)."""
    if isinstance(node, ast.BinOp) and isinstance(node.op, ast.Mult):
        c = const_of(node.left)
        if c is not None:
            return c, node.right
    return None, None

def strip_outer(node):
    """Peel (c)*(E), (c)*(E)+(c')*(E), (E)*(E) squares.  Return (mult, core_node, is_square)."""
    # square: E*E with identical dumps
    if isinstance(node, ast.BinOp) and isinstance(node.op, ast.Mult):
        ld, rd = ast.dump(node.left), ast.dump(node.right)
        if ld == rd:
            return 1, node.left, True
        c = const_of(node.left)
        if c is not None:
            m, core, sq = strip_outer(node.right)
            return c*m, core, sq
    if isinstance(node, ast.BinOp) and isinstance(node.op, ast.Add):
        c1, a1 = split_scalar(node.left); c2, a2 = split_scalar(node.right)
        if c1 is not None and c2 is not None and ast.dump(a1) == ast.dump(a2):
            m, core, sq = strip_outer(a1)
            return (c1+c2)*m, core, sq
    return 1, node, False

def build(cache='atoms.pkl'):
    raw = H.load_raw()
    trees = []
    atom_ids = {}
    atom_src = []
    eq_terms = []   # per eq: (outer_mult, is_square, [(coeff, atom_id)])
    for i, s in enumerate(raw):
        t = ast.parse(s, mode='eval').body
        m, core, sq = strip_outer(t)
        terms = unfold(core)
        tl = []
        for c, a in terms:
            key = ast.unparse(a)
            aid = atom_ids.get(key)
            if aid is None:
                aid = len(atom_src); atom_ids[key] = aid; atom_src.append(key)
            tl.append((c, aid))
        eq_terms.append((m, sq, tl))
    return atom_src, eq_terms

if __name__ == '__main__':
    t0 = time.time()
    atom_src, eq_terms = build()
    print(f'atoms={len(atom_src)} eqs={len(eq_terms)} in {time.time()-t0:.1f}s')
    # compile atoms
    VAR_RE = re.compile(r'x_(\d+)')
    codes = [compile(VAR_RE.sub(r'v[\1]', s), '<a>', 'eval') for s in atom_src]
    pickle.dump({'atom_src': atom_src, 'eq_terms': eq_terms}, open('atoms.pkl','wb'))
    # validate against raw equations at the partial
    v = H.load_assignment('../best/new_instance_partial_39022.json')
    ns = {'v': v, '__builtins__': {}}
    av = [eval(c, ns) for c in codes]
    eqc, _ = H.load_equations()
    bad = 0
    for i,(m,sq,tl) in enumerate(eq_terms):
        s = sum(c*av[a] for c,a in tl)
        val = m*(s*s if sq else s)
        if val != H.resid(eqc, v, i):
            bad += 1
            if bad < 5: print('MISMATCH eq', i)
    print('mismatched equations:', bad)
    nz = [i for i,x in enumerate(av) if x != 0]
    print('nonzero atoms at partial:', len(nz), nz[:20])
