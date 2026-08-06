"""Recover the setter's canonical gate orientation: out = first bare-variable top-level term."""
import ast, pickle, collections, re

NVARS = 38748

def top_terms(node, s=1, out=None):
    """Flatten top-level +/- chain into (sign, node) list."""
    if out is None: out = []
    if isinstance(node, ast.BinOp) and isinstance(node.op, (ast.Add, ast.Sub)):
        top_terms(node.left, s, out)
        top_terms(node.right, s*(1 if isinstance(node.op, ast.Add) else -1), out)
    elif isinstance(node, ast.UnaryOp) and isinstance(node.op, ast.USub):
        top_terms(node.operand, -s, out)
    else:
        out.append((s, node))
    return out

def bare_var(node):
    """If node is c*x_i or x_i (no other variable), return (coeff, vid) else None."""
    if isinstance(node, ast.Name): return (1, int(node.id[2:]))
    if isinstance(node, ast.BinOp) and isinstance(node.op, ast.Mult):
        for a, b in ((node.left, node.right), (node.right, node.left)):
            if isinstance(b, ast.Name):
                try: c = ast.literal_eval(a)
                except Exception: continue
                if isinstance(c, int): return (c, int(b.id[2:]))
    return None

def outputs_of(src):
    """Ordered list of (coeff, vid) candidate outputs for an atom source string."""
    node = ast.parse(src, mode='eval').body
    res = []
    for s, t in top_terms(node):
        bv = bare_var(t)
        if bv is not None: res.append((s*bv[0], bv[1]))
    return res

if __name__ == '__main__':
    d = pickle.load(open('atoms.pkl','rb')); src = d['atom_src']
    polys = pickle.load(open('polys.pkl','rb'))
    cand = []
    for i, s in enumerate(src):
        if max((len(m) for m in polys[i]), default=0) >= 3:
            cand.append([])          # squares are checks
            continue
        higher = set(v for m in polys[i] if len(m) >= 2 for v in m)
        cand.append([(c, v) for c, v in outputs_of(s) if v not in higher])
    # assign definers: shortest atoms first, prefer |coeff|=1
    order_atoms = sorted(range(len(src)), key=lambda a: (len(polys[a]), len(src[a])))
    definer = {}; atom_out = {}
    for a in order_atoms:
        for c, v in cand[a]:
            if v not in definer:
                definer[v] = a; atom_out[a] = (c, v); break
    print(f'variables with a definer: {len(definer)}')
    print(f'atoms used as gates: {len(atom_out)}   checks: {len(src)-len(atom_out)}')
    free = [v for v in range(NVARS) if v not in definer]
    print(f'free inputs: {len(free)}')
    for nm, a in [('C1',22229),('C2',22231)]:
        print(nm, 'candidates', cand[a], '-> assigned', atom_out.get(a))
    pickle.dump({'cand':cand,'definer':definer,'atom_out':atom_out,'free':free}, open('gates.pkl','wb'))
