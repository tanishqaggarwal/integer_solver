"""Expand every atom into a sparse polynomial: dict monomial(tuple sorted vids) -> int coeff."""
import ast, pickle, time, sys, collections

def pmul(a, b):
    out = collections.defaultdict(int)
    for m1, c1 in a.items():
        for m2, c2 in b.items():
            out[tuple(sorted(m1 + m2))] += c1 * c2
    return {m: c for m, c in out.items() if c}

def padd(a, b, s=1):
    out = dict(a)
    for m, c in b.items():
        out[m] = out.get(m, 0) + s * c
        if out[m] == 0: del out[m]
    return out

def topoly(node):
    if isinstance(node, ast.Constant): return {(): node.value} if node.value else {}
    if isinstance(node, ast.Name): return {(int(node.id[2:]),): 1}
    if isinstance(node, ast.UnaryOp):
        p = topoly(node.operand)
        return {m: -c for m, c in p.items()} if isinstance(node.op, ast.USub) else p
    if isinstance(node, ast.BinOp):
        l = topoly(node.left); r = topoly(node.right)
        if isinstance(node.op, ast.Add): return padd(l, r, 1)
        if isinstance(node.op, ast.Sub): return padd(l, r, -1)
        if isinstance(node.op, ast.Mult): return pmul(l, r)
    raise ValueError(ast.dump(node))

def degree(p):
    return max((len(m) for m in p), default=0)

if __name__ == '__main__':
    d = pickle.load(open('atoms.pkl','rb')); src = d['atom_src']
    t0 = time.time(); polys = []
    for i, s in enumerate(src):
        polys.append(topoly(ast.parse(s, mode='eval').body))
    print(f'expanded {len(polys)} atoms in {time.time()-t0:.1f}s')
    deg = collections.Counter(degree(p) for p in polys)
    print('degree histogram:', dict(sorted(deg.items())))
    pickle.dump(polys, open('polys.pkl','wb'))
