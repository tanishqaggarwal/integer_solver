"""Replace square check atoms E^2 by their degree-2 root E; rebuild residual model."""
import ast, pickle, collections
import poly as PY

d = pickle.load(open('atoms.pkl','rb')); src = d['atom_src']
polys = pickle.load(open('polys.pkl','rb'))
g = pickle.load(open('gates.pkl','rb')); atom_out = g['atom_out']

def root_of(s):
    """If source is (E)*(E) return E's source, else None."""
    n = ast.parse(s, mode='eval').body
    if isinstance(n, ast.BinOp) and isinstance(n.op, ast.Mult):
        if ast.dump(n.left) == ast.dump(n.right):
            return ast.unparse(n.left)
    return None

if __name__ == '__main__':
    roots = {}
    nsq = 0
    for a, s in enumerate(src):
        if a in atom_out: continue
        r = root_of(s)
        if r is not None:
            roots[a] = PY.topoly(ast.parse(r, mode='eval').body); nsq += 1
    print(f'square check atoms with extracted root: {nsq}')
    deg = collections.Counter(max((len(m) for m in P), default=0) for P in roots.values())
    print('root degree histogram:', dict(deg))
    # any remaining degree>=3 checks that are NOT clean squares?
    left = [a for a,P in enumerate(polys) if a not in atom_out and a not in roots
            and max((len(m) for m in P), default=0) >= 3]
    print(f'non-square high-degree checks left: {len(left)} {left[:10]}')
    pickle.dump(roots, open('roots.pkl','wb'))
