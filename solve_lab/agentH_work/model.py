"""Agent H independent model: parse EQUATIONS.txt into atoms + equations, build gate DAG."""
import ast, re, json, pickle, sys, time, os

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
EQ_PATH = os.path.join(ROOT, 'EQUATIONS.txt')
NVARS = 38748
VAR_RE = re.compile(r'x_(\d+)')

def load_raw():
    out = []
    with open(EQ_PATH) as f:
        for line in f:
            line = line.strip()
            if not line: continue
            out.append(line.rsplit('=',1)[0])
    return out

def const_of(node):
    if isinstance(node, ast.Constant) and isinstance(node.value, int): return node.value
    if isinstance(node, ast.UnaryOp) and isinstance(node.op, ast.USub):
        c = const_of(node.operand); return None if c is None else -c
    return None

def split_scalar(node):
    if isinstance(node, ast.BinOp) and isinstance(node.op, ast.Mult):
        c = const_of(node.left)
        if c is not None: return c, node.right
    return None, None

def strip_outer(node):
    if isinstance(node, ast.BinOp) and isinstance(node.op, ast.Mult):
        ld, rd = ast.dump(node.left), ast.dump(node.right)
        if ld == rd: return 1, node.left, True
        c = const_of(node.left)
        if c is not None:
            m, core, sq = strip_outer(node.right); return c*m, core, sq
    if isinstance(node, ast.BinOp) and isinstance(node.op, ast.Add):
        c1,a1 = split_scalar(node.left); c2,a2 = split_scalar(node.right)
        if c1 is not None and c2 is not None and ast.dump(a1)==ast.dump(a2):
            m, core, sq = strip_outer(a1); return (c1+c2)*m, core, sq
    return 1, node, False

def unfold(node):
    out = []
    while isinstance(node, ast.BinOp) and isinstance(node.op, ast.Add):
        c,a = split_scalar(node.right)
        if c is None: break
        out.append((c,a)); node = node.left
    out.append((1,node)); out.reverse()
    return out

def build():
    raw = load_raw()
    atom_ids = {}; atom_src = []; eq_terms = []
    for s in raw:
        t = ast.parse(s, mode='eval').body
        m, core, sq = strip_outer(t)
        tl = []
        for c,a in unfold(core):
            key = ast.unparse(a)
            aid = atom_ids.get(key)
            if aid is None:
                aid = len(atom_src); atom_ids[key] = aid; atom_src.append(key)
            tl.append((c,aid))
        eq_terms.append((m, sq, tl))
    return atom_src, eq_terms

CACHE = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'model.pkl')

def get():
    if os.path.exists(CACHE):
        return pickle.load(open(CACHE,'rb'))
    atom_src, eq_terms = build()
    atom_vars = [tuple(sorted(set(int(m) for m in VAR_RE.findall(s)))) for s in atom_src]
    d = {'atom_src': atom_src, 'eq_terms': eq_terms, 'atom_vars': atom_vars}
    pickle.dump(d, open(CACHE,'wb'))
    return d

if __name__ == '__main__':
    t0=time.time(); d = get()
    print('atoms', len(d['atom_src']), 'eqs', len(d['eq_terms']), 'in %.1fs'%(time.time()-t0))
    from collections import Counter
    print('eq term-count histogram:', Counter(len(t[2]) for t in d['eq_terms']).most_common(12))
    print('squares:', sum(1 for t in d['eq_terms'] if t[1]))
    print('atom var-count hist:', Counter(len(v) for v in d['atom_vars']).most_common(12))
