"""Agent D: build my own atom/gate/poly caches into agentD_work/cache/.

Independent rebuild (does not touch s9/ or any shared file).
"""
import ast, re, json, pickle, sys, time, os, collections, heapq

HERE = os.path.dirname(os.path.abspath(__file__))
LAB = os.path.dirname(HERE)
EQ_PATH = os.path.join(LAB, '..', 'EQUATIONS.txt')
CACHE = os.path.join(HERE, 'cache')
os.makedirs(CACHE, exist_ok=True)
NVARS = 38748
VAR_RE = re.compile(r'x_(\d+)')


def load_raw():
    out = []
    with open(EQ_PATH) as f:
        for line in f:
            line = line.strip()
            if line:
                out.append(line.rsplit('=', 1)[0])
    return out


# ---------- atomize ----------
def const_of(node):
    if isinstance(node, ast.Constant) and isinstance(node.value, int):
        return node.value
    if isinstance(node, ast.UnaryOp) and isinstance(node.op, ast.USub):
        c = const_of(node.operand)
        return None if c is None else -c
    return None


def split_scalar(node):
    if isinstance(node, ast.BinOp) and isinstance(node.op, ast.Mult):
        c = const_of(node.left)
        if c is not None:
            return c, node.right
    return None, None


def strip_outer(node):
    if isinstance(node, ast.BinOp) and isinstance(node.op, ast.Mult):
        ld, rd = ast.dump(node.left), ast.dump(node.right)
        if ld == rd:
            return 1, node.left, True
        c = const_of(node.left)
        if c is not None:
            m, core, sq = strip_outer(node.right)
            return c * m, core, sq
    if isinstance(node, ast.BinOp) and isinstance(node.op, ast.Add):
        c1, a1 = split_scalar(node.left)
        c2, a2 = split_scalar(node.right)
        if c1 is not None and c2 is not None and ast.dump(a1) == ast.dump(a2):
            m, core, sq = strip_outer(a1)
            return (c1 + c2) * m, core, sq
    return 1, node, False


def unfold(node):
    out = []
    while isinstance(node, ast.BinOp) and isinstance(node.op, ast.Add):
        c, a = split_scalar(node.right)
        if c is None:
            break
        out.append((c, a))
        node = node.left
    out.append((1, node))
    out.reverse()
    return out


# ---------- poly ----------
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
        if out[m] == 0:
            del out[m]
    return out


def topoly(node):
    if isinstance(node, ast.Constant):
        return {(): node.value} if node.value else {}
    if isinstance(node, ast.Name):
        return {(int(node.id[2:]),): 1}
    if isinstance(node, ast.UnaryOp):
        p = topoly(node.operand)
        return {m: -c for m, c in p.items()} if isinstance(node.op, ast.USub) else p
    if isinstance(node, ast.BinOp):
        l = topoly(node.left)
        r = topoly(node.right)
        if isinstance(node.op, ast.Add):
            return padd(l, r, 1)
        if isinstance(node.op, ast.Sub):
            return padd(l, r, -1)
        if isinstance(node.op, ast.Mult):
            return pmul(l, r)
    raise ValueError(ast.dump(node))


# ---------- gates ----------
def top_terms(node, s=1, out=None):
    if out is None:
        out = []
    if isinstance(node, ast.BinOp) and isinstance(node.op, (ast.Add, ast.Sub)):
        top_terms(node.left, s, out)
        top_terms(node.right, s * (1 if isinstance(node.op, ast.Add) else -1), out)
    elif isinstance(node, ast.UnaryOp) and isinstance(node.op, ast.USub):
        top_terms(node.operand, -s, out)
    else:
        out.append((s, node))
    return out


def bare_var(node):
    if isinstance(node, ast.Name):
        return (1, int(node.id[2:]))
    if isinstance(node, ast.BinOp) and isinstance(node.op, ast.Mult):
        for a, b in ((node.left, node.right), (node.right, node.left)):
            if isinstance(b, ast.Name):
                try:
                    c = ast.literal_eval(a)
                except Exception:
                    continue
                if isinstance(c, int):
                    return (c, int(b.id[2:]))
    return None


def outputs_of(node):
    res = []
    for s, t in top_terms(node):
        bv = bare_var(t)
        if bv is not None:
            res.append((s * bv[0], bv[1]))
    return res


def main():
    t0 = time.time()
    raw = load_raw()
    print(f'raw equations {len(raw)}  ({time.time()-t0:.1f}s)')
    atom_ids = {}
    atom_src = []
    atom_node = []
    eq_terms = []
    for i, s in enumerate(raw):
        t = ast.parse(s, mode='eval').body
        m, core, sq = strip_outer(t)
        tl = []
        for c, a in unfold(core):
            key = ast.unparse(a)
            aid = atom_ids.get(key)
            if aid is None:
                aid = len(atom_src)
                atom_ids[key] = aid
                atom_src.append(key)
                atom_node.append(a)
            tl.append((c, aid))
        eq_terms.append((m, sq, tl))
    print(f'atoms={len(atom_src)} eqs={len(eq_terms)}  ({time.time()-t0:.1f}s)')

    polys = [topoly(n) for n in atom_node]
    print(f'polys expanded ({time.time()-t0:.1f}s)')
    degh = collections.Counter(max((len(m) for m in p), default=0) for p in polys)
    print('degree hist:', dict(sorted(degh.items())))

    cand = []
    for i in range(len(atom_src)):
        if max((len(m) for m in polys[i]), default=0) >= 3:
            cand.append([])
            continue
        higher = set(v for m in polys[i] if len(m) >= 2 for v in m)
        cand.append([(c, v) for c, v in outputs_of(atom_node[i]) if v not in higher])
    order_atoms = sorted(range(len(atom_src)), key=lambda a: (len(polys[a]), len(atom_src[a])))
    definer = {}
    atom_out = {}
    for a in order_atoms:
        for c, v in cand[a]:
            if v not in definer:
                definer[v] = a
                atom_out[a] = (c, v)
                break
    free = [v for v in range(NVARS) if v not in definer]
    print(f'definers={len(definer)} gates={len(atom_out)} checks={len(atom_src)-len(atom_out)} free={len(free)}')

    avars = [set(v for m in P for v in m) for P in polys]
    indeg = {}
    users = collections.defaultdict(list)
    for v, a in definer.items():
        deps = [u for u in avars[a] if u != v and u in definer]
        indeg[v] = len(deps)
        for u in deps:
            users[u].append(v)
    q = collections.deque([v for v in definer if indeg[v] == 0])
    topo = []
    while q:
        v = q.popleft()
        topo.append(v)
        for w in users[v]:
            indeg[w] -= 1
            if indeg[w] == 0:
                q.append(w)
    print(f'topo covers {len(topo)}/{len(definer)} ({len(definer)-len(topo)} in cycles)')

    pickle.dump({'atom_src': atom_src, 'eq_terms': eq_terms}, open(os.path.join(CACHE, 'atoms.pkl'), 'wb'))
    pickle.dump(polys, open(os.path.join(CACHE, 'polys.pkl'), 'wb'))
    pickle.dump({'cand': cand, 'definer': definer, 'atom_out': atom_out, 'free': free},
                open(os.path.join(CACHE, 'gates.pkl'), 'wb'))
    pickle.dump({'topo': topo}, open(os.path.join(CACHE, 'topo.pkl'), 'wb'))
    print(f'done {time.time()-t0:.1f}s')


if __name__ == '__main__':
    main()
