"""Agent D core library: loads my own caches, exact integer atom/equation evaluation."""
import pickle, collections, os, sys, json, heapq

HERE = os.path.dirname(os.path.abspath(__file__))
LAB = os.path.dirname(HERE)
CACHE = os.path.join(HERE, 'cache')
NVARS = 38748
P = 2**256 - 2**32 - 977

_a = pickle.load(open(os.path.join(CACHE, 'atoms.pkl'), 'rb'))
atom_src = _a['atom_src']
eq_terms = _a['eq_terms']
polys = pickle.load(open(os.path.join(CACHE, 'polys.pkl'), 'rb'))
_g = pickle.load(open(os.path.join(CACHE, 'gates.pkl'), 'rb'))
definer, atom_out, free = _g['definer'], _g['atom_out'], _g['free']
topo = pickle.load(open(os.path.join(CACHE, 'topo.pkl'), 'rb'))['topo']

NA = len(polys)
NEQ = len(eq_terms)
avars = [set(v for m in Pp for v in m) for Pp in polys]
var_atoms = collections.defaultdict(list)
for _i, _s in enumerate(avars):
    for _u in _s:
        var_atoms[_u].append(_i)
tidx = {v: i for i, v in enumerate(topo)}
BIG = len(topo) + 10
freeset = set(free)

atom2eq = collections.defaultdict(dict)
eq_atoms = []
for _i, (_m, _sq, _tl) in enumerate(eq_terms):
    co = collections.defaultdict(int)
    for _c, _a in _tl:
        co[_a] += _c
    co = {k: c for k, c in co.items() if c}
    eq_atoms.append((_m, _sq, co))
    for _a, _c in co.items():
        atom2eq[_a][_i] = _c
atom2eq = dict(atom2eq)
checks = [a for a in range(NA) if a not in atom_out]


def evalpoly(Pp, v):
    s = 0
    for m, c in Pp.items():
        t = c
        for u in m:
            t *= v[u]
        s += t
    return s


def all_atom_values(v):
    return [evalpoly(Pp, v) for Pp in polys]


def eq_value(i, av):
    m, sq, co = eq_atoms[i]
    s = 0
    for a, c in co.items():
        s += c * av[a]
    return m * (s * s if sq else s)


def failing_eqs(av):
    out = []
    for i in range(NEQ):
        m, sq, co = eq_atoms[i]
        s = 0
        for a, c in co.items():
            s += c * av[a]
        if s:
            out.append(i)
    return out


def score(v):
    av = all_atom_values(v)
    return NEQ - len(failing_eqs(av)), av


def solve_for(a, t, v):
    Pp = polys[a]
    c = 0
    for m, cc in Pp.items():
        if len(m) == 1 and m[0] == t:
            c += cc
        elif t in m:
            return None
    if c == 0:
        return None
    old = v[t]
    v[t] = 0
    rest = evalpoly(Pp, v)
    v[t] = old
    if rest % c:
        return None
    return -rest // c


def ripple(v, seeds, maxsteps=500000, block=()):
    changed = {}
    pq = []
    block = set(block)
    for u, val in seeds.items():
        if v[u] == val:
            continue
        v[u] = val
        changed[u] = val
        heapq.heappush(pq, (tidx.get(u, BIG), u))
    steps = 0
    seedset = set(seeds)
    while pq:
        _, u = heapq.heappop(pq)
        steps += 1
        if steps > maxsteps:
            break
        for a in var_atoms[u]:
            if a not in atom_out or a in block:
                continue
            c, t = atom_out[a]
            if t == u or t in seedset:
                continue
            nv = solve_for(a, t, v)
            if nv is None or nv == v[t]:
                continue
            v[t] = nv
            changed[t] = nv
            heapq.heappush(pq, (tidx.get(t, BIG), t))
    return changed, steps


def load(path):
    with open(path) as f:
        d = json.load(f)
    v = [0] * NVARS
    for k, val in d.items():
        idx = int(k[2:]) if k.startswith('x_') else int(k)
        v[idx] = int(val)
    return v


def save(v, path):
    json.dump({f'x_{i}': v[i] for i in range(NVARS) if v[i] != 0}, open(path, 'w'))
