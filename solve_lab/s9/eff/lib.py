"""Shared library for the effective-footprint knob census (agent `eff`).

Loads the s9 caches read-only, provides:
  * base state loading (39,022 / 39,024 partials)
  * exact atom evaluation
  * gate-preserving ripple (copied from repair.py, extended to report broken gates)
  * effective footprint of a variable move
"""
import pickle, collections, heapq, os, sys, json, time, random

HERE = os.path.dirname(os.path.abspath(__file__))
S9 = os.path.dirname(HERE)
LAB = os.path.dirname(S9)
if S9 not in sys.path:
    sys.path.insert(0, S9)
if HERE not in sys.path:
    sys.path.insert(0, HERE)
os.chdir(S9)                      # caches use relative paths
import harness as H               # noqa: E402

NVARS = 38748
P = 2**256 - 2**32 - 977

_g = pickle.load(open('gates.pkl', 'rb'))
polys = pickle.load(open('polys.pkl', 'rb'))
topo = pickle.load(open('topo.pkl', 'rb'))['topo']
_d = pickle.load(open('atoms.pkl', 'rb'))
atom_src = _d['atom_src']
eq_terms = _d['eq_terms']
definer, atom_out = _g['definer'], _g['atom_out']
NA = len(polys)
NEQ = len(eq_terms)

avars = [set(u for m in Pp for u in m) for Pp in polys]
var_atoms = collections.defaultdict(list)
for _a, _s in enumerate(avars):
    for _u in _s:
        var_atoms[_u].append(_a)
tidx = {v: i for i, v in enumerate(topo)}
BIG = len(topo) + 10

# atom -> [(eq, coeff)] ; eq -> (mult, is_square, {atom: coeff})
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

# variable -> set of equations mentioning it (SYNTACTIC footprint)
var_eqs = collections.defaultdict(set)
for _a, _s in enumerate(avars):
    for _i in atom2eq.get(_a, ()):
        for _u in _s:
            var_eqs[_u].add(_i)


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
    s = sum(c * av[a] for a, c in co.items())
    return m * (s * s if sq else s)


def failing_eqs(av):
    return [i for i in range(NEQ) if eq_value(i, av) != 0]


def solve_for(a, t, v):
    """Value of var t making atom a vanish (t must occur linearly), else None."""
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
    """Forward gate ripple.  `block` = atoms that must NOT be repaired (we control
    their value explicitly).  Returns (changed_vars_dict, n_steps)."""
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


def touched_atoms(v, base_av, changed):
    """Exact set of atoms whose value differs from base_av, given the changed vars."""
    cand = set()
    for u in changed:
        cand.update(var_atoms[u])
    out = {}
    for a in cand:
        nv = evalpoly(polys[a], v)
        if nv != base_av[a]:
            out[a] = nv
    return out


def eqs_of_atoms(atoms):
    s = set()
    for a in atoms:
        s.update(atom2eq.get(a, ()))
    return s


def load(path):
    return H.load_assignment(path)


def save(v, path):
    H.save_assignment(v, path)


BEST24 = os.path.join(LAB, 'best', 'new_instance_partial_39024.json')
BEST22 = os.path.join(LAB, 'best', 'new_instance_partial_39022.json')
