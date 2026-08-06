"""Effective-footprint model:  sacrifice set S -> confined atoms -> EFFECTIVE knobs
(gate-preserving ripple, with the confined atoms blocked from repair) -> exact integer
lattice -> max number of equations in S that can be simultaneously zeroed.

Square atoms are represented by their degree-2 ROOT (roots.pkl).
"""
import pickle, random, sys, itertools, collections
import lib as L
from isolve import solve_int2 as solve_int

roots = pickle.load(open('roots.pkl', 'rb'))
prim = [roots[a] if a in roots else L.polys[a] for a in range(L.NA)]
IS_SQ = [a in roots for a in range(L.NA)]
prim_vars = [set(u for m in Pp for u in m) for Pp in prim]
prim_var_atoms = collections.defaultdict(list)
for _a, _s in enumerate(prim_vars):
    for _u in _s:
        prim_var_atoms[_u].append(_a)

BASEP = None
_C2 = None


def prim_val(a, v):
    return L.evalpoly(prim[a], v)


def load_census(which='24'):
    global _C2
    if _C2 is None:
        _C2 = pickle.load(open(f'census2_{which}.pkl', 'rb'))
    return _C2


def confined_atoms(S):
    S = set(S)
    return sorted(a for a in L.atom2eq if set(L.atom2eq[a]) <= S)


def move(v0, seeds, block):
    """Apply seeds and the blocked ripple to a COPY; return (v, touched_prim_atoms)."""
    v = list(v0)
    ch, st = L.ripple(v, seeds, block=block)
    cand = set()
    for u in ch:
        cand.update(prim_var_atoms[u])
    tou = {a: prim_val(a, v) for a in cand}
    tou = {a: x for a, x in tou.items() if x != BASEP[a]}
    return v, tou


def build(S, v0, which='24', verbose=True, extra_cand=()):
    """Build the affine model.  Returns dict or None."""
    global BASEP
    if BASEP is None:
        BASEP = [prim_val(a, v0) for a in range(L.NA)]
    C2 = load_census(which)
    rev = C2['rev']
    A = confined_atoms(S)
    Aset = set(A)
    cand = set(extra_cand)
    for a in A:
        cand |= rev.get(a, set())
    knobs = []
    M = {a: {} for a in A}
    nonlin = []
    for x in sorted(cand):
        v1, t1 = move(v0, {x: v0[x] + 1}, Aset)
        if not t1 or not set(t1) <= Aset:
            continue
        v5, t5 = move(v0, {x: v0[x] + 5}, Aset)
        if not set(t5) <= Aset:
            continue
        ok = True
        col = {}
        for a in Aset:
            d1 = t1.get(a, BASEP[a]) - BASEP[a]
            d5 = t5.get(a, BASEP[a]) - BASEP[a]
            if d5 != 5 * d1:
                ok = False
                break
            if d1:
                col[a] = d1
        if not ok:
            nonlin.append(x)
            continue
        knobs.append(x)
        for a, c in col.items():
            M[a][x] = c
    # joint-linearity verification
    random.seed(11)
    joint_ok = True
    for _ in range(3):
        dv = {x: random.randint(-10**5, 10**5) for x in knobs}
        vj, tj = move(v0, {x: v0[x] + dv[x] for x in knobs}, Aset)
        if not set(tj) <= Aset:
            joint_ok = False
            break
        for a in A:
            pred = BASEP[a] + sum(M[a].get(x, 0) * dv[x] for x in knobs)
            if tj.get(a, BASEP[a]) != pred:
                joint_ok = False
                break
        if not joint_ok:
            break
    if verbose and not joint_ok:
        print('  WARNING: joint linearity failed')
    conds = []
    for i in sorted(S):
        m, sq, co = L.eq_atoms[i]
        lin = [(a, c) for a, c in co.items() if not IS_SQ[a]]
        sqs = [(a, c) for a, c in co.items() if IS_SQ[a]]
        rows = []

        def linrow(terms):
            row = {x: 0 for x in knobs}
            rhs = 0
            for a, c in terms:
                rhs -= c * BASEP[a]
                for x, mc in M.get(a, {}).items():
                    row[x] += c * mc
            return row, rhs
        if not sqs:
            rows.append(linrow(lin))
        elif not lin and len(sqs) == 1:
            a, c = sqs[0]
            rows.append(({x: M.get(a, {}).get(x, 0) for x in knobs}, -BASEP[a]))
        else:
            for a, c in sqs:
                rows.append(({x: M.get(a, {}).get(x, 0) for x in knobs}, -BASEP[a]))
            rows.append(linrow(lin))
        conds.append((i, rows))
    return {'knobs': knobs, 'A': A, 'M': M, 'conds': conds, 'nonlin': nonlin,
            'joint_ok': joint_ok, 'S': sorted(S)}


def solvable(mod, T):
    knobs = mod['knobs']
    rows = []
    rhs = []
    for i, rr in mod['conds']:
        if i not in T:
            continue
        for row, b in rr:
            rows.append([row.get(x, 0) for x in knobs])
            rhs.append(b)
    if not rows:
        return {}
    z = solve_int(rows, rhs)
    if z is None:
        return None
    return {x: z[j] for j, x in enumerate(knobs)}


def maximise(mod, S, lo=0, verbose=False):
    """Largest-first search for the biggest simultaneously-zeroable subset of S."""
    S = sorted(S)
    n = len(S)
    for size in range(n, lo - 1, -1):
        for T in itertools.combinations(S, size):
            z = solvable(mod, set(T))
            if z is not None:
                return size, T, z
    return 0, (), {}


def realise(v0, mod, z, verify=True):
    """Apply the delta vector; return (v, ok, failing_eqs)."""
    Aset = set(mod['A'])
    seeds = {x: v0[x] + d for x, d in z.items() if d}
    v, tou = move(v0, seeds, Aset)
    if not verify:
        return v, None, None
    av = L.all_atom_values(v)
    fails = L.failing_eqs(av)
    return v, len(fails), fails
