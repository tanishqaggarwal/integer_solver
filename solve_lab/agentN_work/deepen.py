"""Deepen the frame: more knobs at the SAME state.

`Frame(POOL)` makes every pool variable free, so `make(D)` changes only the row targets b, never
the knob set -- which is why the whole detach lattice collapsed to 16 states.  The untested lever is
the frame's DEPTH: `pool.py` stops two levels above the region atoms.  Detaching further up adds
free inputs without changing the witness assignment (the extra variables are re-attached to their
gate values, which equal their witness values), so it is a strict enlargement of the knob set at a
fixed state.

Re-orientation is subsumed by this: a check atom `x_v - rest` can be turned into a definition of
x_v only if x_v is free or freeable, and once x_v is free we may simply CHOOSE the value that zeroes
that atom.  Making a variable free is exactly detaching it.  So "which atom defines x_v" is a value
choice inside a deep enough frame, and depth is the thing that actually varies.

For each depth: rebuild the frame, rebuild the region, take every free input that touches it, and
compute the exact maximum number of integrally zeroable region rows.
"""
import json, time, sys, os, ast, re
from collections import defaultdict
import ev, model
import frameB as FB
from frameB import Frame, State
import zsolve

HERE = os.path.dirname(os.path.abspath(__file__))
VAR_RE = re.compile(r'x_(\d+)')
d = model.get()
atom_src = d['atom_src']
atom_vars = d['atom_vars']
eq_terms = d['eq_terms']
NV = 38748
FREE0 = set(ev.F['free0'])
definer = ev.F['definer']

atom_eqs = defaultdict(set)
for i, (m, sq, tl) in enumerate(eq_terms):
    for c, a in tl:
        atom_eqs[a].add(i)

REGION_ATOMS = [22229, 22230, 22231, 35758, 35759, 35760, 35761, 35762, 37887]
WIT = [642, 28730, 29854, 31864]

W = json.load(open(os.path.join(HERE, '..', 'best', 'new_instance_partial_39026.json')))
wv = [0] * NV
for k, val in W.items():
    wv[int(k[2:]) if k.startswith('x_') else int(k)] = int(val)


def build_pool(depth):
    """variables of the region atoms, then `depth-1` further levels of their definers' inputs"""
    pool = set()
    for a in REGION_ATOMS:
        for v in atom_vars[a]:
            if v not in FREE0:
                pool.add(v)
    frontier = set(pool)
    for _ in range(depth - 1):
        nxt = set()
        for v in frontier:
            a = definer[v]
            if a >= 0:
                for u in atom_vars[a]:
                    if u != v and u not in FREE0 and u not in pool:
                        nxt.add(u)
        pool |= nxt
        frontier = nxt
        if not nxt:
            break
    for v in WIT:
        pool.add(v)
    return sorted(pool)


def analyse(depth):
    t0 = time.time()
    POOL = build_pool(depth)
    fr = Frame(POOL)
    FREE = set(fr.free)
    FR0 = fr.free
    BASEFV = {u: wv[u] for u in fr.free if wv[u] != 0}
    st = State(fr, BASEFV)
    sc = st.score()
    NZ = [q for q in st.av if st.av[q] != 0]
    R = set()
    for q in NZ:
        R |= atom_eqs[q]
    Rl = sorted(R)
    tbuild = time.time() - t0

    def bits(x):
        o = []
        while x:
            q = x & -x
            o.append(q.bit_length() - 1)
            x ^= q
        return o

    atoms_R = set()
    for e in Rl:
        for c, a in eq_terms[e][2]:
            atoms_R.add(a)
    cands = set()
    for q in atoms_R:
        if q in fr.csup:
            cands.update(FR0[b] for b in bits(fr.csup[q]))
    cands = sorted(y for y in cands if y in FREE)

    def inner(s, e):
        m, sq, tl = eq_terms[e]
        t = 0
        for c, a in tl:
            x = s.av.get(a)
            if x:
                t += c * x
        return t

    Rset = set(Rl)
    b = [inner(st, e) for e in Rl]
    knobs, cols, outside_touch = [], [], set()
    for Y in cands:
        h = st.clone().set_free({Y: st.fv.get(Y, 0) + 1})
        moved = [q for q in h.av if h.av[q] != st.av[q]]
        if not moved:
            continue
        eqs = set()
        for q in moved:
            eqs |= atom_eqs[q]
        if not (eqs & Rset):
            continue
        knobs.append(Y)
        cols.append([inner(h, e) - b[i] for i, e in enumerate(Rl)])
        outside_touch |= (eqs - Rset)
    k = len(knobs)
    M = [[cols[j][i] for j in range(k)] for i in range(len(Rl))]
    # zero-collateral subset of the knobs, for comparability with the shallow frame
    narrow = []
    for j, Y in enumerate(knobs):
        h = st.clone().set_free({Y: st.fv.get(Y, 0) + 1})
        eqs = set()
        for q in [q for q in h.av if h.av[q] != st.av[q]]:
            eqs |= atom_eqs[q]
        if eqs <= Rset:
            narrow.append(j)
    opt_w, rows_w, exh_w, _ = zsolve.max_zero_rows(M, b, k, len(Rl))
    Mn = [[M[i][j] for j in narrow] for i in range(len(Rl))]
    opt_n, rows_n, exh_n, _ = zsolve.max_zero_rows(Mn, b, len(narrow), len(Rl)) if narrow else (0, [], True, 0)
    print('depth %d: pool=%-5d frame free=%-5d score=%d  |R|=%-3d nonzero atoms=%-3d '
          'knobs wide=%-3d narrow=%-3d  OPT wide=%-3d narrow=%-3d  failing=%d  (%.0fs)'
          % (depth, len(POOL), len(fr.free), sc, len(Rl), len(NZ), k, len(narrow),
             opt_w, opt_n, len(Rl) - opt_n, time.time() - t0), flush=True)
    return dict(depth=depth, pool=len(POOL), free=len(fr.free), score=sc, R=len(Rl),
                nz=len(NZ), kwide=k, knarrow=len(narrow), opt_wide=opt_w, opt_narrow=opt_n,
                exh_wide=str(exh_w), exh_narrow=str(exh_n), outside=len(outside_touch),
                best_failing=len(Rl) - opt_n)


if __name__ == '__main__':
    dmax = int(sys.argv[1]) if len(sys.argv) > 1 else 5
    res = []
    for depth in range(2, dmax + 1):
        try:
            res.append(analyse(depth))
        except Exception as e:
            print('depth %d ERROR %s' % (depth, e), flush=True)
            break
        json.dump(res, open(os.path.join(HERE, 'runs', 'deepen.json'), 'w'), indent=1)
