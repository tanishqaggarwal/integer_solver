"""S11 step 61: frame-space search that actually USES the new freedom.

The first attempt detached variables and re-evaluated -- and every frame scored
exactly 39,026, because a freed variable keeps a value that already satisfies its
atom.  A frame only matters if the freedom it creates is exercised.  So after
building the frame, run a short enriched repair (two-level handle moves) using the
newly free parameters, and record the best score reached.

Chunked/resumable:  fsearch2.py START END [SEED]
"""
import os, sys, random
HERE = '/home/user/integer_solver/solve_lab/s10'
LAB = '/home/user/integer_solver/solve_lab'
sys.path.insert(0, os.path.join(LAB, 's9', 'eff'))
sys.path.insert(0, HERE)
import lib as L, tools as T, ad
from chunk import sweep, load

P = ad.P
SEVEN = [22229, 22230, 35758, 35759, 35760, 35761, 35762]
BASE_DET = {7068: 22229, 28730: 22230, 29854: 35758, 31864: 35761, 642: 35762}
FORBID = {2081, 4287}
base = L.load(os.path.join(LAB, 'best', 'new_instance_partial_39026.json'))

cone, st, seen = set(), list(SEVEN), set()
while st:
    a = st.pop()
    if a in seen:
        continue
    seen.add(a)
    for w in L.avars[a]:
        cone.add(w)
        for b in L.var_atoms[w]:
            if b not in seen and len(seen) < 250:
                st.append(b)
CAND = sorted(w for w in cone if w in L.definer)


def evaluate(extra):
    DET = dict(BASE_DET)
    for w in extra:
        a = L.definer.get(w)
        if a is not None:
            DET[w] = a
    definer = {t: a for t, a in L.definer.items() if t not in DET}
    ORDER = [t for t in ad.ORDER if t not in DET]
    FREE = set(t for t in range(L.NVARS) if t not in definer)

    def fwd(v, r=4):
        for _ in range(r):
            for u in ORDER:
                nv = T.solve_lin(definer[u], u, v)
                if nv is not None:
                    v[u] = nv
        return v

    def pot(v):
        av = L.all_atom_values(v)
        nz = [a for a in range(L.NA) if av[a]]
        return (L.NEQ - len(L.failing_eqs(av)), -len(nz)), av, nz

    v = list(base)
    fwd(v)
    cur, av, nz = pot(v)
    best = cur[0]
    for _ in range(4):
        got = None
        for a in nz:
            cands = []
            for w in sorted(set(L.avars[a])):
                if w in FORBID:
                    continue
                tgt = T.solve_lin(a, w, v)
                if tgt is None or tgt == v[w]:
                    continue
                if w in FREE:
                    cands.append((w, tgt))
                else:
                    d = definer.get(w)
                    if d is None:
                        continue
                    vv = list(v)
                    vv[w] = tgt
                    for u in sorted(set(L.avars[d])):
                        if u == w or u not in FREE or u in FORBID:
                            continue
                        nv2 = T.solve_lin(d, u, vv)
                        if nv2 is not None:
                            cands.append((u, nv2))
            for u, nv2 in cands[:6]:
                tr = list(v)
                tr[u] = nv2
                fwd(tr)
                p2, av2, nz2 = pot(tr)
                if p2 > cur:
                    got = (p2, tr, av2, nz2)
                    break
            if got:
                break
        if not got:
            break
        cur, v, av, nz = got
        best = max(best, cur[0])
    if best > 39026:
        T.save(v, os.path.join(HERE, 'FS2_%d.json' % best))
    return {'score': best, 'nz': len(nz), 'det': list(extra)}


start, end = int(sys.argv[1]), int(sys.argv[2])
seed = int(sys.argv[3]) if len(sys.argv) > 3 else 0
random.seed(seed * 100003 + start)
specs = [tuple(sorted(random.sample(CAND, random.randint(1, 4))))
         for _ in range(end - start)]
sweep('fs2_s%d' % seed, specs, evaluate, 0, len(specs),
      keyfn=lambda s: ','.join(map(str, s)), budget=480)
rs = load('fs2_s%d' % seed)
print('seed %d: %d frames, best %d'
      % (seed, len(rs), max((r['score'] for r in rs), default=0)))
