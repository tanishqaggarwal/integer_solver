"""jm step 17: COMPLETE second-move search from each cheapest relaxed state.

A second move can only lower out12 if it changes an atom occurring in one of the
still-failing equations.  So the complete candidate set is
    F = free inputs in the cones of the variables of the atoms of those equations
which is far smaller than 957 and provably exhaustive for one further move.
Sweep F x deltas, measure exactly.

Also reports the linear COMPENSATION kernel of those equations: which atom
vectors would cancel the residual atom inside its own equations.

usage: python3 jm_17_comp.py <state> START END      (chunked, resumable)
       python3 jm_17_comp.py <state> kern           (kernel report only)
states: C1_6418 (13), C2_28730 (11), BOTH_sel01 (14)
"""
import os, sys, json, time, collections
sys.path.insert(0, '/home/user/integer_solver/solve_lab/s10')
import jm_lib as J
import lib as L, tools as T, ad
import jm_05_engine as EN
from jm_14_sel import build, relax

P = J.P
W = J.base_state()
LOG = '/home/user/integer_solver/solve_lab/s10/jm_comp.jsonl'
BEST = '/home/user/integer_solver/solve_lab/s10/jm_best.json'


def state(name):
    if name == 'C1_6418':
        v = list(W); v[6418] += 1000003; J.fwd2(v, 2)
    elif name == 'C2_28730':
        v = list(W); v[28730] += 1000003; J.fwd2(v, 2)
        v[24548] += v[25442] - W[25442]; J.fwd2(v, 2)
    elif name == 'BOTH_sel01':
        v = build(0, 1, {28730: 1000003})
        v[12553] = T.solve_lin(3582, 12553, v); J.fwd2(v, 2)
        v[6418] = T.solve_lin(3580, 6418, v); J.fwd2(v, 2)
    elif name.startswith('AZERO'):
        import jm_19_azero as Z
        s1, s2 = int(name[5]), int(name[6])
        v = Z.azero(build(s1, s2), True)
        v, _o, _s, _n = EN.repair(list(v), verbose=False, maxit=30)
        for _ in range(3):
            v = Z.azero(v, True)
            v, _o, _s, _n = EN.repair(list(v), verbose=False, maxit=20)
    else:
        raise SystemExit('unknown')
    return v


def candidates(v):
    o, s, nz, av = EN.state(v)
    f = set(L.failing_eqs(av)) - J.E12
    atoms = set()
    for e in f:
        atoms |= set(L.eq_atoms[e][2])
    vs = set()
    for a in atoms:
        vs |= set(L.avars[a])
    F = set()
    for u in vs:
        F |= set(x for x in J.cone(u) if x in J.FREESET)
    return sorted(F), sorted(f), sorted(atoms), (o, s, nz)


def done_set():
    s = set()
    if os.path.exists(LOG):
        for ln in open(LOG):
            try:
                r = json.loads(ln)
                s.add((r['state'], r['u'], r['d']))
            except Exception:
                pass
    return s


if __name__ == '__main__':
    name = sys.argv[1]
    v0 = state(name)
    F, feqs, atoms, (o0, s0, nz0) = candidates(v0)
    r0 = relax(v0)
    print(f'{name}: out12={o0} score={s0} broken={nz0} relax={r0}')
    print(f'  {len(feqs)} failing eqs outside the twelve, {len(atoms)} atoms in '
          f'them, {len(F)} candidate free inputs', flush=True)

    if len(sys.argv) > 2 and sys.argv[2] == 'kern':
        # linear compensation kernel of the failing equations over their atoms
        al = sorted(atoms)
        ai = {a: i for i, a in enumerate(al)}
        M = []
        for e in feqs:
            m, sq, co = L.eq_atoms[e]
            M.append({ai[a]: c for a, c in co.items()})
        piv = {}
        for r in M:
            r = dict(r)
            while True:
                h = next((c for c in r if c in piv), None)
                if h is None:
                    break
                fq = r[h]
                for k, vv in piv[h].items():
                    nv = r.get(k, 0) - fq * vv
                    if nv:
                        r[k] = nv
                    else:
                        r.pop(k, None)
            if not r:
                continue
            c0 = min(r)
            piv[c0] = {k: v / r[c0] for k, v in r.items()}
        rank = len(piv)
        print(f'  compensation system {len(feqs)} eqs x {len(al)} atoms, '
              f'rank {rank}, kernel dim {len(al)-rank}')
        seedcol = [ai[a] for a in nz0 if a in ai]
        print(f'  residual atoms among columns: {[a for a in nz0 if a in ai]}')
        print(f'  pivots on residual cols: '
              f'{[c for c in seedcol if c in piv]} (free -> compensable)')
        sys.exit(0)

    a, b = int(sys.argv[2]), int(sys.argv[3])
    have = done_set()
    f = open(LOG, 'a')
    t0 = time.time()
    best = o0
    bests = s0
    DS = [1, 1000003, P, -1]
    for i in range(a, min(b, len(F))):
        u = F[i]
        ds = list(DS)
        for at in nz0:
            if u in L.avars[at]:
                t = T.solve_lin(at, u, v0)
                if t is not None:
                    ds.append(t - v0[u])
        for d in ds:
            if (name, u, str(d)) in have:
                continue
            v = list(v0)
            v[u] = v[u] + d
            J.fwd2(v, 2)
            if relax(v) != r0:
                continue
            o, s, nz, av = EN.state(v)
            f.write(json.dumps({'state': name, 'u': u, 'd': str(d), 'out12': o,
                                'score': s, 'broken': nz}) + '\n')
            if s > bests:
                bests = s
                best = min(best, o)
                print(f'  * x_{u} += {str(d)[:14]}: out12={o} SCORE={s} '
                      f'broken={nz[:10]}', flush=True)
                if s > 39026:
                    json.dump({f'x_{k}': v[k] for k in range(L.NVARS)
                               if v[k] != 0}, open(BEST, 'w'))
                    print('   >>> saved jm_best.json', flush=True)
        f.flush()
        if (i - a) % 60 == 0:
            print(f'  ..{i}/{len(F)} ({time.time()-t0:.0f}s)', flush=True)
    f.close()
    print(f'batch done, best out12={best} best score={bests} '
          f'({time.time()-t0:.0f}s)', flush=True)
