"""Lazy-constraint simultaneous solve: knobs = FULL cone closure (no affinity filter),
   rows activated on demand from the exact re-propagation.  Answers whether E's closure
   was knob-starved or the pin repairs are genuinely unreachable."""
import sys, os, time, json
sys.path.insert(0, '/home/user/integer_solver/solve_lab/agentO_work')
import simO
import engine as E, sparse

OD = '/home/user/integer_solver/solve_lab/agentO_work'


def run(seed, frozen, maxr=5, maxv=4000, iters=25, maxcore=600, maxcorebits=8_000_000,
        log=print):
    v0 = E.forward(seed); bad0 = E.badatoms(v0)
    if not bad0:
        return 0, seed, [], v0, []
    S, cols, nonlin, rounds = simO.closure(v0, bad0, frozen, maxr, maxv, verbose=False)
    allrows = set(bad0)
    for f in S:
        allrows |= set(cols[f])
    rowmap = {a: {} for a in allrows}
    for f in S:
        for a, c in cols[f].items():
            rowmap[a][f] = c
    log(f'  closure: {len(S)} knobs, {len(allrows)} reachable rows, bad0={sorted(bad0)}')
    A = set(bad0)
    best = (len(E.eqfails(bad0)), dict(seed), sorted(bad0), v0)
    hist = []
    for it in range(iters):
        use = sorted(A)
        rows = [rowmap.get(a, {}) for a in use]
        rhs = [-bad0.get(a, 0) for a in use]
        t0 = time.time()
        sol, msg, _ = sparse.solve_sparse(rows, rhs, names=use, verbose=False,
                                          maxcore=maxcore, maxcorebits=maxcorebits)
        if sol is None:
            log(f'  it{it}: |A|={len(A)} UNSAT {msg[:90]} ({time.time()-t0:.0f}s)')
            hist.append((it, len(A), 'unsat', msg))
            return best + (hist,)
        ns = dict(seed)
        for f, d in sol.items():
            if d:
                ns[f] = v0[f] + d
        v = E.forward(ns); av = E.badatoms(v); nf = len(E.eqfails(av))
        moved = sum(1 for d in sol.values() if d)
        log(f'  it{it}: |A|={len(A)} moved={moved} -> fails={nf} score={39033-nf} '
            f'nbad={len(av)} bad={sorted(av)[:10]} ({time.time()-t0:.0f}s)')
        hist.append((it, len(A), nf, sorted(av)))
        if nf < best[0]:
            best = (nf, dict(ns), sorted(av), v)
        if nf == 0:
            return 0, ns, [], v, hist
        newr = set(av) - A
        if not newr:
            log('  no new rows -> stall')
            break
        A |= newr
    return best + (hist,)


if __name__ == '__main__':
    bits = [] if sys.argv[1] == 'empty' else [int(x) for x in sys.argv[1].split(',')]
    maxr = int(sys.argv[2]) if len(sys.argv) > 2 else 5
    maxv = int(sys.argv[3]) if len(sys.argv) > 3 else 4000
    s = dict(simO.C.base)
    for b in bits:
        s[b] = 1
    print(f'=== lazy bits={bits} maxr={maxr} maxv={maxv}', flush=True)
    t0 = time.time()
    r = run(s, set(bits), maxr=maxr, maxv=maxv)
    print(f'  BEST fails={r[0]} score={39033-r[0]} bad={r[2]} ({time.time()-t0:.0f}s)', flush=True)
    if r[0] < 28:
        tag = "_".join(map(str, bits)) or "empty"
        json.dump({f"x_{i}": str(int(r[3][i])) for i in range(E.NV) if r[3][i] != 0},
                  open(f'{OD}/lazy_{tag}_{39033-r[0]}.json', 'w'))
