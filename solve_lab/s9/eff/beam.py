"""Beam search over sacrifice sets S minimising  f(S) = |S| - D(S)
(failing = f(S) + C, C = number of surviving mod-P congruences, empirically 2).

Candidate extensions: atoms that occur in an already-sacrificed equation (only those can
add a column to the equations we are trying to zero), priced at their extra equations.
"""
import pickle, sys, time, collections
import lib as L, model as MD, opt

v0 = opt.init()
S13 = frozenset([2554, 6816, 8124, 8680, 9123, 9421, 12231, 12270, 12350, 14584, 18673, 22044, 29125])


def candidates(S, maxextra=12):
    S = frozenset(S)
    A = set(MD.confined_atoms(S))
    ats = set()
    for i in S:
        ats |= set(L.eq_atoms[i][2])
    out = []
    for a in ats - A:
        ex = frozenset(L.atom2eq[a]) - S
        if 0 < len(ex) <= maxextra:
            out.append((len(ex), a, ex))
    out.sort()
    return out


def f_of(S):
    mod = MD.build(S, v0, verbose=False)
    D = opt.rank_of(mod)
    return len(S) - D, D, mod


if __name__ == '__main__':
    width = int(sys.argv[1]) if len(sys.argv) > 1 else 8
    depth = int(sys.argv[2]) if len(sys.argv) > 2 else 4
    f0, D0, _ = f_of(S13)
    print(f'START |S|=13 D={D0} f={f0}')
    beam = [(f0, S13)]
    seen = {S13}
    best = (f0, S13)
    for d in range(depth):
        nxt = []
        t0 = time.time()
        for f, S in beam:
            for nex, a, ex in candidates(S):
                S2 = frozenset(S) | ex
                if S2 in seen:
                    continue
                seen.add(S2)
                f2, D2, _ = f_of(S2)
                nxt.append((f2, S2, a, nex, D2))
        nxt.sort(key=lambda t: (t[0], len(t[1])))
        print(f'--- depth {d+1}: {len(nxt)} candidates, {time.time()-t0:.0f}s; best 12:')
        for f2, S2, a, nex, D2 in nxt[:12]:
            print(f'    f={f2}  |S|={len(S2)} D={D2}  (+atom {a}, +{nex} eqs)')
        if not nxt:
            break
        if nxt[0][0] < best[0]:
            best = (nxt[0][0], nxt[0][1])
        beam = [(f2, S2) for f2, S2, a, nex, D2 in nxt[:width]]
    print('BEST f =', best[0], ' |S| =', len(best[1]))
    pickle.dump(best, open('beam_best.pkl', 'wb'))
