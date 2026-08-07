"""Fast advice-congruence solver: u := w exactly, handle := 0, sweep to fixed point."""
import collections, sys, time, json
import dlib as L
import engine2 as E
P = L.P

occ = collections.Counter()
for a in range(L.NA):
    for u in L.avars[a]:
        occ[u] += 1
solo = {u for u in L.freeset if occ[u] == 1}


def cone_free(target):
    seen, stack, out = set(), [target], set()
    while stack:
        u = stack.pop()
        if u in seen:
            continue
        seen.add(u)
        a = L.definer.get(u)
        if a is None:
            out.add(u)
            continue
        for w in L.avars[a]:
            if w != u:
                stack.append(w)
    return out


def build():
    rows = []
    for a in range(L.NA):
        if a in L.atom_out:
            continue
        p = L.polys[a]
        if max((len(m) for m in p), default=0) != 1 or len(p) != 3 or p.get(()):
            continue
        lin = {m[0]: c for m, c in p.items()}
        ones = [u for u, c in lin.items() if abs(c) == 1]
        ks = [u for u, c in lin.items() if abs(c) != 1]
        if len(ones) != 1 or len(ks) != 2:
            continue
        hv = ones[0]
        frees = [u for u in ks if u in L.freeset]
        if len(frees) != 1:
            continue
        u = frees[0]
        w = ks[0] if ks[1] == u else ks[1]
        hs = [z for z in cone_free(hv) if z in solo]
        rows.append((a, u, w, lin[u], lin[w], lin[hv], hv, hs))
    return rows


ADV = build()


def sweep(st, rounds=12, verbose=False):
    for r in range(rounds):
        nchg = 0
        for a, u, w, ku, kw, kh, hv, hs in ADV:
            if st.av[a] == 0:
                continue
            seeds = {}
            for h in hs:
                seeds[h] = 0
            # after handle -> 0, need ku*u + kw*w == 0
            target = -kw * st.v[w]
            if target % ku:
                continue
            seeds[u] = target // ku
            r0 = st.apply(seeds)
            if st.av[a] != 0:
                st.revert(r0)
                continue
            nchg += 1
        if verbose:
            print(f'  sweep {r}: changed {nchg} score={st.score} nz={len(st.nz())}')
        if nchg == 0:
            break
    return st


if __name__ == '__main__':
    print('advice rows', len(ADV), ' with handle:', sum(1 for r in ADV if r[-1]))
    st = E.St(L.load(sys.argv[1] if len(sys.argv) > 1 else 'D_state1.json'))
    print('start', st.score, st.nz())
    t0 = time.time()
    sweep(st, verbose=True)
    print('after sweep', st.score, st.nz(), f'{time.time()-t0:.1f}s')
    out = sys.argv[2] if len(sys.argv) > 2 else 'D_adv.json'
    L.save(st.v, out)
