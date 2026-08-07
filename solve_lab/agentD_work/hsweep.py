"""Handle sweep: for every nonzero check atom, try to zero it with a solo free
handle (a free input occurring in exactly one atom).  Exact integer division only."""
import collections, sys, time
import dlib as L
import engine2 as E
import rad
P = L.P

occ = collections.Counter()
for a in range(L.NA):
    for u in L.avars[a]:
        occ[u] += 1
SOLO = {u for u in L.freeset if occ[u] == 1}


def sweep(st, rounds=6, verbose=False, extra=()):
    for r in range(rounds):
        nchg = 0
        for c in st.nz():
            c0 = st.av[c]
            if c0 == 0:
                continue
            kn = rad.free_knobs(c, st.v)
            cands = [u for u in kn if u in SOLO] + [u for u in extra if u in kn]
            for u in cands:
                b = st.v[u]
                rr = st.apply({u: b + 1})
                c1 = st.av[c]
                st.revert(rr)
                s = c1 - c0
                if s == 0 or c0 % s:
                    continue
                rr = st.apply({u: b - c0 // s})
                if st.av[c] != 0:
                    st.revert(rr)
                    continue
                nchg += 1
                break
        if verbose:
            print(f'  hsweep {r}: fixed {nchg} score={st.score} nz={len(st.nz())}', flush=True)
        if nchg == 0:
            break
    return st


if __name__ == '__main__':
    st = E.St(L.load(sys.argv[1]))
    print('start', st.score, st.nz())
    sweep(st, verbose=True)
    print('end', st.score, st.nz())
    out = sys.argv[2] if len(sys.argv) > 2 else 'D_hs.json'
    L.save(st.v, out)
