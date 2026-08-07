"""Advice congruence machinery: the 193 checks  K*(u_free - w) - hvar = 0."""
import collections, sys, random, time
import dlib as L
import engine2 as E
P = L.P

ADV = []      # (atom, u_free, w, K, hvar, handle_free)
occ = collections.Counter()
for a in range(L.NA):
    for u in L.avars[a]:
        occ[u] += 1
solo = {u for u in L.freeset if occ[u] == 1}

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
    K = abs(lin[ks[0]])
    frees = [u for u in ks if u in L.freeset]
    if len(frees) != 1:
        continue
    u = frees[0]
    w = ks[0] if ks[1] == u else ks[1]
    # handle: hv should be defined as wire*handle
    hd = L.definer.get(hv)
    hf = None
    if hd is not None:
        cand = [z for z in L.avars[hd] if z in solo]
        if len(cand) == 1:
            hf = cand[0]
    ADV.append((a, u, w, K, lin[u], lin[w], lin[hv], hv, hf))

if __name__ == '__main__':
    print('advice checks parsed:', len(ADV))
    nohandle = [r for r in ADV if r[-1] is None]
    print('  without identified solo handle:', len(nohandle))
    st = E.St(L.load(sys.argv[1] if len(sys.argv) > 1 else 'D_state1.json'))
    print('score', st.score, 'nz', st.nz())
    bad = [r for r in ADV if st.av[r[0]] != 0]
    print('violated advice checks:', [r[0] for r in bad])
