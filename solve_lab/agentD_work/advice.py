"""Find all 'advice' checks: K*(u - w) - h  or  u - C - h, with h a p-handle,
and report the dependency structure among them."""
import collections, sys
import dlib as L
import engine2 as E
P = L.P

occ = collections.Counter()
for a in range(L.NA):
    for u in L.avars[a]:
        occ[u] += 1
solo = {u for u in L.freeset if occ[u] == 1}


def cone_free(target, stopat=None):
    seen = set()
    stack = [target]
    out = set()
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


if __name__ == '__main__':
    st = E.St(L.load(sys.argv[1] if len(sys.argv) > 1 else 'D_state1.json'))
    print('score', st.score, 'nz', st.nz())
    # every check atom whose value is currently zero but whose structure is
    # K*(u-w) - h : find them by pattern.
    adv = []
    for a in range(L.NA):
        if a in L.atom_out:
            continue
        p = L.polys[a]
        if max((len(m) for m in p), default=0) != 2:
            continue
        lin = {m[0]: c for m, c in p.items() if len(m) == 1}
        quad = {m: c for m, c in p.items() if len(m) == 2}
        if len(quad) != 1 or len(lin) != 2 or p.get(()):
            continue
        (qm, qc), = quad.items()
        # handle term: qm = (wire, handle) with handle solo-free
        hs = [u for u in qm if u in solo]
        if not hs:
            continue
        us = sorted(lin)
        cs = [lin[u] for u in us]
        if abs(cs[0]) != abs(cs[1]):
            continue
        adv.append((a, us[0], us[1], cs[0], qm))
    print('advice-shaped checks K*(u-w) - wire*handle :', len(adv))
    freelhs = [(a, u, w) for a, u, w, c, q in adv if u in L.freeset or w in L.freeset]
    print('  with a FREE side:', len(freelhs))
    for a, u, w in freelhs:
        fu = u if u in L.freeset else w
        fw = w if u in L.freeset else u
        cf = cone_free(fw)
        print(f'  a{a:<6} eqs={len(L.atom2eq.get(a,{})):<3} free x_{fu}  ==  x_{fw} (mod p) ; '
              f'cone_free(rhs)={len(cf)}  val_ok={st.av[a]==0}')
