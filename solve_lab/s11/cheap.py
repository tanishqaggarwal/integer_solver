"""The cheap branch: U=V=1, a*b = c*d = 0.  Only the arithmetic cluster + 2 bits' load pins."""
import sys, os, json, time
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 's9', 'eff'))
import lib as L, fw, deep
P = L.P
C0 = L.polys[688][()]
MM = 8863713
G0 = (-C0 * pow(MM, -1, P)) % P
C0B = L.polys[1618][()]
NAT = {u: len(L.var_atoms[u]) for u in range(L.NVARS)}


def free_cands(a, locked):
    out = []
    for u in L.avars[a]:
        if L.definer.get(u) is not None or u in locked:
            continue
        if any(mm.count(u) > 1 for mm in L.polys[a]):
            continue
        out.append(u)
    out.sort(key=lambda u: (NAT[u], u))
    return out


def build(BITS, verbose=True):
    v = [0] * L.NVARS
    for b in BITS:
        v[b] = 1
    fw.forward(v)
    assert v[15298] == 1 and v[38170] == 0 and v[3896] == 0
    # n = m = 0 exactly ; gaps are already 0 so a29539 / a26731 come along
    v[14853] = v[12186]
    v[16742] = v[24908]
    v[8386] = 0
    v[21868] = 0
    fw.forward(v)
    # arithmetic cluster, exact
    v[30213] = G0
    v[22820] = 0
    v[7497] = (C0 + MM * G0) // P
    v[22162] = -C0B
    v[14393] = 0
    v[11436] = 0
    fw.forward(v)
    if verbose:
        print("  a688=%s a1618=%s a40608=%s  n=%d m=%d" %
              (fw.evalpoly(L.polys[688], v), fw.evalpoly(L.polys[1618], v),
               fw.evalpoly(L.polys[40608], v), v[29322], v[3558]))
    LOCK = set(BITS) | {14853, 16742, 30213, 22820, 7497, 22162, 14393, 11436, 8386, 21868, 12186, 24908}
    # close remaining, protecting the structure
    for rnd in range(30):
        bad = fw.bad_checks(v)
        if not bad:
            break
        prog = False
        for a in sorted(bad, key=lambda a: (len(L.atom2eq.get(a, {})), a)):
            if fw.evalpoly(L.polys[a], v) == 0:
                continue
            cands = [(t, None) for t in free_cands(a, LOCK)]
            try:
                hs, base = deep.handles(v, a, locked=LOCK)
                cands += [(t, d) for t, d in sorted(hs, key=lambda kv: (NAT[kv[0]], kv[0]))]
            except Exception:
                pass
            for t, d in cands:
                old = v[t]
                if d is None:
                    x = fw.solve_lin(a, t, v)
                    if x is None or x == old:
                        continue
                else:
                    bs = fw.evalpoly(L.polys[a], v)
                    if not d or bs % d:
                        continue
                    x = old - bs // d
                v[t] = x
                fw.forward(v)
                if fw.evalpoly(L.polys[a], v) == 0:
                    prog = True
                    break
                v[t] = old
                fw.forward(v)
        nb = fw.bad_checks(v)
        if verbose:
            print(f"    rnd{rnd}: bad={len(nb)} {nb[:12]}", flush=True)
        if not prog or set(nb) == set(bad):
            break
    fw.forward(v)
    return v


if __name__ == '__main__':
    best = None
    for BITS in [(542, 438), (47, 438), (542, 91), (47, 91), (24601, 438), (2081, 438)]:
        try:
            v = build(BITS)
        except AssertionError:
            print(f"BITS={BITS}: wrong channel, skipped")
            continue
        b = fw.bad_checks(v)
        av = L.all_atom_values(v)
        f = L.failing_eqs(av)
        print(f"BITS={BITS}: bad={len(b)} failing={len(f)} score={L.NEQ-len(f)}  {b}", flush=True)
        if best is None or len(f) < best[0]:
            best = (len(f), [x for x in v], BITS)
    print(f"\nBEST: failing={best[0]} score={L.NEQ-best[0]} bits={best[2]}")
    json.dump({str(i): best[1][i] for i in range(L.NVARS)},
              open(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'data', 'cheap_best.json'), 'w'))
