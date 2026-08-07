"""Full construction in channel U=1,V=0 (a*b=1, c=d=0)."""
import sys, os, json, time
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 's9', 'eff'))
import lib as L, fw, deep
from gfp import gauss_solve
P = L.P
HERE = os.path.dirname(os.path.abspath(__file__))
LD = json.load(open(os.path.join(HERE, 'data', 'loads.json')))['loads']
C0 = L.polys[688][()]
MM = 8863713
G0 = (-C0 * pow(MM, -1, P)) % P
C0B = L.polys[1618][()]
NAT = {u: len(L.var_atoms[u]) for u in range(L.NVARS)}
LINK = [(15030, 13222), (21853, 30709), (30067, 7418), (37735, 38667)]


def pins_of(bit):
    return [(rec[0], rec[1]) for rec in LD.get(str(bit), [])]


def solve_for(v, a, t):
    x = fw.solve_lin(a, t, v)
    if x is not None and x != v[t]:
        v[t] = x
        fw.forward(v)
        return True
    return False


def build(BITS, mirror_ctrl=(14681, 28486), verbose=True):
    v = [0] * L.NVARS
    for b in BITS:
        v[b] = 1
    fw.forward(v)
    for b in BITS:
        for a, X in pins_of(b):
            solve_for(v, a, X)
    fw.forward(v)
    # channel invariants
    v[22162] = 0
    v[30213] = 0
    v[8386] = 0
    v[21868] = 0
    fw.forward(v)
    # arithmetic: x37892 = x16742 (x34606 = 1)
    v[16742] = G0
    v[22820] = 0
    fw.forward(v)
    num = C0 + MM * v[37892]
    if num % P:
        if verbose:
            print("   a688 residue mismatch")
    v[7497] = num // P
    fw.forward(v)
    # a1618: x13682 = x12186 ; drive x12186 to -C0B mod P via x5096 (slope +1)
    tgt = (-C0B) % P
    v[5096] = (v[5096] + (tgt - v[12186])) % P
    fw.forward(v)
    d = v[13682] + C0B
    v[14393] = 0
    v[11436] = d // P if d % P == 0 else 0
    fw.forward(v)
    # a26731 : x19083 == x16742  (x21589 slope)
    for _ in range(3):
        r = (v[16742] - v[19083]) % P
        if r == 0:
            break
        old = v[21589]
        v[21589] = old + 1
        fw.forward(v)
        s = ((v[16742] - v[19083]) % P - r) % P
        v[21589] = old
        fw.forward(v)
        if s == 0:
            break
        v[21589] = (old + (-r) * pow(s, -1, P)) % P
        fw.forward(v)
    # a29539 : x14853 == x1308
    v[14853] = v[1308]
    fw.forward(v)
    # linking checks
    for a, t in LINK:
        solve_for(v, a, t)
    fw.forward(v)
    # group-2 mirror: x25614 = x34220 = 0 via mirror_ctrl (2x2 Newton mod p)
    for it in range(12):
        r = [v[25614] % P, v[34220] % P]
        if not any(r):
            break
        J = [[0, 0], [0, 0]]
        for j, c in enumerate(mirror_ctrl):
            old = v[c]
            v[c] = old + 1
            fw.forward(v)
            J[0][j] = (v[25614] % P - r[0]) % P
            J[1][j] = (v[34220] % P - r[1]) % P
            v[c] = old
            fw.forward(v)
        dd = gauss_solve(J, [(-x) % P for x in r], P)
        if dd is None:
            if verbose:
                print("   mirror system singular")
            break
        for j, c in enumerate(mirror_ctrl):
            v[c] = (v[c] + dd[j]) % P
        fw.forward(v)
    if verbose:
        print(f"   x25614%p==0:{v[25614]%P==0} x34220%p==0:{v[34220]%P==0} "
              f"a688={fw.evalpoly(L.polys[688],v)==0} a1618={fw.evalpoly(L.polys[1618],v)==0}")
    return v


def cands(v, a, locked):
    out = []
    for u in L.avars[a]:
        if L.definer.get(u) is not None or u in locked:
            continue
        if any(mm.count(u) > 1 for mm in L.polys[a]):
            continue
        out.append((u, None))
    out.sort(key=lambda kv: (NAT[kv[0]], kv[0]))
    try:
        hs, base = deep.handles(v, a, locked=locked)
        out += [(t, d) for t, d in sorted(hs, key=lambda kv: (NAT[kv[0]], kv[0]))]
    except Exception:
        pass
    return out


def finish(v, locked, rounds=20, verbose=True):
    for rnd in range(rounds):
        bad = fw.bad_checks(v)
        if not bad:
            break
        prog = False
        for a in sorted(bad, key=lambda a: (len(L.atom2eq.get(a, {})), a)):
            if fw.evalpoly(L.polys[a], v) == 0:
                continue
            for t, d in cands(v, a, locked):
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
            print(f"    finish{rnd}: bad={len(nb)} {nb[:12]}", flush=True)
        if not prog or set(nb) == set(bad):
            break
    return v


if __name__ == '__main__':
    best = None
    for BITS in [(542, 47), (1685, 47), (542, 1502), (853, 112), (1357, 1438)]:
        print(f"=== BITS={BITS}", flush=True)
        try:
            v = build(BITS)
        except Exception as e:
            print("   build failed:", e)
            continue
        LOCK = set(BITS) | {16742, 22820, 7497, 5096, 14393, 11436, 21589, 14853,
                            22162, 30213, 8386, 21868, 14681, 28486} | {t for _, t in LINK}
        v = finish(v, LOCK)
        b = fw.bad_checks(v)
        av = L.all_atom_values(v)
        f = L.failing_eqs(av)
        print(f"   RESULT bad={len(b)} failing={len(f)} score={L.NEQ-len(f)}  {b}", flush=True)
        if best is None or len(f) < best[0]:
            best = (len(f), [x for x in v], BITS)
    print(f"\nBEST failing={best[0]} score={L.NEQ-best[0]} bits={best[2]}")
    json.dump({str(i): best[1][i] for i in range(L.NVARS)},
              open(os.path.join(HERE, 'data', 'uv_best.json'), 'w'))
