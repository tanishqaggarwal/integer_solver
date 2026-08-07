"""Full solve attempt in channel U=0,V=1 (c*d=1), bits (490,91).
   Fixed-point: close linking checks -> re-solve mirror -> re-drive arithmetic -> repeat."""
import sys, os, json, time, random
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 's9', 'eff'))
import lib as L, fw, deep, uv01, polyroot
from gfp import gauss_solve
P = L.P
HERE = os.path.dirname(os.path.abspath(__file__))
C0 = L.polys[688][()]
MM = 8863713
G0 = (-C0 * pow(MM, -1, P)) % P
C0B = L.polys[1618][()]
NAT = {u: len(L.var_atoms[u]) for u in range(L.NVARS)}
LINK = [(7881, 2751), (21050, 16441), (26839, 18751), (40065, 28955)]
MIRROR_CTRL = [22917, 31339, 33708]


def drive(v, ctrl, get, target=0, tries=4):
    for _ in range(tries):
        r = (get(v) - target) % P
        if r == 0:
            return True
        old = v[ctrl]
        v[ctrl] = old + 1
        fw.forward(v)
        s = ((get(v) - target) % P - r) % P
        v[ctrl] = old
        fw.forward(v)
        if s == 0:
            return False
        v[ctrl] = (old + (-r) * pow(s, -1, P)) % P
        fw.forward(v)
    return (get(v) - target) % P == 0


def structural(v):
    v[22162] = 0
    v[30213] = 0
    v[8386] = 0
    v[21868] = 0
    fw.forward(v)
    drive(v, 19750, lambda vv: vv[37892], G0)
    num = C0 + MM * v[37892]
    v[7497] = num // P if num % P == 0 else 0
    v[22820] = 0
    fw.forward(v)
    v[14853] = (-C0B) % P
    fw.forward(v)
    d = v[13682] + C0B
    v[14393] = 0
    v[11436] = (d // P) if d % P == 0 else 0
    fw.forward(v)
    drive(v, 14515, lambda vv: vv[1308] - vv[14853], 0)
    v[16742] = v[19083]
    fw.forward(v)
    return v


def mirror(v, rng, iters=10):
    """drive x3719 = x25118 = 0 using MIRROR_CTRL; Newton, then cubic fallback on x33708"""
    for it in range(iters):
        r = [v[3719] % P, v[25118] % P]
        if not any(r):
            return True
        J = [[0] * len(MIRROR_CTRL) for _ in range(2)]
        for j, c in enumerate(MIRROR_CTRL):
            old = v[c]
            v[c] = old + 1
            fw.forward(v)
            J[0][j] = (v[3719] % P - r[0]) % P
            J[1][j] = (v[25118] % P - r[1]) % P
            v[c] = old
            fw.forward(v)
        d = gauss_solve(J, [(-x) % P for x in r], P)
        if d is None:
            break
        for j, c in enumerate(MIRROR_CTRL):
            v[c] = (v[c] + d[j]) % P
        fw.forward(v)
    if not (v[3719] % P or v[25118] % P):
        return True
    # cubic fallback: G(t) = x23776*x2401^2 - x26196^2  as a function of t = x33708
    pts = []
    old = v[33708]
    for t in range(5):
        v[33708] = t
        fw.forward(v)
        pts.append((t, (v[23776] * v[2401] * v[2401] - v[26196] * v[26196]) % P))
    v[33708] = old
    fw.forward(v)
    n = len(pts)
    coeff = [0] * n
    for i, (xi, yi) in enumerate(pts):
        num, den = [1], 1
        for j, (xj, _) in enumerate(pts):
            if i == j:
                continue
            num = polyroot.pmul(num, [(-xj) % P, 1], P)
            den = den * ((xi - xj) % P) % P
        f = yi * pow(den, -1, P) % P
        for k, c in enumerate(num):
            coeff[k] = (coeff[k] + c * f) % P
    fpoly = polyroot.norm(coeff, P)
    if len(fpoly) < 2:
        return False
    for r0 in polyroot.roots(fpoly, P, rng):
        v[33708] = r0
        fw.forward(v)
        if drive(v, 31339, lambda vv: vv[25118], 0) and v[3719] % P == 0:
            return True
    v[33708] = old
    fw.forward(v)
    return False


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


def close_one(v, a, locked):
    if fw.evalpoly(L.polys[a], v) == 0:
        return True
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
            return True
        v[t] = old
        fw.forward(v)
    return False


if __name__ == '__main__':
    rng = random.Random(11)
    BITS = (490, 91)
    v = uv01.state(BITS)
    structural(v)
    LOCK = set(BITS) | {19750, 7497, 22820, 14853, 14393, 11436, 14515, 16742,
                        22162, 30213, 8386, 21868} | set(MIRROR_CTRL) | {t for _, t in LINK}
    t0 = time.time()
    best = (len(L.failing_eqs(L.all_atom_values(v))), [x for x in v])
    for rnd in range(25):
        for a, t in LINK:
            fw.solve_lin(a, t, v)
            x = fw.solve_lin(a, t, v)
            if x is not None:
                v[t] = x
                fw.forward(v)
        m = mirror(v, rng)
        structural(v)
        for a in fw.bad_checks(v):
            if a not in [x for x, _ in LINK]:
                close_one(v, a, LOCK)
        bad = fw.bad_checks(v)
        av = L.all_atom_values(v)
        f = L.failing_eqs(av)
        if len(f) < best[0]:
            best = (len(f), [x for x in v])
        print(f"  rnd{rnd}: mirror={m} bad={len(bad)} failing={len(f)} score={L.NEQ-len(f)} "
              f"({time.time()-t0:.0f}s) {bad[:12]}", flush=True)
        if not bad:
            break
    v[:] = best[1]
    fw.forward(v)
    b = fw.bad_checks(v)
    f = L.failing_eqs(L.all_atom_values(v))
    print(f"BEST failing={len(f)} score={L.NEQ-len(f)} bad={len(b)} {b}")
    json.dump({str(i): v[i] for i in range(L.NVARS)}, open(os.path.join(HERE, 'data', 'uv01_full.json'), 'w'))
