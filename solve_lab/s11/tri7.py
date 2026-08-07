"""Correct triangular solve for channel U=0,V=1, bits (490,91).
   Links in THIS channel:  a7881<-x2751, a21050<-x16441, a26839<-x18751, a40065<-x28955
   Mirror controls:        x31339, x33708
   Then close every check exactly with its own cone handle.
"""
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
MIRROR = [31339, 33708]
LINKS = [(21050, 16441, lambda v: v[16441] - v[4920]),
         (40065, 28955, lambda v: v[28955] - v[11408]),
         (7881,  2751,  lambda v: v[2751] - v[1085]),
         (26839, 18751, lambda v: v[18751] - v[33091])]


def drive(v, ctrl, get, target=0, tries=5):
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


def arithmetic(v):
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
    dd = v[13682] + C0B
    v[14393] = 0
    v[11436] = (dd // P) if dd % P == 0 else 0
    fw.forward(v)
    drive(v, 14515, lambda vv: vv[1308] - vv[14853], 0)
    v[16742] = v[19083]
    fw.forward(v)


def mirror(v, rng, iters=15):
    for it in range(iters):
        r = [v[3719] % P, v[25118] % P]
        if not any(r):
            return True
        J = [[0] * len(MIRROR) for _ in range(2)]
        for j, c in enumerate(MIRROR):
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
        for j, c in enumerate(MIRROR):
            v[c] = (v[c] + d[j]) % P
        fw.forward(v)
    if not (v[3719] % P or v[25118] % P):
        return True
    old = v[33708]
    pts = []
    for t in range(5):
        v[33708] = t
        fw.forward(v)
        pts.append((t, (v[23776] * v[2401] * v[2401] - v[26196] * v[26196]) % P))
    v[33708] = old
    fw.forward(v)
    coeff = [0] * len(pts)
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
    fp = polyroot.norm(coeff, P)
    if len(fp) >= 2:
        for r0 in polyroot.roots(fp, P, rng):
            v[33708] = r0
            fw.forward(v)
            if drive(v, 31339, lambda vv: vv[25118], 0) and v[3719] % P == 0:
                return True
    v[33708] = old
    fw.forward(v)
    return False


def close_all(v, locked, rounds=12, verbose=False):
    for rnd in range(rounds):
        bad = fw.bad_checks(v)
        if not bad:
            break
        prog = False
        for a in sorted(bad, key=lambda a: (len(L.atom2eq.get(a, {})), a)):
            if fw.evalpoly(L.polys[a], v) == 0:
                continue
            cs = []
            for u in L.avars[a]:
                if L.definer.get(u) is None and u not in locked and \
                        not any(mm.count(u) > 1 for mm in L.polys[a]):
                    cs.append((u, None))
            cs.sort(key=lambda kv: (NAT[kv[0]], kv[0]))
            try:
                hs, base = deep.handles(v, a, locked=locked)
                cs += [(t, d) for t, d in sorted(hs, key=lambda kv: (NAT[kv[0]], kv[0]))]
            except Exception:
                pass
            for t, d in cs:
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
            print(f"     close{rnd}: bad={len(nb)} {nb[:12]}", flush=True)
        if not prog or set(nb) == set(bad):
            break
    return v


def run(seed):
    rng = random.Random(seed)
    v = uv01.state((490, 91))
    for c in MIRROR:
        v[c] = rng.randrange(0, 1 << 60)
    fw.forward(v)
    arithmetic(v)
    for a, ctrl, get in LINKS:
        drive(v, ctrl, get, 0)
    m = mirror(v, rng)
    for a, ctrl, get in LINKS:
        drive(v, ctrl, get, 0)
    arithmetic(v)
    res = [v[3719] % P == 0, v[25118] % P == 0] + [g(v) % P == 0 for _, _, g in LINKS]
    LOCK = {490, 91, 19750, 7497, 22820, 14853, 14393, 11436, 14515, 16742,
            22162, 30213, 8386, 21868} | set(MIRROR) | {c for _, c, _ in LINKS}
    close_all(v, LOCK)
    return v, res


if __name__ == '__main__':
    t0 = time.time()
    best = None
    for seed in range(30):
        v, res = run(seed)
        bad = fw.bad_checks(v)
        f = L.failing_eqs(L.all_atom_values(v))
        if best is None or len(f) < best[0]:
            best = (len(f), [x for x in v], seed)
        print(f"seed{seed}: modp={res} bad={len(bad)} failing={len(f)} score={L.NEQ-len(f)} "
              f"{bad[:12]} ({time.time()-t0:.0f}s)", flush=True)
        if not bad:
            print("FULL SOLVE")
            break
    print(f"\nBEST failing={best[0]} score={L.NEQ-best[0]} seed={best[2]}")
    json.dump({str(i): best[1][i] for i in range(L.NVARS)}, open(os.path.join(HERE, 'data', 'tri7_best.json'), 'w'))
