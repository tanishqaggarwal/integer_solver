"""Triangular solve in channel U=0,V=1 (bits 490,91):
     arithmetic drives (x19750, x14853, x14515, x16742)
  -> mirror  x3719 = x25118 = 0   with {16441, 28955, 31339, 33708}  (Newton + cubic fallback)
  -> a34580 <- x33129 ;  a14445 <- x18751 ;  a33796 <- x37088 ;  a27139 <- x2751
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
MIRROR = [16441, 28955, 31339, 33708]
LINKS = [(34580, 33129, lambda v: v[33708] - v[10170]),
         (14445, 18751, lambda v: v[33129] - v[3757]),
         (33796, 37088, lambda v: v[31339] - v[6858]),
         (27139, 2751,  lambda v: v[37088] - v[13585])]


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
    # cubic fallback in x33708:  G(t) = x23776*x2401^2 - x26196^2
    old = v[33708]
    pts = []
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


def run(seed, verbose=True):
    rng = random.Random(seed)
    v = uv01.state((490, 91))
    for c in MIRROR:
        v[c] = rng.randrange(0, 1 << 60)
    fw.forward(v)
    arithmetic(v)
    m = mirror(v, rng)
    for a, ctrl, get in LINKS:
        drive(v, ctrl, get, 0)
    arithmetic(v)
    ok = (v[3719] % P == 0 and v[25118] % P == 0 and
          all(g(v) % P == 0 for _, _, g in LINKS))
    if verbose:
        print(f"  seed{seed}: mirror={m} all6={ok} "
              f"[{v[3719]%P==0},{v[25118]%P==0},"
              f"{','.join(str(g(v)%P==0) for _,_,g in LINKS)}]", flush=True)
    return v, ok


if __name__ == '__main__':
    t0 = time.time()
    best = None
    for seed in range(40):
        v, ok = run(seed)
        bad = fw.bad_checks(v)
        f = L.failing_eqs(L.all_atom_values(v))
        if best is None or len(f) < best[0]:
            best = (len(f), [x for x in v], seed)
        print(f"     bad={len(bad)} failing={len(f)} score={L.NEQ-len(f)} {bad[:12]} ({time.time()-t0:.0f}s)", flush=True)
        if ok and not bad:
            break
    print(f"\nBEST failing={best[0]} score={L.NEQ-best[0]} seed={best[2]}")
    json.dump({str(i): best[1][i] for i in range(L.NVARS)}, open(os.path.join(HERE, 'data', 'tri6_best.json'), 'w'))
