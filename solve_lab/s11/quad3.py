"""gamma(k,l) = x_12000/P after x31339 += k*P, x33708 += l*P.
   Determine the true bidegree empirically, interpolate exactly, verify, then solve mod 53*163027."""
import sys, os, json, time
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 's9', 'eff'))
import lib as L, fw, deep
from quad8640431 import quad_roots, load
P = L.P
HERE = os.path.dirname(os.path.abspath(__file__))
M = 8640431
C1, C2 = 31339, 33708


def sample(v, k, l, b1, b2):
    v[C1] = b1 + k * P
    v[C2] = b2 + l * P
    fw.forward(v)
    if v[3719] % P or v[25118] % P:
        return None
    return (v[12000] // P) % M


def interp1(ys, m):
    """Newton forward differences -> coefficients of the polynomial in one variable, mod m."""
    n = len(ys)
    diff = [ys[:]]
    for i in range(1, n):
        diff.append([(diff[-1][j + 1] - diff[-1][j]) % m for j in range(len(diff[-1]) - 1)])
    # p(x) = sum_i diff[i][0] * C(x, i)  -> expand to power basis
    coef = [0] * n
    binom = [[0] * n for _ in range(n)]
    for i in range(n):
        binom[i][0] = 0
    # build falling factorial expansion
    poly = [1] + [0] * (n - 1)   # current product (x)(x-1)...(x-i+1)/i!
    cur = [1] + [0] * (n - 1)
    fact = 1
    for i in range(n):
        if i > 0:
            new = [0] * n
            for d in range(n - 1):
                new[d + 1] = (new[d + 1] + cur[d]) % m
                new[d] = (new[d] - (i - 1) * cur[d]) % m
            cur = new
            fact = fact * i % m
        invf = pow(fact, -1, m)
        for d in range(n):
            coef[d] = (coef[d] + diff[i][0] * cur[d] % m * invf) % m
    return coef


if __name__ == '__main__':
    v = load('three.json')
    b1, b2 = v[C1], v[C2]
    print("bad:", fw.bad_checks(v), "failing:", len(L.failing_eqs(L.all_atom_values(v))))
    # empirical degrees
    ysk = [sample(v, k, 0, b1, b2) for k in range(6)]
    ysl = [sample(v, 0, l, b1, b2) for l in range(7)]
    v[C1], v[C2] = b1, b2
    fw.forward(v)
    def deg_of(ys, m):
        d = ys[:]
        for order in range(len(ys)):
            if all(x % m == 0 for x in d[1:]) and len(d) > 1 and all(x % m == 0 for x in d):
                return order - 1
            nd = [(d[i + 1] - d[i]) % m for i in range(len(d) - 1)]
            if all(x == 0 for x in nd):
                return order
            d = nd
        return len(ys) - 1
    dk = deg_of(ysk, M)
    dl = deg_of(ysl, M)
    print(f"  degree in k = {dk}, degree in l = {dl}")
    DK, DL = dk + 1, dl + 1
    grid = [[sample(v, k, l, b1, b2) for l in range(DL)] for k in range(DK)]
    v[C1], v[C2] = b1, b2
    fw.forward(v)
    if any(x is None for row in grid for x in row):
        print("  a sample broke the mirror"); sys.exit()
    # interpolate in l for each k, then in k for each l-coefficient
    rows = [interp1(grid[k], M) for k in range(DK)]
    coef = [[0] * DL for _ in range(DK)]
    for j in range(DL):
        col = interp1([rows[k][j] for k in range(DK)], M)
        for i in range(DK):
            coef[i][j] = col[i]

    def gam(k, l, m=M):
        s = 0
        for i in range(DK):
            for j in range(DL):
                s += coef[i][j] * pow(k, i, m) % m * pow(l, j, m)
        return s % m
    ok = all(gam(k, l) == grid[k][l] for k in range(DK) for l in range(DL))
    extra = sample(v, DK + 1, DL + 1, b1, b2)
    v[C1], v[C2] = b1, b2
    fw.forward(v)
    ok2 = (extra is not None and gam(DK + 1, DL + 1) == extra)
    print(f"  interpolation exact on grid: {ok} ; on a held-out point: {ok2}")
    if not (ok and ok2):
        print("  model wrong - abort"); sys.exit()
    # solve mod 53 (brute) and mod 163027 (iterate l, solve in k if dk==2)
    s53 = [(k, l) for k in range(53) for l in range(53) if gam(k, l, 53) % 53 == 0]
    print(f"  solutions mod 53: {len(s53)}")
    q = 163027
    s163 = []
    for l in range(q):
        # polynomial in k of degree dk
        cs = [sum(coef[i][j] * pow(l, j, q) for j in range(DL)) % q for i in range(DK)]
        if dk == 2:
            for k in quad_roots(cs[2], cs[1], cs[0], q):
                s163.append((k, l))
        elif dk == 1:
            if cs[1] % q:
                s163.append(((-cs[0]) * pow(cs[1], -1, q) % q, l))
        if len(s163) >= 30:
            break
    print(f"  solutions mod 163027: {len(s163)}")
    if not s53 or not s163:
        print("  none"); sys.exit()
    inv53 = pow(163027, -1, 53)
    inv163 = pow(53, -1, 163027)
    hits = 0
    for (k1, l1) in s53[:6]:
        for (k2, l2) in s163[:6]:
            k = (k1 * 163027 * inv53 + k2 * 53 * inv163) % M
            l = (l1 * 163027 * inv53 + l2 * 53 * inv163) % M
            g = sample(v, k, l, b1, b2)
            bad = fw.bad_checks(v)
            f = L.failing_eqs(L.all_atom_values(v))
            print(f"  k,l: gamma={g} bad={len(bad)} failing={len(f)} score={L.NEQ-len(f)} {bad[:10]}", flush=True)
            if g == 0:
                hits += 1
                json.dump({str(i): v[i] for i in range(L.NVARS)},
                          open(os.path.join(HERE, 'data', 'quad3_hit.json'), 'w'))
                print("    *** gamma == 0 ; saved")
            v[C1], v[C2] = b1, b2
            fw.forward(v)
            if hits:
                break
        if hits:
            break
