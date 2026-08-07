import sys, os, json, time, random
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 's9', 'eff'))
import lib as L, fw, quick, polyroot
P = L.P


def solve1(th, ctrl, idx, tries=4):
    for _ in range(tries):
        r = quick.six(quick.ev(th))[idx]
        if r == 0:
            return True
        t2 = dict(th)
        t2[ctrl] = th.get(ctrl, 0) + 1
        s = (quick.six(quick.ev(t2))[idx] - r) % P
        if s == 0:
            return False
        th[ctrl] = (th.get(ctrl, 0) + (-r) * pow(s, -1, P)) % P
    return quick.six(quick.ev(th))[idx] == 0


def closure(th):
    solve1(th, 5096, 2)     # x25614  (linear)
    solve1(th, 21589, 3)    # x34220  (linear)
    solve1(th, 14515, 4)    # n-gap
    solve1(th, 19750, 5)    # m-gap
    return th


def lagrange(pts, p):
    """interpolate polynomial through (x,y) points over GF(p)"""
    n = len(pts)
    coeff = [0] * n
    for i, (xi, yi) in enumerate(pts):
        num = [1]
        den = 1
        for j, (xj, _) in enumerate(pts):
            if i == j:
                continue
            num = polyroot.pmul(num, [(-xj) % p, 1], p)
            den = den * ((xi - xj) % p) % p
        f = yi * pow(den, -1, p) % p
        for k, c in enumerate(num):
            coeff[k] = (coeff[k] + c * f) % p
    return polyroot.norm(coeff, p)


def blockA_solve(th, rng):
    """set x33708 (cubic root) and x31339 so that x3719 = x25118 = 0"""
    # G(t) = x23776 * x2401^2 - x26196^2   as a function of t = x33708
    pts = []
    for t in range(5):
        t2 = dict(th)
        t2[33708] = t
        v = quick.ev(t2)
        G = (v[23776] * v[2401] * v[2401] - v[26196] * v[26196]) % P
        pts.append((t, G))
    f = lagrange(pts, P)
    if len(f) < 2:
        return None
    rs = polyroot.roots(f, P, rng)
    for r in rs:
        cand = dict(th)
        cand[33708] = r
        # solve x25118 = 0 for x31339 (linear)
        if not solve1(cand, 31339, 1):
            continue
        six = quick.six(quick.ev(cand))
        if six[0] == 0 and six[1] == 0:
            return cand
    return None


def attempt(seed):
    rng = random.Random(seed)
    th = {c: rng.randrange(1, 1 << 80) for c in (16441,)}
    th[33708] = 0
    th[31339] = 0
    closure(th)
    got = blockA_solve(th, rng)
    if got is None:
        return None
    closure(got)
    got2 = blockA_solve(got, rng)
    if got2 is None:
        return None
    closure(got2)
    if not any(quick.six(quick.ev(got2))):
        return got2
    return None


if __name__ == '__main__':
    t0 = time.time()
    for seed in range(300):
        r = attempt(seed)
        if r is not None:
            print(f"seed {seed}: ALL SIX TARGETS ZERO  ({time.time()-t0:.0f}s)")
            print("  theta:", {k: str(x)[:22] for k, x in r.items()})
            json.dump({str(k): x for k, x in r.items()}, open('theta_solveB.json', 'w'))
            break
        if seed % 25 == 0:
            print(f"  seed{seed} no ({time.time()-t0:.0f}s)", flush=True)
    else:
        print("none converged")
