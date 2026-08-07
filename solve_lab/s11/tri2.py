import sys, os, json, time, random
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 's9', 'eff'))
import lib as L, fw, fast
from gfp import gauss_solve
P = L.P
BITS = fast.BITS
A_CTRL = [16441, 22917, 31339, 33708]


def raw(th):
    v = [0] * L.NVARS
    for b in BITS:
        v[b] = 1
    for k, x in th.items():
        v[k] = x
    fw.forward(v)
    return v


def solve1(th, ctrl, get, tries=3):
    """solve get(v)==0 for a single control assumed linear in it."""
    for _ in range(tries):
        v = raw(th)
        r = get(v) % P
        if r == 0:
            return True
        t2 = dict(th)
        t2[ctrl] = th.get(ctrl, 0) + 1
        s = (get(raw(t2)) - get(v)) % P
        if s == 0:
            return False
        th[ctrl] = (th.get(ctrl, 0) + (-r) * pow(s, -1, P)) % P
    return get(raw(th)) % P == 0


def closure(th):
    """triangular: x5096 <- x25614 ; x21589 <- x34220 ; x14515 <- n-gap ; x19750 <- m-gap"""
    ok = True
    ok &= solve1(th, 5096, lambda v: v[25614])
    ok &= solve1(th, 21589, lambda v: v[34220])
    ok &= solve1(th, 14515, lambda v: v[12186] - v[1308])
    ok &= solve1(th, 19750, lambda v: v[24908] - v[19083])
    return ok


def blkA(v): return [v[3719] % P, v[25118] % P]


def full(th):
    th = dict(th)
    closure(th)
    return th, raw(th)


def run(seed, iters=30):
    rnd = random.Random(seed)
    th = {c: rnd.randrange(1, 1 << 80) for c in A_CTRL}
    for it in range(iters):
        th, v = full(th)
        r = blkA(v)
        chk = (v[25614] % P, v[34220] % P, (v[12186]-v[1308]) % P, (v[24908]-v[19083]) % P)
        if not any(r) and not any(chk):
            return th, True
        J = [[0] * len(A_CTRL) for _ in r]
        for j, c in enumerate(A_CTRL):
            t2 = dict(th)
            t2[c] = th[c] + 1
            _, v1 = full(t2)
            r1 = blkA(v1)
            for i in range(len(r)):
                J[i][j] = (r1[i] - r[i]) % P
        d = gauss_solve(J, [(-x) % P for x in r], P)
        if d is None:
            return th, False
        for j, c in enumerate(A_CTRL):
            th[c] = (th[c] + d[j]) % P
    th, v = full(th)
    return th, not any(blkA(v))


if __name__ == '__main__':
    t0 = time.time()
    for seed in range(30):
        th, ok = run(seed)
        th, v = full(th)
        chk = (v[3719] % P, v[25118] % P, v[25614] % P, v[34220] % P,
               (v[12186]-v[1308]) % P, (v[24908]-v[19083]) % P)
        print(f"seed{seed}: ok={ok} residuals_zero={not any(chk)} ({time.time()-t0:.0f}s)", flush=True)
        if not any(chk):
            json.dump({str(k): x for k, x in th.items()}, open('theta_tri2.json', 'w'))
            b = fw.bad_checks(v)
            print(f"  ALL SIX TARGETS ZERO. bad_checks(pre-repair)={len(b)}: {b}")
            break
