import sys, os, json, time, random
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 's9', 'eff'))
import lib as L, fw, quick
from gfp import gauss_solve
P = L.P
A_CTRL = [16441, 22917, 31339, 33708]
NAMES = ['x3719', 'x25118', 'x25614', 'x34220', 'n-gap', 'm-gap']


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
    solve1(th, 5096, 2)     # x25614
    solve1(th, 21589, 3)    # x34220
    solve1(th, 14515, 4)    # n-gap
    solve1(th, 19750, 5)    # m-gap
    return th


def full(th):
    th = dict(th)
    closure(th)
    return th, quick.six(quick.ev(th))


def run(seed, iters=40):
    rnd = random.Random(seed)
    th = {c: rnd.randrange(1, 1 << 80) for c in A_CTRL}
    for it in range(iters):
        th, r = full(th)
        if not any(r):
            return th, True
        J = [[0] * len(A_CTRL) for _ in range(2)]
        for j, c in enumerate(A_CTRL):
            t2 = dict(th)
            t2[c] = th[c] + 1
            _, r1 = full(t2)
            for i in range(2):
                J[i][j] = (r1[i] - r[i]) % P
        d = gauss_solve(J, [(-r[0]) % P, (-r[1]) % P], P)
        if d is None:
            return th, False
        for j, c in enumerate(A_CTRL):
            th[c] = (th[c] + d[j]) % P
    th, r = full(th)
    return th, not any(r)


if __name__ == '__main__':
    t0 = time.time()
    best = None
    for seed in range(400):
        th, ok = run(seed)
        th, r = full(th)
        if not any(r):
            print(f"seed{seed}: ALL SIX TARGETS ZERO ({time.time()-t0:.0f}s)")
            json.dump({str(k): x for k, x in th.items()}, open('theta_tri3.json', 'w'))
            best = th
            break
        if seed % 40 == 0:
            print(f"  seed{seed}: nz={[NAMES[i] for i,x in enumerate(r) if x]} ({time.time()-t0:.0f}s)", flush=True)
    if best is None:
        print("no seed converged")
