import sys, os, json, time, random
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 's9', 'eff'))
import lib as L, fw, fast
from gfp import gauss_solve
P = L.P
CTRL = [14515, 16441, 22917, 31339, 33708, 19750, 5096, 13222, 14681, 28486, 38667, 21589]


def newton(seed, iters=25, verbose=True):
    rnd = random.Random(seed)
    th = {c: rnd.randrange(1, 1 << 60) for c in CTRL}
    for it in range(iters):
        r = fast.targets(fast.light(th))
        if not any(r):
            return th, True
        J = [[0] * len(CTRL) for _ in range(6)]
        for j, c in enumerate(CTRL):
            t2 = dict(th)
            t2[c] = th[c] + 1
            r1 = fast.targets(fast.light(t2))
            for i in range(6):
                J[i][j] = (r1[i] - r[i]) % P
        d = gauss_solve(J, [(-x) % P for x in r], P)
        if d is None:
            if verbose:
                print(f"  seed{seed} it{it}: inconsistent")
            return th, False
        for j, c in enumerate(CTRL):
            th[c] = (th[c] + d[j]) % P
        if verbose:
            print(f"  seed{seed} it{it}: nonzero={[fast.NAMES[i] for i,x in enumerate(r) if x]}", flush=True)
    return th, not any(fast.targets(fast.light(th)))


if __name__ == '__main__':
    t0 = time.time()
    for seed in range(6):
        th, ok = newton(seed)
        print(f"seed {seed}: converged={ok} ({time.time()-t0:.0f}s)")
        if ok:
            json.dump({str(k): v for k, v in th.items()}, open('theta_newton.json', 'w'))
            v = fast.light(th)
            print("  targets:", fast.targets(v))
            b = fw.bad_checks(v)
            print(f"  bad_checks (pre-repair) = {len(b)}")
            break
