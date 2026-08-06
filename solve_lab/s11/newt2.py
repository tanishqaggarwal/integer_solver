import sys, os, json, time, random
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 's9', 'eff'))
import lib as L, fw, fast
from gfp import gauss_solve
P = L.P
BITS = fast.BITS
# free parameters for the four group targets
CTRL = [16441, 22917, 31339, 33708, 13222, 14681, 28486, 38667]
NAMES4 = ['x3719', 'x25118', 'x25614', 'x34220']


def light2(th):
    """set controls, forward, then derive x14515 / x21589 so n-gap = m-gap = 0."""
    v = [0] * L.NVARS
    for b in BITS:
        v[b] = 1
    for k, x in th.items():
        v[k] = x
    fw.forward(v)
    for _ in range(3):
        # n-gap: x12186 - x1308 ; x1308 has slope +1 in x14515
        v[14515] = (v[14515] + (v[12186] - v[1308])) % P
        v[21589] = (v[21589] + (v[24908] - v[19083])) % P
        fw.forward(v)
    return v


def t4(v):
    return [v[3719] % P, v[25118] % P, v[25614] % P, v[34220] % P]


def gaps(v):
    return (v[12186] - v[1308]) % P, (v[24908] - v[19083]) % P


def newton(seed, iters=30, verbose=True):
    rnd = random.Random(seed)
    th = {c: rnd.randrange(1, 1 << 60) for c in CTRL}
    for it in range(iters):
        v = light2(th)
        r = t4(v)
        if not any(r) and not any(gaps(v)):
            return th, True
        J = [[0] * len(CTRL) for _ in range(4)]
        for j, c in enumerate(CTRL):
            t2 = dict(th)
            t2[c] = th[c] + 1
            r1 = t4(light2(t2))
            for i in range(4):
                J[i][j] = (r1[i] - r[i]) % P
        d = gauss_solve(J, [(-x) % P for x in r], P)
        if d is None:
            if verbose:
                print(f"  seed{seed} it{it}: inconsistent", flush=True)
            return th, False
        for j, c in enumerate(CTRL):
            th[c] = (th[c] + d[j]) % P
        if verbose:
            print(f"  seed{seed} it{it}: nz={[NAMES4[i] for i,x in enumerate(r) if x]} gaps={gaps(v)!=(0,0)}", flush=True)
    v = light2(th)
    return th, (not any(t4(v))) and gaps(v) == (0, 0)


if __name__ == '__main__':
    t0 = time.time()
    for seed in range(10):
        th, ok = newton(seed)
        print(f"seed {seed}: converged={ok} ({time.time()-t0:.0f}s)", flush=True)
        if ok:
            v = light2(th)
            th2 = dict(th)
            th2[14515] = v[14515]
            th2[21589] = v[21589]
            json.dump({str(k): x for k, x in th2.items()}, open('theta_n2.json', 'w'))
            print("  targets:", t4(v), "gaps:", gaps(v))
            b = fw.bad_checks(v)
            print(f"  bad_checks(pre-repair)={len(b)}: {b}")
            break
