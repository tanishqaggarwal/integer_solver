import sys, os, json, time, random
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 's9', 'eff'))
import lib as L, fw, sys6, polyroot
from gfp import gauss_solve
P = L.P
HERE = os.path.dirname(os.path.abspath(__file__))
LD = json.load(open(os.path.join(HERE, 'data', 'loads.json')))['loads']
BITSET = set(int(b) for b in LD)
S = json.load(open(os.path.join(HERE, 'data', 'sys6.json')))
NAMES = sys6.NAMES
POOL = sorted({u for nm in NAMES for u in S[nm] if u not in BITSET})


def resid(th):
    return sys6.six(sys6.ev(th))


def jac(th, r, ctrl):
    J = [[0] * len(ctrl) for _ in range(6)]
    for j, c in enumerate(ctrl):
        t2 = dict(th)
        t2[c] = th.get(c, sys6.BASE[c]) + 1
        r1 = resid(t2)
        for i in range(6):
            J[i][j] = (r1[i] - r[i]) % P
    return J


def newton(th, ctrl, iters=30):
    for it in range(iters):
        r = resid(th)
        if not any(r):
            return th, True
        J = jac(th, r, ctrl)
        d = gauss_solve(J, [(-x) % P for x in r], P)
        if d is None:
            return th, False
        for j, c in enumerate(ctrl):
            th[c] = (th.get(c, sys6.BASE[c]) + d[j]) % P
    return th, not any(resid(th))


if __name__ == '__main__':
    t0 = time.time()
    got = None
    for seed in range(200):
        rnd = random.Random(seed)
        th = {c: sys6.BASE[c] + rnd.randrange(0, 1 << 60) for c in POOL}
        th, ok = newton(th, POOL)
        r = resid(th)
        if ok:
            print(f"seed{seed}: ALL SIX RESIDUALS ZERO ({time.time()-t0:.0f}s)")
            got = th
            break
        if seed % 20 == 0:
            print(f"  seed{seed}: nz={[NAMES[i] for i,x in enumerate(r) if x]} ({time.time()-t0:.0f}s)", flush=True)
    if got is None:
        print("no seed converged")
    else:
        json.dump({str(k): x for k, x in got.items()}, open(os.path.join(HERE, 'data', 'sys6b_theta.json'), 'w'))
        v = sys6.BASE[:]
        for k, x in got.items():
            v[k] = x
        fw.forward(v)
        print("  six on full forward:", sys6.six(v))
        bad = fw.bad_checks(v)
        f = L.failing_eqs(L.all_atom_values(v))
        print(f"  full state: bad={len(bad)} failing={len(f)} score={L.NEQ-len(f)}")
        print("  bad:", bad)
        json.dump({str(i): v[i] for i in range(L.NVARS)}, open(os.path.join(HERE, 'data', 'sys6b_state.json'), 'w'))
