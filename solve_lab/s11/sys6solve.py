import sys, os, json, time, random
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 's9', 'eff'))
import lib as L, fw, sys6, uv01full
from gfp import gauss_solve
P = L.P
HERE = os.path.dirname(os.path.abspath(__file__))
LD = json.load(open(os.path.join(HERE, 'data', 'loads.json')))['loads']
BITSET = set(int(b) for b in LD)
S = json.load(open(os.path.join(HERE, 'data', 'sys6.json')))
NAMES = sys6.NAMES

# controls: union of scan hits, excluding message bits; prefer vars appearing in few atoms
NAT = {u: len(L.var_atoms[u]) for u in range(L.NVARS)}
pool = []
for nm in NAMES:
    for u in S[nm]:
        if u not in BITSET and u not in pool:
            pool.append(u)
pool.sort(key=lambda u: (NAT[u], u))
print(f"control pool (non-bit): {len(pool)}")
CTRL = pool[:24]
print("using:", CTRL)


def resid(th):
    return sys6.six(sys6.ev(th))


def newton(th, iters=25, verbose=True):
    for it in range(iters):
        r = resid(th)
        if not any(r):
            return th, True
        J = [[0] * len(CTRL) for _ in range(6)]
        for j, c in enumerate(CTRL):
            t2 = dict(th)
            t2[c] = th.get(c, sys6.BASE[c]) + 1
            r1 = resid(t2)
            for i in range(6):
                J[i][j] = (r1[i] - r[i]) % P
        d = gauss_solve(J, [(-x) % P for x in r], P)
        if d is None:
            if verbose:
                print(f"  it{it}: inconsistent linearisation")
            return th, False
        for j, c in enumerate(CTRL):
            th[c] = (th.get(c, sys6.BASE[c]) + d[j]) % P
        if verbose:
            print(f"  it{it}: nz={[NAMES[i] for i,x in enumerate(r) if x]}", flush=True)
    return th, not any(resid(th))


if __name__ == '__main__':
    t0 = time.time()
    for seed in range(40):
        rnd = random.Random(seed)
        th = {c: sys6.BASE[c] + rnd.randrange(0, 1 << 40) for c in CTRL}
        th, ok = newton(th, verbose=(seed == 0))
        r = resid(th)
        print(f"seed{seed}: converged={ok} nz={[NAMES[i] for i,x in enumerate(r) if x]} ({time.time()-t0:.0f}s)", flush=True)
        if ok:
            json.dump({str(k): x for k, x in th.items()}, open(os.path.join(HERE, 'data', 'sys6_theta.json'), 'w'))
            # apply to the FULL state and measure
            v = sys6.BASE[:]
            for k, x in th.items():
                v[k] = x
            fw.forward(v)
            print("  six residuals on full forward:", sys6.six(v))
            bad = fw.bad_checks(v)
            f = L.failing_eqs(L.all_atom_values(v))
            print(f"  full state: bad={len(bad)} failing={len(f)} score={L.NEQ-len(f)} {bad[:15]}")
            json.dump({str(i): v[i] for i in range(L.NVARS)}, open(os.path.join(HERE, 'data', 'sys6_state.json'), 'w'))
            break
