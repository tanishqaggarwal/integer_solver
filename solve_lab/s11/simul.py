"""Simultaneous exact-integer repair of all bad checks."""
import sys, os, json, time
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 's9', 'eff'))
import lib as L, fw, deep
from zsolve import solve_int
P = L.P
NAT = {u: len(L.var_atoms[u]) for u in range(L.NVARS)}


def simul(v, locked, rounds=12, verbose=True):
    t0 = time.time()
    best = (len(fw.bad_checks(v)), [x for x in v])
    for rnd in range(rounds):
        bad = fw.bad_checks(v)
        if not bad:
            break
        H = {}
        ctrl = []
        for a in bad:
            try:
                hs, base = deep.handles(v, a, locked=locked)
            except Exception:
                hs, base = [], fw.evalpoly(L.polys[a], v)
            H[a] = dict(hs)
            for t, _ in hs:
                if t not in ctrl:
                    ctrl.append(t)
        if not ctrl:
            break
        M = [[H[a].get(c, 0) for c in ctrl] for a in bad]
        rhs = [-fw.evalpoly(L.polys[a], v) for a in bad]
        d = solve_int(M, rhs)
        if verbose:
            print(f"  round{rnd}: bad={len(bad)} ctrl={len(ctrl)} solved={d is not None} ({time.time()-t0:.0f}s)", flush=True)
        if d is None:
            # fall back: solve the largest solvable subset (drop rows greedily)
            order = sorted(range(len(bad)), key=lambda i: -sum(1 for x in M[i] if x))
            keep = list(range(len(bad)))
            for i in order:
                trial = [j for j in keep if j != i]
                if not trial:
                    break
                if solve_int([M[j] for j in trial], [rhs[j] for j in trial]) is not None:
                    keep = trial
                    break
            d = solve_int([M[j] for j in keep], [rhs[j] for j in keep])
            if d is None:
                break
        for j, c in enumerate(ctrl):
            v[c] += d[j]
        fw.forward(v)
        nb = fw.bad_checks(v)
        if len(nb) < best[0]:
            best = (len(nb), [x for x in v])
        if verbose:
            print(f"    -> bad={len(nb)} {nb[:14]}", flush=True)
    if len(fw.bad_checks(v)) > best[0]:
        v[:] = best[1]
        fw.forward(v)
    return v, fw.bad_checks(v)


if __name__ == '__main__':
    BITS = (542, 47, 438, 91)
    C0 = L.polys[688][()]
    MM = 8863713
    G0 = (-C0 * pow(MM, -1, P)) % P
    C0B = L.polys[1618][()]
    TH = {int(k): x for k, x in json.load(open('theta_solveB.json')).items()}
    v = [0] * L.NVARS
    for b in BITS:
        v[b] = 1
    for k, x in TH.items():
        v[k] = x
    fw.forward(v)
    v[14853] = v[12186]
    v[16742] = v[24908]
    v[30213] = G0
    v[22820] = 0
    v[7497] = (C0 + MM * G0) // P
    v[22162] = -C0B
    v[14393] = 0
    v[11436] = 0
    fw.forward(v)
    LOCK = set(BITS) | set(TH) | {14853, 16742, 30213, 22820, 7497, 22162, 14393, 11436}
    v, b = simul(v, LOCK)
    av = L.all_atom_values(v)
    f = L.failing_eqs(av)
    print(f"FINAL bad={len(b)} failing={len(f)} score={L.NEQ-len(f)}")
    print("bad:", b)
    json.dump({str(i): v[i] for i in range(L.NVARS)}, open('simul.json', 'w'))
