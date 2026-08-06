"""IP #7 -- the checkpoint integer program, reduced by locality.

Only 56 variables can move the 7 failing equations.  Only the equations those 56 actually
disturb (the COLLATERAL) need compensating, and only variables touching the collateral can
compensate.  That shrinks the system enormously versus IP #6's 913 x 1463.

    minimise |violated|  s.t.  G_core*k_core + G_comp*k_comp = (-b on FAIL, 0 on COLLATERAL)
"""
import sys, os, json, itertools, time, collections
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 's9', 'eff'))
import lib as L
from zsolve import solve_int
P = L.P
HERE = os.path.dirname(os.path.abspath(__file__))
sys.set_int_max_str_digits(400000)


def atomval(a, v):
    s = 0
    for m, c in L.polys[a].items():
        t = c
        for u in m:
            t *= v[u]
        s += t
    return s


def load_raw(path):
    v = [0] * L.NVARS
    for k, x in json.load(open(path)).items():
        v[int(k[2:]) if k.startswith('x_') else int(k)] = int(x)
    return v


def deltas(v, AV, u, step=1):
    old = v[u]
    v[u] = old + step
    d = collections.defaultdict(int)
    for a in L.var_atoms[u]:
        dv = atomval(a, v) - AV[a]
        if dv:
            for e, c in L.atom2eq.get(a, {}).items():
                d[e] += c * dv
    v[u] = old
    return d


def run(v, maxallow=6):
    AV = [atomval(a, v) for a in range(L.NA)]
    def eqs(e): return sum(c * AV[a] for a, c in L.eq_atoms[e][2].items())
    FAIL = [e for e in range(L.NEQ) if eqs(e) != 0]
    print(f"  failing={len(FAIL)} score={L.NEQ-len(FAIL)}")
    if not FAIL:
        return v, 0
    core = set()
    for e in FAIL:
        for a in L.eq_atoms[e][2]:
            core |= set(L.avars[a])
    core = sorted(core)
    dc, uc = [], []
    collat = set()
    for u in core:
        d1 = deltas(v, AV, u, 1)
        if not d1:
            continue
        d2 = deltas(v, AV, u, 2)
        if not all(d2.get(e, 0) == 2 * d1.get(e, 0) for e in set(d1) | set(d2)):
            continue
        dc.append(d1)
        uc.append(u)
        collat |= {e for e in d1 if e not in set(FAIL)}
    print(f"  core vars with exact linear effect: {len(uc)} ; collateral equations: {len(collat)}")
    comp = set()
    for e in collat:
        for a in L.eq_atoms[e][2]:
            comp |= set(L.avars[a])
    comp -= set(uc)
    comp = sorted(comp)
    print(f"  compensator candidates: {len(comp)}")
    dcomp, ucomp = [], []
    for u in comp:
        d1 = deltas(v, AV, u, 1)
        if not d1 or not (set(d1) & collat):
            continue
        d2 = deltas(v, AV, u, 2)
        if not all(d2.get(e, 0) == 2 * d1.get(e, 0) for e in set(d1) | set(d2)):
            continue
        # a compensator must not disturb anything outside FAIL u COLLATERAL
        if any(e not in collat and e not in set(FAIL) for e in d1):
            continue
        dcomp.append(d1)
        ucomp.append(u)
    print(f"  usable compensators (no new collateral): {len(ucomp)}")
    used = uc + ucomp
    cols = dc + dcomp
    ROWS = list(FAIL) + sorted(collat)
    ridx = {e: i for i, e in enumerate(ROWS)}
    M = [[cols[j].get(e, 0) for j in range(len(cols))] for e in ROWS]
    rhs = [-eqs(e) for e in FAIL] + [0] * len(collat)
    print(f"  system: {len(ROWS)} x {len(cols)}")
    t0 = time.time()
    for allow in range(0, min(maxallow, len(FAIL)) + 1):
        for combo in itertools.combinations(range(len(FAIL)), allow):
            drop = {ridx[FAIL[i]] for i in combo}
            keep = [i for i in range(len(ROWS)) if i not in drop]
            x = solve_int([M[i] for i in keep], [rhs[i] for i in keep])
            if x is None:
                continue
            snap = [v[u] for u in used]
            for j, u in enumerate(used):
                v[u] += x[j]
            AV2 = [atomval(a, v) for a in range(L.NA)]
            f2 = [e for e in range(L.NEQ)
                  if sum(c * AV2[a] for a, c in L.eq_atoms[e][2].items()) != 0]
            print(f"  model allows {allow}; APPLIED -> failing={len(f2)} score={L.NEQ-len(f2)}"
                  f" ({time.time()-t0:.0f}s)", flush=True)
            if len(f2) <= len(FAIL):
                return v, len(f2)
            for j, u in enumerate(used):
                v[u] = snap[j]
    print(f"  no solution with <= {maxallow} allowed failures ({time.time()-t0:.0f}s)")
    return v, len(FAIL)


if __name__ == '__main__':
    LAB = os.path.join(HERE, '..')
    src = sys.argv[1] if len(sys.argv) > 1 else os.path.join(LAB, 'best', 'new_instance_partial_39026.json')
    v = load_raw(src)
    print(f"=== {os.path.basename(src)} (locality-reduced integer program)")
    best = None
    prev = None
    for it in range(8):
        v, f = run(v)
        print(f"  it{it}: failing={f} score={L.NEQ-f}", flush=True)
        if best is None or f < best[0]:
            best = (f, [x for x in v])
        if f == 0 or f == prev:
            break
        prev = f
    print(f"BEST failing={best[0]} score={L.NEQ-best[0]}")
    if best[0] < 7:
        json.dump({('x_%d' % i): best[1][i] for i in range(L.NVARS)},
                  open(os.path.join(HERE, 'data', 'ip7_best_named.json'), 'w'))
        print("saved data/ip7_best_named.json")
