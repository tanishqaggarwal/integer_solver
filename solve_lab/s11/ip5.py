"""IP #5 -- the instance AS AN INTEGER PROGRAM, with no circuit orientation at all.

    38,748 integer variables, 39,033 polynomial equations.
    At a point, perturbing variables gives an exact integer lattice of equation-value changes.
    Objective: minimise the number of violated equations.

No forward evaluation: the state is just the raw assignment, so this applies to the 39,026
checkpoint (whose score depends on five DELIBERATELY BROKEN gates and would be destroyed by
forward-evaluating).
"""
import sys, os, json, itertools, time, collections
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 's9', 'eff'))
import lib as L
from zsolve import solve_int
from ip3 import int_kernel
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


def eqsum(e, v):
    return sum(c * atomval(a, v) for a, c in L.eq_atoms[e][2].items())


def failing(v):
    return [e for e in range(L.NEQ) if eqsum(e, v) != 0]


def load_raw(path):
    v = [0] * L.NVARS
    for k, x in json.load(open(path)).items():
        v[int(k[2:]) if k.startswith('x_') else int(k)] = int(x)
    return v


def solve_at(v, maxallow=6, verbose=True):
    FAIL = failing(v)
    if verbose:
        print(f"  failing = {len(FAIL)}  score = {L.NEQ-len(FAIL)}   {FAIL}")
    if not FAIL:
        return v, 0
    # every variable occurring in any atom of any failing equation
    cands = set()
    for e in FAIL:
        for a in L.eq_atoms[e][2]:
            cands |= set(L.avars[a])
    cands = sorted(cands)
    # equations any candidate can touch
    E = set(FAIL)
    for u in cands:
        for a in L.var_atoms[u]:
            E |= set(L.atom2eq.get(a, {}))
    E = sorted(E)
    if verbose:
        print(f"  candidate variables = {len(cands)} ; touched equations = {len(E)}")
    base = [eqsum(e, v) for e in E]
    idx = {e: i for i, e in enumerate(E)}
    cols, used = [], []
    for u in cands:
        old = v[u]
        v[u] = old + 1
        d1 = [eqsum(e, v) - base[idx[e]] for e in E]
        v[u] = old + 2
        d2 = [eqsum(e, v) - base[idx[e]] for e in E]
        v[u] = old
        if any(d1) and all(d2[i] == 2 * d1[i] for i in range(len(E))):
            cols.append(d1)
            used.append(u)
    if verbose:
        print(f"  variables with an EXACT LINEAR effect: {len(used)}")
    if not used:
        return v, len(FAIL)
    KEEP = [e for e in E if e not in set(FAIL)]
    Gk = [[cols[j][idx[e]] for j in range(len(used))] for e in KEEP]
    ker = int_kernel(Gk) if Gk else [[1 if i == j else 0 for j in range(len(used))]
                                     for i in range(len(used))]
    if verbose:
        print(f"  kernel dim (preserving all {len(KEEP)} satisfied equations) = {len(ker)}")
    if not ker:
        return v, len(FAIL)
    Gf = [[sum(cols[j][idx[e]] * ker[t][j] for j in range(len(used))) for t in range(len(ker))]
          for e in FAIL]
    bf = [base[idx[e]] for e in FAIL]
    for allow in range(0, min(maxallow, len(FAIL)) + 1):
        for combo in itertools.combinations(range(len(FAIL)), allow):
            keep = [i for i in range(len(FAIL)) if i not in combo]
            x = solve_int([Gf[i] for i in keep], [-bf[i] for i in keep])
            if x is None:
                continue
            kk = [sum(ker[t][j] * x[t] for t in range(len(ker))) for j in range(len(used))]
            snap = [v[u] for u in used]
            for j, u in enumerate(used):
                v[u] += kk[j]
            f2 = failing(v)
            if verbose:
                print(f"  model allows {allow}; APPLIED -> failing={len(f2)} score={L.NEQ-len(f2)}")
            if len(f2) <= len(FAIL):
                return v, len(f2)
            for j, u in enumerate(used):
                v[u] = snap[j]
    return v, len(FAIL)


if __name__ == '__main__':
    LAB = os.path.join(HERE, '..')
    src = sys.argv[1] if len(sys.argv) > 1 else os.path.join(LAB, 'best', 'new_instance_partial_39026.json')
    v = load_raw(src)
    print(f"=== {os.path.basename(src)} (raw, no forward-eval)")
    t0 = time.time()
    best = None
    for it in range(12):
        v, f = solve_at(v)
        print(f"  it{it}: failing={f} score={L.NEQ-f} ({time.time()-t0:.0f}s)", flush=True)
        if best is None or f < best[0]:
            best = (f, [x for x in v])
        if f == 0:
            break
        if it and f == prev:
            break
        prev = f
    print(f"BEST failing={best[0]} score={L.NEQ-best[0]}")
    if best[0] < 7:
        json.dump({('x_%d' % i): best[1][i] for i in range(L.NVARS)},
                  open(os.path.join(HERE, 'data', 'ip5_best_named.json'), 'w'))
        print("saved data/ip5_best_named.json")
