"""IP #4 -- iterated minimum-weight coset (integer Newton on the lattice).

IP #3 solves the LINEARISED coset problem; the handles interact, so the exact map is
polynomial.  Iterate: re-linearise at the current point, solve the coset IP, apply, repeat.
Accept a step only if the true failing count does not increase.
"""
import sys, os, json, itertools, time
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 's9', 'eff'))
import lib as L, fw
from zsolve import solve_int
from ip3 import int_kernel, load
P = L.P
HERE = os.path.dirname(os.path.abspath(__file__))
sys.set_int_max_str_digits(400000)


def step(v, maxallow=6, verbose=True):
    av = L.all_atom_values(v)
    FAIL = sorted(L.failing_eqs(av))
    if not FAIL:
        return v, 0, None
    region_atoms = set()
    for e in FAIL:
        region_atoms |= set(L.eq_atoms[e][2])
    cands = sorted({u for a in region_atoms for u in L.avars[a] if L.definer.get(u) is None})
    E = set(FAIL)
    for t in cands:
        for a in L.var_atoms[t]:
            E |= set(L.atom2eq.get(a, {}))
    E = sorted(E)

    def inner(vv):
        a2 = L.all_atom_values(vv)
        return [sum(c * a2[a] for a, c in L.eq_atoms[e][2].items()) for e in E]

    b = inner(v)
    cols, used = [], []
    for t in cands:
        old = v[t]
        v[t] = old + 1
        fw.forward(v)
        e1 = inner(v)
        v[t] = old
        fw.forward(v)
        d1 = [e1[i] - b[i] for i in range(len(E))]
        if any(d1):
            cols.append(d1)
            used.append(t)
    if not used:
        return v, len(FAIL), None
    idx = {e: i for i, e in enumerate(E)}
    KEEP = [e for e in E if e not in set(FAIL)]
    Gk = [[cols[j][idx[e]] for j in range(len(used))] for e in KEEP]
    ker = int_kernel(Gk) if Gk else [[1 if i == j else 0 for j in range(len(used))]
                                     for i in range(len(used))]
    if not ker:
        return v, len(FAIL), None
    Gf = [[sum(cols[j][idx[e]] * ker[t][j] for j in range(len(used))) for t in range(len(ker))]
          for e in FAIL]
    bf = [b[idx[e]] for e in FAIL]
    for allow in range(0, min(maxallow, len(FAIL)) + 1):
        for combo in itertools.combinations(range(len(FAIL)), allow):
            keep = [i for i in range(len(FAIL)) if i not in combo]
            x = solve_int([Gf[i] for i in keep], [-bf[i] for i in keep])
            if x is None:
                continue
            kk = [sum(ker[t][j] * x[t] for t in range(len(ker))) for j in range(len(used))]
            snap = [v[t] for t in used]
            for j, t in enumerate(used):
                v[t] += kk[j]
            fw.forward(v)
            f2 = L.failing_eqs(L.all_atom_values(v))
            if len(f2) <= len(FAIL):
                if verbose:
                    print(f"    model says {allow}; applied -> failing={len(f2)} "
                          f"(region {len(E)} eqs, {len(used)} handles, ker {len(ker)})", flush=True)
                return v, len(f2), allow
            for j, t in enumerate(used):
                v[t] = snap[j]
            fw.forward(v)
    return v, len(FAIL), None


if __name__ == '__main__':
    for nm in ['closehit2.json', 'finish3.json', 'three.json']:
        v = load(nm)
        cur = len(L.failing_eqs(L.all_atom_values(v)))
        print(f"=== {nm}: start failing={cur} score={L.NEQ-cur}", flush=True)
        best = (cur, [x for x in v])
        t0 = time.time()
        for it in range(20):
            v, f, allow = step(v)
            print(f"  it{it}: failing={f} score={L.NEQ-f} ({time.time()-t0:.0f}s)", flush=True)
            if f < best[0]:
                best = (f, [x for x in v])
            if f == 0:
                break
            if allow is None or f == cur:
                break
            cur = f
        print(f"  BEST failing={best[0]} score={L.NEQ-best[0]}")
        if best[0] < 15:
            json.dump({('x_%d' % i): best[1][i] for i in range(L.NVARS)},
                      open(os.path.join(HERE, 'data', f'ip4_{nm}'), 'w'))
            print(f"  saved data/ip4_{nm}")
