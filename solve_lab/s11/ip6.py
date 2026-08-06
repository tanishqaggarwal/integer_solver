"""IP #6 -- the strong form of the integer program at the checkpoint.

Same objective (minimise violated equations) but with a much larger variable set: every
variable touching ANY equation in the constrained region, not just the failing ones.  Extra
variables cannot change the failing equations directly, but they enlarge the kernel of
'keep the satisfied equations satisfied', which is what buys freedom.

Deltas are computed incrementally (only the atoms containing the perturbed variable change),
so thousands of candidates are affordable.
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


def load_raw(path):
    v = [0] * L.NVARS
    for k, x in json.load(open(path)).items():
        v[int(k[2:]) if k.startswith('x_') else int(k)] = int(x)
    return v


def run(v, maxallow=5, cap=6000, verbose=True):
    AV = [atomval(a, v) for a in range(L.NA)]
    FAIL = [e for e in range(L.NEQ)
            if sum(c * AV[a] for a, c in L.eq_atoms[e][2].items()) != 0]
    print(f"  failing={len(FAIL)} score={L.NEQ-len(FAIL)}  {FAIL}")
    if not FAIL:
        return v, 0
    # region: equations reachable from the failing ones through shared atoms/variables
    core_vars = set()
    for e in FAIL:
        for a in L.eq_atoms[e][2]:
            core_vars |= set(L.avars[a])
    E = set(FAIL)
    for u in core_vars:
        for a in L.var_atoms[u]:
            E |= set(L.atom2eq.get(a, {}))
    E = sorted(E)
    idx = {e: i for i, e in enumerate(E)}
    # candidates: EVERY variable occurring in any atom of any equation of E
    cands = set()
    for e in E:
        for a in L.eq_atoms[e][2]:
            cands |= set(L.avars[a])
    cands = sorted(cands)[:cap]
    if verbose:
        print(f"  region {len(E)} equations ; {len(cands)} candidate variables")

    def delta(u, step):
        old = v[u]
        v[u] = old + step
        d = collections.defaultdict(int)
        for a in L.var_atoms[u]:
            dv = atomval(a, v) - AV[a]
            if dv:
                for e, c in L.atom2eq.get(a, {}).items():
                    if e in idx:
                        d[idx[e]] += c * dv
        v[u] = old
        return d

    cols, used = [], []
    for u in cands:
        d1 = delta(u, 1)
        if not d1:
            continue
        d2 = delta(u, 2)
        if all(d2.get(i, 0) == 2 * d1.get(i, 0) for i in set(d1) | set(d2)):
            cols.append(d1)
            used.append(u)
    if verbose:
        print(f"  exact-linear variables: {len(used)}")
    if not used:
        return v, len(FAIL)
    FS = set(FAIL)
    KEEP = [e for e in E if e not in FS]
    Gk = [[cols[j].get(idx[e], 0) for j in range(len(used))] for e in KEEP]
    t0 = time.time()
    ker = int_kernel(Gk)
    if verbose:
        print(f"  kernel dim over {len(KEEP)} kept equations = {len(ker)}  ({time.time()-t0:.0f}s)")
    if not ker:
        return v, len(FAIL)
    Gf = [[sum(cols[j].get(idx[e], 0) * ker[t][j] for j in range(len(used)))
           for t in range(len(ker))] for e in FAIL]
    bf = [sum(c * AV[a] for a, c in L.eq_atoms[e][2].items()) for e in FAIL]
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
            AV2 = [atomval(a, v) for a in range(L.NA)]
            f2 = [e for e in range(L.NEQ)
                  if sum(c * AV2[a] for a, c in L.eq_atoms[e][2].items()) != 0]
            print(f"  model allows {allow}; APPLIED -> failing={len(f2)} score={L.NEQ-len(f2)}",
                  flush=True)
            if len(f2) <= len(FAIL):
                return v, len(f2)
            for j, u in enumerate(used):
                v[u] = snap[j]
    return v, len(FAIL)


if __name__ == '__main__':
    LAB = os.path.join(HERE, '..')
    src = sys.argv[1] if len(sys.argv) > 1 else os.path.join(LAB, 'best', 'new_instance_partial_39026.json')
    v = load_raw(src)
    print(f"=== {os.path.basename(src)}  (raw integer program, no circuit orientation)")
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
                  open(os.path.join(HERE, 'data', 'ip6_best_named.json'), 'w'))
        print("saved data/ip6_best_named.json")
