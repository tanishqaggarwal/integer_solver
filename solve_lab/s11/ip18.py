"""IP #18 -- two-stage solve on the LOCALITY-REDUCED system (small enough for exact work).

    ip8.build gives  M (rows = FAIL then COLLATERAL)  x  (cols = core vars + safe compensators)
    Stage A : k in ker_Z(M_collateral)  with  b_fail + M_fail k = 0  (mod p)
    Stage B : exact integer solve of the now p-divisible residual
"""
import sys, os, json, time
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 's9', 'eff'))
import lib as L
from zsolve import solve_int
from ip3 import int_kernel
from ip7 import atomval, load_raw
from ip8 import build
from ip14 import gf_solve
P = L.P
HERE = os.path.dirname(os.path.abspath(__file__))
sys.set_int_max_str_digits(400000)


def run(path, name):
    v = load_raw(path)
    print(f"=== {name}")
    v, FAIL, used, M, rhs, nf = build(v)
    nF = len(FAIL)
    Mf = M[:nF]
    Mc = M[nF:]
    bf = [-r for r in rhs[:nF]]          # bf = current failing values
    t0 = time.time()
    ker = int_kernel(Mc) if Mc else [[1 if i == j else 0 for j in range(len(M[0]))]
                                     for i in range(len(M[0]))]
    print(f"  kernel of the collateral block: dim {len(ker)} of {len(M[0])} "
          f"({time.time()-t0:.0f}s)", flush=True)
    if not ker:
        return
    Gf = [[sum(Mf[i][j] * ker[t][j] for j in range(len(M[0]))) % P for t in range(len(ker))]
          for i in range(nF)]
    y = gf_solve(Gf, [(-b) % P for b in bf], P)
    print(f"  Stage A (failing = 0 mod p inside the kernel): "
          f"{'SOLVABLE' if y is not None else 'NOT solvable'}")
    if y is None:
        # how much of it is reachable?
        for drop in range(1, min(6, nF) + 1):
            import itertools
            hit = None
            for combo in itertools.combinations(range(nF), drop):
                keep = [i for i in range(nF) if i not in combo]
                if gf_solve([Gf[i] for i in keep], [(-bf[i]) % P for i in keep], P) is not None:
                    hit = combo
                    break
            if hit:
                print(f"    ({nF-drop} of {nF} failing values can be made p-divisible)")
                break
        return
    kk = [sum(ker[t][j] * y[t] for t in range(len(ker))) for j in range(len(M[0]))]
    snap = [v[u] for u in used]
    for j, u in enumerate(used):
        v[u] += kk[j]
    AV = [atomval(a, v) for a in range(L.NA)]
    def eqs(e): return sum(c * AV[a] for a, c in L.eq_atoms[e][2].items())
    F2 = [e for e in range(L.NEQ) if eqs(e) != 0]
    pd = sum(1 for e in F2 if eqs(e) % P == 0)
    print(f"  after Stage A: failing={len(F2)} score={L.NEQ-len(F2)} ; "
          f"p-divisible {pd} of {len(F2)}", flush=True)
    if len(F2) > 4 * nF:
        print("  reverting (kernel step not exact)")
        for j, u in enumerate(used):
            v[u] = snap[j]
        return
    v2, F2b, used2, M2, rhs2, nf2 = build(v)
    t1 = time.time()
    x = solve_int(M2, rhs2)
    print(f"  Stage B exact integer solve: {'FOUND' if x else 'none'} ({time.time()-t1:.0f}s)")
    if x:
        for j, u in enumerate(used2):
            v[u] += x[j]
        AV3 = [atomval(a, v) for a in range(L.NA)]
        F3 = [e for e in range(L.NEQ)
              if sum(c * AV3[a] for a, c in L.eq_atoms[e][2].items()) != 0]
        print(f"  FINAL failing={len(F3)} score={L.NEQ-len(F3)}")
        json.dump({('x_%d' % i): v[i] for i in range(L.NVARS)},
                  open(os.path.join(HERE, 'data', f'ip18_{name}.json'), 'w'))


if __name__ == '__main__':
    LAB = os.path.join(HERE, '..')
    for rel, nm in [(os.path.join(LAB, 'best', 'new_instance_partial_39026.json'), 'checkpoint39026'),
                    (os.path.join(HERE, 'data', 'finish3_named.json'), 's11best39018')]:
        try:
            run(rel, nm)
        except Exception as e:
            print(nm, 'error', repr(e))
