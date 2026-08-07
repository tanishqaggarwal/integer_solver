"""Joint 6x6 Newton over the six residuals with the six controls freed by breaking
   a41332 [1 eq] and a36244 [4 eqs].  If it closes, only those 5 equations fail -> 39,028."""
import sys, os, json, time, random
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 's9', 'eff'))
import lib as L, fw, deep, uvskip3
from gfp import gauss_solve
P = L.P
HERE = os.path.dirname(os.path.abspath(__file__))
NAT = {u: len(L.var_atoms[u]) for u in range(L.NVARS)}
sys.set_int_max_str_digits(300000)
fwdskip = uvskip3.fwdskip
CT = [3432, 24453, 33708, 31339, 33129, 37088]
GET = [lambda v: v[25118], lambda v: v[3719],
       lambda v: v[33708] - v[10170], lambda v: v[31339] - v[6858],
       lambda v: v[33129] - v[3757], lambda v: v[37088] - v[13585]]
NM = ['x25118', 'x3719', 'a34580', 'a33796', 'a14445', 'a27139']


def res(v):
    return [g(v) % P for g in GET]


def newton(v, iters=30, verbose=True):
    for it in range(iters):
        r = res(v)
        if not any(r):
            return True
        J = [[0] * len(CT) for _ in range(6)]
        for j, c in enumerate(CT):
            old = v[c]
            v[c] = old + 1
            fwdskip(v)
            r1 = res(v)
            v[c] = old
            fwdskip(v)
            for i in range(6):
                J[i][j] = (r1[i] - r[i]) % P
        d = gauss_solve(J, [(-x) % P for x in r], P)
        if d is None:
            if verbose:
                print(f"    it{it}: singular/inconsistent")
            return False
        for j, c in enumerate(CT):
            v[c] = (v[c] + d[j]) % P
        fwdskip(v)
        if verbose and it % 5 == 0:
            print(f"    it{it}: nz={[NM[i] for i,x in enumerate(r) if x]}", flush=True)
    return not any(res(v))


base = [0] * L.NVARS
for k, x in json.load(open(os.path.join(HERE, 'data', 'closehit2.json'))).items():
    base[int(k)] = int(x)

t0 = time.time()
best = None
for seed in range(12):
    v = base[:]
    if seed:
        rnd = random.Random(seed)
        for c in CT:
            v[c] = rnd.randrange(0, 1 << 60)
    fwdskip(v)
    ok = newton(v, verbose=(seed == 0))
    f = L.failing_eqs(L.all_atom_values(v))
    print(f"seed{seed}: all6={ok} failing={len(f)} score={L.NEQ-len(f)} ({time.time()-t0:.0f}s)", flush=True)
    if ok:
        LOCK = {490, 91, 19750, 7497, 22820, 14853, 14393, 11436, 14515, 16742,
                22162, 30213, 8386, 21868, 16441, 28955, 2751, 18751} | set(CT)
        for rnd2 in range(10):
            bad = [a for a in range(L.NA) if L.atom_out.get(a) is None
                   and fw.evalpoly(L.polys[a], v) != 0]
            if not bad:
                break
            prog = False
            for a in sorted(bad, key=lambda a: len(L.atom2eq.get(a, {}))):
                if fw.evalpoly(L.polys[a], v) == 0:
                    continue
                cs = [(u, None) for u in L.avars[a]
                      if L.definer.get(u) is None and u not in LOCK
                      and not any(mm.count(u) > 1 for mm in L.polys[a])]
                cs.sort(key=lambda kv: (NAT[kv[0]], kv[0]))
                try:
                    hs, bb = deep.handles(v, a, locked=LOCK)
                    cs += [(t, d) for t, d in sorted(hs, key=lambda kv: (NAT[kv[0]], kv[0]))]
                except Exception:
                    pass
                for t, d in cs:
                    old = v[t]
                    if d is None:
                        x = fw.solve_lin(a, t, v)
                        if x is None or x == old:
                            continue
                    else:
                        bs = fw.evalpoly(L.polys[a], v)
                        if not d or bs % d:
                            continue
                        x = old - bs // d
                    v[t] = x
                    fwdskip(v)
                    if fw.evalpoly(L.polys[a], v) == 0:
                        prog = True
                        break
                    v[t] = old
                    fwdskip(v)
            f = L.failing_eqs(L.all_atom_values(v))
            print(f"    close{rnd2}: failing={len(f)} score={L.NEQ-len(f)}", flush=True)
            if best is None or len(f) < best[0]:
                best = (len(f), [x for x in v])
            if not prog:
                break
        break
if best:
    print(f"BEST failing={best[0]} score={L.NEQ-best[0]}")
    json.dump({('x_%d' % i): best[1][i] for i in range(L.NVARS)},
              open(os.path.join(HERE, 'data', 'joint6_named.json'), 'w'))
