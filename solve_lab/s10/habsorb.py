"""S11 step 116: absorb, do not zero -- and solve the handle lattice jointly.

handzero set the handles to zero and lost points.  The diagnosis: with A = B = 0 the
numerators are multiples of p, NOT zero:

    x11150 = 11356049836642346573658703...   ≡ 0 (mod p) but not 0 over Z
    a19297 = x11150*x15298 + x4007 ,  x4007 = p*x30317
           = p * (x11150/p + x30317)

so the handle must ABSORB -- x30317 <- -x11150/p -- not be set to zero.  The generic
lift missed it because it demands a strictly-improving move and absorbing a19297
disturbs the bundles a36185 and a40812 that contain it.

The right object is the joint one: the handles shift their checks by exact multiples
of p, the equations are linear in the check values, so choose the handles TOGETHER to
leave the fewest equations broken.  Enumerate the handle set, take the exact integer
coefficient of each handle on each residual check (intad.jacZ), and search.

Usage: habsorb.py [state.json]
"""
import os, sys, itertools, time
HERE = '/home/user/integer_solver/solve_lab/s10'
LAB = '/home/user/integer_solver/solve_lab'
sys.path.insert(0, os.path.join(LAB, 's9', 'eff'))
sys.path.insert(0, HERE)
import lib as L, tools as T, ad
from intad import jacZ
import suppfree
P = ad.P
src = sys.argv[1] if len(sys.argv) > 1 else 'PF_best_39015.json'
base = L.load(src if os.path.isabs(src) else os.path.join(HERE, src))
ad.fwd(base, rounds=6)


def report(v, tag):
    av = L.all_atom_values(v)
    f = L.failing_eqs(av)
    s = L.NEQ - len(f)
    nz = [a for a in range(L.NA) if a not in L.atom_out and av[a]]
    print('%-46s score %-6d failing %-3d checks %s' % (tag, s, len(f), nz),
          flush=True)
    return s, av, nz


BASE, av0, NZ = report(base, 'start (%s)' % os.path.basename(src))
_, freelist, SVS = suppfree.build(base, modp=None)
# every free input with an exact integer coefficient on some residual check
CAND = {}
for c in NZ:
    m = suppfree.atom_supp(c, base, SVS, modp=None)
    for i in range(len(freelist)):
        if not ((m >> i) & 1):
            continue
        u = freelist[i]
        CAND.setdefault(u, {})
t0 = time.time()
COEF = {}
for u in sorted(CAND):
    col = jacZ(u, base, NZ)
    if col:
        COEF[u] = col
print('%d free inputs with an exact integer effect on the residual (%.0fs)'
      % (len(COEF), time.time() - t0), flush=True)
ABS = []
for u, col in COEF.items():
    hits = [c for c in NZ if col.get(c) and av0[c] % col[c] == 0]
    if hits:
        ABS.append((u, col, hits))
print('%d of them can EXACTLY absorb at least one check:' % len(ABS), flush=True)
for u, col, hits in ABS[:20]:
    print('   x%-6d absorbs %s  (coefficients %s)'
          % (u, hits, {c: ('%s' % col[c])[:14] for c in hits}), flush=True)

best, bestv = BASE, list(base)
print('\n--- single absorptions ---', flush=True)
single = []
for u, col, hits in ABS:
    for c in hits:
        v = list(base)
        v[u] = v[u] - av0[c] // col[c]
        ad.fwd(v, rounds=6)
        aw = L.all_atom_values(v)
        s = L.NEQ - len(L.failing_eqs(aw))
        if aw[c] == 0:
            single.append((s, u, c))
            if s >= BASE - 2:
                print('   x%-6d absorbs a%-6d -> score %d' % (u, c, s), flush=True)
            if s > best:
                best, bestv = s, list(v)
                T.save(v, os.path.join(HERE, 'HA_%d.json' % s))
                print('      *** NEW BEST %d' % s, flush=True)
single.sort(reverse=True)
print('\n--- pairs and triples of the best absorptions ---', flush=True)
top = [(u, c) for _, u, c in single[:12]]
for k in (2, 3):
    for S in itertools.combinations(top, k):
        if len({u for u, _ in S}) < k or len({c for _, c in S}) < k:
            continue
        v = list(base)
        ok = True
        for u, c in S:
            av = L.all_atom_values(v)
            col = jacZ(u, v, [c]).get(c, 0)
            if not col or av[c] % col:
                ok = False
                break
            v[u] = v[u] - av[c] // col
            ad.fwd(v, rounds=6)
        if not ok:
            continue
        aw = L.all_atom_values(v)
        s = L.NEQ - len(L.failing_eqs(aw))
        if s >= BASE:
            print('   %s -> score %d' % (['x%d/a%d' % t for t in S], s), flush=True)
        if s > best:
            best, bestv = s, list(v)
            T.save(v, os.path.join(HERE, 'HA_%d.json' % s))
            print('      *** NEW BEST %d' % s, flush=True)
report(bestv, 'BEST')
print('\nbest %d (was %d)' % (best, BASE))
