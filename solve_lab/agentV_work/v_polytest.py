#!/usr/bin/env python3
"""agent V -- self-test of v_poly against brute force.  A root finder that silently misses roots
turns a NEGATIVE result into a wrong one, so it is tested before it is trusted, at exactly the
prime sizes the |S|=17 residue actually uses (39703, 46273, 4787, 3449, ...)."""
import random, itertools, sys
import v_poly as P

rnd = random.Random(20260807)
fails = 0

# ---- 1. univariate root finding mod a prime, against brute force
print('--- roots_pp vs brute force, e=1 ---')
for q in [11, 19, 43, 89, 199, 463, 3449, 4787, 39703, 46273]:
    for trial in range(6):
        deg = rnd.choice([1, 2, 3, 4])
        mo = [rnd.randrange(q) for _ in range(deg+1)]
        if mo[-1] == 0:
            mo[-1] = 1
        # to Newton form by fitting
        u = P.nd_fit(lambda idx: P.mono_eval_exact(mo, idx[0]), 1, deg)
        un = [u[(i,)] for i in range(deg+1)]
        got = sorted(P.roots_pp(un, q, 1))
        exp = sorted(t for t in range(q) if P.mono_eval(mo, t, q) == 0)
        if got != exp:
            print('  FAIL q=%d deg=%d got=%s exp=%s' % (q, deg, got[:8], exp[:8])); fails += 1
    print('  q=%-8d ok' % q)

# ---- 2. prime powers
print('--- roots_pp vs brute force, e=2,3 ---')
for q, e in [(3, 3), (7, 2), (11, 2), (19, 2), (43, 2), (2, 5), (5, 4)]:
    m = q**e
    for trial in range(6):
        deg = rnd.choice([1, 2, 3])
        mo = [rnd.randrange(m) for _ in range(deg+1)]
        u = P.nd_fit(lambda idx: P.mono_eval_exact(mo, idx[0]), 1, deg)
        un = [u[(i,)] for i in range(deg+1)]
        got = sorted(P.roots_pp(un, q, e))
        exp = sorted(t for t in range(m) if P.mono_eval(mo, t, m) == 0)
        if got != exp:
            print('  FAIL q^e=%d deg=%d got=%s exp=%s' % (m, deg, got[:8], exp[:8])); fails += 1
    print('  q^e=%-8d ok' % m)

# ---- 3. multivariate Newton fit exactness, k = 1..4, degree <= 3 per variable
print('--- nd_fit exactness, k=1..4 ---')
for k in [1, 2, 3, 4]:
    D = 3
    coef = {idx: rnd.randrange(-10**6, 10**6) for idx in itertools.product(range(D+1), repeat=k)}

    def f(ts):
        tot = 0
        for idx, c in coef.items():
            v = c
            for j, i in enumerate(idx):
                v *= pow(ts[j], i)
            tot += v
        return tot
    cf = P.nd_fit(lambda idx: f(idx), k, D)
    ok = True
    for _ in range(40):
        ts = tuple(rnd.randrange(-10**5, 10**5) for _ in range(k))
        if P.nd_eval(cf, ts) != f(ts):
            ok = False
    print('  k=%d  exact at 40 random points: %s   degrees %s' % (k, ok, P.nd_degrees(cf)))
    if not ok:
        fails += 1

print('\nSELF-TEST %s (%d failures)' % ('PASS' if fails == 0 else 'FAIL', fails))
sys.exit(1 if fails else 0)
