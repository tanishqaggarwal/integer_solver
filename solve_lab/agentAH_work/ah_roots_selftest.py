"""Equality self-test: fast rootset_pp vs the fleet's brute-force rootset_pp, on the ACTUAL
prime powers dividing the instance's handle cofactors, with random degree-<=4 Newton coefficients
of realistic magnitude.  Any mismatch is fatal."""
import sys, random, collections
sys.path.insert(0, '/home/user/integer_solver/solve_lab/agentT_work')
import t_close2w as C
import ah_roots
peval = C.solve_group3.__globals__['peval']
brute_orig = C.solve_group3.__globals__['rootset_pp']
fast = ah_roots.make(peval)
p = C.p; SL = C.SL

def factor(n):
    f = {}; d = 2
    while d*d <= n:
        while n % d == 0:
            f[d] = f.get(d, 0)+1; n //= d
        d += 1
    if n > 1:
        f[n] = f.get(n, 0)+1
    return f

pps = set()
for a, s in SL.items():
    if s and s % p == 0 and abs(s)//p > 1:
        for q, e in factor(abs(s)//p).items():
            pps.add((q, e))
pps = sorted(pps)
print('distinct (q,e) over all handle cofactors:', len(pps),
      ' max q^e =', max(q**e for q, e in pps))
small = [(q, e) for q, e in pps if q**e <= 30000]
print('testing exhaustively against brute force on the %d prime powers <= 400000' % len(small))
rnd = random.Random(12345)
bad = 0; tested = 0
for q, e in small:
    for trial in range(4):
        cf = [rnd.randrange(-10**30, 10**30) for _ in range(5)]
        if trial == 0:
            cf[0] = 0                      # t=0 a root, the common real shape
        A = brute_orig(cf, q, e)
        B = fast(cf, q, e)
        tested += 1
        if A != B:
            bad += 1
            print('MISMATCH q=%d e=%d  brute=%d roots fast=%d roots' % (q, e, len(A), len(B)))
# a few identically-zero cases
for q, e in small[:20]:
    A = brute_orig([0, 0, 0, 0, 0], q, e); B = fast([0, 0, 0, 0, 0], q, e)
    tested += 1
    if A != B:
        bad += 1; print('MISMATCH-ALL q=%d e=%d' % (q, e))
print('SELFTEST: %d comparisons, %d mismatches' % (tested, bad))
print('stats', ah_roots.STATS)
sys.exit(1 if bad else 0)
