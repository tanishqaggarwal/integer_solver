"""S11 step 117: the last lock is ONE divisibility, and the k*p freedom can move it.

With A = x35389 and B = x6671 both ≡ 0 (mod p), the three primitives are exactly

    x11150 = 8646263*A  + 1073965*B     a19297 = x11150  + p*x30317
    x25739 = 10159099*A + 6926539*B     a19299 = x25739  - 6672769*p*x5146
    x37758 = 8272701*A  + 5921311*B     a30984 = 537773*x37758 - p*x2936

so a19297 and a30984 are absorbed by their handles the moment A and B are multiples
of p -- and habsorb confirms x30317 absorbs a19297 exactly.  a19299 is the exception:
its handle enters with coefficient 6672769*p, so it needs

    6672769  |  (10159099*a + 6926539*b),        A = p*a,  B = p*b

one extra divisibility.  And a and b are NOT fixed: every advice value is k*p + r, and
bumping the k of x22162 or x30213 shifts A and B by exact multiples of p while leaving
every residue -- hence every congruence -- untouched.  So a and b move on a lattice,
and the lock is a linear congruence mod 6672769 that can simply be solved.

Usage: divlock.py [state.json] [RANGE]
"""
import os, sys, math, time
HERE = '/home/user/integer_solver/solve_lab/s10'
LAB = '/home/user/integer_solver/solve_lab'
sys.path.insert(0, os.path.join(LAB, 's9', 'eff'))
sys.path.insert(0, HERE)
import lib as L, tools as T, ad
from intad import jacZ
import suppfree
P = ad.P
src = sys.argv[1] if len(sys.argv) > 1 else 'PF_best_39015.json'
RANGE = int(sys.argv[2]) if len(sys.argv) > 2 else 40
base = L.load(src if os.path.isabs(src) else os.path.join(HERE, src))
ad.fwd(base, rounds=6)


def report(v, tag):
    av = L.all_atom_values(v)
    f = L.failing_eqs(av)
    s = L.NEQ - len(f)
    nz = [a for a in range(L.NA) if a not in L.atom_out and av[a]]
    print('%-40s score %-6d failing %-3d checks %s' % (tag, s, len(f), nz),
          flush=True)
    return s, av, nz


BASE, av0, NZ = report(base, 'start (%s)' % os.path.basename(src))
A, B = base[35389], base[6671]
print('   A = %s  (A mod p = %d)' % (str(A)[:26], A % P), flush=True)
print('   B = %s  (B mod p = %d)' % (str(B)[:26], B % P), flush=True)
if A % P or B % P:
    print('A or B not ≡ 0 mod p: the primitives cannot be absorbed here')
    sys.exit()
a, b = A // P, B // P
Q = 6672769
print('   a = A/p = %s' % str(a)[:30], flush=True)
print('   b = B/p = %s' % str(b)[:30], flush=True)
print('   6672769 factors:', flush=True)
n, f = Q, []
d = 2
while d * d <= n:
    while n % d == 0:
        f.append(d)
        n //= d
    d += 1
if n > 1:
    f.append(n)
print('      %s = %s' % (Q, ' * '.join(map(str, f))), flush=True)
tgt = (10159099 * a + 6926539 * b) % Q
print('   (10159099*a + 6926539*b) mod 6672769 = %d  %s'
      % (tgt, 'ALREADY ZERO' if tgt == 0 else 'must be driven to 0'), flush=True)

# how do a and b move when we bump the k*p part of x22162 (w5) and x30213 (w6)?
print('\n--- effect of bumping the k*p part of the two free advice values ---',
      flush=True)
EFF = {}
for u in (22162, 30213):
    v = list(base)
    v[u] = v[u] + P
    ad.fwd(v, rounds=6)
    da = (v[35389] - A)
    db = (v[6671] - B)
    ok = (da % P == 0 and db % P == 0)
    EFF[u] = (da // P if da % P == 0 else None, db // P if db % P == 0 else None)
    print('   x%-6d += p  ->  dA/p = %s, dB/p = %s   (exact multiples: %s)'
          % (u, str(EFF[u][0])[:24], str(EFF[u][1])[:24], ok), flush=True)

c1 = (10159099 * (EFF[22162][0] or 0) + 6926539 * (EFF[22162][1] or 0)) % Q
c2 = (10159099 * (EFF[30213][0] or 0) + 6926539 * (EFF[30213][1] or 0)) % Q
print('\n   step of the target per unit k:  x22162 -> %d,  x30213 -> %d  (mod %d)'
      % (c1, c2, Q), flush=True)
g = math.gcd(math.gcd(c1, c2), Q)
print('   gcd(steps, %d) = %d ;  target %d ;  solvable: %s'
      % (Q, g, tgt, tgt % g == 0), flush=True)

sols = []
if tgt % g == 0 and (c1 or c2):
    for n1 in range(-RANGE, RANGE + 1):
        r = (tgt + c1 * n1) % Q
        if c2:
            gg = math.gcd(c2, Q)
            if r % gg == 0:
                n2 = (-r // gg) * pow(c2 // gg, -1, Q // gg) % (Q // gg)
                sols.append((n1, n2))
        elif r == 0:
            sols.append((n1, 0))
print('   %d candidate (k1, k2) shifts found' % len(sols), flush=True)

best, bestv = BASE, list(base)
t0 = time.time()
_, freelist, SVS = suppfree.build(base, modp=None)


def absorb_all(v):
    for _ in range(12):
        av = L.all_atom_values(v)
        cur = L.NEQ - len(L.failing_eqs(av))
        moved = False
        for c in [x for x in range(L.NA) if x not in L.atom_out and av[x]]:
            m = suppfree.atom_supp(c, v, SVS, modp=None)
            for i in range(len(freelist)):
                if not ((m >> i) & 1):
                    continue
                u = freelist[i]
                col = jacZ(u, v, [c]).get(c, 0)
                if not col or av[c] % col:
                    continue
                w = list(v)
                w[u] = w[u] - av[c] // col
                ad.fwd(w, rounds=6)
                aw = L.all_atom_values(w)
                if aw[c] == 0 and L.NEQ - len(L.failing_eqs(aw)) >= cur:
                    v, moved = w, True
                    break
            if moved:
                break
        if not moved:
            break
    return v


for n1, n2 in sols[:60]:
    if time.time() - t0 > 1500:
        break
    v = list(base)
    v[22162] += n1 * P
    v[30213] += n2 * P
    ad.fwd(v, rounds=6)
    A2, B2 = v[35389], v[6671]
    if A2 % P or B2 % P:
        continue
    chk = (10159099 * (A2 // P) + 6926539 * (B2 // P)) % Q
    v = absorb_all(v)
    aw = L.all_atom_values(v)
    s = L.NEQ - len(L.failing_eqs(aw))
    nz = [x for x in range(L.NA) if x not in L.atom_out and aw[x]]
    if chk == 0 or s > BASE - 2:
        print('   k1=%-5d k2=%-22s  lock %d  score %d  checks %s'
              % (n1, str(n2)[:20], chk, s, nz), flush=True)
    if s > best:
        best, bestv = s, list(v)
        T.save(v, os.path.join(HERE, 'DL_%d.json' % s))
        print('      *** NEW BEST %d -- saved DL_%d.json' % (s, s), flush=True)
print('\nbest %d (was %d)' % (best, BASE))
