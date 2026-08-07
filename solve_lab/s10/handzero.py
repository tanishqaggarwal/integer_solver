"""S11 step 115: with A = B = 0 the three primitives are PURE HANDLE, so set the handle.

Breaking the two pins is far cheaper than assumed: moving x3 off a1618 costs +2
equations and moving y3 off a688 costs +4, and both together cost +4 (they overlap).
So closing the addition should be nearly free -- yet PF_best_39015 still shows six
nonzero checks, three of which should have vanished:

    a19297 = x11150*x15298 + x4007 ,   x4007 = x5101*x30317 = p*x30317
    a19299 = x15298*x25739 - 6672769*x29804 ,  x29804 = x5146*x32017 = p*x5146
    a30984 = 537773*x15298*x37758 - x35605 ,   x35605 = x2936*x26789 = p*x2936

With A = B = 0 the three numerators x11150, x25739 and x37758 are all zero, so each of
these is exactly `p times its handle` -- and the handles x30317, x5146 and x2936 are
free.  Setting them to zero kills the checks outright.  The atom-level lift missed
them because it only fires on checks that are ≡ 0 (mod p) AND whose repair does not
lower the score, and it oscillated instead.

Usage: handzero.py [state.json]
"""
import os, sys
HERE = '/home/user/integer_solver/solve_lab/s10'
LAB = '/home/user/integer_solver/solve_lab'
sys.path.insert(0, os.path.join(LAB, 's9', 'eff'))
sys.path.insert(0, HERE)
import lib as L, tools as T, ad
from intad import jacZ
import suppfree
P = ad.P
src = sys.argv[1] if len(sys.argv) > 1 else 'PF_best_39015.json'
v = L.load(src if os.path.isabs(src) else os.path.join(HERE, src))
ad.fwd(v, rounds=6)


def report(v, tag):
    av = L.all_atom_values(v)
    f = L.failing_eqs(av)
    s = L.NEQ - len(f)
    nz = [a for a in range(L.NA) if a not in L.atom_out and av[a]]
    print('%-40s score %-6d failing %-3d checks %s' % (tag, s, len(f), nz),
          flush=True)
    return s, av


report(v, 'start (%s)' % os.path.basename(src))
print('   A = x35389 = %s ; B = x6671 = %s'
      % (v[35389] % P == 0, v[6671] % P == 0), flush=True)
print('   handles: x30317 = %s, x5146 = %s, x2936 = %s'
      % (str(v[30317])[:18], str(v[5146])[:18], str(v[2936])[:18]), flush=True)

best, bestv = report(v, 'baseline')[0], list(v)
# zero the three handles, together and in every combination
import itertools
H = [30317, 5146, 2936]
for k in range(1, 4):
    for S in itertools.combinations(H, k):
        w = list(v)
        for h in S:
            w[h] = 0
        ad.fwd(w, rounds=6)
        s, av = report(w, '   zero %s' % ','.join('x%d' % h for h in S))
        if s > best:
            best, bestv = s, list(w)
            T.save(w, os.path.join(HERE, 'HZ_%d.json' % s))
            print('      *** NEW BEST %d -- saved HZ_%d.json' % (s, s), flush=True)

# then a general handle sweep: every nonzero check that is p times a free handle
v = list(bestv)
_, freelist, SVS = suppfree.build(v, modp=None)
for rnd in range(25):
    av = L.all_atom_values(v)
    NZ = [a for a in range(L.NA) if a not in L.atom_out and av[a]]
    cur = L.NEQ - len(L.failing_eqs(av))
    moved = False
    for c in NZ:
        m = suppfree.atom_supp(c, v, SVS, modp=None)
        for i in range(len(freelist)):
            if not ((m >> i) & 1):
                continue
            u = freelist[i]
            g = jacZ(u, v, [c]).get(c, 0)
            if not g or av[c] % g:
                continue
            w = list(v)
            w[u] = w[u] - av[c] // g
            ad.fwd(w, rounds=6)
            aw = L.all_atom_values(w)
            s2 = L.NEQ - len(L.failing_eqs(aw))
            if aw[c] == 0 and s2 > cur:
                v, cur, moved = w, s2, True
                print('   check a%d zeroed by x%d -> score %d' % (c, u, s2),
                      flush=True)
                break
        if moved:
            break
    if not moved:
        break
s, av = report(v, 'FINAL')
if s > best:
    best, bestv = s, list(v)
    T.save(v, os.path.join(HERE, 'HZ_%d.json' % s))
    print('*** saved HZ_%d.json' % s)
print('\nbest %d' % best)
