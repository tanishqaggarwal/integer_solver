"""S11 step 97: THE QUADRANT IS THE BARRIER, and x2081 / x24601 are its two switches.

The structural free-input supports settle what the selectors actually are:

    s = x7715  is 1 only because of x24601      (its 178-bit support has one nonzero)
    t = x34554 is 1 only because of x2081       (its  78-bit support has one nonzero)

and the three quadrant indicators are exactly

    x15298 = s*t        x34606 = s*(1-t)        x5647 = (1-s)*t

So the instance is a point-addition MULTIPLEXER, and reading a1618 and a688 without
assuming a branch gives

    x24468 ≡ x5647*x2 + x34606*x1 + x15298*x3   ≡ C2
    x18956 ≡ x5647*y2 + x34606*y1 + x15298*y3   ≡ C1/8863713

  (s,t) = (1,1)   x15298 = 1   ->  P3 = P1 + P2 must close, and x3, y3 are pinned
  (s,t) = (1,0)   x34606 = 1   ->  the addition check VANISHES; instead x1 ≡ C2, y1 ≡ C1'
  (s,t) = (0,1)   x5647  = 1   ->  the addition check VANISHES; instead x2 ≡ C2, y2 ≡ C1'
  (s,t) = (0,0)                ->  x24468 ≡ 0 ≢ C2: dead

i.e. "P3 = P1 + P2, or P2 is infinity and P3 = P1, or P1 is infinity and P3 = P2".
Every session of this lab has worked inside (1,1) -- x2081 and x24548 have been on the
FORBID list since Session 9 -- which is why the addition looked like a wall.  In the
other two quadrants the addition is not checked at all.

This module goes to a quadrant, re-solves the advice DAG there, and then solves the
checks that the quadrant switches ON, exactly: each is linear in a released or advice
value, so one unit probe per knob gives the exact column.

Usage: quadrant.py [q01|q10] [state.json]
"""
import os, sys, time
HERE = '/home/user/integer_solver/solve_lab/s10'
LAB = '/home/user/integer_solver/solve_lab'
sys.path.insert(0, os.path.join(LAB, 's9', 'eff'))
sys.path.insert(0, HERE)
import lib as L, tools as T, ad
from intad import jacZ
import suppfree
P = ad.P
q = sys.argv[1] if len(sys.argv) > 1 else 'q01'
src = sys.argv[2] if len(sys.argv) > 2 else ('AG_38993.json' if q == 'q01'
                                             else 'AG_38977.json')
v = L.load(src if os.path.isabs(src) else os.path.join(HERE, src))
ad.fwd(v, rounds=6)


def report(v, tag):
    av = L.all_atom_values(v)
    s = L.NEQ - len(L.failing_eqs(av))
    nz = [a for a in range(L.NA) if a not in L.atom_out and av[a]]
    print('%-40s score %-6d checks %s' % (tag, s, nz), flush=True)
    return s, av, nz


s0, av, NZ = report(v, '%s start (%s)' % (q, src))
print('   quadrant: s=%d t=%d  ->  x15298=%d x34606=%d x5647=%d'
      % (v[7715], v[34554], v[15298], v[34606], v[5647]), flush=True)

_, freelist, SVS = suppfree.build(v, modp=None)
FREE = set(freelist)
# knobs: every free input that structurally reaches a failing check, ranked by how
# few checks it touches (the released constants and advice values come out on top)
KN = set()
for c in NZ:
    m = suppfree.atom_supp(c, v, SVS, modp=None)
    KN |= {freelist[i] for i in range(len(freelist)) if (m >> i) & 1}
PREF = [22152, 33462, 6418, 12553, 22162, 30213, 14853, 16742, 8778, 22649,
        31339, 14623, 24548, 12186, 24908]
KN = [u for u in PREF if u in KN] + sorted(KN - set(PREF))
print('   %d knobs reach the residual; preferred head %s'
      % (len(KN), KN[:10]), flush=True)


def val(v, c):
    return L.all_atom_values(v)[c] % P


def solve(v, targets, knobs, rounds=3):
    """Exact linear solve: one unit probe per knob gives the exact column."""
    for _ in range(rounds):
        av = L.all_atom_values(v)
        b = [(-av[c]) % P for c in targets]
        if not any(b):
            return v, True
        cols = []
        for u in knobs:
            w = list(v)
            w[u] = w[u] + 1
            ad.fwd(w, rounds=6)
            aw = L.all_atom_values(w)
            cols.append([(aw[c] - av[c]) % P for c in targets])
        n, m = len(targets), len(knobs)
        A = [[cols[j][i] for j in range(m)] + [b[i]] for i in range(n)]
        piv, r_ = [], 0
        for j in range(m):
            k = next((i for i in range(r_, n) if A[i][j]), None)
            if k is None:
                continue
            A[r_], A[k] = A[k], A[r_]
            inv = pow(A[r_][j], -1, P)
            A[r_] = [x * inv % P for x in A[r_]]
            for i in range(n):
                if i != r_ and A[i][j]:
                    f = A[i][j]
                    A[i] = [(x - f * z) % P for x, z in zip(A[i], A[r_])]
            piv.append(j)
            r_ += 1
        if any(A[i][m] for i in range(r_, n)):
            return v, False
        d = [0] * m
        for i, j in enumerate(piv):
            d[j] = A[i][m]
        w = list(v)
        for j, u in enumerate(knobs):
            if d[j]:
                w[u] = w[u] + d[j]
        ad.fwd(w, rounds=6)
        aw = L.all_atom_values(w)
        if all(aw[c] % P == 0 for c in targets):
            return w, True
        v = w
    return v, False


print('\n--- solving the checks this quadrant switches on ---', flush=True)
t0 = time.time()
for k in range(1, min(len(NZ), 8) + 1):
    for size in (len(KN[:10]),):
        tg = NZ[:k]
        w, ok = solve(v, tg, KN[:10])
        s, aw, nz = report(w, '  targets %s' % tg)
        if ok and s >= s0:
            v, s0, NZ = w, s, nz
            print('     *** accepted (%.0fs)' % (time.time() - t0), flush=True)
            T.save(v, os.path.join(HERE, 'Q_%s_%d.json' % (q, s)))
            break

print('\n--- integer lift ---', flush=True)
_, fl, S = suppfree.build(v, modp=None)
for _ in range(20):
    av = L.all_atom_values(v)
    todo = [a for a in range(L.NA) if a not in L.atom_out and av[a] and av[a] % P == 0]
    cur = L.NEQ - len(L.failing_eqs(av))
    moved = False
    for c in todo:
        m = suppfree.atom_supp(c, v, S, modp=None)
        for i in range(len(fl)):
            if not ((m >> i) & 1):
                continue
            u = fl[i]
            g = jacZ(u, v, [c]).get(c, 0)
            if not g or g % P or av[c] % g:
                continue
            w = list(v)
            w[u] = w[u] - av[c] // g
            ad.fwd(w, rounds=6)
            a2 = L.all_atom_values(w)
            if a2[c] == 0 and L.NEQ - len(L.failing_eqs(a2)) >= cur:
                v, av, moved = w, a2, True
                break
        if moved:
            break
    if not moved:
        break
s, av, nz = report(v, 'FINAL %s' % q)
T.save(v, os.path.join(HERE, 'Q_%s_%d.json' % (q, s)))
print('saved Q_%s_%d.json' % (q, s))
