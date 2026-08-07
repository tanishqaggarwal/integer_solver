"""S11 step 119: at 39,017 the residual is two numbers -- what can their RATIO save?

FIN_39017 (checker-verified) leaves exactly three nonzero checks:

    a688   in 11 equations      a1618  in 14 equations      a40608 in 1

and 16 distinct failing equations, so ten of them contain BOTH a688 and a1618.  In
each of those the combination is `c1*a688 + c2*a1618` with every other atom zero, so
the equation holds as soon as the two values sit in the right ratio -- exactly the
cancellation mechanism that makes the 39,026 deliverable a coding optimum (§152).

The residues of a688 and a1618 mod p are fixed (the primitive pins w5 and w6 mod
p), but their handles add arbitrary multiples of p.  So write

    a688 = r688 + p*s ,   a1618 = r1618 + p*t ,    s, t free integers

and an equation is satisfiable iff  c1*r688 + c2*r1618 = 0 (mod p)  -- a fixed test --
and then  c1*s + c2*t = q  for the resulting quotient.  Two free integers, so the
question is how many of those linear conditions can hold at once.

Usage: ratio.py [state.json]
"""
import os, sys, itertools
HERE = '/home/user/integer_solver/solve_lab/s10'
LAB = '/home/user/integer_solver/solve_lab'
sys.path.insert(0, os.path.join(LAB, 's9', 'eff'))
sys.path.insert(0, HERE)
import lib as L, tools as T, ad
from intad import jacZ
import suppfree
P = ad.P
src = sys.argv[1] if len(sys.argv) > 1 else 'FIN_39017.json'
v = L.load(src if os.path.isabs(src) else os.path.join(HERE, src))
ad.fwd(v, rounds=6)
av = L.all_atom_values(v)
BASE = L.NEQ - len(L.failing_eqs(av))
R = [a for a in range(L.NA) if a not in L.atom_out and av[a]]
F = sorted(L.failing_eqs(av))
print('%s: score %d; residual %s; %d failing' % (src, BASE, R, len(F)), flush=True)
r688, r1618 = av[688] % P, av[1618] % P
print('   a688 ≡ %d\n   a1618 ≡ %d' % (r688, r1618), flush=True)

rows = []
for e in F:
    co = L.eq_atoms[e][2]
    c1, c2, c3 = co.get(688, 0), co.get(1618, 0), co.get(40608, 0)
    other = [a for a in co if a not in (688, 1618, 40608) and av[a]]
    val = sum(c * av[a] for a, c in co.items() if av[a])
    ok = (c1 * r688 + c2 * r1618 + c3 * (av[40608] % P)) % P == 0
    rows.append((e, c1, c2, c3, val, ok, other))
print('\n  eq      c(a688) c(a1618) c(a40608)  passes the mod-p test?')
for e, c1, c2, c3, val, ok, other in rows:
    print('  %-7d %-8d %-8d %-8d  %s%s'
          % (e, c1, c2, c3, ok, '   other nonzero atoms: %s' % other if other else ''),
          flush=True)
good = [r for r in rows if r[5]]
print('\n%d of %d failing equations pass the mod-p ratio test' % (len(good), len(F)),
      flush=True)
if not good:
    print('none: no choice of the two handle multiples can save any equation')
    sys.exit()

# for those, the integer condition is c1*s + c2*t = q with q = -val/p
cond = []
for e, c1, c2, c3, val, ok, other in good:
    if val % P:
        continue
    cond.append((e, c1, c2, -(val // P)))
print('conditions c1*s + c2*t = q :', flush=True)
for e, c1, c2, q in cond:
    print('   eq %-7d %d*s %+d*t = %s' % (e, c1, c2, str(q)[:26]), flush=True)


def solve2(sub):
    """Integer solution of the listed 2-variable conditions, or None."""
    import math
    S = T2 = None
    # reduce pairwise
    A = [[c1, c2, q] for _, c1, c2, q in sub]
    # gaussian over Z on 2 columns
    r = 0
    for j in range(2):
        piv = None
        for i in range(r, len(A)):
            if A[i][j]:
                piv = i
                break
        if piv is None:
            continue
        A[r], A[piv] = A[piv], A[r]
        for i in range(len(A)):
            if i != r and A[i][j]:
                g = math.gcd(A[r][j], A[i][j])
                m1, m2 = A[i][j] // g, A[r][j] // g
                A[i] = [A[i][k] * m2 - A[r][k] * m1 for k in range(3)]
        r += 1
    for i in range(r, len(A)):
        if A[i][2]:
            return None
    st = [0, 0]
    for j in range(min(r, 2) - 1, -1, -1):
        row = A[j]
        rest = row[2] - sum(row[k] * st[k] for k in range(2) if k != j)
        if row[j] == 0:
            if rest:
                return None
            continue
        if rest % row[j]:
            return None
        st[j] = rest // row[j]
    for c1, c2, q in [(x[0], x[1], x[2]) for x in A]:
        pass
    for _, c1, c2, q in sub:
        if c1 * st[0] + c2 * st[1] != q:
            return None
    return st


best = (0, None)
for k in range(len(cond), 0, -1):
    got = False
    for sub in itertools.combinations(cond, k):
        st = solve2(list(sub))
        if st is not None:
            print('\n*** %d conditions solvable together: eqs %s  with (s,t) = %s'
                  % (k, [x[0] for x in sub], [str(x)[:20] for x in st]), flush=True)
            best = (k, (sub, st))
            got = True
            break
    if got:
        break
print('\nceiling from the ratio: %d + %d = %d'
      % (BASE, best[0], BASE + best[0]))
