"""S11 step 14: solve for the RIGHT-HAND SIDE that makes 6 of 12 satisfiable.

Achievable atom vectors A (frame 2, detached x_7068/x_28730/x_29854/x_31864/x_642):
    A2 = a35758 = x_29854 - p*x_1329        A3 = a35759 = -x_29854 + 5113045*x_9118
    A4 = a35760 = x_31864 - p*x_10903       A5 = a35761 =  x_31864 + x_8731
  =>  A2 + A3 == 5113045*x_9118  (mod p)  =: R1
      A5 - A4 == x_8731          (mod p)  =: R2
and A0 (a22229), A1 (a22230), A6 (a35762) are unconstrained.

lattice3.py tested every 6-subset against the FIXED R1,R2 and found none solvable.
But R1 and R2 are set by two free inputs.  Solve for the ratio R1:R2 that MAKES a
6-subset solvable, then hit it with one Newton move.
"""
import os, sys, itertools, json
from fractions import Fraction
HERE = os.path.dirname(os.path.abspath(__file__))
LAB = os.path.dirname(HERE)
sys.path.insert(0, os.path.join(LAB, 's9', 'eff'))
sys.path.insert(0, HERE)
import lib as L, ad
P = ad.P
SEVEN = [22229, 22230, 35758, 35759, 35760, 35761, 35762]
IDX = {a: i for i, a in enumerate(SEVEN)}
E = sorted(set().union(*[set(L.atom2eq[a]) for a in SEVEN]))
print(f'{len(E)} equations touched by the seven')
rows = []
for e in E:
    m, sq, co = L.eq_atoms[e]
    rows.append([co.get(a, 0) for a in SEVEN])
for e, r in zip(E, rows): print(f'  eq {e:>6}: {r}')

v = L.load(os.path.join(LAB, 'best', 'new_instance_partial_39026.json'))
av = L.all_atom_values(v)
A_now = [av[a] for a in SEVEN]
R1 = (5113045 * v[7075] * v[9118]) % P
R2 = v[8731] % P
print(f'\ncheck the congruences at the delivered witness:')
print(f'  A2+A3 - R1 mod p = {(A_now[2] + A_now[3] - R1) % P}')
print(f'  A5-A4 - R2 mod p = {(A_now[5] - A_now[4] - R2) % P}')

def kernel_q(M):
    """rational kernel basis of the integer matrix M (list of rows)."""
    n = len(M); m = len(M[0])
    A = [[Fraction(x) for x in r] for r in M]
    piv, r_ = [], 0
    for j in range(m):
        k = next((i for i in range(r_, n) if A[i][j] != 0), None)
        if k is None: continue
        A[r_], A[k] = A[k], A[r_]
        pv = A[r_][j]; A[r_] = [x / pv for x in A[r_]]
        for i in range(n):
            if i != r_ and A[i][j] != 0:
                f = A[i][j]
                A[i] = [x - f * y for x, y in zip(A[i], A[r_])]
        piv.append(j); r_ += 1
    ker = []
    ps = set(piv)
    for fc in range(m):
        if fc in ps: continue
        w = [Fraction(0)] * m; w[fc] = Fraction(1)
        for i, pj in enumerate(piv): w[pj] = -A[i][fc]
        den = 1
        for x in w: den = den * x.denominator // __import__('math').gcd(den, x.denominator)
        ker.append([int(x * den) for x in w])
    return ker

print('\n=== 6-subsets of the 12 equations with a nonzero kernel ===')
found = []
for S in itertools.combinations(range(len(E)), 6):
    M = [rows[i] for i in S]
    K = kernel_q(M)
    if not K: continue
    for w in K:
        al = w[2] + w[3]          # coefficient of R1
        be = w[5] - w[4]          # coefficient of R2
        found.append((S, w, al, be))
print(f'{len(found)} (subset, kernel-vector) pairs')
ok = []
for S, w, al, be in found:
    if al == 0 and be == 0:
        # A = lambda*w needs R1 == R2 == 0
        ok.append(('needs R1=R2=0', S, w)); continue
    if al == 0:
        ok.append((f'needs R1 == 0, lambda = R2/{be}', S, w)); continue
    if be == 0:
        ok.append((f'needs R2 == 0, lambda = R1/{al}', S, w)); continue
    ok.append((f'needs R1*{be} == R2*{al} (mod p)', S, w))
seen = set()
for cond, S, w in ok:
    key = cond
    if key in seen: continue
    seen.add(key)
    print(f'  {cond:<44} eqs {[E[i] for i in S]}  w={w}')
    if len(seen) > 14: break
json.dump({'E': E, 'rows': rows, 'R1': str(R1), 'R2': str(R2),
           'found': [[list(S), w, al, be] for S, w, al, be in found]},
          open(os.path.join(HERE, 'rhs.json'), 'w'))
print(f'\nsaved rhs.json ({len(found)} candidates)')
