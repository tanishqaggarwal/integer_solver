"""S11 step 26: the optimum with a22231 as an eighth free value.

Detach x_4432 from a22231.  Then
    a22231 = (x_4432 - x_19964) - x_28730 = K - x_28730 ,  K fixed
    a22230 = x_28730 - p*x_9413
 => A1 + B == K (mod p),  otherwise both free.
And detaching x_4432 severs x_28730 -> a7930, so A1 is no longer pinned.
Cost: a22231 has ZERO equations outside the twelve; the only casualty is a37887,
a single-equation check.
Maximise the number of the 12 equations satisfied over this admissible set.
"""
import os, sys, itertools
from fractions import Fraction
HERE = os.path.dirname(os.path.abspath(__file__))
LAB = os.path.dirname(HERE)
sys.path.insert(0, os.path.join(LAB, 's9', 'eff'))
sys.path.insert(0, HERE)
import lib as L, ad
P = ad.P
SEVEN = [22229, 22230, 35758, 35759, 35760, 35761, 35762]
ATOMS = SEVEN + [22231]
E = sorted(set().union(*[set(L.atom2eq[a]) for a in SEVEN]))
rows = []
for e in E:
    m, sq, co = L.eq_atoms[e]
    rows.append([co.get(a, 0) for a in ATOMS])
print('12 x 8 coefficient rows (A0..A6, B=a22231):')
for e, r in zip(E, rows): print(f'  eq {e:>6}: {r}')

def solvable(sel):
    """is there an admissible (A,B) killing exactly the equations in sel?
    Admissible: A2,A3,A4,A5 free; A0 + 7376877*A6 == C0 (mod p);
                A1 + B == K (mod p).  Both congruences have FREE right-hand
                sides here (C0 shifts by multiples of p, K is fixed but A1 and B
                trade off), so over Q the only constraint is the linear system."""
    M = [[Fraction(rows[i][k]) for k in range(8)] for i in sel]
    n = len(M); r_ = 0; piv = []
    for j in range(8):
        k = next((i for i in range(r_, n) if M[i][j] != 0), None)
        if k is None: continue
        M[r_], M[k] = M[k], M[r_]
        pv = M[r_][j]; M[r_] = [x / pv for x in M[r_]]
        for i in range(n):
            if i != r_ and M[i][j] != 0:
                f = M[i][j]
                M[i] = [x - f * y for x, y in zip(M[i], M[r_])]
        piv.append(j); r_ += 1
    return 8 - r_          # kernel dimension

best = None
for k in range(12, 4, -1):
    hit = []
    for sel in itertools.combinations(range(12), k):
        d = solvable(sel)
        if d > 0: hit.append((sel, d))
    print(f'subsets of size {k}: {len(hit)} have a nonzero kernel')
    if hit:
        best = (k, hit[:4])
        break
if best:
    k, hit = best
    print(f'\nMAX SATISFIABLE (over Q, ignoring the mod-p congruences): {k} of 12')
    print(f'  -> {12-k} of the twelve fail, plus a37887 = {12-k+1} equations total')
    print(f'  -> score {L.NEQ - (12-k+1)}')
    for sel, d in hit[:3]:
        print(f'   eqs {[E[i] for i in sel]}  kernel dim {d}')
