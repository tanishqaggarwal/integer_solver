"""S11 step 36: (a) do any of the 792 kernels already satisfy the compatibility
condition?  (b) is there a cheap repair variable for a37887?"""
import os, sys, itertools, collections
from fractions import Fraction
HERE = os.path.dirname(os.path.abspath(__file__))
LAB = os.path.dirname(HERE)
sys.path.insert(0, os.path.join(LAB, 's9', 'eff'))
sys.path.insert(0, HERE)
import lib as L, tools as T, ad
from frame3 import DETACH, definer, ORDER, FREE, CHECKS, fwd, score, SSET
P = ad.P
SEVEN = [22229, 22230, 35758, 35759, 35760, 35761, 35762]
ATOMS = SEVEN + [22231]
E = sorted(set().union(*[set(L.atom2eq[a]) for a in SEVEN]))
rows = []
for e in E:
    m, sq, co = L.eq_atoms[e]
    rows.append([co.get(a, 0) for a in ATOMS])
base = L.load(os.path.join(LAB, 'best', 'new_instance_partial_39026.json'))
fwd(base)
C0 = (base[7068] - base[2099]) % P
K = (base[4432] - base[19964]) % P
print(f'C0 = {str(C0)[:34]}\nK  = {str(K)[:34]}')

def kernel(sel):
    M = [[Fraction(rows[i][k]) for k in range(8)] for i in sel]
    n = len(M); piv, r_ = [], 0
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
    out, ps = [], set(piv)
    for fc in range(8):
        if fc in ps: continue
        w = [Fraction(0)] * 8; w[fc] = Fraction(1)
        for i, pj in enumerate(piv): w[pj] = -M[i][fc]
        den = 1
        import math
        for x in w: den = den * x.denominator // math.gcd(den, x.denominator)
        out.append([int(x * den) for x in w])
    return out

hits = 0; tested = 0
for sel in itertools.combinations(range(12), 7):
    for w in kernel(sel):
        al = (w[0] + 7376877 * w[6]) % P
        be = (w[1] + w[7]) % P
        if al == 0 or be == 0: continue
        tested += 1
        if (C0 * be - K * al) % P == 0:
            hits += 1
            print(f'  *** COMPATIBLE: eqs {[E[i] for i in sel]}  w={w}')
print(f'tested {tested} kernels; compatible: {hits}')

print('\n=== a37887: coefficient of each variable (as a multiple of x_4432) ===')
lin = collections.defaultdict(int)
other = 0
for m, c in L.polys[37887].items():
    r = [z for z in m if z != 4432]
    if len(m) - len(r) == 1 and len(r) == 1:
        lin[r[0]] += c
    else:
        other += 1
print(f'{len(lin)} variables enter as x_4432*w; {other} monomials are other shapes')
v = list(base)
cheap = sorted(lin.items(), key=lambda kv: (abs(kv[1]), len(L.var_atoms[kv[0]])))
for w, c in cheap[:14]:
    print(f'  coef {c:>12}  x_{w:<7} {"FREE" if w in FREE else "gate":<5} '
          f'atoms {sorted(L.var_atoms[w])[:6]}')
