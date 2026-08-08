#!/usr/bin/env python3
"""
reduce.py -- collapse the 39,031-equation instance to its minimal decision core,
using ONLY data extracted from EQUATIONS.txt (via pinrec.json / huge_consts.json)
and arithmetic verified here from scratch.

Output: anneal/core.json
"""
import json, collections, os

HERE = os.path.dirname(os.path.abspath(__file__))
LAB  = os.path.dirname(HERE)
p = 2**256 - 2**32 - 977          # the modulus that every pin residue lives in

# ---------- 1. the 256 gated constant pairs ----------
P = json.load(open(os.path.join(LAB, 'pinrec.json')))
byb = collections.defaultdict(dict)
for _f0, b, t, C, _coef, _f5 in P:
    byb[b][t] = C % p
assert len(byb) == 256 and all(len(v) == 2 for v in byb.values())

# ---------- 2. fit a single cubic  y^2 = x^3 + a2 x^2 + a4 x + a6 ----------
def fit(pts):
    A = [[x*x % p, x % p, 1, (y*y - x*x*x) % p] for x, y in pts]
    for i in range(3):
        pi = next((r for r in range(i, 3) if A[r][i] % p), None)
        if pi is None: return None
        A[i], A[pi] = A[pi], A[i]
        inv = pow(A[i][i], -1, p); A[i] = [v*inv % p for v in A[i]]
        for r in range(3):
            if r != i and A[r][i] % p:
                f = A[r][i]; A[r] = [(A[r][c]-f*A[i][c]) % p for c in range(4)]
    return A[0][3], A[1][3], A[2][3]

bits = sorted(byb)
raw  = [tuple(v for _, v in sorted(byb[b].items())) for b in bits]

curve = None
for seed in range(0, 40):
    for mask in range(8):
        trial = [(raw[seed+j][mask >> j & 1], raw[seed+j][1 - (mask >> j & 1)]) for j in range(3)]
        c = fit(trial)
        if not c: continue
        a2, a4, a6 = c
        good = sum(1 for u, v in raw
                   if (v*v-u*u*u-a2*u*u-a4*u-a6) % p == 0 or (u*u-v*v*v-a2*v*v-a4*v-a6) % p == 0)
        if good == 256: curve = c; break
    if curve: break
assert curve, "no single cubic fits all 256 pairs"
a2, a4, a6 = curve

def on(x, y): return (y*y - x*x*x - a2*x*x - a4*x - a6) % p == 0
pts = {}
for b, (u, v) in zip(bits, raw):
    pts[b] = (u, v) if on(u, v) else (v, u)
assert all(on(*q) for q in pts.values())

# ---------- 3. the target pair ----------
H = json.load(open(os.path.join(LAB, 'huge_consts.json')))
T = (int(H['C1']) % p, int(H['C2']) % p)
assert on(*T), "target is not on the curve"

# ---------- 4. depress the cubic: x -> x - a2/3 gives y^2 = x^3 + A x + B ----------
i3 = pow(3, -1, p)
A = (a4 - a2*a2*i3) % p
B = (a6 - a2*a4*i3 + 2*pow(a2, 3, p)*pow(i3, 3, p)) % p
sh = a2 * i3 % p
dep = lambda q: ((q[0] + sh) % p, q[1] % p)

json.dump({
    'p': str(p), 'a2': str(a2), 'a4': str(a4), 'a6': str(a6),
    'A': str(A), 'B': str(B),
    'bit_vars': bits,
    'points': {str(b): [str(x) for x in pts[b]] for b in bits},
    'points_dep': {str(b): [str(x) for x in dep(pts[b])] for b in bits},
    'target': [str(x) for x in T],
    'target_dep': [str(x) for x in dep(T)],
}, open(os.path.join(HERE, 'core.json'), 'w'), indent=1)

print(f"p  = {p}")
print(f"curve  y^2 = x^3 + a2 x^2 + a4 x + a6")
print(f"  a2 = {a2}\n  a4 = {a4}\n  a6 = {a6}")
print(f"depressed:  y^2 = x^3 + {A} x + {B}")
print(f"  A == 0 ?  {A == 0}")
print(f"256 gated pairs: all on the curve.  target on the curve: True")
print("wrote anneal/core.json")
