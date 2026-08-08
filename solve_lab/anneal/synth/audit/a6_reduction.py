#!/usr/bin/env python3
"""ATTACK 3: reduction soundness, recomputed INDEPENDENTLY from raw data.

Trust nothing from reduce.py/core.json. Reload pinrec.json + huge_consts.json,
re-fit the cubic from scratch with an independent linear solve, re-verify every
gated pair lies on it, re-verify the target, re-derive the depressed curve, and
re-check the doubling chain P_i == 2^i G and the prime order.  Then compare to
what instance.py / core.json actually feed the encoder, and probe:
  * does the fitted (p,B,G,T) match the raw constants bit-for-bit?
  * is the doubling chain EXACTLY 2^i G, in the order the encoder uses?
  * solution multiplicity: how many k in [0,2^256) satisfy k*G==T  (k vs k+n)?
"""
import json, os, sys
LAB = '/home/user/integer_solver/solve_lab'
AN = os.path.join(LAB, 'anneal')
sys.path.insert(0, AN)
p = 2**256 - 2**32 - 977

P = json.load(open(os.path.join(LAB, 'pinrec.json')))
H = json.load(open(os.path.join(LAB, 'huge_consts.json')))
core = json.load(open(os.path.join(AN, 'core.json')))

import collections
byb = collections.defaultdict(dict)
for rec in P:
    _f0, b, t, C, _coef, _f5 = rec
    byb[b][t] = C % p
print(f"pinrec: {len(byb)} selector bits, each gating {set(len(v) for v in byb.values())} constants")
assert len(byb) == 256 and all(len(v) == 2 for v in byb.values())

bits = sorted(byb)
raw = [tuple(v for _, v in sorted(byb[b].items())) for b in bits]

# --- independent cubic fit: pick 3 pairs, solve for (a2,a4,a6), then verify all 256
def solve3(rows):
    # rows: 3 x [c0,c1,c2 | rhs] over F_p, Gaussian elim
    M = [r[:] for r in rows]
    for i in range(3):
        piv = next((r for r in range(i, 3) if M[r][i] % p), None)
        if piv is None: return None
        M[i], M[piv] = M[piv], M[i]
        inv = pow(M[i][i], -1, p)
        M[i] = [v * inv % p for v in M[i]]
        for r in range(3):
            if r != i and M[r][i] % p:
                f = M[r][i]; M[r] = [(M[r][c] - f * M[i][c]) % p for c in range(4)]
    return M[0][3], M[1][3], M[2][3]

# use points (x,y): y^2 - x^3 = a2 x^2 + a4 x + a6
curve = None
import itertools
for combo in itertools.combinations(range(min(20, len(raw))), 3):
    for orient in itertools.product((0, 1), repeat=3):
        pts = [(raw[combo[j]][orient[j]], raw[combo[j]][1 - orient[j]]) for j in range(3)]
        rows = [[x*x % p, x % p, 1, (y*y - x*x*x) % p] for x, y in pts]
        c = solve3(rows)
        if not c: continue
        a2, a4, a6 = c
        good = sum(1 for u, v in raw
                   if (v*v - u*u*u - a2*u*u - a4*u - a6) % p == 0
                   or (u*u - v*v*v - a2*v*v - a4*v - a6) % p == 0)
        if good == 256:
            curve = c; break
    if curve: break
assert curve, "independent fit found no cubic on all 256 pairs"
a2, a4, a6 = curve
print(f"independent cubic fit: a2,a4,a6 recovered; all 256 pairs on curve = True")
print(f"  matches core.json a2/a4/a6: {str(a2)==core['a2'] and str(a4)==core['a4'] and str(a6)==core['a6']}")

# orient each gated pair onto the curve
def on_aff(x, y): return (y*y - x*x*x - a2*x*x - a4*x - a6) % p == 0
pts = {}
for b, (u, v) in zip(bits, raw):
    pts[b] = (u, v) if on_aff(u, v) else (v, u)
assert all(on_aff(*q) for q in pts.values())

# target
T_raw = (int(H['C1']) % p, int(H['C2']) % p)
print(f"target on fitted curve: {on_aff(*T_raw)}   matches core target: "
      f"{[str(T_raw[0]),str(T_raw[1])]==core['target']}")

# depress x -> x - a2/3 : y^2 = x^3 + A x + B
i3 = pow(3, -1, p)
A = (a4 - a2*a2*i3) % p
Bc = (a6 - a2*a4*i3 + 2*pow(a2, 3, p)*pow(i3, 3, p)) % p
sh = a2 * i3 % p
dep = lambda q: ((q[0] + sh) % p, q[1] % p)
print(f"depressed: A == 0 ? {A == 0}   B matches core: {str(Bc)==core['B']}")

# group law on depressed curve
def add(Pt, Qt):
    if Pt is None: return Qt
    if Qt is None: return Pt
    x1,y1=Pt; x2,y2=Qt
    if x1==x2 and (y1+y2)%p==0: return None
    lam = (3*x1*x1%p*pow(2*y1,-1,p) if Pt==Qt else (y2-y1)*pow(x2-x1,-1,p))%p
    x3=(lam*lam-x1-x2)%p
    return (x3,(lam*(x1-x3)-y1)%p)
def mul(k,Pt):
    R=None
    while k:
        if k&1: R=add(R,Pt)
        Pt=add(Pt,Pt); k>>=1
    return R

n = 0xFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEBAAEDCE6AF48A03BBFD25E8CD0364141
Tdep = dep(T_raw)
depPts = {b: dep(pts[b]) for b in bits}

# doubling chain: head -> ... , verify P_i == 2^i G in the encoder's order
idx = {}
for b,q in depPts.items(): idx.setdefault(q,b)
dbl = {b: add(q,q) for b,q in depPts.items()}
child = {b: idx[d] for b,d in dbl.items() if d in idx}
heads = [b for b in depPts if b not in set(child.values())]
chain=[]; cur=heads[0]
while True:
    chain.append(cur)
    if cur not in child: break
    cur=child[cur]
G = depPts[chain[0]]
chain_ok = all(depPts[chain[i]]==mul(1<<i,G) for i in range(len(chain)))
print(f"doubling chain: heads={len(heads)} len={len(chain)}  P_i==2^i G for all i: {chain_ok}")

# does the encoder's PTS order (instance.py) equal 2^i G ?
import instance
enc_ok = all(instance.PTS[i]==mul(1<<i, instance.G) for i in range(256)) and instance.G==G and instance.T==Tdep
print(f"instance.py PTS[i]==2^i G, G and T match independent recompute: {enc_ok}")
print(f"  n prime & n*G==O & n*T==O: {instance.mul(n,G) is None and instance.mul(n,Tdep) is None}")

# --- solution multiplicity: k vs k+n aliasing ---
gap = (1 << 256) - n
print(f"\nSOLUTION MULTIPLICITY:")
print(f"  2^256 - n = {gap}  (~2^{gap.bit_length()-1})")
print(f"  For a residue r=dlog(T) mod n, #k in [0,2^256) with k==r (mod n) is 2 if r < 2^256-n else 1.")
print(f"  => IF the instance's dlog < {gap} (~2^128), EQUATIONS.txt has TWO valid bit-vectors")
print(f"     (k and k+n), both decoding to T; else exactly one. dlog is unknown (hard),")
print(f"     so this is an instance-dependent non-uniqueness, NOT an encoder soundness break:")
print(f"     both bit-vectors satisfy sum b_i P_i == T faithfully.")
