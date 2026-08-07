#!/usr/bin/env python3
"""agent AF, step 20: re-derive the block law shape, and how the two free output wires enter it."""
import sys, os, pickle, itertools
from math import gcd
from collections import Counter, defaultdict
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from af6_expand import atoms, defc, val, find, pp, expand, varsof, shape_of, Pval
from af1_parse import is_const
C = pickle.load(open(os.path.join(HERE, 'af_cond.pkl'), 'rb'))
M = pickle.load(open(os.path.join(HERE, 'af_map.pkl'), 'rb'))
LAM = pickle.load(open(os.path.join(HERE, 'af_lam.pkl'), 'rb'))
conds = C['conds']; info = M['info']; pure = LAM['pure']

def raw(n):
    if n[0] == 'c':
        return str(n[1]) if abs(n[1]) < 10**14 else ('P' if n[1] == Pval else 'BIG')
    if n[0] == 'v':
        return 'x%d' % find(n[1])
    if n[0] == 'neg':
        return '-' + raw(n[1])
    return '(%s %s %s)' % (raw(n[1]), n[0], raw(n[2]))

def term(t):
    if t[0] == '*':
        k = is_const(t[1])
        if k is not None:
            return (k, t[2])
        k = is_const(t[2])
        if k is not None:
            return (k, t[1])
    if t[0] == 'neg':
        k, x = term(t[1]); return (-k, x)
    return (1, t)

# ---------- 1. multivariate polynomial expansion of the residual in terms of leaves ----------
# sparse polynomial: dict{ monomial (sorted tuple of var ids with repetition) : coeff }
def pmul(A, Bp):
    R = defaultdict(int)
    for ma, ca in A.items():
        for mb, cb in Bp.items():
            R[tuple(sorted(ma + mb))] += ca * cb
    return {k: v for k, v in R.items() if v}

def padd(A, Bp, s=1):
    R = defaultdict(int, A)
    for m, c in Bp.items():
        R[m] += s * c
    return {k: v for k, v in R.items() if v}

def poly(n, stop, depth=0, memo=None):
    if memo is None:
        memo = {}
    t = n[0]
    if t == 'c':
        return {(): n[1]}
    if t == 'v':
        r = find(n[1])
        if r in val:
            return {(): val[r]}
        if r in stop or depth > 60:
            return {(r,): 1}
        dl = defc.get(r)
        if not dl or len(dl) != 1:
            return {(r,): 1}
        if r in memo:
            return memo[r]
        memo[r] = {(r,): 1}
        res = poly(dl[0][1], stop, depth + 1, memo)
        memo[r] = res
        return res
    if t == 'neg':
        return {m: -c for m, c in poly(n[1], stop, depth, memo).items()}
    a = poly(n[1], stop, depth, memo); b = poly(n[2], stop, depth, memo)
    if t == '*':
        return pmul(a, b)
    return padd(a, b, 1 if t == '+' else -1)

# pick a small merge block whose two children are leaves
cand = sorted(pure.items(), key=lambda kv: len(kv[1][0]) + len(kv[1][1]))
for g, (I, J) in cand[:1]:
    print('=== block gate x%d  I=%s J=%s' % (g, sorted(I), sorted(J)))
    ops = [(info[i]['other'], conds[i][1]) for i in range(len(conds))
           if info[i]['cls'] == 'offpin' and info[i]['gate'] == g]
    print('  free chord-output wires (off-pin wires):', ops)
    outs = set(w for w, c in ops)
    for i in range(len(conds)):
        if info[i]['cls'] != 'cong' or info[i]['gate'] != g:
            continue
        c = conds[i][1]
        rhs = defc[info[i]['other']][0][1]
        (a1, D1), (a2, D2) = term(rhs[1]), term(rhs[2])
        print('  row c=%-12d  alpha=%d beta=%d' % (c, a1, a2))
        for nm, D in (('N1?', D1), ('N2?', D2)):
            p = poly(D, outs)
            deg = Counter(len(m) for m in p)
            # degree in each free output wire
            per = {}
            for w in outs:
                per[w] = max((m.count(w) for m in p), default=0)
            print('     %s : %d monomials, total-degree profile %s, degree in free outs %s'
                  % (nm, len(p), dict(sorted(deg.items())), {('x%d' % k): v for k, v in per.items()}))
            # the coefficient of each free output wire (as a polynomial in the rest)
            for w in outs:
                lin = {m: cc for m, cc in p.items() if m.count(w) == 1}
                if lin:
                    red = {tuple(x for x in m if x != w): cc for m, cc in lin.items()}
                    print('        d/d x%d  = %d monomials, degrees %s'
                          % (w, len(red), sorted(set(len(m) for m in red))))
        break
    break

# ---------- 2. the (alpha,beta) triples per block ----------
byblk = defaultdict(list)
for i in range(len(conds)):
    if info[i]['cls'] != 'cong':
        continue
    rhs = defc[info[i]['other']][0][1]
    (a1, D1), (a2, D2) = term(rhs[1]), term(rhs[2])
    byblk[info[i]['gate']].append((conds[i][1], a1, a2, D1, D2))
print('\n(alpha,beta) triples, 3 example blocks:')
for g in list(byblk)[:3]:
    print('  block x%d:' % g, [(a1, a2) for (c, a1, a2, D1, D2) in byblk[g]])
    for (c, a1, a2, D1, D2) in byblk[g]:
        print('     c=%-12d  %d*[%s]  +  %d*[%s]' % (c, a1, raw(D1)[:60], a2, raw(D2)[:60]))
