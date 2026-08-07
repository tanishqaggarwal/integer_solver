#!/usr/bin/env python3
"""agent AF, step 8: the exact condition list  c*P | Expr, and validation vs the deliverable."""
import sys, os, pickle, json
from collections import Counter, defaultdict
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from af6_expand import (atoms, defs, defc, val, lift, Pval, find, pp, expand, varsof, shape_of)

real = [t for t in lift if t[4] == Pval]
Rof = {}                       # slack class -> (defining atom, u class)
for (aid, Rw, mv, uv, M) in real:
    Rof[find(Rw)] = (aid, find(uv))

v2a = defaultdict(list)
for aid, a in enumerate(atoms):
    for v in varsof(a, set()):
        v2a[find(v)].append(aid)

def coef_of_R(a, R):
    """If atom `a` is  E  -/+  c*R  (R appearing exactly once, linearly), return (c, E_ast, sign).
    Representation: a == E - c*R  after normalisation. Returns None if not that shape."""
    # a = ('-', X, Y) or ('+', X, Y)
    if a[0] not in ('+', '-'):
        return None
    op, X, Y = a
    def islinR(n):
        # returns c if n == c*R or R, else None
        if n == ('v', R):
            return 1
        if n[0] == 'v' and find(n[1]) == R:
            return 1
        if n[0] == '*':
            for u, w in ((n[1], n[2]), (n[2], n[1])):
                if w[0] == 'v' and find(w[1]) == R:
                    from af1_parse import is_const
                    k = is_const(u)
                    if k is not None:
                        return k
        return None
    cy = islinR(Y)
    if cy is not None and R not in varsof(X, set()) | set(find(z) for z in varsof(X, set())):
        pass
    if cy is not None:
        Xv = set(find(z) for z in varsof(X, set()))
        if R not in Xv:
            return (cy if op == '-' else -cy, X)
    cx = islinR(X)
    if cx is not None:
        Yv = set(find(z) for z in varsof(Y, set()))
        if R not in Yv:
            # a = c*R (op) Y  ->  -(Y (op') c*R) ... normalise to E - c'*R  with E = -/+Y
            if op == '-':      # cR - Y = -(Y - cR)
                return (cx, ('neg', Y))
            else:              # cR + Y  ->  Y - (-c)R
                return (-cx, Y)
    return None

conds = []      # (Rclass, c, Expr_ast, use_atom_id, u_class)
unresolved = []
for R, (daid, uc) in Rof.items():
    cands = []
    for aid in v2a[R]:
        if aid == daid:
            continue
        a = atoms[aid]
        if a[0] == 'v':
            continue
        r = coef_of_R(a, R)
        if r is not None:
            cands.append((aid, r[0], r[1]))
    if not cands:
        unresolved.append(R); continue
    # prefer the smallest atom (the local consumer), not the giant composite equation bodies
    cands.sort(key=lambda t: len(str(atoms[t[0]])))
    aid, c, Ex = cands[0]
    conds.append((R, abs(c), Ex, aid, uc))

print('slack wires with an identified consumer: %d   unresolved: %d' % (len(conds), len(unresolved)))
cc = Counter(c for (_, c, _, _, _) in conds)
n1 = cc.get(1, 0)
print('  c == 1 : %d      c > 1 : %d' % (n1, len(conds) - n1))
print('  distinct c>1 values: %d   min %d  max %d' % (
    len([c for c in cc if c != 1]),
    min([c for c in cc if c != 1]), max([c for c in cc if c != 1])))
shc = Counter(shape_of(atoms[aid]) for (_, c, _, aid, _) in conds if c > 1)
print('\nshapes of the 927-family consumer atoms:')
for s, k in shc.most_common(12):
    print('   %5d  %s' % (k, s))

# ---------------- validation against the deliverable ----------------
ASG = os.path.join(HERE, '..', 'best', 'new_instance_partial_39026.json')
d = json.load(open(ASG))
asg = {}
for k, v in d.items():
    asg[int(k[2:]) if k.startswith('x_') else int(k)] = int(v)
print('\ndeliverable assigns %d variables' % len(asg))

# P aliases
Pcls = find(26064)
pal = [v for v in range(38748) if find(v) == Pcls]
got = [asg[v] for v in pal if v in asg]
print('P-alias classes: %d vars ; assigned in deliverable: %d ; all == P : %s' %
      (len(pal), len(got), all(g == Pval for g in got)))

def evalast(n, A):
    t = n[0]
    if t == 'c':
        return n[1]
    if t == 'v':
        r = find(n[1])
        if r in val:
            return val[r]
        return A.get(n[1], A.get(r, 0))
    if t == 'neg':
        return -evalast(n[1], A)
    a = evalast(n[1], A); b = evalast(n[2], A)
    return a + b if t == '+' else (a - b if t == '-' else a * b)

A = defaultdict(int)
for v, x in asg.items():
    A[v] = x
    A[find(v)] = x
nz = 0
for (aid, Rw, mv, uv, M) in real:
    if evalast(atoms[aid], A) != 0:
        nz += 1
print('lift-defining atoms  R - P*u  nonzero on the deliverable: %d / %d' % (nz, len(real)))
nz2 = 0
for (R, c, Ex, aid, uc) in conds:
    if evalast(atoms[aid], A) != 0:
        nz2 += 1
print('consumer atoms nonzero on the deliverable: %d / %d' % (nz2, len(conds)))

pickle.dump({'conds': conds, 'Rof': Rof, 'unresolved': unresolved},
            open(os.path.join(HERE, 'af_cond.pkl'), 'wb'), 2)
