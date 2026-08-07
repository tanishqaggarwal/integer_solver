#!/usr/bin/env python3
"""agent AF, step 5: constant propagation; classify the 3707 `W - M*u` lift atoms."""
import sys, os, pickle
from collections import Counter, defaultdict
sys.setrecursionlimit(100000)
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from af2_atoms import shape_of, varsof

def main():
    D = pickle.load(open(os.path.join(HERE, 'af_atoms.pkl'), 'rb'))
    E = pickle.load(open(os.path.join(HERE, 'af_defs.pkl'), 'rb'))
    F = pickle.load(open(os.path.join(HERE, 'af_P.pkl'), 'rb'))
    atoms = D['atoms']; defs = E['defs']; Pval = F['Pval']
    parent = F['parent']; cval = F['cval']
    def find(x):
        parent.setdefault(x, x)
        r = x
        while parent[r] != r:
            r = parent[r]
        while parent[x] != r:
            parent[x], x = r, parent[x]
        return r

    # ---------- constant propagation over the definition DAG ----------
    val = {}
    for v, c in cval.items():
        val[v] = c
    def ev(n, depth=0):
        """return int value if the AST node is constant under `val`, else None"""
        if depth > 40:
            return None
        t = n[0]
        if t == 'c':
            return n[1]
        if t == 'v':
            return val.get(find(n[1]))
        if t == 'neg':
            a = ev(n[1], depth+1); return None if a is None else -a
        a = ev(n[1], depth+1)
        if a is None:
            return None
        b = ev(n[2], depth+1)
        if b is None:
            return None
        return a + b if t == '+' else (a - b if t == '-' else a * b)

    for _ in range(12):
        changed = 0
        for v, lst in defs.items():
            r = find(v)
            if r in val:
                continue
            for aid, rhs in lst:
                x = ev(rhs)
                if x is not None:
                    val[r] = x; changed += 1; break
        if not changed:
            break
    print('constant-valued variable classes: %d' % len(val))

    # ---------- the lift atoms ----------
    lift = []   # (atom_id, Rwire, Mvar, uvar, M)
    other = []
    for aid, a in enumerate(atoms):
        if not (a[0] == '-' and a[1][0] == 'v' and a[2][0] == '*'):
            continue
        L, Rr = a[2][1], a[2][2]
        if L[0] != 'v' or Rr[0] != 'v':
            continue
        va, vb = L[1], Rr[1]
        ca, cb = val.get(find(va)), val.get(find(vb))
        cand = None
        if ca is not None and ca % Pval == 0:
            cand = (ca, va, vb)
        elif cb is not None and cb % Pval == 0:
            cand = (cb, vb, va)
        if cand is None:
            continue
        M, mv, uv = cand
        lift.append((aid, a[1][1], mv, uv, M))
    print('lift atoms  W - (M*u) with P | M : %d' % len(lift))
    cc = Counter(M // Pval for (_, _, _, _, M) in lift)
    print('distinct multipliers c = M/P : %d' % len(cc))
    print('  c == 1 : %d      c > 1 : %d' % (cc.get(1, 0), sum(k for c, k in cc.items() if c != 1)))
    print('  top c values:', cc.most_common(12))
    cs = sorted(c for c in cc if c != 1)
    print('  c range: %d .. %d   (%d distinct)' % (cs[0], cs[-1], len(cs)))

    pickle.dump({'val': val, 'lift': lift}, open(os.path.join(HERE, 'af_lift.pkl'), 'wb'), 2)

if __name__ == '__main__':
    main()
