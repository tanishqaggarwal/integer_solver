#!/usr/bin/env python3
"""Look for an EXACT global certificate:  sum_e y_e * body_e  ==  c * alpha  identically.
E' is closed under (eq -> its atoms -> equations containing them)?  NO -- we only need
that every atom of every equation in E' is a COLUMN of the system, which is automatic
if columns = union of atoms of E'.  Then y^T C = c*e_alpha is a genuine identity over
ALL atoms, hence alpha = 0 in every solution."""
import pickle, sys
from fractions import Fraction
from collections import defaultdict

D = pickle.load(open('_om_parsed2.pkl', 'rb'))
eqatoms = D['eqatoms']
ainc = defaultdict(list)
for e, d in enumerate(eqatoms):
    for k in d: ainc[k].append(e)
ALPHA = '((x7068-x2099)-(7376877*x642))'
BETA = '((x4432-x19964)-x28730)'
GAM = '(x28730-(x17499*x9413))'
DEL = '(x642-(x28599*x17325))'

def solve_cert(Eset, target):
    """solve y^T C = e_target over columns = union atoms(E').  Returns y dict or None."""
    E = sorted(Eset)
    cols = set()
    for e in E: cols |= set(eqatoms[e])
    cols = sorted(cols)
    if target not in cols: return None
    ci = {k: i for i, k in enumerate(cols)}
    # system: for each column a:  sum_e y_e C[e,a] = [a==target]
    rows = []
    rhs = []
    for a in cols:
        rows.append([Fraction(eqatoms[e].get(a, 0)) for e in E])
        rhs.append(Fraction(1 if a == target else 0))
    nr = len(rows); nc = len(E)
    piv = []
    r = 0
    for c in range(nc):
        s = None
        for i in range(r, nr):
            if rows[i][c] != 0: s = i; break
        if s is None: continue
        rows[r], rows[s] = rows[s], rows[r]; rhs[r], rhs[s] = rhs[s], rhs[r]
        pv = rows[r][c]
        rows[r] = [x / pv for x in rows[r]]; rhs[r] = rhs[r] / pv
        for i in range(nr):
            if i != r and rows[i][c] != 0:
                f = rows[i][c]
                rows[i] = [a2 - f * b for a2, b in zip(rows[i], rows[r])]
                rhs[i] = rhs[i] - f * rhs[r]
        piv.append(c); r += 1
        if r == nr: break
    for i in range(r, nr):
        if rhs[i] != 0: return None   # inconsistent
    y = [Fraction(0)] * nc
    for i, c in enumerate(piv): y[c] = rhs[i]
    return {E[i]: y[i] for i in range(nc) if y[i] != 0}

def grow(seed_eqs, steps):
    E = set(seed_eqs)
    for _ in range(steps):
        A = set()
        for e in E: A |= set(eqatoms[e])
        newE = set(E)
        for a in A: newE |= set(ainc[a])
        if newE == E: break
        E = newE
    return E

if __name__ == '__main__':
    for tname, target in [('ALPHA', ALPHA), ('BETA', BETA), ('GAMMA', GAM), ('DELTA', DEL)]:
        found = None
        E = set(ainc[target])
        for step in range(6):
            y = solve_cert(E, target)
            if y is not None:
                found = (step, y); break
            E = grow(E, 1)
            if len(E) > 4000: break
        if found:
            step, y = found
            print('%-6s: CERTIFIED zero.  |E|=%d  support=%d equations' % (tname, len(E), len(y)))
            items = sorted(y.items())
            print('        y =', {k: str(v) for k, v in items[:12]})
        else:
            print('%-6s: no certificate found up to |E|=%d' % (tname, len(E)))
