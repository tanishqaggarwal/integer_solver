"""W: WHY is 9123 essential to the lattice but absent from U's route discount?
For each of the six essential rows, extract the cocircuit direction (the knob delta that breaks
ONLY that row among all 198 SAT rows) and characterise it."""
import sys, os, json
from fractions import Fraction
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
sys.set_int_max_str_digits(20_000_000)
import w_setup2 as S, frameB
from collections import Counter
KN = S.KNOB; N = len(KN)
ESS = [2554, 6816, 8124, 9123, 9421, 'S']
V = {e: [Fraction(S.rows[e].get(u, 0)) for u in KN] for e in S.names}

def nullspace(rows):
    """basis of {x : rows . x = 0} in Q^N"""
    M = [r[:] for r in rows]; piv = []; r = 0
    for c in range(N):
        p = next((i for i in range(r, len(M)) if M[i][c] != 0), None)
        if p is None: continue
        M[r], M[p] = M[p], M[r]
        d = M[r][c]; M[r] = [x/d for x in M[r]]
        for i in range(len(M)):
            if i != r and M[i][c] != 0:
                f = M[i][c]; M[i] = [a-f*b for a, b in zip(M[i], M[r])]
        piv.append(c); r += 1
        if r == len(M): break
    free = [c for c in range(N) if c not in piv]
    basis = []
    for f in free:
        x = [Fraction(0)]*N; x[f] = Fraction(1)
        for i, c in enumerate(piv): x[c] = -M[i][f]
        basis.append(x)
    return basis

OFFPIN_EQ = {6816, 8124, 9123, 9421, 12231, 12270, 12350, 14584, 18673}
U_FIVE = {2554, 6816, 8124, 8680, 9421}
print('block-7181 off-pin equations :', sorted(OFFPIN_EQ))
print("U's five discount equations  :", sorted(U_FIVE))
print('my six essential rows        :', ESS, '  (S = eq8680)')
print('  in BOTH            :', sorted({6816, 8124, 9421}))
print('  essential, not U   :', 9123, '  <-- the target')
print('  U, not off-pin site:', sorted({2554, 8680}))
print()
for e in ESS:
    others = [x for x in S.SAT if x != e]
    ns = nullspace([V[x] for x in others])
    # a direction in ns that actually breaks e
    d = next((x for x in ns if sum(a*b for a, b in zip(V[e], x)) != 0), None)
    if d is None:
        print('  %-6s : NO direction breaks it alone (not essential?)' % str(e)); continue
    sup = sorted(KN[i] for i in range(N) if d[i] != 0)
    ee = 8680 if e == 'S' else e
    natoms = len(frameB.eq_terms[ee][2])
    print('  row %-6s : cocircuit direction support = %2d knobs %s'
          % (str(e), len(sup), sup))
    print('             eq has %d atoms ; in the 7181 off-pin set: %s ; in U\'s five: %s'
          % (natoms, ee in OFFPIN_EQ, ee in U_FIVE))
