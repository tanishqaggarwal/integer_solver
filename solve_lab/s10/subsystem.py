"""S10 step 2: the residual sub-system, exactly.

The 7 nonzero atoms appear in exactly 12 equations. Ask:
  (a) what is the rank of the 12x7 coefficient matrix?  (does a=0 get forced?)
  (b) which variables feed those atoms, are they free inputs, and where else do they occur?
"""
import os, sys, json, collections
from fractions import Fraction

HERE = os.path.dirname(os.path.abspath(__file__))
LAB = os.path.dirname(HERE)
sys.path.insert(0, os.path.join(LAB, 's9', 'eff'))
import lib as L

P = 2**256 - 2**32 - 977
BEST = os.path.join(LAB, 'best', 'new_instance_partial_39026.json')
v = L.load(BEST)
av = L.all_atom_values(v)
NZ = [a for a in range(L.NA) if av[a] != 0]
EQS = sorted(L.eqs_of_atoms(NZ))
fail = set(L.failing_eqs(av))
print(f'nonzero atoms {NZ}')
print(f'equations touching them: {len(EQS)} {EQS}')
print(f'  of which failing: {sorted(set(EQS)&fail)}')
print(f'  of which satisfied: {sorted(set(EQS)-fail)}')
print()

# --- (a) the 12 x 7 matrix over Q -------------------------------------------
M = [[L.eq_atoms[i][2].get(a, 0) for a in NZ] for i in EQS]
print('coefficient matrix (rows=eqs, cols=atoms %s):' % NZ)
for i, r in enumerate(M):
    print(f'  eq {EQS[i]:>6} sq={int(L.eq_atoms[EQS[i]][1])} : {r}')

def rank_q(mat):
    m = [[Fraction(x) for x in row] for row in mat]
    rows, cols = len(m), len(m[0]); r = 0
    for c in range(cols):
        piv = next((i for i in range(r, rows) if m[i][c] != 0), None)
        if piv is None: continue
        m[r], m[piv] = m[piv], m[r]
        pv = m[r][c]
        m[r] = [x / pv for x in m[r]]
        for i in range(rows):
            if i != r and m[i][c] != 0:
                f = m[i][c]
                m[i] = [m[i][j] - f * m[r][j] for j in range(cols)]
        r += 1
        if r == rows: break
    return r

print(f'\nrank over Q = {rank_q(M)} of {len(NZ)} atom-columns')
print('  -> if rank == #atoms, the ONLY way to satisfy all 12 with no other atom')
print('     becoming nonzero is to zero all 7 atoms simultaneously.')

# --- (b) variables of the nonzero atoms -------------------------------------
print('\n=== variables feeding the nonzero atoms ===')
gvars = set()
for a in NZ:
    gvars |= L.avars[a]
definer = L.definer
for u in sorted(gvars):
    kind = 'FREE' if u not in definer else f'gate(def by atom {definer[u]})'
    # all atoms mentioning u, and whether they are currently zero
    ats = L.var_atoms[u]
    nzats = [x for x in ats if av[x] != 0]
    # equations effectively reached (syntactic)
    eqs = L.var_eqs[u]
    print(f'x_{u:<6} val={str(v[u])[:28]:<30} {kind:<26} '
          f'#atoms={len(ats):<4} #nz_atoms={len(nzats):<3} #eqs={len(eqs)}')
    if v[u] == P: print('        ^^ equals p')

# --- (c) which of them are p, which are free handles ------------------------
print('\n=== atom sources ===')
for a in NZ:
    print(f'  a{a}: {L.atom_src[a]}   = {av[a] % P} (mod p)')
