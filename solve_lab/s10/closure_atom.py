"""S11 step 44: RECURSIVE COMPENSATION CLOSURE on the equation-atom bipartite graph.

An equation with several atoms is satisfied iff its COMBINATION vanishes, so an
atom forced nonzero can be paid for by another atom in the same equation -- and
that atom's own equations can in turn be paid for.  Propagate:

  ACTIVE  = atoms allowed nonzero (seed: the seven residual atoms)
  OBLIG   = equations touching ACTIVE  (all must still vanish)
  need    rank(M[OBLIG, ACTIVE]) < |ACTIVE|, i.e. a kernel with support on the seed

If the closure ever reaches a state with a kernel vector nonzero on the seed, ALL
equations hold and the instance is SOLVED.  If it terminates without one, the
minimum failing count is |OBLIG| - rank, and we learn exactly where it stops.
"""
import os, sys, collections, time
from fractions import Fraction
HERE='/home/user/integer_solver/solve_lab/s10'; LAB='/home/user/integer_solver/solve_lab'
sys.path.insert(0, os.path.join(LAB,'s9','eff')); sys.path.insert(0, HERE)
import lib as L, ad
P = ad.P
SEED = [22229, 22230, 35758, 35759, 35760, 35761, 35762]

# single-atom equations force their atom to zero
forced = set()
for e in range(L.NEQ):
    m, sq, co = L.eq_atoms[e]
    nz = [a for a, c in co.items() if c]
    if len(nz) == 1: forced.add(nz[0])
print(f'atoms forced to zero by a single-atom equation: {len(forced)}')
print(f'  any of the seven forced? {[a for a in SEED if a in forced]}')

ACTIVE = set(SEED)
t0 = time.time()
for rnd in range(12):
    OBLIG = set()
    for a in ACTIVE: OBLIG |= set(L.atom2eq[a])
    # atoms appearing in OBLIG equations that are NOT forced zero -> candidate compensators
    cand = set()
    for e in OBLIG:
        m, sq, co = L.eq_atoms[e]
        for a, c in co.items():
            if c and a not in forced: cand.add(a)
    new = cand - ACTIVE
    print(f'round {rnd}: |ACTIVE| {len(ACTIVE)}  |OBLIG| {len(OBLIG)}  '
          f'new compensators available {len(new)}  ({time.time()-t0:.0f}s)', flush=True)
    if not new:
        print('  closure complete (no further compensators)')
        break
    ACTIVE |= new
    if len(ACTIVE) > 4000:
        print('  ACTIVE exceeded 4000 -- stopping expansion')
        break

OBLIG = set()
for a in ACTIVE: OBLIG |= set(L.atom2eq[a])
rows = sorted(OBLIG); cols = sorted(ACTIVE)
print(f'\nfinal system: {len(rows)} equations x {len(cols)} atoms')
print(f'  rank deficit if full rank: {len(rows) - len(cols)}')
print(f'  a kernel with support on the seed needs rank < {len(cols)}')
# how many of the OBLIG equations are single-atom (unpayable)?
bad = 0
for e in rows:
    m, sq, co = L.eq_atoms[e]
    nz = [a for a, c in co.items() if c]
    if len(nz) == 1 and nz[0] in ACTIVE: bad += 1
print(f'  single-atom equations inside OBLIG whose atom is ACTIVE: {bad}'
      f'   (each forces that atom to zero)')
