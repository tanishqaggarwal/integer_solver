"""S11 step 52: the COMBINED closure.

Compensation needs atoms that share an EQUATION; realisability needs atoms that
share a detached VARIABLE.  Close over both and see whether the fixed point still
has a kernel touching the seed -- and how big it gets.
"""
import os, sys, collections, time
HERE='/home/user/integer_solver/solve_lab/s10'; LAB='/home/user/integer_solver/solve_lab'
sys.path.insert(0, os.path.join(LAB,'s9','eff')); sys.path.insert(0, HERE)
import lib as L, ad
P = ad.P
SEED = [22229, 22230, 35758, 35759, 35760, 35761, 35762]
forced = set()
for e in range(L.NEQ):
    m, sq, co = L.eq_atoms[e]
    nz = [a for a, c in co.items() if c]
    if len(nz) == 1: forced.add(nz[0])
print(f'atoms forced to zero by single-atom equations: {len(forced)}')

ACTIVE = set(SEED)
t0 = time.time()
for rnd in range(30):
    prev = len(ACTIVE)
    # (a) compensation: atoms sharing an equation with an ACTIVE atom
    OB = set()
    for a in ACTIVE: OB |= set(L.atom2eq[a])
    for e in OB:
        m, sq, co = L.eq_atoms[e]
        for a, c in co.items():
            if c and a not in forced: ACTIVE.add(a)
    # (b) realisability: detach each ACTIVE gate atom's output variable; any atom
    #     sharing that variable must also be allowed nonzero
    for a in list(ACTIVE):
        ov = L.atom_out.get(a)
        if ov is None: continue
        for b in L.var_atoms[ov[1]]:
            if b not in forced: ACTIVE.add(b)
    if len(ACTIVE) == prev:
        print(f'round {rnd}: FIXED POINT at {len(ACTIVE)} atoms '
              f'({time.time()-t0:.0f}s)'); break
    if rnd % 3 == 0 or len(ACTIVE) > 8000:
        print(f'round {rnd}: |ACTIVE| {len(ACTIVE)}  ({time.time()-t0:.0f}s)',
              flush=True)
    if len(ACTIVE) > 12000:
        print('  exceeded 12000 -- the combined closure blows up'); break
OB = set()
for a in ACTIVE: OB |= set(L.atom2eq[a])
print(f'\ncombined closure: {len(ACTIVE)} atoms, {len(OB)} equations')
print(f'  atoms - equations = {len(ACTIVE) - len(OB)}  '
      f'(positive means a kernel is guaranteed)')
print(f'  as a fraction of the instance: atoms {100*len(ACTIVE)/L.NA:.1f}%, '
      f'equations {100*len(OB)/L.NEQ:.1f}%')
forced_in = [a for a in ACTIVE if a in forced]
print(f'  forced-zero atoms inside ACTIVE: {len(forced_in)}')
