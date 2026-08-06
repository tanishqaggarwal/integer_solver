"""S10 step 14: which knob is NOT free?  Move each one alone (with the residual
atoms blocked) and report the collateral atoms it breaks."""
import os, sys
HERE = os.path.dirname(os.path.abspath(__file__))
LAB = os.path.dirname(HERE)
sys.path.insert(0, os.path.join(LAB, 's9', 'eff'))
sys.path.insert(0, HERE)
import lib as L

P = 2**256 - 2**32 - 977
NZ = [22229, 22230, 35758, 35759, 35760, 35761, 35762]
BLOCK = set(NZ) | {22231}
BEST = os.path.join(LAB, 'best', 'new_instance_partial_39026.json')
base = L.load(BEST)
base_av = L.all_atom_values(base)
BASE_NZ = set(a for a in range(L.NA) if base_av[a])

KNOBS = [642, 17325, 9413, 1329, 10903, 29854, 31864, 9118, 8731, 28730, 4432]
for u in KNOBS:
    for tag, nv in (('+1', base[u] + 1), ('=0', 0)):
        v = list(base)
        seeds = {u: nv}
        if u == 28730:                      # keep a22231 = 0 by hand
            seeds[4432] = base[19964] + nv
        ch, _ = L.ripple(v, seeds, block=BLOCK)
        av = L.all_atom_values(v)
        nz = set(a for a in range(L.NA) if av[a])
        extra = sorted(nz - BASE_NZ)
        fail = L.failing_eqs(av)
        print(f'x_{u:<6} {tag:<4} changed={len(ch):<4} new_nonzero_atoms={extra} '
              f'failing={len(fail)}')
