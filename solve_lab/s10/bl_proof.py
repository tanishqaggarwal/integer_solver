"""bl_proof: empirical check of the structural cap.

Claim: in the witness frame the seven residual atoms depend only on 39 variables
whose free inputs contain no boolean except x_2081/x_4287.  So ANY set of
non-MUX boolean flips leaves all seven atom VALUES bit-identical.
Test with mass flips (10, 50, 200, all neutral, all 1154).
"""
import os, sys, json, time, random
HERE='/home/user/integer_solver/solve_lab/s10'; LAB='/home/user/integer_solver/solve_lab'
sys.path.insert(0, os.path.join(LAB,'s9','eff')); sys.path.insert(0, HERE)
import lib as L, tools as T, ad
from bl_core import BOOL, CANON, F2, pot, cheap, FORBID
random.seed(7)
SEVEN=[22229,22230,35758,35759,35760,35761,35762]
w = L.load(os.path.join(LAB,'best','new_instance_partial_39026.json')); F2.fwd(w)
av0 = L.all_atom_values(w)
ref = [av0[a] for a in SEVEN]
BF = sorted((BOOL & CANON.FREE) - FORBID)
S = json.load(open(os.path.join(HERE,'bl_single_f2.json')))
neutral = sorted(u for sc, nz, u, ic in S if sc == 39026)
print(f'free booleans {len(BF)}; score-neutral in the witness frame: {len(neutral)}', flush=True)
for tag, flips in (('rand10', random.sample(BF, 10)), ('rand50', random.sample(BF, 50)),
                   ('rand200', random.sample(BF, 200)), ('neutral-all', neutral),
                   ('ALL-1154', BF),
                   ('neutral-rand-half', random.sample(neutral, len(neutral)//2))):
    v = list(w)
    for u in flips: v[u] = 1 - v[u]
    F2.fwd(v, rounds=8)
    av = L.all_atom_values(v)
    same = all(av[a] == r for a, r in zip(SEVEN, ref))
    p, _, nz = pot(v)
    print(f'  {tag:<18} ({len(flips):>4} flips): score {p[0]}  nz {len(nz)}  '
          f'seven-atom values identical: {same}', flush=True)
