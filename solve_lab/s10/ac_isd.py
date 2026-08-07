"""S12 step 18: information-set decoding -- how far can the closure's coset
leader actually be pushed?  The score after any linearised move inside the
closure is 39033 - |eqs_of_atoms(D)|, D = supp(Mx - r).  Randomised weighted
information sets bound the frame's mod-p ceiling.
"""
import os, sys, json, time, random, collections
HERE='/home/user/integer_solver/solve_lab/s10'; LAB='/home/user/integer_solver/solve_lab'
sys.path.insert(0, os.path.join(LAB,'s9','eff')); sys.path.insert(0, HERE)
import lib as L, tools as T, ad
from ac_sacr import closure_of, decode
P = ad.P
v = L.load(os.path.join(HERE,'mod9118_0.json'))
av = L.all_atom_values(v)
BAD = [21617, 29539]
t0 = time.time()
rows, Us, cols = closure_of(v, BAD)
print(f'closure {len(rows)} x {len(Us)} ({time.time()-t0:.0f}s)', flush=True)
w = {c: len(L.atom2eq.get(c, {})) for c in rows}
random.seed(1234)
best = (10**9, None)
t0 = time.time(); N = 0
while time.time() - t0 < 900:
    N += 1
    T_ = random.choice([0.0, 0.3, 1.0, 3.0])
    o = sorted(rows, key=lambda c: -(w[c] + random.gauss(0, T_ * 4)))
    x, D, rank, nbad = decode(rows, Us, cols, av, o)
    e = len(L.eqs_of_atoms(D))
    if e < best[0]:
        best = (e, D)
        print(f'  trial {N} T={T_}: |D|={len(D)} cost {e} equations -> '
              f'score {L.NEQ-e}  D={D}  ({time.time()-t0:.0f}s)', flush=True)
print(f'\n{N} information sets tried; BEST equation cost {best[0]} '
      f'-> mod-p ceiling for this frame: {L.NEQ-best[0]}')
print(f'  best violated set D = {best[1]}')
json.dump({'cost': best[0], 'D': best[1], 'trials': N},
          open(os.path.join(HERE,'ac_isd.json'),'w'))
