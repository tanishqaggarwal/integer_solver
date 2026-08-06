import os, sys, json, time
HERE='/home/user/integer_solver/solve_lab/s10'; LAB='/home/user/integer_solver/solve_lab'
sys.path.insert(0, os.path.join(LAB,'s9','eff')); sys.path.insert(0, HERE)
import lib as L, tools as T, ad
BEST = os.path.join(LAB,'best','new_instance_partial_39026.json')
v = L.load(BEST)
t=time.time(); av=L.all_atom_values(v); print('all_atom_values', time.time()-t)
t=time.time(); f=L.failing_eqs(av); print('failing_eqs', time.time()-t, len(f))
w=list(v); t=time.time(); ad.fwd(w,1); print('fwd 1 round', time.time()-t)
av2=L.all_atom_values(w); print('after fwd score', L.NEQ-len(L.failing_eqs(av2)))
