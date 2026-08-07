import os, sys, collections, json
HERE='/home/user/integer_solver/solve_lab/s10'; LAB='/home/user/integer_solver/solve_lab'
sys.path.insert(0, os.path.join(LAB,'s9','eff')); sys.path.insert(0, HERE)
import lib as L
h=collections.Counter(len(co) for m,sq,co in L.eq_atoms)
print('eq size hist:', sorted(h.items()))
print('total eqs', sum(h.values()))
S=[22229, 22230, 35758, 35759, 35760, 35761, 35762]
VARS=set()
for a in S: VARS|=set(L.avars[a])
print('\nvars of S:', sorted(VARS))
v = L.load(os.path.join(LAB,'best','new_instance_partial_39026.json'))
P=2**256-2**32-977
for u in sorted(VARS):
    d=L.definer.get(u)
    others=[a for a in L.var_atoms[u] if a not in S]
    print(f'x_{u:<6} val={str(v[u])[:26]:<28} definer=a{d} ({L.atom_src[d][:42] if d is not None else "FREE"})')
    print(f'        other atoms: {[(a,len(L.atom2eq[a])) for a in others]}')
