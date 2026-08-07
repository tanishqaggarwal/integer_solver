import os, sys, collections
HERE='/home/user/integer_solver/solve_lab/s10'; LAB='/home/user/integer_solver/solve_lab'
sys.path.insert(0, os.path.join(LAB,'s9','eff')); sys.path.insert(0, HERE)
import lib as L, tools as T, ad
P = 2**256-2**32-977
v=L.load(os.path.join(LAB,'best','new_instance_partial_39026.json'))
for a in (3575, 3577, 3579, 3581, 36602):
    print(f'a{a}: out={L.atom_out.get(a)} neq={len(L.atom2eq.get(a,{}))}')
    print('  SRC:', L.atom_src[a][:400])
    print('  polys:', dict(list(L.polys[a].items())[:10]))
    print()
for u in (26777, 13458, 36602):
    pass
# what defines the inputs of a3575 / a3577
import math
print('gcd(15804267, p) =', math.gcd(15804267, P))
print('15804267 factor?', [d for d in range(2,2000) if 15804267%d==0][:10])
